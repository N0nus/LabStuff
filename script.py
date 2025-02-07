import pandas as pd
import os
from typing import Dict, Literal, Tuple
from enum import Enum
from rdkit import Chem
from rdkit.Chem import AllChem
import requests
import time
import csv

class OutputFormat(Enum):
    GAUSSIAN = "gaussian"
    AVOGADRO = "avogadro"

class InputType(Enum):
    FROM_CSV = "csv"
    DIRECT_SMILES = "smiles"

def fetch_smiles_from_pubchem(chemical_name: str) -> Tuple[bool, str]:
    """
    Fetch SMILES notation from PubChem API for a given chemical name.
    Returns (success, smiles_or_error_message)
    """
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chemical_name}/property/CanonicalSMILES/TXT"
        response = requests.get(url)
        
        if response.status_code == 200:
            return True, response.text.strip()
        else:
            return False, f"PubChem API returned status code: {response.status_code}"
            
    except Exception as e:
        return False, str(e)

def add_to_pubchemicals_csv(chemical_name: str, smiles: str):
    """Add a new chemical and its SMILES notation to PubChemicals.csv"""
    file_exists = os.path.isfile('PubChemicals.csv')
    
    with open('PubChemicals.csv', mode='a', newline='') as file:
        writer = csv.writer(file)
        if not file_exists:
            writer.writerow(['Chemical Name', 'SMILES Notation'])
        writer.writerow([chemical_name, smiles])

def load_chemicals(csv_path: str) -> Dict[str, str]:
    """
    Load chemicals from CSV file. If a chemical isn't found,
    try to fetch it from PubChem and add it to the CSV.
    """
    # Create CSV if it doesn't exist
    if not os.path.exists(csv_path):
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Chemical Name', 'SMILES Notation'])
    
    df = pd.read_csv(csv_path)
    chemicals = {row['Chemical Name']: row['SMILES'] for _, row in df.iterrows()}
    return chemicals

def get_chemical_smiles(chemical_name: str, csv_path: str = "PubChemicals.csv") -> str:
    """
    Get SMILES for a chemical, trying CSV first then PubChem if not found.
    Raises ValueError if chemical cannot be found.
    """
    chemicals = load_chemicals(csv_path)
    
    if chemical_name in chemicals:
        return chemicals[chemical_name]
    
    # Try PubChem
    print(f"Chemical '{chemical_name}' not found in CSV, checking PubChem...")
    success, result = fetch_smiles_from_pubchem(chemical_name)
    
    if success:
        print(f"Found '{chemical_name}' in PubChem, adding to CSV...")
        add_to_pubchemicals_csv(chemical_name, result)
        return result
    else:
        raise ValueError(f"Could not find chemical '{chemical_name}' in CSV or PubChem: {result}")

def build_polymer_smiles(repeat_unit: str, n_units: int, end_groups: Dict[str, str], repeat_units: Dict[str, str]) -> str:
    """
    Build polymer SMILES with more robust connection and representation.
    
    Args:
        repeat_unit (str): Name of the repeat unit
        n_units (int): Number of repeat units
        end_groups (dict): Start and end groups for the polymer
        repeat_units (dict): Dictionary of known repeat units with connection points
    
    Returns:
        str: SMILES representation of the polymer
    """
    if repeat_unit not in repeat_units:
        raise ValueError(f"Unknown repeat unit: {repeat_unit}")
    
    # Get the base repeat unit SMILES, ensuring connection points are preserved
    base_unit = repeat_units[repeat_unit]
    
    # Identify connection points (marked with *)
    connection_points = base_unit.count('*')
    if connection_points != 2:
        raise ValueError(f"Repeat unit must have exactly 2 connection points, found {connection_points}")
    
    # Remove connection point markers for base unit
    base_cleaned = base_unit.replace('*', '')
    
    # Start group handling
    start_group = end_groups.get('start', '')
    end_group = end_groups.get('end', '')
    
    # Polymer construction
    if n_units == 1:
        # Single unit case
        polymer = (f"{start_group or ''}{base_cleaned}{end_group or ''}")
    else:
        # Multiple units with linear connection
        # Use a linking pattern that connects the units
        polymer_units = [base_cleaned] * n_units
        polymer = start_group + ''.join(polymer_units) + end_group
    
    return polymer.strip()

def smiles_to_3d_coords(smiles: str, add_hydrogens: bool = True) -> str:
    """Convert SMILES to 3D coordinates using RDKit."""
    try:
        # Create RDKit molecule object
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Failed to parse SMILES string")

        # Add hydrogens if requested
        if add_hydrogens:
            mol = Chem.AddHs(mol)

        # Generate 3D conformation
        success = AllChem.EmbedMolecule(mol, randomSeed=42)
        if success == -1:
            # Try with different parameters if initial embedding fails
            success = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
            if success == -1:
                raise ValueError("3D embedding failed")

        # Optimize the structure using MMFF94
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)

        # Convert to XYZ format
        num_atoms = mol.GetNumAtoms()
        xyz_lines = [str(num_atoms), "Generated by RDKit"]
        
        conf = mol.GetConformer()
        for i in range(num_atoms):
            pos = conf.GetAtomPosition(i)
            atom_symbol = mol.GetAtomWithIdx(i).GetSymbol()
            xyz_lines.append(f"{atom_symbol:<2} {pos.x:>10.4f} {pos.y:>10.4f} {pos.z:>10.4f}")

        return "\n".join(xyz_lines)

    except Exception as e:
        print(f"Primary 3D conversion failed: {str(e)}")
        print("Attempting emergency structure recovery...")

        try:
            # Alternative embedding parameters
            mol = Chem.MolFromSmiles(smiles)
            if add_hydrogens:
                mol = Chem.AddHs(mol)
            
            # Try different conformer generation parameters
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            success = AllChem.EmbedMolecule(mol, params)
            
            if success == -1:
                raise ValueError("Emergency structure recovery failed")

            # Optimize with more iterations
            AllChem.MMFFOptimizeMolecule(mol, maxIters=5000)

            # Convert to XYZ format
            num_atoms = mol.GetNumAtoms()
            xyz_lines = [str(num_atoms), "Generated by RDKit (emergency recovery)"]
            
            conf = mol.GetConformer()
            for i in range(num_atoms):
                pos = conf.GetAtomPosition(i)
                atom_symbol = mol.GetAtomWithIdx(i).GetSymbol()
                xyz_lines.append(f"{atom_symbol:<2} {pos.x:>10.4f} {pos.y:>10.4f} {pos.z:>10.4f}")

            return "\n".join(xyz_lines)

        except Exception as fallback_error:
            raise RuntimeError(f"All conversion attempts failed: {str(fallback_error)}")

def generate_output_file(molecule_name: str, config: Dict, xyz_coords: str, output_format: OutputFormat):
    """Generate output file in specified format."""
    os.makedirs("outputs2", exist_ok=True)
    
    if output_format == OutputFormat.GAUSSIAN:
        filename = f"outputs2/{molecule_name}.com"
        input_content = [
            f"%chk={molecule_name}.chk",
            f"%mem={config.get('memory', '4GB')}",
            f"%nproc={config.get('nproc', 2)}",
            f"#P {config['method']}/{config['basis_set']} "
            f"{'opt' if config['optimization'] else ''} "
            f"{'freq' if config['frequency'] else ''} "
            f"{config.get('additional_keywords', '')}",
            "",
            config["title"],
            "",
            f"{config['charge']} {config['multiplicity']}",
            *xyz_coords.split('\n')[2:],  # Skip first two lines of XYZ
            ""
        ]
        content = '\n'.join(input_content)
    
    elif output_format == OutputFormat.AVOGADRO:
        filename = f"outputs2/{molecule_name}.xyz"
        content = xyz_coords
    
    with open(filename, 'w') as f:
        f.write(content)
    
    print(f"Generated {filename}")
    return filename

# Configuration template
MOLECULE_CONFIG = {
    "name": "Example",
    "input_type": InputType.FROM_CSV,
    "input_source": "",
    "repeat_unit": "Custom",
    "n_units": 3,
    "end_groups": {
        "start": "",
        "end": ""
    },
    "config": {
        "title": "Example Molecule",
        "method": "B3LYP",
        "basis_set": "3-21G",
        "charge": 0,
        "multiplicity": 1,
        "optimization": True,
        "frequency": True,
        "memory": "4GB",
        "nproc": 2,
        "additional_keywords": "",
        "output_format": OutputFormat.AVOGADRO
    }
}

# Example configurations
WATER_DIRECT = {
    **MOLECULE_CONFIG,
    "name": "Ethane-1,2-diol;terephthalic acid",
    "input_type": InputType.DIRECT_SMILES,
    "input_source": "C1=CC(=CC=C1C(=O)O)C(=O)O.C(CO)O"
}

# Example configurations
CUSTO = {
    **MOLECULE_CONFIG,
    "name": "NewChemical",
    "input_type": InputType.FROM_CSV,
    "input_source": "Tetrafluoroethylene"
}

# Modified polymer configurations

POLYETHYLENE = {
    **MOLECULE_CONFIG,
    "name": "Polyethylene",
    "input_type": InputType.DIRECT_SMILES,
    "input_source": "*CC*",  # Ethylene repeat unit
    "repeat_unit": "ethylene",
    "n_units": 6,
    "end_groups": {"start": "H", "end": "H"},
}

POLYPROPYLENE = {
    **MOLECULE_CONFIG,
    "name": "Polypropylene",
    "input_type": InputType.DIRECT_SMILES,
    "input_source": "*CC(C)*",  # Propylene repeat unit
    "repeat_unit": "propylene",
    "n_units": 6,
    "end_groups": {"start": "H", "end": "H"},
}

POLYSTYRENE = {
    **MOLECULE_CONFIG,
    "name": "Polystyrene",
    "input_type": InputType.DIRECT_SMILES,
    "input_source": "*CC(c1ccccc1)*",  # Styrene repeat unit
    "repeat_unit": "styrene",
    "n_units": 6,
    "end_groups": {"start": "H", "end": "H"},
}

POLYTETRAFLUOROETHYLENE = {
    **MOLECULE_CONFIG,
    "name": "Polytetrafluoroethylene",
    "input_type": InputType.DIRECT_SMILES,
    "input_source": "*C(F)(F)C(F)(F)*",  # Tetrafluoroethylene repeat unit
    "repeat_unit": "tetrafluoroethylene",
    "n_units": 6,
    "end_groups": {"start": "H", "end": "H"},
}

POLYMETHYL_METHACRYLATE = {
    **MOLECULE_CONFIG,
    "name": "PolymethylMethacrylate",
    "input_type": InputType.DIRECT_SMILES,
    "input_source": "*CC(C(=O)OC)*",  # Methyl methacrylate repeat unit
    "repeat_unit": "methyl methacrylate",
    "n_units": 6,
    "end_groups": {"start": "H", "end": "H"},
}

# List of polymer configurations
POLYMER_CONFIGS = [
    POLYETHYLENE,
    POLYPROPYLENE,
    POLYSTYRENE,
    POLYTETRAFLUOROETHYLENE,
    POLYMETHYL_METHACRYLATE,
]

def generate_polymer_3d_coords(repeat_unit_smiles: str, n_units: int) -> str:
    """
    Generate 3D coordinates for a polymer with proper connectivity
    """
    # Remove all asterisks from the SMILES
    clean_smiles = repeat_unit_smiles.replace('*', '')
    
    # Create multiple repeat units
    polymer_smiles = clean_smiles * n_units
    
    # Create molecule from SMILES
    mol = Chem.MolFromSmiles(polymer_smiles)
    if not mol:
        raise ValueError("Invalid polymer SMILES structure")
    
    # Add Hydrogens
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates with improved parameters
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.maxIterations = 500
    params.optimizeConstrainedGroups = True
    
    try:
        # First attempt at embedding
        status = AllChem.EmbedMolecule(mol, params)
        if status == -1:
            # Fallback to random coords
            status = AllChem.EmbedMolecule(mol, useRandomCoords=True)
        
        if status == -1:
            raise ValueError("3D embedding failed")
        
        # Improved optimization
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
        
        # Convert to XYZ format
        xyz_lines = []
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            xyz_lines.append(f"{mol.GetAtomWithIdx(i).GetSymbol():<2} {pos.x:>10.4f} {pos.y:>10.4f} {pos.z:>10.4f}")
        
        return f"{len(xyz_lines)}\nPolymer Generated by RDKit\n" + "\n".join(xyz_lines)
    
    except Exception as e:
        raise RuntimeError(f"Failed to generate 3D structure: {str(e)}")
    
    
# Modify process_molecule function to use this for polymers
def process_molecule(config: Dict):
    """Process a single molecule configuration."""
    try:
        # Get SMILES string based on input type
        if config["input_type"] == InputType.FROM_CSV:
            smiles = get_chemical_smiles(config["input_source"])
        else:  # InputType.DIRECT_SMILES
            smiles = config["input_source"]
        
        # Check if it's a polymer configuration
        if config.get('n_units', 1) > 1:
            # Use polymer generation for multiple units
            xyz = generate_polymer_3d_coords(smiles, config['n_units'])
        else:
            # Use existing 3D conversion for single molecules
            xyz = smiles_to_3d_coords(smiles)
        
        # Generate output file
        generate_output_file(
            config["name"],
            config["config"],
            xyz,
            config["config"]["output_format"]
        )
        
    except Exception as e:
        print(f"Error processing {config['name']}: {str(e)}")

def main():
    # Process all polymer configurations
    for config in POLYMER_CONFIGS:
        process_molecule(config)

if __name__ == "__main__":
    main()