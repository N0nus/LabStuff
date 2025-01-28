import pandas as pd
import os
from openbabel import pybel
from typing import Dict

openbabel_data_path = r"C:\Users\djloc\AppData\Roaming\OpenBabel-3.1.1\data"
os.environ['OBABEL_DATA'] = openbabel_data_path

# Load validated chemicals from the CSV
def load_chemicals(csv_path: str) -> Dict[str, str]:
    df = pd.read_csv(csv_path)
    chemicals = {row['Chemical Name']: row['SMILES Notation'] for _, row in df.iterrows()}
    return chemicals

# Build polymer SMILES with explicit connection points
def build_polymer_smiles(repeat_unit: str, n_units: int, end_groups: Dict[str, str], repeat_units: Dict[str, str]) -> str:
    """Build polymer SMILES with explicit connection points."""
    if repeat_unit not in repeat_units:
        raise ValueError(f"Unknown repeat unit: {repeat_unit}")
    
    base = repeat_units[repeat_unit].strip('*')
    polymer = []
    
    # Start group
    if "start" in end_groups:
        polymer.append(f"{end_groups['start']}")
    
    # Repeat units (connect with single bonds)
    polymer.append(f"({base})" * n_units)
    
    # End group
    if "end" in end_groups:
        polymer.append(f"{end_groups['end']}")
    
    return "".join(polymer).replace('*', '')  # Remove connection points after joining

# Convert SMILES to 3D coordinates
def smiles_to_3d_coords(smiles: str) -> str:
    """Convert SMILES to 3D coordinates using Open Babel."""
    try:
        mol = pybel.readstring("smi", smiles)
        mol.addh()
        mol.make3D(forcefield="mmff94", steps=500)
        return mol.write("xyz")
    except Exception as e:
        raise RuntimeError(f"3D conversion failed: {str(e)}")

# Generate Gaussian input file
def generate_gaussian_input(polymer_name: str, config: Dict, xyz_coords: str):
    """Generate Gaussian input file."""
    os.makedirs("outputs2", exist_ok=True)
    filename = f"outputs2/{polymer_name}.com"
    
    input_content = [
        f"%chk={polymer_name}.chk",
        "%mem=4GB",
        "%nproc=2",
        f"#P {config['method']}/{config['basis_set']} {'opt' if config['optimization'] else ''} {'freq' if config['frequency'] else ''}",
        "",
        config["title"],
        "",
        f"{config['charge']} {config['multiplicity']}",
        *xyz_coords.split('\n')[2:],  # Skip first two lines of XYZ
        ""
    ]
    
    with open(filename, 'w') as f:
        f.write('\n'.join(input_content))
    
    print(f"Generated {filename}")
    return filename

# Test Case
TEST_CASE = {
    "name": "Poly(ethylene glycol)",
    "repeat_unit": "Custom",  # SMILES for the repeating unit of PEG
    "n_units": 3,  # Number of repeating units
    "end_groups": {
        "start": "[H]O",  # Hydroxyl group as the starting end group
        "end": "O[H]"  # Hydroxyl group as the ending end group
    },
    "config": {
        "title": "Polyethylene Glycol Chain",
        "method": "B3LYP",  # Computational method
        "basis_set": "3-21G",  # Basis set
        "charge": 0,  # Neutral molecule
        "multiplicity": 1,  # Singlet state
        "optimization": True,  # Perform geometry optimization
        "frequency": True  # Perform frequency analysis
    }
}

def main():
    try:
        repeat_units = load_chemicals("chemicals.csv")  # Load the chemicals from CSV
        case = TEST_CASE
        

        # Build the polymer SMILES
        try:
            smiles = build_polymer_smiles(
                case["repeat_unit"],
                case["n_units"],
                case["end_groups"],
                repeat_units  # Pass the repeat_units dictionary
            )
            
            # Convert SMILES to 3D coordinates
            xyz = smiles_to_3d_coords(smiles)
            
            # Generate Gaussian input file for each chemical
            generate_gaussian_input(case["name"], case["config"], xyz)
        
        except Exception as e:
            # Only print errors if something goes wrong
            print(f"Error processing {case["name"]}: {str(e)}")
    
    except Exception as e:
        # Catch any global exceptions and print the error
        print(f"Error: {str(e)}")



if __name__ == "__main__":
    main()
