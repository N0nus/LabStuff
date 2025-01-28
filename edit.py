import os
from typing import List

# Configuration Dictionary (Edit this directly in the script)
config = {
    "file_name": "Valine.com",
    "charge_multiplicity": {
        "charge": 0,
        "multiplicity": 1
    },
    "method_basis_set": {
        "method": "B3LYP",
        "basis_set": "6-31G"
    },
    "optimization_frequency": {
        "optimization": True,
        "frequency": True
    },
    "molecule_name": "Polymer_Chain",
    "shift_coordinates": {
        "x_shift": 0.5,
        "y_shift": 0.0,
        "z_shift": 0.0
    }
}

def check_file_exists(file_path: str) -> bool:
    """Check if the specified file exists."""
    if os.path.isfile(file_path):
        return True
    else:
        print(f"File not found: {file_path}")
        return False

def read_gaussian_input(file_path: str) -> List[str]:
    """Read the contents of the Gaussian input file."""
    with open(file_path, 'r') as f:
        return f.readlines()

def write_gaussian_input(file_path: str, lines: List[str]):
    """Write the modified content back to the Gaussian input file."""
    with open(file_path, 'w') as f:
        f.writelines(lines)

def modify_charge_multiplicity(lines: List[str], charge: int, multiplicity: int) -> List[str]:
    """Modify the charge and multiplicity of the molecule."""
    lines[6] = f"{charge} {multiplicity}\n"
    return lines

def modify_method_basis_set(lines: List[str], method: str, basis_set: str) -> List[str]:
    """Modify the method and basis set in the Gaussian input."""
    lines[3] = f"#P {method}/{basis_set} opt freq\n"
    return lines

def modify_optimization_frequency(lines: List[str], optimization: bool, frequency: bool) -> List[str]:
    """Modify optimization and frequency calculation settings."""
    if optimization:
        lines[3] = lines[3].replace("opt", "opt")
    else:
        lines[3] = lines[3].replace("opt", "")
    
    if frequency:
        lines[3] = lines[3].replace("freq", "freq")
    else:
        lines[3] = lines[3].replace("freq", "")
    
    return lines

def modify_molecule_name(lines: List[str], molecule_name: str) -> List[str]:
    """Change the molecule name in the Gaussian input file."""
    lines[5] = f"{molecule_name}\n"
    return lines

def shift_coordinates(lines: List[str], x_shift: float, y_shift: float, z_shift: float) -> List[str]:
    """Shift atomic coordinates in the Gaussian input file."""
    coords_start_index = 7  # Starting line for coordinates in Gaussian input
    for i in range(coords_start_index, len(lines)):
        if lines[i].strip() == "":
            break  # End of coordinates section
        
        coords = lines[i].split()
        if len(coords) < 4:
            print(f"Skipping invalid coordinate line: {lines[i].strip()}")
            continue  # Skip lines that don't have at least 4 elements (atom + x, y, z)
        
        try:
            # Shift the coordinates (x, y, z)
            coords[1] = str(float(coords[1]) + x_shift)  # Shift x
            coords[2] = str(float(coords[2]) + y_shift)  # Shift y
            coords[3] = str(float(coords[3]) + z_shift)  # Shift z
        except ValueError:
            print(f"Skipping line due to invalid number format: {lines[i].strip()}")
            continue  # Skip lines where conversion to float fails
        
        lines[i] = "  ".join(coords) + "\n"
    
    return lines


def modify_input_file(file_path: str, config: dict):
    """Modify the Gaussian input file based on the configuration."""
    if not check_file_exists(file_path):
        return

    lines = read_gaussian_input(file_path)

    # Apply modifications from the config dictionary
    lines = modify_charge_multiplicity(lines, config["charge_multiplicity"]["charge"], config["charge_multiplicity"]["multiplicity"])
    lines = modify_method_basis_set(lines, config["method_basis_set"]["method"], config["method_basis_set"]["basis_set"])
    lines = modify_optimization_frequency(lines, config["optimization_frequency"]["optimization"], config["optimization_frequency"]["frequency"])
    lines = modify_molecule_name(lines, config["molecule_name"])
    lines = shift_coordinates(lines, config["shift_coordinates"]["x_shift"], config["shift_coordinates"]["y_shift"], config["shift_coordinates"]["z_shift"])

    # Write the modified content back to the Gaussian input file
    write_gaussian_input(file_path, lines)
    print(f"Gaussian input file {file_path} has been modified.")

def main():
    """Main function to run the input file modification."""
    file_name = config["file_name"]

    # Call the modify_input_file function to modify the Gaussian input file
    modify_input_file(file_name, config)

if __name__ == "__main__":
    main()
