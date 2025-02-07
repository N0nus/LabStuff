Converts from SMILES string representation to 3d xqz coordinates using RDkit library.
Also uses PubChem API






Notes:

-Information about file can be specified in configuration. e.g, 

# Configuration template
MOLECULE_CONFIG = {
    "name": "Example",                       Name of what your making
    "input_type": InputType.FROM_CSV,        Where the input is coming from, direct string or csv/pubchem
    "input_source": "",                      SMILES string or name of chemical
    "repeat_unit": "Custom",                 Name of the unit repeating
    "n_units": 3,                            Number if times it repeats (I think this isnt working at makes too many, idk tho)
    "end_groups": {
        "start": "",                         Molecule on start and end 
        "end": ""
    },
    "config": {
        "title": "Example Molecule",         Title of the file
        "method": "B3LYP",                   Idk what this stuff does  |
        "basis_set": "3-21G",                                          |
        "charge": 0,                                                   | 
        "multiplicity": 1,                                             |
        "optimization": True,                                          |
        "frequency": True,                                             |
        "memory": "4GB",                                               |
        "nproc": 2,                                                    |
        "additional_keywords": "",           __________________________|
        "output_format": OutputFormat.AVOGADRO           Whether output is in avagadro .xyz file or gaussian .com file.
    }
}

-Uses either name of chemical or directly from smiles to generate

-Can fetch SMILES string for chemical by name from either PubChemicals.csv or https://pubchem.ncbi.nlm.nih.gov/. If the chemical name was not found in csv file it will look for it on the pubchem website and, if it finds it, add it to the csv file

-Use *'s to specify connection points when making polymers. like so:  "input_source": "*CC*",  # Ethylene repeat unit

-Important: I have no idea what im doing and I used help of AI to make this, please don't sue me.
