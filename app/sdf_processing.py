from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel
from logger import log_to_general_log

def process_sdf(files, callback):
    print("Starting process_sdf")
    log_to_general_log(f"Found {len(files)} files to process")
    log_to_general_log(f"Starting to process SDF files: {files}")
    molecules = []
    fragments = []

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("sdf", "mol")
    print("OpenBabel conversion set for sdf to mol")

    obConversion_smiles = openbabel.OBConversion()
    obConversion_smiles.SetOutFormat("smiles")
    print("OpenBabel conversion set for mol to smiles")

    for sdf_file in files:
        print(f"Processing file: {sdf_file}")
        log_to_general_log(f"Starting to process: {sdf_file}")
        callback(f"Processing file: {sdf_file}")

        # Load the SDF file using RDKit to extract all properties
        supplier = Chem.SDMolSupplier(sdf_file)
        if not supplier:
            print(f"RDKit could not read file: {sdf_file}")
            continue

        for rdkit_mol in supplier:
            if rdkit_mol is None:
                print(f"Failed to parse molecule in file: {sdf_file}")
                continue

            try:
                # Convert RDKit molecule to OpenBabel molecule
                mol_block = Chem.MolToMolBlock(rdkit_mol)
                obMol = openbabel.OBMol()
                obConversion.ReadString(obMol, mol_block)

                # Further processing...

            except Exception as e:
                print(f"Error processing molecule in file {sdf_file}: {str(e)}")
                callback(f"Error processing molecule in file {sdf_file}: {str(e)}")

    print(f"Total molecules extracted: {len(molecules)}")
    callback(f"Total molecules extracted: {len(molecules)}")
    return molecules, fragments

# Other SDF related functions such as log_duplicate, etc.
