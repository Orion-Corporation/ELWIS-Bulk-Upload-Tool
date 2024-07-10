from rdkit import Chem

import cdxml_converter

def process_sdf(files, callback):
    molecules = []
    for sdf_file in files:
        callback(f"Processing file: {sdf_file}")
        try:
            supplier = Chem.SDMolSupplier(sdf_file)
            for mol in supplier:
                if mol is not None:
                    molecule_data = {}
                    for prop_name in mol.GetPropNames():
                        molecule_data[prop_name] = mol.GetProp(prop_name)

                    # Extract atomic coordinates and bonds
                    atom_block = []
                    for atom in mol.GetAtoms():
                        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
                        atom_info = {
                            "symbol": atom.GetSymbol(),
                            "x": pos.x,
                            "y": pos.y,
                            "z": pos.z
                        }
                        atom_block.append(atom_info)

                    bond_block = []
                    for bond in mol.GetBonds():
                        bond_info = {
                            "begin_atom_idx": bond.GetBeginAtomIdx(),
                            "end_atom_idx": bond.GetEndAtomIdx(),
                            "bond_type": bond.GetBondTypeAsDouble()
                        }
                        bond_block.append(bond_info)

                    molecule_data["atom_block"] = atom_block
                    molecule_data["bond_block"] = bond_block
                    molecule_data["cdxml"] = cdxml_converter.convert_mol_to_cdxml(molecule_data)

                    molecules.append(molecule_data)
                else:
                    callback(f"Error: Molecule in file {sdf_file} could not be parsed and will be skipped.")
        except Exception as e:
            callback(f"Error processing file {sdf_file}: {str(e)}")

    callback(f"Total molecules extracted: {len(molecules)}")
    return molecules
