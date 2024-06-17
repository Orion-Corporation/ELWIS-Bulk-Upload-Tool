from rdkit import Chem

def process_sdf(files, params, callback):
    molecules = []
    for sdf_file in files:
        callback(f"Processing file: {sdf_file}")
        supplier = Chem.SDMolSupplier(sdf_file)
        print(f"Supplier: {supplier}")
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

                molecules.append(molecule_data)
            else:
                callback(f"Warning: Molecule in file {sdf_file} could not be parsed and will be skipped.")

    callback(f"Total molecules extracted: {len(molecules)}")
    return molecules
