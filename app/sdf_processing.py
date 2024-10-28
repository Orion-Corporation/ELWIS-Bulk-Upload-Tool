# sdf_processing.py

from rdkit import Chem
from openbabel import openbabel
import requests
import json
from logger import log_to_general_log
from config import BATCH_FIELDS_CONFIG, API_ENDPOINTS, SDF_PROPERTIES_CONFIG, VIABLE_SUPPLIERS
from config import load_viable_suppliers

# SDF processing functions using OpenBabel
def process_sdf(files, callback, project_value):
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

                # Separate fragments using OpenBabel
                separated_fragments = obMol.Separate()
                
                if len(separated_fragments) == 0:
                    # No fragments to separate
                    main_molecule = obMol
                    fragment = None
                elif len(separated_fragments) == 1:
                    # Only one molecule, main molecule
                    main_molecule = separated_fragments[0]
                    fragment = None
                elif len(separated_fragments) >= 2:
                    # Multiple fragments, assume the first is the salt and the second is the main molecule ??? works so far, will this be a problem?
                    main_molecule = separated_fragments[1]
                    fragment = separated_fragments[0]
                else:
                    print("Error: Unexpected number of fragments.")
                    continue

                # Extract molecular data using RDKit
                smiles = obConversion_smiles.WriteString(main_molecule).strip().upper()
                print(f"SMILES: {smiles}")

                # Extract properties using RDKit
                chemical_name = get_property(rdkit_mol, SDF_PROPERTIES_CONFIG.get("Chemical name", ["Chemical name", "Systematic name", "IUPAC"]), '')
                supplier_code = get_property(rdkit_mol, SDF_PROPERTIES_CONFIG.get("Supplier code", ["ID", "Query Mcule ID", "MOLPORTID"]), '')
                supplier_name = get_normalized_property(rdkit_mol, SDF_PROPERTIES_CONFIG.get("Supplier name", ["Supplier name", "SUPPLIER NAME"]), 'Unknown')
                amount_mg = get_property(rdkit_mol, SDF_PROPERTIES_CONFIG.get("Amount_mg", ["Amount_mg", "Amount (mg)", "QUANTITY"]), 0)
                compound_id = get_property(rdkit_mol, SDF_PROPERTIES_CONFIG.get("Compound ID", ["ID", "Delivered Mcule ID", "MOLPORTID"]), '')
                formula = get_property(rdkit_mol, SDF_PROPERTIES_CONFIG.get("Formula", ["Formula", "MOL FORMULA"]), '')
                purity = get_property(rdkit_mol, SDF_PROPERTIES_CONFIG.get("Purity", ["Purity", "Guaranteed purity (%)", "PURITY CLASSIFIED"]), '')
                po = get_property(rdkit_mol, SDF_PROPERTIES_CONFIG.get("PO", ["PO", "Customer PO", "PO NUMBER FROM CLIENT"]), '')
                plate_id = get_property(rdkit_mol, SDF_PROPERTIES_CONFIG.get("Plate_ID", ["Plate_ID", "Multi container ID", "BOX_NAME"]), '')
                well = get_property(rdkit_mol, SDF_PROPERTIES_CONFIG.get("Well", ["Well", "Single container position", "BOX_ROW"]), '')
                barcode = get_property(rdkit_mol, SDF_PROPERTIES_CONFIG.get("Barcode", ["Barcode", "VIAL_BARCODE"]), '')
                stereochemistry = get_property(rdkit_mol, SDF_PROPERTIES_CONFIG.get("Stereochemistry", ["Stereochem.data", "STEREOCHEMISTRY"]), 'No stereochemistry')

                molecule_data = {
                    "Chemical name": chemical_name,
                    "Supplier code": supplier_code,
                    "Supplier name": supplier_name,
                    "MolecularFormula": main_molecule.GetFormula().upper(),
                    "MW": main_molecule.GetMolWt(),
                    "Smile": smiles,
                    "Amount_mg": amount_mg,
                    "ID": compound_id,
                    "Formula": formula,
                    "Purity": purity,
                    "PO": po,
                    "Plate_ID": plate_id,
                    "Well": well,
                    "Barcode": barcode,
                    "Stereochemistry": stereochemistry
                }

                print(f"Extracted molecule data: {molecule_data}")

                # Extract fragment properties using RDKit
                fragment_data = {}
                if fragment is not None:
                    fragment_salt_name = get_property(rdkit_mol, ["Salt_Name", "Salt name"], '')
                    mw_salt = get_property(rdkit_mol, ["MW_salt", "MOL WEIGHT TOTAL"], fragment.GetMolWt())
                    mf_salt = get_property(rdkit_mol, ["Salt smiles", "MOL FORMULA"], fragment.GetFormula().upper())
                    
                    fragment_data = {
                        "Salt_name": fragment_salt_name,
                        "MolecularFormula": mf_salt,
                        "MW_salt": mw_salt
                    }
                    print(f"Extracted fragment data: {fragment_data}")

                # Extract atom and bond blocks using OpenBabel
                atom_block = []
                for atom in openbabel.OBMolAtomIter(main_molecule):
                    atomic_num = atom.GetAtomicNum()
                    if atomic_num == 0:
                        print(f"Invalid element symbol: {atom.GetType()}")
                        continue
                    element_symbol = openbabel.GetSymbol(atomic_num)
                    atom_block.append({
                        "symbol": element_symbol,
                        "x": atom.GetX(),
                        "y": atom.GetY(),
                        "z": atom.GetZ()
                    })
                print(f"Extracted atom block: {atom_block}")

                bond_block = [{
                    "begin_atom_idx": bond.GetBeginAtomIdx() - 1,
                    "end_atom_idx": bond.GetEndAtomIdx() - 1,
                    "bond_type": bond.GetBondOrder()
                } for bond in openbabel.OBMolBondIter(main_molecule)]
                print(f"Extracted bond block: {bond_block}")

                molecule_data.update({"atom_block": atom_block, "bond_block": bond_block})
                molecule_data["cdxml"] = convert_mol_to_cdxml(molecule_data)
                print(f"Converted to CDXML format")

                molecules.append(molecule_data)
                fragments.append(fragment_data)

            except Exception as e:
                print(f"Error processing molecule in file {sdf_file}: {str(e)}")
                callback(f"Error processing molecule in file {sdf_file}: {str(e)}")

    print(f"Total molecules extracted: {len(molecules)}")
    callback(f"Total molecules extracted: {len(molecules)}")
    return molecules, fragments

def get_property(mol, prop_names, default_value):
    for prop_name in prop_names:
        if mol.HasProp(prop_name):
            return mol.GetProp(prop_name)
    return default_value

import re

def normalize_name(name):
    # Remove non-alphabetic characters and convert to lowercase
    normalized_name = re.sub(r'[^a-zA-Z]', '', name).lower()
    # Remove common suffixes like "inc", "co", "ltd", etc., if they appear at the end of the string
    normalized_name = re.sub(r'(inc|co|ltd|corp|llc|gmbh|bv|plc)$', '', normalized_name)
    # Return the cleaned name with leading/trailing whitespace stripped
    return normalized_name.strip()

def get_normalized_property(mol, prop_names, default_value):
    VIABLE_SUPPLIERS = load_viable_suppliers()
    for prop_name in prop_names:
        if mol.HasProp(prop_name):
            # Get the raw supplier name from the SDF and normalize it
            sdf_supplier_name = mol.GetProp(prop_name).strip()
            normalized_sdf_name = normalize_name(sdf_supplier_name)
            print(f"Normalized SDF name: {normalized_sdf_name}")
            print(f"Loaded VIABLE_SUPPLIERS: {VIABLE_SUPPLIERS}")

            # Search for an exact match in normalized form with VIABLE_SUPPLIERS
            for viable_supplier in VIABLE_SUPPLIERS:
                normalized_viable_name = normalize_name(viable_supplier)
                print(f"Normalized Viable Supplier: {normalized_viable_name}")

                # If a match is found, return the exact form from VIABLE_SUPPLIERS
                if normalized_sdf_name == normalized_viable_name:
                    print(f"Matched SDF name '{sdf_supplier_name}' with viable supplier '{viable_supplier}'")
                    return viable_supplier

    # Return default if no properties match
    return default_value


def convert_mol_to_cdxml(molecule_data):
    obConversion = openbabel.OBConversion()
    if not obConversion.SetInAndOutFormats("mol", "cdxml"):
        print("Error: Could not set input/output formats to mol/cdxml")
        return None

    mol = openbabel.OBMol()
    
    # Add atoms and bonds without altering hydrogen count
    for atom_info in molecule_data["atom_block"]:
        atom = mol.NewAtom()
        atomic_num = openbabel.GetAtomicNum(atom_info["symbol"])
        if atomic_num == 0:
            print(f"Invalid element symbol: {atom_info['symbol']}")
            continue
        atom.SetAtomicNum(atomic_num)
        atom.SetVector(atom_info["x"], atom_info["y"], atom_info["z"])

    for bond_info in molecule_data["bond_block"]:
        mol.AddBond(bond_info["begin_atom_idx"] + 1, bond_info["end_atom_idx"] + 1, int(bond_info["bond_type"]))

    output_cdxml = obConversion.WriteString(mol)
    if not output_cdxml:
        print("Error: Could not write the CDXML content")
        return None

    print("Successfully converted molecule to CDXML format")
    return output_cdxml

def construct_payload(molecule_data, salt_id, fragment_data, project_value, library_id):
    # Ensure name is a string
    salt_name = fragment_data.get('Salt_name', '') or fragment_data.get('Salt_Name', '')
    if isinstance(salt_name, (list, tuple)):
        salt_name = salt_name[0] if salt_name else ''
    else:
        salt_name = str(salt_name)

    solvate_name = fragment_data.get('Solvate_name', '') or fragment_data.get('Solvate_Name', '')
    if isinstance(solvate_name, (list, tuple)):
        solvate_name = solvate_name[0] if solvate_name else ''
    else:
        solvate_name = str(solvate_name)
        
    source_value = BATCH_FIELDS_CONFIG.get("source", {}).get("value", "Acquired")
    # project_value = BATCH_FIELDS_CONFIG.get("project", {}).get("value", "Unspecified") USES DROPDOWN VALUE
    synthesis_datetime_value = BATCH_FIELDS_CONFIG.get("synthesis_datetime", {}).get("value", "2011-10-10T14:48:00Z")
    chemist_value = BATCH_FIELDS_CONFIG.get("chemist", {}).get("value", "TestUser MCChemist")
    batch_purpose_value = BATCH_FIELDS_CONFIG.get("batch_purpose", {}).get("value", "Dummy compound")
    batch_type_value = BATCH_FIELDS_CONFIG.get("batch_type", {}).get("value", "Discovery")
    
    print("Retrieved source value:", source_value)
    data = {
        "data": {
            "type": "asset",
            "attributes": {
                "synonyms": [
                    molecule_data.get('Smile', ''),
                    molecule_data.get('MolecularFormula', '')
                ],
                "fields": [
                    {
                        "id": "5d6e0287ee35880008c18d6d",
                        "value": molecule_data.get("cdxml", "")
                    },
                    {
                        "id": "62f9fe5b74770f14d1de43a8",
                        "value": molecule_data.get("Stereochemistry", "No stereochemistry")
                    },
                    # Chemical Name assigned by ELWIS automatically
                    # {
                    #     "id": "5d6e0287ee35880008c18db6",
                    #     "value": molecule_data.get("Chemical name", "")
                    # },
                    {
                        "id": "5d6e0287ee35880008c18db7",
                        "value": {
                            "rawValue": str(molecule_data.get("MW", "")),
                            "displayValue": f"{molecule_data.get('MW', '')} g/mol"
                        }
                    },
                    {
                        "id": "5d6e0287ee35880008c18db8",
                        "value": float(molecule_data.get("Amount_mg", 0))
                    },
                    {
                        "id": "5d6e0287ee35880008c18db9",
                        "value": {
                            "rawValue": molecule_data.get("MolecularFormula", ""),
                            "displayValue": molecule_data.get("MolecularFormula", "")
                        }
                    }
                ]
            },
            "relationships": {
                "batch": {
                    "data": {
                        "type": "batch",
                        "attributes": {
                            "fields": [
                                {
                                "id": "62fcceeb19660304d1e5beee",
                                "value": source_value     
                                },
                                {
                                "id": "63469c69ed8a726a31923537",
                                "value": project_value
                                },
                                {
                                "id": "62fcceeb19660304d1e5bef1",
                                "value": synthesis_datetime_value
                                },
                                {
                                "id": "6384a1270d28381d21deaca7",
                                "value": chemist_value
                                },
                                {
                                "id": "62fa096d19660304d1e5b2db",
                                "value": batch_purpose_value
                                },
                                {
                                "id": "62fa096d19660304d1e5b2da",
                                "value": batch_type_value
                                },
                                {
                                "id": "62fcceeb19660304d1e5bef2",
                                "value": library_id
                                },
                                {
                                "id": "62fa0b5b19660304d1e5b2de", 
                                "value": molecule_data.get("Supplier name", "Unknown")
                                }
                            ]
                        }
                    }
                }
            }
        }
    }



    if salt_id:
        fragments = {}
        if salt_id:
            fragments["salts"] = [{
                "type": "SALT",
                "name": salt_name,
                "mf": fragment_data.get('MolecularFormula', ''),
                "mw": {
                    "rawValue": str(fragment_data.get('MW_salt', '')),
                    "displayValue": f"{fragment_data.get('MW_salt', '')} g/mol"
                },
                "id": f"SALT:{salt_id}",
                "coefficient": 1
            }]
        data["data"]["relationships"]["batch"]["data"]["attributes"]["fragments"] = fragments

    return data

def check_uniqueness(molecule_data, api_key, project_value, library_id):
    headers = {
        'Content-Type': 'application/vnd.api+json',
        'accept': 'application/vnd.api+json',
        'x-api-key': api_key
    }

    # Dynamically load values from BATCH_FIELDS_CONFIG
    source_value = BATCH_FIELDS_CONFIG.get("source", {}).get("value", "Acquired")
    # project_value = BATCH_FIELDS_CONFIG.get("project", {}).get("value", "Unspecified") USES DROPDOWN VALUE
    synthesis_datetime_value = BATCH_FIELDS_CONFIG.get("synthesis_datetime", {}).get("value", "2011-10-10T14:48:00Z")
    chemist_value = BATCH_FIELDS_CONFIG.get("chemist", {}).get("value", "TestUser MCChemist")
    batch_purpose_value = BATCH_FIELDS_CONFIG.get("batch_purpose", {}).get("value", "Dummy compound")
    batch_type_value = BATCH_FIELDS_CONFIG.get("batch_type", {}).get("value", "Discovery")

    data = {
        "data": {
            "type": "asset",
            "attributes": {
                "synonyms": [
                    molecule_data.get('Smile', ''),
                    molecule_data.get('MolecularFormula', '')
                ],
                "fields": [
                    {
                        "id": "5d6e0287ee35880008c18d6d",
                        "value": molecule_data.get("cdxml", "")
                    },
                    {
                        "id": "62f9fe5b74770f14d1de43a8",
                        "value": molecule_data.get("Stereochemistry", "No stereochemistry")
                    },
                    # Chemical Name assigned by ELWIS automatically
                    # {
                    #     "id": "5d6e0287ee35880008c18db6",
                    #     "value": molecule_data.get("Chemical name", "")
                    # },
                    {
                        "id": "5d6e0287ee35880008c18db7",
                        "value": {
                            "rawValue": str(molecule_data.get("MW", "")),
                            "displayValue": f"{molecule_data.get('MW', '')} g/mol"
                        }
                    },
                    {
                        "id": "5d6e0287ee35880008c18db8",
                        "value": float(molecule_data.get("Amount_mg", 0))
                    },
                    {
                        "id": "5d6e0287ee35880008c18db9",
                        "value": {
                            "rawValue": molecule_data.get("MolecularFormula", ""),
                            "displayValue": molecule_data.get("MolecularFormula", "")
                        }
                    }
                ]
            },
            "relationships": {
                "batch": {
                    "data": {
                        "type": "batch",
                        "attributes": {
                            "fields": [
                                {
                                    "id": "62fcceeb19660304d1e5beee",  # Source
                                    "value": source_value
                                },
                                {
                                    "id": "63469c69ed8a726a31923537",  # Project
                                    "value": project_value
                                },
                                {
                                    "id": "62fcceeb19660304d1e5bef1",  # Synthesis datetime
                                    "value": synthesis_datetime_value
                                },
                                {
                                    "id": "6384a1270d28381d21deaca7",  # Chemist
                                    "value": chemist_value
                                },
                                {
                                    "id": "62fa096d19660304d1e5b2db",  # Batch Purpose
                                    "value": batch_purpose_value
                                },
                                {
                                    "id": "62fa096d19660304d1e5b2da",  # Batch Type
                                    "value": batch_type_value
                                },
                                {
                                    "id": "62fcceeb19660304d1e5bef2",
                                    "value": library_id
                                },
                                {
                                "id": "62fa0b5b19660304d1e5b2de", 
                                "value": molecule_data.get("Supplier name", "Unknown")
                                }
                            ]
                        }
                    }
                }
            }
        }
    }

    uniqueness_endpoint = f"{API_ENDPOINTS['Compound Endpoint']}/uniquenessCheck"
    response = requests.post(uniqueness_endpoint, headers=headers, data=json.dumps(data))
    print(f"Uniqueness response: {response.text}")
    if response.status_code in [200, 201]:
        return response.json()
    else:
        print(f"Error checking uniqueness: {response.status_code} - {response.text}")
        return None
