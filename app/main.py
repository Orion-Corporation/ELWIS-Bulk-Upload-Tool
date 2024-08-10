from kivy.config import Config
Config.set('kivy', 'window_icon', 'burana.ico')

import kivy
import logging
import datetime
from dotenv import load_dotenv
import os
from kivy.app import App
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.button import Button
from kivy.uix.label import Label
from kivy.uix.filechooser import FileChooserIconView
from kivy.uix.popup import Popup
from kivy.uix.scrollview import ScrollView
from kivy.uix.textinput import TextInput
from kivy.core.window import Window
from kivy.properties import ListProperty, BooleanProperty
from kivy.lang import Builder
from kivy.clock import Clock
import json
import shutil
from openbabel import openbabel
import requests

# Set the window size
Window.size = (1200, 600)

# Load environment variables
load_dotenv()
api_key = os.getenv('API_KEY')
kivy.logger.Logger.setLevel(logging.ERROR)

# Load Kivy file
try:
    Builder.load_file('styles.kv')
except Exception as e:
    print(f"Error loading Kv file: {e}")

# Load configuration files
def load_config(file_path, default):
    try:
        with open(file_path, 'r') as file:
            return json.load(file)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error loading config from {file_path}: {e}")
        return default

OUTPUT_PATHS = load_config('config/output_paths.json', {
    "success_log": "Successfull/success.txt",
    "duplicate_log": "Failed/duplicates.txt",
    "failed_log": "Failed/failed.txt",
    "general_log": "logs.txt",
    "download_folder": "downloads"
})

API_ENDPOINTS = load_config('config/api_endpoints.json', {
    "Fragment Endpoint": "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/fragments",
    "Compound Endpoint": "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/materials/Compounds/assets"
})

BATCH_FIELDS_CONFIG = load_config('config/batch_fields_config.json', {})

# API functions
def upload_fragment(fragment_data, fragment_type, callback, headers):
    fragment_endpoint = f"{API_ENDPOINTS['Fragment Endpoint']}/{fragment_type}"
    response = requests.post(fragment_endpoint, data=json.dumps(fragment_data), headers=headers)
    if response.status_code in [200, 201]:
        fragment_id = response.json()['data']['id']
        return fragment_id
    callback(f"Failed to upload {fragment_type}. Status code: {response.status_code}, response: {response.text}")
    return None

def get_existing_fragment_details(fragment_name, fragment_type, api_key):
    headers = {
        'accept': 'application/vnd.api+json',
        'x-api-key': api_key
    }
    params = {
        'filter[name]': fragment_name
    }
    fragment_endpoint = f"{API_ENDPOINTS['Fragment Endpoint']}/{fragment_type}"
    response = requests.get(fragment_endpoint, headers=headers, params=params)

    if response.status_code in [200, 201]:
        results = response.json().get('data', [])
        for fragment_data in results:
            # Ensure the name matches exactly
            if fragment_data['attributes'].get('name', '').lower() == fragment_name.lower():
                print(f"Exact match found: {fragment_data}")
                return {
                    'id': fragment_data['id'],
                    'name': fragment_data['attributes'].get('name', ''),
                    'mf': fragment_data['attributes'].get('mf', ''),
                    'mw': fragment_data['attributes'].get('mw', '')
                }
        print(f"No exact match found for salt name: {fragment_name}")
    return None

def upload_fragments(fragment_data, callback, headers, api_key):
    salt_details = None

    # Ensure the correct salt_name field is checked and used
    salt_name = fragment_data.get('Salt_name', '') or fragment_data.get('Salt_Name', '')

    # Ensure salt_name is a string
    if isinstance(salt_name, (list, tuple)):
        salt_name = salt_name[0] if salt_name else ''
    else:
        salt_name = str(salt_name)

    # Debug log for salt name
    print(f"Debug: Retrieved salt name from fragment_data: '{salt_name}'")

    if salt_name:
        callback(f"Salt detected: {salt_name}")
        print(f"Salt detected in fragment data: {salt_name}")  # Debug print for salt name

        salt_details = get_existing_fragment_details(salt_name, "salts", api_key)
        print(f"Debug: Retrieved salt details from API: '{salt_details}' for salt_name: '{salt_name}'")

        if salt_details and salt_details.get('id'):
            callback(f"Salt '{salt_name}' already exists with ID: {salt_details['id']}, not uploading duplicate")
        else:
            salt_data = {
                "data": {
                    "type": "salt",
                    "attributes": {
                        "name": salt_name,
                        "mf": fragment_data.get('MolecularFormula', ''),
                        "mw": f"{fragment_data.get('MW_salt', '')} g/mol"
                    }
                }
            }
            salt_id = upload_fragment(salt_data, "salts", callback, headers)
            if salt_id:
                callback(f"Successfully uploaded salt: {salt_name}")
                salt_details = {'id': salt_id, 'name': salt_name, 'mf': salt_data['data']['attributes']['mf'], 'mw': salt_data['data']['attributes']['mw']}
            else:
                callback(f"Failed to upload salt: {salt_name}")
    else:
        print("No valid salt name found, skipping fragment upload.")
    
    return salt_details

# SDF processing functions using OpenBabel
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import AllChem

def process_sdf(files, callback):
    print("Starting process_sdf")
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

                # Separate fragments using OpenBabel
                separated_fragments = obMol.Separate()
                
                if len(separated_fragments) == 0:
                    # No fragments to separate
                    main_molecule = obMol
                    fragment = None
                elif len(separated_fragments) == 1:
                    # Only one molecule, assume it's the main molecule
                    main_molecule = separated_fragments[0]
                    fragment = None
                elif len(separated_fragments) >= 2:
                    # Multiple fragments, assume the first is the salt/solvent and the second is the main molecule
                    main_molecule = separated_fragments[1]
                    fragment = separated_fragments[0]
                else:
                    print("Error: Unexpected number of fragments.")
                    continue

                # Extract molecular data using RDKit
                smiles = obConversion_smiles.WriteString(main_molecule).strip().upper()
                print(f"SMILES: {smiles}")

                # Extract properties using RDKit
                chemical_name = rdkit_mol.GetProp("Chemical name") if rdkit_mol.HasProp("Chemical name") else ''
                # Remove everything after the first whitespace in the chemical name
                if chemical_name:
                    chemical_name = chemical_name.split(' ')[0]
                    
                amount_mg = rdkit_mol.GetProp("Amount_mg") if rdkit_mol.HasProp("Amount_mg") else 0
                compound_id = rdkit_mol.GetProp("ID") if rdkit_mol.HasProp("ID") else ''
                formula = rdkit_mol.GetProp("Formula") if rdkit_mol.HasProp("Formula") else ''
                purity = rdkit_mol.GetProp("Purity") if rdkit_mol.HasProp("Purity") else ''
                po = rdkit_mol.GetProp("PO") if rdkit_mol.HasProp("PO") else ''
                plate_id = rdkit_mol.GetProp("Plate_ID") if rdkit_mol.HasProp("Plate_ID") else ''
                well = rdkit_mol.GetProp("Well") if rdkit_mol.HasProp("Well") else ''
                barcode = rdkit_mol.GetProp("Barcode") if rdkit_mol.HasProp("Barcode") else ''
                # stereochemistry
                stereochemistry = rdkit_mol.GetProp("Stereochem.data") if rdkit_mol.HasProp("Stereochem.data") else 'No stereochemistry'

                molecule_data = {
                    "Chemical name": chemical_name,
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
                    if rdkit_mol.HasProp("Salt_Name"):
                        fragment_salt_name = rdkit_mol.GetProp("Salt_Name")
                    elif rdkit_mol.HasProp("Salt_name"):
                        fragment_salt_name = rdkit_mol.GetProp("Salt_name")
                    else:
                        fragment_salt_name = ''
                    
                    mw_salt = rdkit_mol.GetProp("MW_salt") if rdkit_mol.HasProp("MW_salt") else fragment.GetMolWt()
                    mf_salt = (rdkit_mol.GetProp("Salt smiles").strip('[]') if rdkit_mol.HasProp("Salt smiles") else fragment.GetFormula().upper())
                    
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

def construct_payload(molecule_data, salt_id, fragment_data):
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
    project_value = BATCH_FIELDS_CONFIG.get("project", {}).get("value", "Unspecified")
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
                    {
                        "id": "5d6e0287ee35880008c18db6",
                        "value": molecule_data.get("Chemical name", "")
                    },
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

def send_request(data, file, callback, endpoint, headers, OUTPUT_PATHS):
    data_json = json.dumps(data)
    try:
        response = requests.post(endpoint, data=data_json, headers=headers)
        # write response to json file
        with open('response.json', 'w') as f:
            json.dump(response.json(), f, indent=4)
        if response.status_code in [200, 201]:
            handle_success(file, data, OUTPUT_PATHS, callback)
            return True
        handle_failure(file, data, response, OUTPUT_PATHS, callback)
        return False
    except requests.exceptions.RequestException as e:
        callback(f"Network error during request: {str(e)}")
        return False

def check_uniqueness(molecule_data, api_key):
    headers = {
        'Content-Type': 'application/vnd.api+json',
        'accept': 'application/vnd.api+json',
        'x-api-key': api_key
    }

    # Dynamically load values from BATCH_FIELDS_CONFIG
    source_value = BATCH_FIELDS_CONFIG.get("source", {}).get("value", "Acquired")
    project_value = BATCH_FIELDS_CONFIG.get("project", {}).get("value", "Unspecified")
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
                        "value": "No stereochemistry"
                    },
                    {
                        "id": "5d6e0287ee35880008c18db6",
                        "value": molecule_data.get("Chemical name", "")
                    },
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

def post_to_api(molecule_data, fragment_data, file, callback, api_key, OUTPUT_PATHS):
    # Check uniqueness of the compound
    uniqueness_result = check_uniqueness(molecule_data, api_key)
    
    if uniqueness_result is None:
        callback(f"Could not verify uniqueness for {file}. Aborting upload.")
        return False
    
    if uniqueness_result.get("data"):  # If the result is not empty, it's a duplicate
        orm_code = uniqueness_result["data"][0]["attributes"].get("name", "Unknown")
        log_duplicate(file, molecule_data, callback, orm_code)
        return False
    
    # Prepare headers
    headers = {
        'Content-Type': 'application/vnd.api+json',
        'accept': 'application/vnd.api+json',
        'x-api-key': api_key
    }
    
    # Upload fragments (salts/solvates)
    salt_id = upload_fragments(fragment_data, callback, headers, api_key)
    
    # Construct the payload
    payload = construct_payload(molecule_data, salt_id, fragment_data)
    
    # Send the request
    success = send_request(payload, file, callback, API_ENDPOINTS['Compound Endpoint'], headers, OUTPUT_PATHS)
    print(f"Payload sent for file {file}: {payload}")  # Debug print for payload sent
    # write payload to json
    with open('payload.json', 'w') as f:
        json.dump(payload, f, indent=4)
        # also write the response in full
    return success

def handle_success(file, data, OUTPUT_PATHS, callback):
    # Log the success information using the molecular formula instead of SMILES
    molecular_formula = data['data']['attributes']['synonyms'][1]  # Assuming the second synonym is the molecular formula
    with open(OUTPUT_PATHS['success_log'], 'a') as success_log:
        success_log.write(f"File: {file}, Molecule: {molecular_formula}, API Response Status Code: 200\n")
    callback(f"Success: {molecular_formula} logged successfully.")


def handle_failure(file, data, response, OUTPUT_PATHS, callback):
    log_file = OUTPUT_PATHS['failed_log']
    
    with open(log_file, 'a') as log:
        log.write(f"File: {file}, Molecule: {data['data']['attributes']['synonyms'][0]}, API Response Status Code: {response.status_code}, response text: {response.text}\n")
    
    callback(f"Failed: {data['data']['attributes']['synonyms'][0]} with status code {response.status_code}. Logged failure.")

    try:
        response_json = response.json()
        with open(f"{log_file}_response.json", 'w') as f:
            json.dump(response_json, f, indent=4)
    except ValueError:
        with open(f"{log_file}_response.json", 'w') as f:
            f.write(response.text)

def log_duplicate(file, molecule_data, callback, orm_code):
    # Log the duplicate information
    with open(OUTPUT_PATHS['duplicate_log'], 'a') as duplicate_log:
        duplicate_log.write(f"File: {file}, Molecule: {molecule_data.get('MolecularFormula', '')}, ORM Code: {orm_code}\n")
    
    callback(f"Duplicate: {molecule_data.get('MolecularFormula', '')}, ORM Code: {orm_code}. Logged duplicate.")

# Kivy app classes
class CustomFileChooserIconView(FileChooserIconView):
    selection = ListProperty([])

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        Window.bind(on_key_down=self._on_key_down)

    def _on_key_down(self, instance, keyboard, keycode, text, modifiers):
        if 'ctrl' in modifiers and text == 'a':
            self.select_all_files()

    def select_all_files(self):
        all_files = [f for f in self.files if f.upper().endswith('.SDF') or f.lower().endswith('.sdf')]
        self.selection = all_files
        self._trigger_update_selection()

    def _trigger_update_selection(self):
        self._update_item_selection()

class FileChooserPopup(Popup):
    def __init__(self, select_callback, **kwargs):
        super().__init__(**kwargs)
        self.select_callback = select_callback
        self.title = "Select SDF Files"
        self.size_hint = (0.9, 0.9)
        
        layout = BoxLayout(orientation='vertical', padding=10, spacing=10)
        
        self.filechooser = CustomFileChooserIconView(multiselect=True, filters=['*.sdf'])
        layout.add_widget(self.filechooser)
        
        button_layout = BoxLayout(size_hint_y=None, height='48dp', spacing=10)
        
        select_button = Button(text="Select", background_normal='', background_color=[0, 0.14, 0.29, 1], color=[1, 1, 1, 1])
        select_button.bind(on_release=self.on_select)
        button_layout.add_widget(select_button)
        
        close_button = Button(text="Close", background_normal='', background_color=[0.5, 0, 0, 1], color=[1, 1, 1, 1])
        close_button.bind(on_release=self.on_close)
        button_layout.add_widget(close_button)
        
        layout.add_widget(button_layout)
        
        self.add_widget(layout)
    
    def on_select(self, *args):
        selected_files = self.filechooser.selection
        if selected_files:
            self.select_callback(selected_files)
            self.dismiss()

    def on_close(self, *args):
        self.dismiss()

class ReadmePopup(Popup):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.title = "README"
        self.size_hint = (0.8, 0.8)
        layout = BoxLayout(orientation='vertical', padding=10, spacing=10)
        scrollview = ScrollView()
        self.text_input = TextInput(readonly=True, font_size='14sp', size_hint_y=None)
        self.text_input.bind(minimum_height=self.text_input.setter('height'))
        scrollview.add_widget(self.text_input)
        layout.add_widget(scrollview)
        close_button = Button(text="Close", size_hint_y=None, height='48dp', background_normal='', background_color=[0.5, 0, 0, 1], color=[1, 1, 1, 1])
        close_button.bind(on_release=self.dismiss)
        layout.add_widget(close_button)
        self.add_widget(layout)
        self.load_readme()

    def load_readme(self):
        try:
            with open('README.md', 'r') as file:
                self.text_input.text = file.read()
        except Exception as e:
            self.text_input.text = f"Failed to load README file:\n{str(e)}"

class MyApp(App):
    filechooser_popup_open = False
    API_ENDPOINTS = API_ENDPOINTS
    OUTPUT_PATHS = OUTPUT_PATHS
    upload_in_progress = BooleanProperty(False)

    def build(self):
        self.title = "Structure-Data Format (SDF) File Processor"
        self.root = Builder.load_file('styles.kv')
        self.label = self.root.ids.label
        self.button_select = self.root.ids.button_select
        self.button_upload = self.root.ids.button_upload
        self.button_stop_upload = self.root.ids.button_stop_upload
        self.button_readme = self.root.ids.button_readme
        self.button_clear_folders = self.root.ids.button_clear_folders
        self.terminal_output = self.root.ids.terminal_output

        self.button_select.bind(on_release=self.show_filechooser)
        self.button_upload.bind(on_release=self.upload_files)
        self.button_stop_upload.bind(on_release=self.stop_upload)
        self.button_readme.bind(on_release=self.show_readme)
        self.button_clear_folders.bind(on_release=self.clear_output_folders)
        
        self.upload_in_progress = False  

        if not api_key:
            self.print_terminal("Error: API key not found. Please ensure '.env' file contains a valid API key.")
            self.button_upload.disabled = True

        return self.root
    
    def show_filechooser(self, instance=None):
        if MyApp.filechooser_popup_open:
            self.print_terminal("File chooser is already open, not opening another one.")
            return
        MyApp.filechooser_popup_open = True
        self.filechooser_popup = FileChooserPopup(select_callback=self.process_files)
        self.filechooser_popup.bind(on_dismiss=self.on_filechooser_dismiss)
        self.filechooser_popup.open()
    
    def on_filechooser_dismiss(self, instance):
        MyApp.filechooser_popup_open = False

    def show_readme(self, instance=None):
        self.readme_popup = ReadmePopup()
        self.readme_popup.open()

    def clear_output_folders(self, instance=None):
        # Use the correct keys from OUTPUT_PATHS
        log_files = [
            OUTPUT_PATHS['success_log'],
            OUTPUT_PATHS['duplicate_log'],
            OUTPUT_PATHS['failed_log'],
            OUTPUT_PATHS['general_log']
        ]

        for log_file in log_files:
            try:
                if os.path.exists(log_file):
                    os.remove(log_file)
            except Exception as e:
                self.print_terminal(f"Error deleting log file {log_file}: {str(e)}")
        
        self.print_terminal("Output logs cleared.")

    def process_files(self, files):
        self.selected_files = files
        self.label.text = f'Selected {len(files)} files.'
        file_list = self.root.ids.file_list
        file_list.clear_widgets()
        for file in files:
            file_list.add_widget(Label(text=file, size_hint_y=None, height='30dp', font_size='12sp', color=[0, 0, 0, 1]))
        self.button_upload.background_color = [0.04, 0.33, 0.64, 1]
        self.button_upload.disabled = False
        self.print_terminal(f'Selected files: {files}')
        
        if hasattr(self, 'filechooser_popup'):
            self.filechooser_popup.filechooser.selection = []

    def upload_files(self, instance=None):
        if self.upload_in_progress:
            self.print_terminal("Upload already in progress, skipping redundant call")
            return

        self.upload_in_progress = True
        self.button_stop_upload.disabled = False

        def update_progress_bar(dt):
            pass
        
        self.schedule_upload(0, update_progress_bar)

    def schedule_upload(self, file_index, update_progress_bar):
        if not self.upload_in_progress or file_index >= len(self.selected_files):
            self.upload_in_progress = False
            self.button_stop_upload.disabled = True
            return
        
        file = self.selected_files[file_index]
        molecules, fragments = process_sdf([file], self.print_terminal)
        if not molecules:
            self.print_terminal(f"No valid molecules found in file: {file}")
            self.schedule_next_file(file_index, update_progress_bar)
            return

        def process_molecule(molecule_index):
            if not self.upload_in_progress or molecule_index >= len(molecules):
                Clock.schedule_once(update_progress_bar)
                self.schedule_next_file(file_index, update_progress_bar)
                return

            molecule_data = molecules[molecule_index]
            fragment_data = fragments[molecule_index] if molecule_index < len(fragments) else {}
            success = post_to_api(molecule_data, fragment_data, file, self.print_terminal, api_key, OUTPUT_PATHS)
            if success:
                self.print_terminal(f'Compound successfully uploaded: {molecule_data.get("MolecularFormula", "")}')
            else:
                self.print_terminal(f'Failed to upload compound: {molecule_data.get("MolecularFormula", "")}')
            
            Clock.schedule_once(lambda dt: process_molecule(molecule_index + 1))

        process_molecule(0)

    def schedule_next_file(self, file_index, update_progress_bar):
        if file_index + 1 >= len(self.selected_files):
            self.upload_in_progress = False
            self.button_stop_upload.disabled = True
            self.print_terminal("End of Upload")
        else:
            Clock.schedule_once(lambda dt: self.schedule_upload(file_index + 1, update_progress_bar))

    def stop_upload(self, instance=None):
        self.upload_in_progress = False
        self.print_terminal("Upload stopped.")

    def print_terminal(self, message):
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_message = f"{timestamp} - {message}"
        
        with open(OUTPUT_PATHS['general_log'], 'a') as f:
            f.write(log_message + '\n')
            f.flush()
        self.terminal_output.text += log_message + '\n'

if __name__ == '__main__':
    MyApp().run()
