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
    "successful_files": "Successfully uploaded files",
    "failed_files": "Failed files",
    "duplicate_files": "Failed files/Duplicate compounds",
    "success_log": "Successfully uploaded files/success.txt",
    "duplicate_log": "Failed files/Duplicate compounds/duplicates.txt",
    "general_log": "logs.txt"
})

API_ENDPOINTS = load_config('config/api_endpoints.json', {
    "Fragment Endpoint": "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/fragments",
    "Compound Endpoint": "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/materials/Compounds/assets"
})

# Helper functions
def log_duplicate(file, molecule_data, callback, orm_code):
    if not os.path.exists(OUTPUT_PATHS['duplicate_files']):
        os.makedirs(OUTPUT_PATHS['duplicate_files'])
    shutil.copy(file, OUTPUT_PATHS['duplicate_files'])

    with open(OUTPUT_PATHS['duplicate_log'], 'a') as duplicate_log:
        duplicate_log.write(f"File: {file}, Molecule: {molecule_data.get('MolecularFormula', '')}, ORM Code: {orm_code}\n")

    callback(f"Duplicate: {molecule_data.get('MolecularFormula', '')}, ORM Code: {orm_code}. Copying file to '{OUTPUT_PATHS['duplicate_files']}' folder.")

# API functions
def upload_fragment(fragment_data, fragment_type, callback, headers):
    fragment_endpoint = f"{API_ENDPOINTS['Fragment Endpoint']}/{fragment_type}"
    response = requests.post(fragment_endpoint, data=json.dumps(fragment_data), headers=headers)
    if response.status_code in [200, 201]:
        fragment_id = response.json()['data']['id']
        return fragment_id
    callback(f"Failed to upload {fragment_type}. Status code: {response.status_code}, response: {response.text}")
    return None

def get_existing_fragment_id(fragment_name, fragment_type, api_key):
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
        if results:
            return results[0]['id']
    return None

def upload_fragments(molecule_data, callback, headers, api_key):
    salt_id = solvate_id = None

    # Check for salt presence
    if 'Salt_name' in molecule_data or 'Salt_Name' in molecule_data:
        salt_name = molecule_data.get('Salt_name', molecule_data.get('Salt_Name', ''))
        compound_name = molecule_data.get('MolecularFormula', '')
        callback(f"Salt for compound {compound_name} detected: {salt_name}")

        salt_id = get_existing_fragment_id(salt_name, "salts", api_key)
        if salt_id:
            callback(f"Salt '{salt_name}' already exists with ID: {salt_id}")
        else:
            salt_data = {
                "data": {
                    "type": "salt",
                    "attributes": {
                        "name": salt_name,
                        "mf": molecule_data.get('Formula', ''),
                        "mw": f"{molecule_data.get('MW_salt', '')} g/mol"
                    }
                }
            }
            salt_id = upload_fragment(salt_data, "salts", callback, headers)
            if salt_id:
                callback(f"Successfully uploaded salt: {salt_name}")
            else:
                callback(f"Failed to upload salt: {salt_name}")

    # Check for solvate presence
    if 'Solvate_name' in molecule_data or 'Solvate_Name' in molecule_data:
        solvate_name = molecule_data.get('Solvate_name', molecule_data.get('Solvate_Name', ''))
        compound_name = molecule_data.get('MolecularFormula', '')
        callback(f"Solvate for compound {compound_name} detected: {solvate_name}")
        
        solvate_id = get_existing_fragment_id(solvate_name, "solvates", api_key)
        if solvate_id:
            callback(f"Solvate '{solvate_name}' already exists with ID: {solvate_id}")
        else:
            solvate_data = {
                "data": {
                    "type": "solvate",
                    "attributes": {
                        "name": solvate_name,
                        "mf": molecule_data.get('Formula', ''),
                        "mw": f"{molecule_data.get('MW', '')} g/mol"
                    }
                }
            }
            solvate_id = upload_fragment(solvate_data, "solvates", callback, headers)
            if solvate_id:
                callback(f"Successfully uploaded solvate: {solvate_name}")
            else:
                callback(f"Failed to upload solvate: {solvate_name}")

    return salt_id, solvate_id

# SDF processing functions using OpenBabel
def process_sdf(files, callback):
    molecules = []
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("sdf", "mol")

    obConversion_smiles = openbabel.OBConversion()
    obConversion_smiles.SetOutFormat("smiles")

    for sdf_file in files:
        callback(f"Processing file: {sdf_file}")
        obMol = openbabel.OBMol()
        
        not_at_end = obConversion.ReadFile(obMol, sdf_file)
        while not_at_end:
            try:
                # Extract molecular data
                smiles = obConversion_smiles.WriteString(obMol).strip().upper()
                molecule_data = {
                    "Chemical name": obMol.GetData("Chemical name").GetValue() if obMol.HasData("Chemical name") else '',
                    "MolecularFormula": obMol.GetFormula().upper(),
                    "MW": obMol.GetMolWt(),
                    "Smile": smiles,
                    "Amount_mg": obMol.GetData("Amount_mg").GetValue() if obMol.HasData("Amount_mg") else 0,
                    "ID": obMol.GetData("ID").GetValue() if obMol.HasData("ID") else '',
                    "Formula": obMol.GetData("Formula").GetValue() if obMol.HasData("Formula") else '',
                    "Purity": obMol.GetData("Purity").GetValue() if obMol.HasData("Purity") else '',
                    "PO": obMol.GetData("PO").GetValue() if obMol.HasData("PO") else '',
                    "Salt_name": obMol.GetData("Salt_name").GetValue() if obMol.HasData("Salt_name") else obMol.GetData("Salt_Name").GetValue() if obMol.HasData("Salt_Name") else '',                    "Salt_ratio": obMol.GetData("Salt_ratio").GetValue() if obMol.HasData("Salt_ratio") else '',
                    "MW_salt": obMol.GetData("MW_salt").GetValue() if obMol.HasData("MW_salt") else '',
                    "Plate_ID": obMol.GetData("Plate_ID").GetValue() if obMol.HasData("Plate_ID") else '',
                    "Well": obMol.GetData("Well").GetValue() if obMol.HasData("Well") else '',
                    "Barcode": obMol.GetData("Barcode").GetValue() if obMol.HasData("Barcode") else ''
                }
                
                # Extract atom and bond blocks
                atom_block = []
                for atom in openbabel.OBMolAtomIter(obMol):
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
                
                bond_block = [{
                    "begin_atom_idx": bond.GetBeginAtomIdx() - 1,
                    "end_atom_idx": bond.GetEndAtomIdx() - 1,
                    "bond_type": bond.GetBondOrder()
                } for bond in openbabel.OBMolBondIter(obMol)]
                
                molecule_data.update({"atom_block": atom_block, "bond_block": bond_block})
                molecule_data["cdxml"] = convert_mol_to_cdxml(molecule_data)
                
                molecules.append(molecule_data)
            except Exception as e:
                callback(f"Error processing molecule in file {sdf_file}: {str(e)}")
            
            not_at_end = obConversion.Read(obMol)
    
    callback(f"Total molecules extracted: {len(molecules)}")
    return molecules

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

def construct_payload(molecule_data, salt_id, solvate_id):
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
                            "rawValue": molecule_data.get("MW", ""),
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
                                "value": "Acquired"
                                },
                                {
                                "id": "63469c69ed8a726a31923537",
                                "value": "Unspecified"
                                },
                                {
                                "id": "62fcceeb19660304d1e5bef1",
                                "value": "2011-10-10T14:48:00Z"
                                },
                                {
                                "id": "6384a1270d28381d21deaca7",
                                "value": "TestUser MCChemist"
                                },
                                {
                                "id": "62fa096d19660304d1e5b2db",
                                "value": "Dummy compound"
                                },
                                {
                                "id": "62fa096d19660304d1e5b2da",
                                "value": "Discovery"
                                }
                            ]
                        }
                    }
                }
            }
        }
    }

    if salt_id or solvate_id:
        fragments = {}
        if salt_id:
            fragments["salts"] = [{
                "type": "SALT",
                "id": salt_id,
                "name": molecule_data.get('Salt_name', ''),
                "mf": molecule_data.get('Formula', ''),
                "mw": {
                    "rawValue": molecule_data.get('MW_salt', ''),
                    "displayValue": f"{molecule_data.get('MW_salt', '')} g/mol"
                },
                "coefficient": 1
            }]
        if solvate_id:
            fragments["solvates"] = [{
                "type": "SOLVATE",
                "id": solvate_id,
                "name": molecule_data.get('Solvate_name', ''),
                "mf": molecule_data.get('Formula', ''),
                "mw": {
                    "rawValue": molecule_data.get('MW', ''),
                    "displayValue": f"{molecule_data.get('MW', '')} g/mol"
                },
                "coefficient": 1
            }]
        data["data"]["relationships"]["batch"]["data"]["attributes"]["fragments"] = fragments

    return data

def send_request(data, file, callback, endpoint, headers, OUTPUT_PATHS):
    data_json = json.dumps(data)
    try:
        response = requests.post(endpoint, data=data_json, headers=headers)
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
                            "rawValue": molecule_data.get("MW", ""),
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
                                    "value": "Acquired"
                                },
                                {
                                    "id": "63469c69ed8a726a31923537",
                                    "value": "Unspecified"
                                },
                                {
                                    "id": "62fcceeb19660304d1e5bef1",
                                    "value": "2011-10-10T14:48:00Z"
                                },
                                {
                                    "id": "6384a1270d28381d21deaca7",
                                    "value": "TestUser MCChemist"
                                },
                                {
                                    "id": "62fa096d19660304d1e5b2db",
                                    "value": "Dummy compound"
                                },
                                {
                                    "id": "62fa096d19660304d1e5b2da",
                                    "value": "Discovery"
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

def post_to_api(molecule_data, file, callback, api_key, OUTPUT_PATHS):
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
    salt_id, solvate_id = upload_fragments(molecule_data, callback, headers, api_key)
    
    # Construct the payload
    payload = construct_payload(molecule_data, salt_id, solvate_id)
    
    # Send the request
    success = send_request(payload, file, callback, API_ENDPOINTS['Compound Endpoint'], headers, OUTPUT_PATHS)
    return success

def handle_success(file, data, OUTPUT_PATHS, callback):
    if not os.path.exists(OUTPUT_PATHS['successful_files']):
        os.makedirs(OUTPUT_PATHS['successful_files'])
    shutil.copy(file, OUTPUT_PATHS['successful_files'])
    
    with open(OUTPUT_PATHS['success_log'], 'a') as success_log:
        success_log.write(f"File: {file}, Molecule: {data['data']['attributes']['synonyms'][0]}, API Response Status Code: 200\n")
    
def handle_failure(file, data, response, OUTPUT_PATHS, callback):
    if "duplicate" in response.text.lower():
        log_folder = OUTPUT_PATHS['duplicate_files']
        log_file = OUTPUT_PATHS['duplicate_log']
    else:
        log_folder = OUTPUT_PATHS['failed_files']
        log_file = OUTPUT_PATHS['general_log']

    if not os.path.exists(log_folder):
        os.makedirs(log_folder)
    shutil.copy(file, log_folder)
    
    with open(log_file, 'a') as log:
        log.write(f"File: {file}, Molecule: {data['data']['attributes']['synonyms'][0]}, API Response Status Code: {response.status_code}, response text: {response.text}\n")
    
    callback(f"Failed: {data['data']['attributes']['synonyms'][0]} with status code {response.status_code}. Copying file to '{log_folder}' folder.")
    try:
        response_json = response.json()
        with open(f"{log_folder}/response.json", 'w') as f:
            json.dump(response_json, f, indent=4)
    except ValueError:
        with open(f"{log_folder}/response.json", 'w') as f:
            f.write(response.text)

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
        folders = [OUTPUT_PATHS['successful_files'], OUTPUT_PATHS['failed_files'], OUTPUT_PATHS['duplicate_files']]
        for folder in folders:
            if os.path.exists(folder):
                for filename in os.listdir(folder):
                    file_path = os.path.join(folder, filename)
                    try:
                        if os.path.isfile(file_path):
                            os.unlink(file_path)
                    except Exception as e:
                        self.print_terminal(f"Error deleting file {file_path}: {str(e)}")
        with open(OUTPUT_PATHS["general_log"], 'w') as f:
            f.write('')
        self.print_terminal("Output folders & logs cleared.")

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
        molecules = process_sdf([file], self.print_terminal)
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
            success = post_to_api(molecule_data, file, self.print_terminal, api_key, OUTPUT_PATHS)
            if success:
                self.print_terminal(f'Success: {molecule_data.get("MolecularFormula", "")}')
            else:
                self.print_terminal(f'Failed: {molecule_data.get("MolecularFormula", "")}')
            
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
