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
from rdkit import Chem
# from cdxml_converter import convert_mol_to_cdxml

Window.size = (1200, 600)

load_dotenv()
api_key = os.getenv('API_KEY')
kivy.logger.Logger.setLevel(logging.ERROR)
try:
    Builder.load_file('styles.kv')
except Exception as e:
    print(f"Error loading Kv file: {e}")

# Load configuration files
def load_config(file_path, default):
    try:
        with open(file_path, 'r') as file:
            return json.load(file)
    except FileNotFoundError:
        return default
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from {file_path}: {e}")
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

class CustomFileChooserIconView(FileChooserIconView):
    selection = ListProperty([])

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        Window.bind(on_key_down=self._on_key_down)

    def _on_key_down(self, instance, keyboard, keycode, text, modifiers):
        if 'ctrl' in modifiers and text == 'a':
            self.select_all_files()

    def select_all_files(self):
        all_files = [f for f in self.files if f.lower().endswith('.sdf')]
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

        # Ensure buttons are bound only once
        self.button_select.bind(on_release=self.show_filechooser)
        self.button_upload.bind(on_release=self.upload_files)
        self.button_stop_upload.bind(on_release=self.stop_upload)
        self.button_readme.bind(on_release=self.show_readme)
        self.button_clear_folders.bind(on_release=self.clear_output_folders)
        
        self.upload_in_progress = False  

        # Check if the API key is set
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
        # clear output folders
        for folder in folders:
            if os.path.exists(folder):
                for filename in os.listdir(folder):
                    file_path = os.path.join(folder, filename)
                    try:
                        if os.path.isfile(file_path):
                            os.unlink(file_path)
                    except Exception as e:
                        self.print_terminal(f"Error deleting file {file_path}: {str(e)}")
        # clear logs.txt
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

        # self.print_terminal("Upload button pressed")

        def update_progress_bar(dt):
            pass

        # self.print_terminal("Uploading files to API...")
        
        # Schedule the upload process to allow frequent checking of the stop button condition
        self.schedule_upload(0, update_progress_bar)
        
    def schedule_upload(self, file_index, update_progress_bar):
        if not self.upload_in_progress or file_index >= len(self.selected_files):
            self.upload_in_progress = False
            self.button_stop_upload.disabled = True
            return
        
        file = self.selected_files[file_index]
        molecules = process_sdf([file], self.print_terminal)
        if len(molecules) == 0:
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
                self.print_terminal(f'Successfully uploaded compound from file: {file}')
            else:
                self.print_terminal(f'Failed to upload compound from file: {file}')
            
            Clock.schedule_once(lambda dt: process_molecule(molecule_index + 1))

        process_molecule(0)

    def schedule_next_file(self, file_index, update_progress_bar):
        Clock.schedule_once(lambda dt: self.schedule_upload(file_index + 1, update_progress_bar))

    def stop_upload(self, instance=None):
        self.upload_in_progress = False
        self.print_terminal("Upload stopped.")

    def print_terminal(self, message):
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_message = f"{timestamp} - {message}"
        
        # Write to log file specified in OUTPUT_PATHS JSON file
        with open(OUTPUT_PATHS['general_log'], 'a') as f:
            f.write(log_message + '\n')
            f.flush()
        self.terminal_output.text += log_message + '\n'

import requests

def post_to_api(molecule_data, file, callback, api_key, OUTPUT_PATHS):
    if not api_key:
        callback("Error: API key not found. Cannot proceed with the upload.")
        return False

    headers = {
        'accept': 'application/vnd.api+json',
        'x-api-key': api_key,
        'Content-Type': 'application/vnd.api+json',
    }

    try:
        # Always attempt to upload salts and solvates
        salt_id, solvate_id = upload_fragments(molecule_data, callback, headers)
    except Exception as e:
        callback(f"Error uploading fragments: {str(e)}")
        return False

    try:
        data = construct_payload(molecule_data, salt_id, solvate_id)
    except Exception as e:
        callback(f"Error constructing payload: {str(e)}")
        return False

    if not data:
        callback("Invalid payload data")
        return False

    try:
        return send_request(data, file, callback, API_ENDPOINTS['Compound Endpoint'], headers, OUTPUT_PATHS)
    except Exception as e:
        callback(f"Error sending request: {str(e)}")
        return False


def upload_fragments(molecule_data, callback, headers):
    salt_id = None
    solvate_id = None

    try:
        if 'Salt_name' in molecule_data:
            salt_data = {
                "data": {
                    "type": "salt",
                    "attributes": {
                        "name": molecule_data.get('Salt_name', ''),
                        "mf": molecule_data.get('Formula', ''),
                        "mw": f"{molecule_data.get('MW_salt', '')} g/mol"
                    }
                }
            }
            salt_id = upload_fragment(salt_data, "salts", callback, headers)
            if salt_id:
                callback(f"Successfully uploaded salt: {molecule_data.get('Salt_name', '')}")
            else:
                callback(f"Failed to upload salt: {molecule_data.get('Salt_name', '')}")
    except Exception as e:
        callback(f"No Salt_name found: {str(e)}")
    try:
        if 'Solvate_name' in molecule_data:
            solvate_data = {
                "data": {
                    "type": "solvate",
                    "attributes": {
                        "name": molecule_data.get('Solvate_name', ''),
                        "mf": molecule_data.get('Formula', ''),
                        "mw": f"{molecule_data.get('MW', '')} g/mol"
                    }
                }
            }
            solvate_id = upload_fragment(solvate_data, "solvates", callback, headers)
            if solvate_id:
                callback(f"Successfully uploaded solvate: {molecule_data.get('Solvate_name', '')}")
            else:
                callback(f"Failed to upload solvate: {molecule_data.get('Solvate_name', '')}")
    except Exception as e:
        callback(f"No Solvate_name found: {str(e)}")

    return salt_id, solvate_id

def upload_fragment(fragment_data, fragment_type, callback, headers):
    fragment_endpoint = f"{API_ENDPOINTS['Fragment Endpoint']}/{fragment_type}"
    response = requests.post(fragment_endpoint, data=json.dumps(fragment_data), headers=headers)
    if response.status_code == 200:
        fragment_id = response.json()['data']['id']
        return fragment_id
    else:
        callback(f"Failed to upload {fragment_type}. Status code: {response.status_code}, response: {response.text}")
        return None
    
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
                    molecule_data["cdxml"] = convert_mol_to_cdxml(molecule_data)

                    callback(f"Molecule Data: {molecule_data}")  # Log molecule data for debugging
                    molecules.append(molecule_data)
                else:
                    callback(f"Error: Molecule in file {sdf_file} could not be parsed and will be skipped.")
        except Exception as e:
            callback(f"Error processing file {sdf_file}: {str(e)}")

    callback(f"Total molecules extracted: {len(molecules)}")
    return molecules

def construct_payload(molecule_data, salt_id, solvate_id):
    data = {
        "data": {
            "type": "asset",
            "attributes": {
                "synonyms": [
                    molecule_data.get('Chemical name', ''),
                    molecule_data.get('Smile', ''),
                    molecule_data.get('Formula', '')
                ],
                "fields": [
                    {
                        "id": "5d6e0287ee35880008c18d6d", # Chemical Structure, CDXML
                        "value": molecule_data.get("cdxml", "")
                    },
                    {
                        "id": "62f9fe5b74770f14d1de43a8", # Stereochemistry, not specified in ENAMINE sdf files
                        "value": "No stereochemistry"
                    },
                    {
                        "id": "5d6e0287ee35880008c18db6", # Chemical Name
                        "value": molecule_data.get("Chemical name", "")
                    },
                    {
                        "id": "5d6e0287ee35880008c18db7", # Molecular Weight
                        "value": {
                            "rawValue": molecule_data.get("MW", ""),
                            "displayValue": f"{molecule_data.get('MW', '')} g/mol"
                        }
                    },
                    {
                        "id": "5d6e0287ee35880008c18db8", # Exact Mass
                        "value": float(molecule_data.get("Amount_mg", 0))
                    },
                    {
                        "id": "5d6e0287ee35880008c18db9", # Molecular Formula
                        "value": {
                            "rawValue": molecule_data.get("Formula", ""),
                            "displayValue": molecule_data.get("Formula", "")
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
                                "id": "62fcceeb19660304d1e5beee", # Source: Internal, Collaboration, Acquired
                                "value": "Acquired"
                                },
                                {
                                "id": "63469c69ed8a726a31923537", # Project
                                "value": "Unspecified"
                                },
                                {
                                "id": "62fcceeb19660304d1e5bef1", # Synthesis Date, ISO 8601 format: 2011-10-10T14:48:00Z
                                "value": "2011-10-10T14:48:00Z"
                                },
                                {
                                "id": "6384a1270d28381d21deaca7", # Chemist
                                "value": "TestUser MCChemist"
                                },
                                {
                                "id": "62fa096d19660304d1e5b2db", # Batch Purpose, 
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

    # Only add fragments if they exist
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
                "coefficient": 1 # Where does this value come from?
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
                "coefficient": 1 # Where does this value come from?
            }]
        data["data"]["relationships"]["batch"]["data"]["attributes"]["fragments"] = fragments

    print(data)  # Log the constructed payload for debugging
    return data


def send_request(data, file, callback, endpoint, headers, OUTPUT_PATHS):
    # Write data to json file
    with open('data.json', 'w') as f:
        f.write(json.dumps(data, indent=4))
    data_json = json.dumps(data)
    # callback(f"Sending POST request for molecule: {data['data']['attributes']['synonyms'][0]} to endpoint: {endpoint}")
    # callback(f"POST payload: {data_json}")

    try:
        response = requests.post(endpoint, data=data_json, headers=headers)
        # callback(f"Response status code: {response.status_code}, response text: {response.text}")

        if response.status_code == 200:
            handle_success(file, data, OUTPUT_PATHS, callback)
            return True
        else:
            handle_failure(file, data, response, OUTPUT_PATHS, callback)
            return False
    except requests.exceptions.RequestException as e:
        callback(f"Network error during request: {str(e)}")
        return False

def handle_success(file, data, OUTPUT_PATHS, callback):
    if not os.path.exists(OUTPUT_PATHS['successful_files']):
        os.makedirs(OUTPUT_PATHS['successful_files'])
    shutil.copy(file, OUTPUT_PATHS['successful_files'])
    
    with open(OUTPUT_PATHS['success_log'], 'a') as success_log:
        success_log.write(f"File: {file}, Molecule: {data['data']['attributes']['synonyms'][0]}, API Response Status Code: 200\n")
    
    callback(f"Successfully uploaded compound from {file}. Copying file to '{OUTPUT_PATHS['successful_files']}' folder.")

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
    
    callback(f"API request failed for {file} with status code {response.status_code}, response: {response.text}. Copying file to '{log_folder}' folder.")
    # Write the response to a file named "response.json" as formatted JSON
    try:
        response_json = response.json()
        with open(f"{log_folder}/response.json", 'w') as f:
            json.dump(response_json, f, indent=4)
    except ValueError:
        # In case the response is not valid JSON, write the raw text
        with open(f"{log_folder}/response.json", 'w') as f:
            f.write(response.text)

def convert_mol_to_cdxml(molecule_data):
    # Template for the CDXML format
    cdxml_template = """<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd" >
<CDXML
 CreationProgram="ChemDraw JS 23.2.0.0"
 BoundingBox="347.74 95.77 502.96 153.53"
 WindowPosition="0 0"
 WindowSize="0 0"
 FractionalWidths="yes"
 InterpretChemically="yes"
 ShowAtomQuery="yes"
 ShowAtomStereo="no"
 ShowAtomEnhancedStereo="yes"
 ShowAtomNumber="no"
 ShowResidueID="no"
 ShowBondQuery="yes"
 ShowBondRxn="yes"
 ShowBondStereo="no"
 ShowTerminalCarbonLabels="no"
 ShowNonTerminalCarbonLabels="no"
 HideImplicitHydrogens="no"
 Magnification="666"
 LabelFont="24"
 LabelSize="10"
 LabelFace="96"
 CaptionFont="24"
 CaptionSize="10"
 HashSpacing="2.50"
 MarginWidth="1.60"
 LineWidth="0.60"
 BoldWidth="2"
 BondLength="14.40"
 BondSpacing="18"
 ChainAngle="120"
 LabelJustification="Auto"
 CaptionJustification="Left"
 AminoAcidTermini="HOH"
 ShowSequenceTermini="yes"
 ShowSequenceBonds="yes"
 ShowSequenceUnlinkedBranches="no"
 ResidueWrapCount="40"
 ResidueBlockCount="10"
 PrintMargins="36 36 36 36"
 MacPrintInfo="0003000000480048000000000318026400000000031802640367052803FC00020000004800480000000003180264000100000064000000010001010100000001270F000100010000000000000000000000000002001901900000000000400000000000000000000100000000000000000000000000000000"
 ChemPropName=""
 ChemPropFormula="Chemical Formula: "
 ChemPropExactMass="Exact Mass: "
 ChemPropMolWt="Molecular Weight: "
 ChemPropMOverZ="m/z: "
 ChemPropAnalysis="Elemental Analysis: "
 ChemPropBoilingPt="Boiling Point: "
 ChemPropMeltingPt="Melting Point: "
 ChemPropCritTemp="Critical Temp: "
 ChemPropCritPres="Critical Pres: "
 ChemPropCritVol="Critical Vol: "
 ChemPropGibbs="Gibbs Energy: "
 ChemPropLogP="Log P: "
 ChemPropMR="MR: "
 ChemPropHenry="Henry&apos;s Law: "
 ChemPropEForm="Heat of Form: "
 ChemProptPSA="tPSA: "
 ChemPropID=""
 ChemPropFragmentLabel=""
 color="0"
 bgcolor="1"
 RxnAutonumberStart="1"
 RxnAutonumberConditions="no"
 RxnAutonumberStyle="Roman"
 RxnAutonumberFormat="(#)"
 MonomerRenderingStyle="graphic"
><colortable>
<color r="1" g="1" b="1"/>
<color r="0" g="0" b="0"/>
<color r="1" g="0" b="0"/>
<color r="1" g="1" b="0"/>
<color r="0" g="1" b="0"/>
<color r="0" g="1" b="1"/>
<color r="0" g="0" b="1"/>
<color r="1" g="0" b="1"/>
</colortable><fonttable>
<font id="24" charset="utf-8" name="Arial"/>
</fonttable><page
 id="477"
 BoundingBox="0 0 629.33 196"
 Width="629.33"
 Height="196"
 HeaderPosition="36"
 FooterPosition="36"
 PageOverlap="0"
 PrintTrimMarks="yes"
 HeightPages="1"
 WidthPages="2"
 DrawingSpace="poster"
><fragment
 id="1"
 BoundingBox="347.74 95.77 502.96 153.53"
 Z="1"
>{atoms}
{bonds}
</fragment></page></CDXML>
"""

    atoms = []
    bonds = []

    # Add atoms
    atom_template = '<n\n id="{id}"\n p="{x} {y}"\n Z="{z}"\n AS="{symbol}"\n AtomID="{atom_id}"\n Element="{element}"\n NumHydrogens="{num_hydrogens}"\n NeedsClean="{needs_clean}"\n/>'
    for i, atom in enumerate(molecule_data["atom_block"]):
        element = atom["symbol"]
        num_hydrogens = 0  # You can adjust this if you have hydrogen information
        needs_clean = "yes" if element in ["N", "O", "F", "S"] else "no"
        atoms.append(atom_template.format(id=i + 1, x=atom["x"], y=atom["y"], z=1, symbol=element, atom_id=i + 1, element=element, num_hydrogens=num_hydrogens, needs_clean=needs_clean))

    # Add bonds
    bond_template = '<b\n id="{id}"\n Z="{z}"\n B="{begin}"\n E="{end}"\n BS="N"\n Order="{order}"\n/>'
    for i, bond in enumerate(molecule_data["bond_block"]):
        bonds.append(bond_template.format(id=i + 1, z=i + 1, begin=bond["begin_atom_idx"] + 1, end=bond["end_atom_idx"] + 1, order=int(bond["bond_type"])))

    # Combine atoms and bonds into the final CDXML content
    cdxml_content = cdxml_template.format(atoms="\n".join(atoms), bonds="\n".join(bonds))

    
    import xml.etree.ElementTree as ET

    # Parse the CDXML content
    root = ET.fromstring(cdxml_content)

    # Count the number of carbon atoms (AS="C")
    carbon_atoms = root.findall(".//n[@AS='C']")
    carbon_count = len(carbon_atoms)
    nitrogen_atoms = root.findall(".//n[@AS='N']")
    nitrogen_count = len(nitrogen_atoms)
    oxygen_atoms = root.findall(".//n[@AS='O']")
    oxygen_count = len(oxygen_atoms)
    sulfur_atoms = root.findall(".//n[@AS='S']")
    sulfur_count = len(sulfur_atoms)
    fluorine_atoms = root.findall(".//n[@AS='F']")
    fluorine_count = len(fluorine_atoms)

    print(f"Number of carbon atoms: {carbon_count}")
    print(f"Number of nitrogen atoms: {nitrogen_count}")
    print(f"Number of oxygen atoms: {oxygen_count}")
    print(f"Number of sulfur atoms: {sulfur_count}")
    print(f"Number of fluorine atoms: {fluorine_count}")

    return cdxml_content

if __name__ == '__main__':
    molecule_data = {
        "atom_block": [
            {"symbol": "C", "x": 0.7145, "y": -0.4125, "z": 0.0000},
            {"symbol": "N", "x": 0.7145, "y": -1.2375, "z": 0.0000},
            {"symbol": "C", "x": 1.3819, "y": -1.7224, "z": 0.0000},
            {"symbol": "S", "x": 2.1665, "y": -1.4675, "z": 0.0000},
            {"symbol": "C", "x": 2.7796, "y": -2.0195, "z": 0.0000},
            {"symbol": "C", "x": 3.5642, "y": -1.7646, "z": 0.0000},
            {"symbol": "O", "x": 3.7358, "y": -0.9576, "z": 0.0000},
            {"symbol": "N", "x": 4.1773, "y": -2.3166, "z": 0.0000},
            {"symbol": "C", "x": 4.9620, "y": -2.0617, "z": 0.0000},
            {"symbol": "C", "x": 5.5751, "y": -2.6137, "z": 0.0000},
            {"symbol": "O", "x": 5.4035, "y": -3.4207, "z": 0.0000},
            {"symbol": "C", "x": 4.6189, "y": -3.6756, "z": 0.0000},
            {"symbol": "C", "x": 4.0058, "y": -3.1236, "z": 0.0000},
            {"symbol": "N", "x": 1.1270, "y": -2.5070, "z": 0.0000},
            {"symbol": "N", "x": 0.3020, "y": -2.5070, "z": 0.0000},
            {"symbol": "C", "x": 0.0470, "y": -1.7224, "z": 0.0000},
            {"symbol": "C", "x": -0.7376, "y": -1.4675, "z": 0.0000},
            {"symbol": "C", "x": -1.3507, "y": -2.0195, "z": 0.0000},
            {"symbol": "C", "x": -2.1353, "y": -1.7646, "z": 0.0000},
            {"symbol": "C", "x": -2.3068, "y": -0.9576, "z": 0.0000},
            {"symbol": "F", "x": -3.0915, "y": -0.7027, "z": 0.0000},
            {"symbol": "C", "x": -1.6937, "y": -0.4056, "z": 0.0000},
            {"symbol": "C", "x": -0.9091, "y": -0.6605, "z": 0.0000}
        ],
        "bond_block": [
            {"begin_atom_idx": 0, "end_atom_idx": 1, "bond_type": 1},
            {"begin_atom_idx": 1, "end_atom_idx": 2, "bond_type": 1},
            {"begin_atom_idx": 2, "end_atom_idx": 3, "bond_type": 1},
            {"begin_atom_idx": 3, "end_atom_idx": 4, "bond_type": 1},
            {"begin_atom_idx": 4, "end_atom_idx": 5, "bond_type": 1},
            {"begin_atom_idx": 5, "end_atom_idx": 6, "bond_type": 2},
            {"begin_atom_idx": 5, "end_atom_idx": 7, "bond_type": 1},
            {"begin_atom_idx": 7, "end_atom_idx": 8, "bond_type": 1},
            {"begin_atom_idx": 8, "end_atom_idx": 9, "bond_type": 1},
            {"begin_atom_idx": 9, "end_atom_idx": 10, "bond_type": 1},
            {"begin_atom_idx": 10, "end_atom_idx": 11, "bond_type": 1},
            {"begin_atom_idx": 11, "end_atom_idx": 12, "bond_type": 1},
            {"begin_atom_idx": 7, "end_atom_idx": 12, "bond_type": 1},
            {"begin_atom_idx": 2, "end_atom_idx": 13, "bond_type": 2},
            {"begin_atom_idx": 13, "end_atom_idx": 14, "bond_type": 1},
            {"begin_atom_idx": 14, "end_atom_idx": 15, "bond_type": 2},
            {"begin_atom_idx": 1, "end_atom_idx": 15, "bond_type": 1},
            {"begin_atom_idx": 15, "end_atom_idx": 16, "bond_type": 1},
            {"begin_atom_idx": 16, "end_atom_idx": 17, "bond_type": 1},
            {"begin_atom_idx": 17, "end_atom_idx": 18, "bond_type": 2},
            {"begin_atom_idx": 18, "end_atom_idx": 19, "bond_type": 1},
            {"begin_atom_idx": 19, "end_atom_idx": 20, "bond_type": 1},
            {"begin_atom_idx": 19, "end_atom_idx": 21, "bond_type": 2},
            {"begin_atom_idx": 21, "end_atom_idx": 22, "bond_type": 1},
            {"begin_atom_idx": 16, "end_atom_idx": 22, "bond_type": 2}
        ]
    }

    # Generate CDXML
    cdxml_content = convert_mol_to_cdxml(molecule_data)
    print(cdxml_content)
    # MyApp().run()


