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
from kivy.uix.gridlayout import GridLayout
from kivy.uix.textinput import TextInput
from kivy.uix.spinner import Spinner
from kivy.core.window import Window
from kivy.properties import ListProperty, BooleanProperty
from kivy.lang import Builder
from kivy.clock import Clock
import json
import shutil

from api_utils import post_to_api
from sdf_utils import process_sdf

Window.size = (1200, 600)

load_dotenv()
api_key = os.getenv('API_KEY')
kivy.logger.Logger.setLevel(logging.ERROR)
try:
    Builder.load_file('styles.kv')
except Exception as e:
    print(f"Error loading Kv file: {e}")

# Load parameter sets from JSON
def load_parameter_sets(file_path='config/parameters.json'):
    try:
        with open(file_path, 'r') as file:
            return json.load(file)
    except FileNotFoundError: 
        return {
            'Salt': ['Salt_name', 'Formula', 'MW_salt']
        }

# Load API endpoints from JSON 
def load_api_endpoints(file_path='config/api_endpoints.json'):
    try:
        with open(file_path, 'r') as file:
            return json.load(file)
    except FileNotFoundError:
        return {
            "Default Endpoint": "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/fragments/salts"
        }
def load_output_paths(file_path='config/output_paths.json'):
    try:
        with open(file_path, 'r') as file:
            return json.load(file)
    except FileNotFoundError:
        return {
            "successful_files": "Successfully uploaded files",
            "failed_files": "Failed files",
            "duplicate_files": "Failed files/Duplicate compounds",
            "success_log": "Successfully uploaded files/success.txt",
            "duplicate_log": "Failed files/Duplicate compounds/duplicates.txt",
            "general_log": "logs.txt"
        }

OUTPUT_PATHS = load_output_paths()
PARAMETER_SETS = load_parameter_sets()
API_ENDPOINTS = load_api_endpoints()

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
    PARAMETER_SETS = PARAMETER_SETS
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
        self.param_spinner = self.root.ids.param_spinner

        # Ensure buttons are bound only once
        self.button_select.bind(on_release=self.show_filechooser)
        self.button_upload.bind(on_release=self.upload_files)
        self.button_stop_upload.bind(on_release=self.stop_upload)
        self.button_readme.bind(on_release=self.show_readme)
        self.button_clear_folders.bind(on_release=self.clear_output_folders)
        
        self.upload_in_progress = False  

        # Set the values for the spinners
        self.param_spinner.values = list(PARAMETER_SETS.keys())

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
        
        # Reset the file chooser selection state
        if hasattr(self, 'filechooser_popup'):
            self.filechooser_popup.filechooser.selection = []

    def upload_files(self, instance=None):
        if self.upload_in_progress:
            self.print_terminal("Upload already in progress, skipping redundant call")
            return

        self.upload_in_progress = True
        self.button_stop_upload.disabled = False

        selected_param_set = self.param_spinner.text

        # Automatically select endpoint based on parameter set
        endpoint = None
        if selected_param_set == "Create Salt":
            endpoint = "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/fragments/salts"
        elif selected_param_set == "Create Solvate":
            endpoint = "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/fragments/solvates"
        elif selected_param_set == "Create Compound (Does not work)":
            endpoint = "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/materials/Compounds/assets"
        else:
            self.print_terminal("Invalid parameter set selected")
            return

        params = PARAMETER_SETS[selected_param_set]
        self.print_terminal("Upload button pressed")

        def update_progress_bar(dt):
            pass

        self.print_terminal("Uploading files to API...")
        
        # Schedule the upload process to allow frequent checking of the stop condition
        # and allow multithreading in the future, if the API can keep up
        self.schedule_upload(0, params, endpoint, update_progress_bar)
        
    def schedule_upload(self, file_index, params, endpoint, update_progress_bar):
        if not self.upload_in_progress or file_index >= len(self.selected_files):
            self.upload_in_progress = False
            self.button_stop_upload.disabled = True
            return
        
        file = self.selected_files[file_index]
        molecules = process_sdf([file], params, self.print_terminal)
        print(molecules)
        if len(molecules) == 0:
            self.print_terminal(f"No valid molecules found in file: {file}")
            self.schedule_next_file(file_index, params, endpoint, update_progress_bar)
            return

        def process_molecule(molecule_index):
            if not self.upload_in_progress or molecule_index >= len(molecules):
                Clock.schedule_once(update_progress_bar)
                self.schedule_next_file(file_index, params, endpoint, update_progress_bar)
                return

            molecule_data = molecules[molecule_index]
            success = post_to_api(molecule_data, file, self.print_terminal, endpoint, self.param_spinner.text, api_key, OUTPUT_PATHS)
            if success:
                self.print_terminal(f'Successfully uploaded file: {file}')
            else:
                self.print_terminal(f'Failed to upload file: {file}')
            
            Clock.schedule_once(lambda dt: process_molecule(molecule_index + 1))

        process_molecule(0)

    def schedule_next_file(self, file_index, params, endpoint, update_progress_bar):
        Clock.schedule_once(lambda dt: self.schedule_upload(file_index + 1, params, endpoint, update_progress_bar))

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

import os
import requests
import shutil
import json

def post_to_api(molecule_data, file, callback, endpoint, param_set, api_key, OUTPUT_PATHS):
    print(f"PARAMS: {param_set}")
    if not api_key:
        callback("Error: API key not found. Cannot proceed with the upload.")
        return False

    headers = {
        'accept': 'application/vnd.api+json',
        'x-api-key': api_key,
        'Content-Type': 'application/vnd.api+json',
    }

    # Construct the API request payload based on the parameter set
    # Todo: load the structure from a json file
    if param_set == "Create Compound (Does not work)":
        data = {
            "data": {
                "type": "asset",
                "attributes": {
                    "synonyms": [
                        molecule_data.get('Chemical name', ''),
                        molecule_data.get('Smile', '')
                    ],
                    "fields": [
                        {
                            "id": "5eaa8792b5ed583958f447cb",
                            "value": molecule_data.get('Chemical name', '')
                        },
                        {
                            "id": "5eaa8792b5ed583958f447cc",
                            "value": {
                                "rawValue": molecule_data.get('MW', ''),
                                "displayValue": f"{molecule_data.get('MW', '')} g/mol"
                            }
                        },
                        {
                            "id": "5eaa8792b5ed583958f447cd",
                            "value": molecule_data.get('MW_salt', '0')
                        },
                        {
                            "id": "5eaa8792b5ed583958f447ce",
                            "value": {
                                "rawValue": molecule_data.get('Formula', ''),
                                "displayValue": molecule_data.get('Formula', '')
                            }
                        },
                        {
                            "id": "5eaa8792b5ed583958f44778",
                            "value": molecule_data.get('CDXML', '')
                        },
                        {
                            "id": "5eaa8792b5ed583958f447de",
                            "value": {
                                "filename": f"{molecule_data.get('Chemical name', '')}.png",
                                "base64": molecule_data.get('Base64', '')
                            }
                        }
                    ]
                },
                "relationships": {
                    "batch": {
                        "data": {
                            "type": "batch",
                            "id": "00011",
                            "attributes": {
                                "fields": [
                                    {
                                        "id": "5eaa8792b5ed583958f447d0",
                                        "value": {
                                            "rawValue": molecule_data.get('Amount_mg', ''),
                                            "displayValue": f"{molecule_data.get('Amount_mg', '')} mg"
                                        }
                                    },
                                    {
                                        "id": "5eaa8792b5ed583958f447d2",
                                        "value": {
                                            "rawValue": molecule_data.get('Purity', ''),
                                            "displayValue": f"{molecule_data.get('Purity', '')} %"
                                        }
                                    },
                                    {
                                        "id": "5eaa8792b5ed583958f447d3",
                                        "value": molecule_data.get('Chemical name', '')
                                    },
                                    {
                                        "id": "5eaa8792b5ed583958f447d4",
                                        "value": {
                                            "rawValue": molecule_data.get('Formula', ''),
                                            "displayValue": molecule_data.get('Formula', '')
                                        }
                                    },
                                    {
                                        "id": "5f02947e9a3a272654dfeb28",
                                        "value": {
                                            "rawValue": molecule_data.get('MW', ''),
                                            "displayValue": f"{molecule_data.get('MW', '')} g/mol"
                                        }
                                    }
                                ]
                            }
                        }
                    }
                }
            }
        }

        fragments = {}
        if 'Salt_name' in molecule_data:
            fragments["salts"] = [
                {
                    "type": "SALT",
                    "name": molecule_data.get('Salt_name', ''),
                    "mf": molecule_data.get('Salt smiles', ''),
                    "mw": {
                        "rawValue": molecule_data.get('MW_salt', ''),
                        "displayValue": f"{molecule_data.get('MW_salt', '')} g/mol"
                    },
                    "id": "SALT:17",
                    "coefficient": 1
                }
            ]
        if 'Solvate_name' in molecule_data:
            fragments["solvates"] = [
                {
                    "type": "SOLVATE",
                    "name": molecule_data.get('Solvate_name', ''),
                    "mf": molecule_data.get('Solvate smiles', ''),
                    "mw": {
                        "rawValue": molecule_data.get('MW_solvate', ''),
                        "displayValue": f"{molecule_data.get('MW_solvate', '')} g/mol"
                    },
                    "id": "Solvate:17",
                    "coefficient": 1
                }
            ]
        
        if fragments:
            data["data"]["relationships"]["batch"]["data"]["attributes"]["fragments"] = fragments

    elif param_set == "Create Salt":
        data = {
            "data": {
                "type": "salt",
                "attributes": {
                    "name": molecule_data.get('Salt_name', ''),
                    "mf": molecule_data.get('Formula', ''),
                    "mw": f"{molecule_data.get('MW_salt', '')} g/mol"
                }
            }
        }
    elif param_set == "Create Solvate":
        data = {
            "data": {
                "type": "solvate",
                "attributes": {
                    "name": molecule_data.get('Solvate_name', ''),
                    "mf": molecule_data.get('Formula', ''),
                    "mw": f"{molecule_data.get('MW', '')} g/mol"
                }
            }
        }
    else:
        callback(f"Invalid parameter in molecule {molecule_data}. Cannot proceed with the upload.")

    # Write data to json file
    with open('data.json', 'w') as f:
        f.write(json.dumps(data, indent=4))
    data_json = json.dumps(data)
    callback(f"Sending POST request for molecule: {molecule_data.get('Chemical name', 'Unknown')} to endpoint: {endpoint}")
    callback(f"POST payload: {data_json}")

    response = requests.post(endpoint, data=data_json, headers=headers)
    callback(f"Response status code: {response.status_code}, response text: {response.text}")

    if response.status_code == 200:
        if not os.path.exists(OUTPUT_PATHS['successful_files']):
            os.makedirs(OUTPUT_PATHS['successful_files'])
        shutil.copy(file, OUTPUT_PATHS['successful_files'])
        
        # Log successes
        with open(OUTPUT_PATHS['success_log'], 'a') as success_log:
            success_log.write(f"File: {file}, Molecule: {molecule_data.get('Chemical name', 'Unknown')}, {molecule_data.get('Formula', 'Unknown')}, API Response Status Code: {response.status_code}, response text: {response.text})\n")
        
        callback(f"Successfully uploaded {file} with molecule data {molecule_data}. Copying file to '{OUTPUT_PATHS['successful_files']}' folder.")
        return True
    else:
        if "duplicate" in response.text.lower():
            duplicate_folder = OUTPUT_PATHS['duplicate_files']
            if not os.path.exists(duplicate_folder):
                os.makedirs(duplicate_folder)
            shutil.copy(file, duplicate_folder)
            
            # Log duplicates
            with open(OUTPUT_PATHS['duplicate_log'], 'a') as duplicate_log:
                duplicate_log.write(f"File: {file}, Molecule: {molecule_data.get('Chemical name', 'Unknown')}, {molecule_data.get('Formula', 'Unknown')}, API Response Status Code: {response.status_code}, response text: {response.text})\n")
            
            callback(f"API request failed for {file} with status code {response.status_code}, response: {response.text}. Copying file to '{OUTPUT_PATHS['duplicate_files']}' folder.")
        else:
            if not os.path.exists(OUTPUT_PATHS['failed_files']):
                os.makedirs(OUTPUT_PATHS['failed_files'])
            shutil.copy(file, OUTPUT_PATHS['failed_files'])
                        # Log duplicates
            with open(OUTPUT_PATHS['failed_log'], 'a') as duplicate_log:
                duplicate_log.write(f"File: {file}, Molecule: {molecule_data.get('Chemical name', 'Unknown')}, {molecule_data.get('Formula', 'Unknown')}, API Response Status Code: {response.status_code}, response text: {response.text})\n")
            
            callback(f"API request failed for {file} with status code {response.status_code}, response: {response.text}. Copying file to '{OUTPUT_PATHS['failed_files']}' folder.")


        return False

from rdkit import Chem

def process_sdf(files, params, callback):
    molecules = []
    for sdf_file in files:
        callback(f"Processing file: {sdf_file}")
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

                molecules.append(molecule_data)
            else:
                callback(f"Error: Molecule in file {sdf_file} could not be parsed and will be skipped.")

    callback(f"Total molecules extracted: {len(molecules)}")
    return molecules


if __name__ == '__main__':
    MyApp().run()
