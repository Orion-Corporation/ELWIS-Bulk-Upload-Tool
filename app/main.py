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
        
        # Reset the file chooser selection state
        if hasattr(self, 'filechooser_popup'):
            self.filechooser_popup.filechooser.selection = []

    def upload_files(self, instance=None):
        if self.upload_in_progress:
            self.print_terminal("Upload already in progress, skipping redundant call")
            return

        self.upload_in_progress = True
        self.button_stop_upload.disabled = False

        self.print_terminal("Upload button pressed")

        def update_progress_bar(dt):
            pass

        self.print_terminal("Uploading files to API...")
        
        # Schedule the upload process to allow frequent checking of the stop condition
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
                self.print_terminal(f'Successfully uploaded file: {file}')
            else:
                self.print_terminal(f'Failed to upload file: {file}')
            
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


import os
import requests
import shutil
import json

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
    except Exception as e:
        callback(f"Error uploading salt: {str(e)}")

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
    except Exception as e:
        callback(f"Error uploading solvate: {str(e)}")

    return salt_id, solvate_id

def upload_fragment(fragment_data, fragment_type, callback, headers):
    fragment_endpoint = f"https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/fragments/{fragment_type}s"
    response = requests.post(fragment_endpoint, data=json.dumps(fragment_data), headers=headers)
    if response.status_code == 200:
        fragment_id = response.json()['data']['id']
        return fragment_id
    else:
        callback(f"Failed to upload {fragment_type}. Status code: {response.status_code}, response: {response.text}")
        return None

def construct_payload(molecule_data, salt_id, solvate_id):
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
    # Write data to json file
    with open('data.json', 'w') as f:
        f.write(json.dumps(data, indent=4))
    data_json = json.dumps(data)
    callback(f"Sending POST request for molecule: {data['data']['attributes']['synonyms'][0]} to endpoint: {endpoint}")
    callback(f"POST payload: {data_json}")

    response = requests.post(endpoint, data=data_json, headers=headers)
    callback(f"Response status code: {response.status_code}, response text: {response.text}")

    if response.status_code == 200:
        handle_success(file, data, OUTPUT_PATHS, callback)
        return True
    else:
        handle_failure(file, data, response, OUTPUT_PATHS, callback)
        return False

def handle_success(file, data, OUTPUT_PATHS, callback):
    if not os.path.exists(OUTPUT_PATHS['successful_files']):
        os.makedirs(OUTPUT_PATHS['successful_files'])
    shutil.copy(file, OUTPUT_PATHS['successful_files'])
    
    with open(OUTPUT_PATHS['success_log'], 'a') as success_log:
        success_log.write(f"File: {file}, Molecule: {data['data']['attributes']['synonyms'][0]}, API Response Status Code: 200\n")
    
    callback(f"Successfully uploaded {file}. Copying file to '{OUTPUT_PATHS['successful_files']}' folder.")

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

                    molecules.append(molecule_data)
                else:
                    callback(f"Error: Molecule in file {sdf_file} could not be parsed and will be skipped.")
        except Exception as e:
            callback(f"Error processing file {sdf_file}: {str(e)}")

    callback(f"Total molecules extracted: {len(molecules)}")
    return molecules


if __name__ == '__main__':
    MyApp().run()
