# main.py
from kivy.config import Config
Config.set('kivy', 'window_icon', 'burana.ico')

import kivy
import logging
import datetime
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
from kivy.properties import ListProperty
from kivy.lang import Builder
from kivy.clock import Clock
import json
from openbabel import openbabel
import requests
from rdkit import Chem
from config import API_ENDPOINTS, OUTPUT_PATHS, api_key
from sdf_processing import process_sdf, check_uniqueness, construct_payload
from api import upload_xlsx_logs, get_existing_fragment_details, upload_fragment, send_request
from logger import log_to_general_log, log_duplicate, log_failed_upload

# Set the window size
Window.size = (1200, 600)
kivy.logger.Logger.setLevel(logging.ERROR)

# Load Kivy file
try:
    Builder.load_file('styles.kv')
except Exception as e:
    print(f"Error loading Kv file: {e}")

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

    def build(self):
        self.title = "Structure-Data Format (SDF) File Processor"
        self.root = Builder.load_file('styles.kv')
        self.label = self.root.ids.label
        self.button_select = self.root.ids.button_select
        self.button_upload = self.root.ids.button_upload
        self.button_readme = self.root.ids.button_readme
        self.button_clear_folders = self.root.ids.button_clear_folders
        self.terminal_output = self.root.ids.terminal_output

        self.button_select.bind(on_release=self.show_filechooser)
        self.button_upload.bind(on_release=self.upload_files)
        self.button_readme.bind(on_release=self.show_readme)
        self.button_clear_folders.bind(on_release=self.clear_output_folders)
        
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
            OUTPUT_PATHS['success_log_excel'],
            OUTPUT_PATHS['duplicate_log'],
            OUTPUT_PATHS['duplicate_log_excel'],
            OUTPUT_PATHS['failed_log'],
            OUTPUT_PATHS['failed_log_excel'],
            # OUTPUT_PATHS['general_log'] # General log is preserved, should only be deleted manually
        ]

        for log_file in log_files:
            try:
                if os.path.exists(log_file):
                    os.remove(log_file)
            except Exception as e:
                self.print_terminal(f"Error deleting log file {log_file}: {str(e)}")
        
        self.print_terminal("Output files & logs cleared.")

    def process_files(self, files):
        if not files:
            self.print_terminal("No files selected.")
            return

        # Avoiding double-processing by ensuring this function is called once
        if hasattr(self, 'selected_files') and self.selected_files == files:
            self.print_terminal("Files are already being processed. Skipping reprocessing.")
            return

        # Proceed with processing the files
        self.selected_files = files
        self.label.text = f'Selected {len(files)} files.'
        
        # Clear the existing file list
        file_list = self.root.ids.file_list
        file_list.clear_widgets()

        for file in files:
            file_list.add_widget(Label(text=file, size_hint_y=None, height='30dp', font_size='12sp', color=[0, 0, 0, 1]))
        
        self.button_upload.background_color = [0.04, 0.33, 0.64, 1]
        self.button_upload.disabled = False
        self.print_terminal(f'Selected files: {files}')
        
        # Ensure filechooser state is reset
        if hasattr(self, 'filechooser_popup'):
            self.filechooser_popup.filechooser.selection = []

    def upload_files(self, instance=None):
        self.print_terminal("Starting upload process...")
        self.print_terminal(f"Processing file(s): {self.selected_files}")

        for file in self.selected_files:
            # Step 1: Process one SDF file, outputs list of molecules and fragments in cdxml format
            molecules, fragments = process_sdf([file], self.print_terminal)

            if not molecules:
                self.print_terminal(f"No valid molecules found in file: {file}")
                continue

            for molecule_data, fragment_data in zip(molecules, fragments):
                # Step 2: Check uniqueness of the molecule
                uniqueness_result = check_uniqueness(molecule_data, api_key)
                if uniqueness_result and uniqueness_result.get("data"):
                    orm_code = uniqueness_result["data"][0]["attributes"].get("name", "Unknown")
                    log_duplicate(file, molecule_data, self.print_terminal, orm_code)
                    self.print_terminal(f'Duplicate compound detected: - ORM Code: {orm_code} - Molecular Formula: {molecule_data.get("MolecularFormula", "")} - From File: {file}')
                    continue

                # Step 3: Check if fragment_data is valid before proceeding
                salt_id = None
                salt_mf = fragment_data.get('MolecularFormula', '').strip().upper()
                if salt_mf:
                    # Check uniqueness of the fragment
                    salt_details = get_existing_fragment_details(salt_mf, "salts", api_key)
                    
                    # Step 4: Upload fragment if not unique
                    if not salt_details:
                        salt_id = upload_fragment(fragment_data, "salts", self.print_terminal, api_key)
                    else:
                        salt_id = salt_details['id']
                        self.print_terminal(f"Fragment with Molecular Formula '{salt_mf}' already exists with ID: {salt_id}")

                # Step 5: Construct payload
                payload = construct_payload(molecule_data, salt_id, fragment_data)

                # Step 6: Send the payload
                success = send_request(payload, file, self.print_terminal, API_ENDPOINTS['Compound Endpoint'], api_key, OUTPUT_PATHS, molecule_data)
                if success:
                    continue
                else:
                    log_failed_upload(file, molecule_data)

        self.print_terminal("All files processed. Now uploading logs to Registration - Bulk registration via API logs Journal.")
        upload_xlsx_logs(api_key)
        self.print_terminal("Log upload complete.")

    def print_terminal(self, message):
        log_to_general_log(message)
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_message = f"{timestamp} - {message}"
        self.terminal_output.text += log_message + '\n'

if __name__ == '__main__':
    MyApp().run()
