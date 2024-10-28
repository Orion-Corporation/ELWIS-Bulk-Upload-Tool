# config.py

import json
import os
from dotenv import load_dotenv

# Load environment variables
load_dotenv()
api_key = os.getenv('API_KEY')

# Load configuration files
def load_config(file_path, default):
    try:
        with open(file_path, 'r') as file:
            return json.load(file)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error loading config from {file_path}: {e}")
        return default

OUTPUT_PATHS = load_config('config/output_paths.json', {
    "success_log": "logs/Successfull/success.txt",
    "success_log_excel": "logs/Successfull/success.xlsx",
    "duplicate_log": "logs/Failed/duplicates.txt",
    "duplicate_log_excel": "logs/Failed/duplicates.xlsx",
    "failed_log": "logs/Failed/failed.txt",
    "failed_log_excel": "logs/Failed/failed.xlsx",
    "general_log": "logs/General Log/logs.txt",
    "general_log_excel": "logs/General Log/logs.xlsx",
    "download_folder": "downloads"
})

API_ENDPOINTS = load_config('config/api_endpoints.json', {
    "Fragment Endpoint": "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/fragments",
    "Compound Endpoint": "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/materials/Compounds/assets"
})

BATCH_FIELDS_CONFIG = load_config('config/batch_fields_config.json', {})

SDF_PROPERTIES_CONFIG = load_config('config/sdf_properties_config.json', {})

def load_viable_suppliers(materials_table_path='config/Materials_table.json'):
    try:
        with open(materials_table_path, 'r') as f:
            data = json.load(f)

            # Traverse through the nested structure to reach "fields"
            for item in data.get("data", []):
                attributes = item.get("attributes", {})
                batches = attributes.get("batches", {})
                fields = batches.get("fields", [])

                # Now search through fields for the target ID
                for field in fields:
                    if field.get("id") == "62fa0b5b19660304d1e5b2de" and field.get("name") == "Supplier Name":
                        return field.get("options", [])

            print("Target ID not found in materials table.")
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error loading materials table: {e}")
    return []

# Load viable suppliers at the start
VIABLE_SUPPLIERS = load_viable_suppliers()

def load_supplier_synonyms(file_path='config/supplier_synonyms.json'):
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        print("Supplier synonyms file not found.")
        return {}

# Load synonyms for supplier names
SUPPLIER_SYNONYMS = load_supplier_synonyms()
