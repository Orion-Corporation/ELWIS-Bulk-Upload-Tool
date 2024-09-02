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
    "duplicate_log": "logs/Failed/duplicates.txt",
    "failed_log": "logs/Failed/failed.txt",
    "general_log": "logs/General Log/logs.txt",
    "download_folder": "downloads"
})

API_ENDPOINTS = load_config('config/api_endpoints.json', {
    "Fragment Endpoint": "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/fragments",
    "Compound Endpoint": "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/materials/Compounds/assets"
})

BATCH_FIELDS_CONFIG = load_config('config/batch_fields_config.json', {})