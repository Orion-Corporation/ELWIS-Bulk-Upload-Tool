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
    "Compound Endpoint": "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/materials/Compounds/assets",
    "Bulk Import Endpoint": "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/materials/Compounds/bulkImport?rule=NO_DUPLICATED&importType=zip"
})

BATCH_FIELDS_CONFIG = load_config('config/batch_fields_config.json', {})

SDF_PROPERTIES_CONFIG = load_config('config/sdf_properties_config.json', {})

import requests

# API endpoint
url = 'https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/materials/Compounds/bulkImport?rule=NO_DUPLICATED&importType=zip'

# Path to the zip file you want to upload
file_path = '/home/robekott/ERAT/examples/compound_test.zip'

# Headers (add authorization if required)
headers = {
    'Content-Type': 'application/zip',
    'accept': 'application/zip',
    'x-api-key': api_key
}

# Read the zip file in binary mode and send it as a file in the request
with open(file_path, 'rb') as file:
    files = {'file': ('file.zip', file, 'application/zip')}
    response = requests.post(url, headers=headers, files=files)

# Check response
if response.status_code == 200:
    print("File uploaded successfully")
else:
    print(f"Failed to upload file. Status code: {response.status_code}, Message: {response.text}")
