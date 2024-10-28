# api.py

import requests
import json
from config import API_ENDPOINTS, OUTPUT_PATHS
from logger import log_to_general_log, handle_success, handle_failure, log_duplicate
from sdf_processing import construct_payload, check_uniqueness
import os
from datetime import datetime

# Function to send the payload to the API
def send_request(data, file, callback, endpoint, api_key, OUTPUT_PATHS, molecule_data):
    headers = {
        'Content-Type': 'application/vnd.api+json',
        'accept': 'application/vnd.api+json',
        'x-api-key': api_key
    }
    data_json = json.dumps(data)
    try:
        response = requests.post(endpoint, data=data_json, headers=headers)
        
        # Write response to JSON for debugging
        with open('response.json', 'w') as f:
            json.dump(response.json(), f, indent=4)
        
        # Extract orm_code from the response
        try:
            orm_code = response.json().get("data", {}).get("attributes", {}).get("name", "Unknown")
        except Exception as e:
            callback(f"ORM code not found in request: {str(e)}")
        
        if response.status_code in [200, 201]:
            handle_success(file, data, orm_code, OUTPUT_PATHS, callback, response, molecule_data)
            return True
        handle_failure(file, data, response, orm_code, OUTPUT_PATHS, callback, response, molecule_data)
        return False
    except requests.exceptions.RequestException as e:
        callback(f"Network error during request: {str(e)}")
        return False
    
# Function to upload fragment (salt)
def upload_fragment(fragment_data, fragment_type, callback, api_key):
    fragment_endpoint = f"{API_ENDPOINTS['Fragment Endpoint']}/{fragment_type}"
    headers = {
        'Content-Type': 'application/vnd.api+json',
        'accept': 'application/vnd.api+json',
        'x-api-key': api_key
    }

    # Ensure fragment data contains valid fields
    fragment_name = fragment_data.get("Salt_name", "").strip()
    fragment_mf = fragment_data.get("MolecularFormula", "").strip()
    fragment_mw = f"{fragment_data.get('MW_salt', '')} g/mol" if fragment_data.get("MW_salt") else ""

    if not fragment_name or not fragment_mf or not fragment_mw:
        callback(f"Fragment data is incomplete: name={fragment_name}, mf={fragment_mf}, mw={fragment_mw}")
        return None

    fragment_payload = {
        "data": {
            "type": fragment_type,
            "attributes": {
                "name": fragment_name,
                "mf": fragment_mf,
                "mw": fragment_mw
            }
        }
    }

    log_to_general_log(f"Attempting to upload fragment to {fragment_endpoint} with fragment data: {fragment_payload}")
    
    response = requests.post(fragment_endpoint, data=json.dumps(fragment_payload), headers=headers)
    log_to_general_log(f"Received response for fragment upload: {response.status_code} - {response.text}")
    
    if response.status_code in [200, 201]:
        fragment_id = response.json()['data']['id']
        return fragment_id
    
    callback(f"Failed to upload {fragment_type}. Status code: {response.status_code}, response: {response.text}")
    return None

# Function to check if a fragment already exists
def get_existing_fragment_details(fragment_mf, fragment_type, api_key):
    headers = {
        'accept': 'application/vnd.api+json',
        'x-api-key': api_key
    }
    params = {
        'filter[mf]': fragment_mf  # Filtering for molecular formula
    }
    fragment_endpoint = f"{API_ENDPOINTS['Fragment Endpoint']}/{fragment_type}"
    response = requests.get(fragment_endpoint, headers=headers, params=params)

    if response.status_code in [200, 201]:
        results = response.json().get('data', [])
        for fragment_data in results:
            # Ensure the molecular formula matches
            if fragment_data['attributes'].get('mf', '').upper() == fragment_mf.upper():
                return {
                    'id': fragment_data['id'],
                    'mf': fragment_data['attributes'].get('mf', ''),
                    'mw': fragment_data['attributes'].get('mw', '')
                }
    return None

def upload_txt_logs(api_key):
    """
    Uploads log files to the specified endpoint as plain text.
    """
    headers = {
        'accept': 'application/vnd.api+json',
        'Content-Type': 'text/plain', # Use plain text for txt files
        'x-api-key': api_key
    }

    # List of log files to upload
    log_files = [
        OUTPUT_PATHS['success_log'],
        OUTPUT_PATHS['duplicate_log'],
        OUTPUT_PATHS['failed_log'],
        OUTPUT_PATHS['general_log']
    ]

    # Base endpoint URL
    base_endpoint = API_ENDPOINTS.get('Log Upload Base Endpoint')

    # Ensure the base endpoint exists
    if not base_endpoint:
        log_to_general_log("Log Upload Base Endpoint is not configured.")
        return False

    for log_file in log_files:
        try:
            # Check if the file exists
            if not os.path.exists(log_file):
                log_to_general_log(f"Log file not found: {log_file}, skipping upload.")
                continue

            # Dynamically create the endpoint URL with the log file name
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            file_name = os.path.basename(log_file)
            file_name_without_extension, file_extension = os.path.splitext(file_name)
            # Append the timestamp to the file name
            new_file_name = f"{file_name_without_extension}_{timestamp}{file_extension}"
            # Use the new file name in the endpoint URL
            endpoint = f"{base_endpoint}{new_file_name}?force=true"

            with open(log_file, 'rb') as f:
                # Read the file content for a plain text upload
                response = requests.post(endpoint, headers=headers, data=f)

                # Log the response status
                if response.status_code in [200, 201]:
                    log_to_general_log(f"Successfully uploaded log file: {log_file}")
                else:
                    log_to_general_log(f"Failed to upload log file: {log_file}. Status code: {response.status_code} - Response: {response.text}")

        except Exception as e:
            log_to_general_log(f"Error uploading log file: {log_file} - {str(e)}")
            continue

    return True

def upload_xlsx_logs(api_key):
    """
    Uploads xlsx log files to the specified endpoint.
    """
    headers = {
        'accept': 'application/vnd.api+json',
        'Content-Type': 'application/octet-stream', # Use binary data for Excel files
        'x-api-key': api_key
    }

    # List of Excel log files to upload
    xlsx_log_files = [
        OUTPUT_PATHS['success_log_excel'],
        OUTPUT_PATHS['duplicate_log_excel'],
        OUTPUT_PATHS['failed_log_excel'],
        OUTPUT_PATHS['general_log_excel']   
    ]

    # Base endpoint URL for uploading logs from the configuration
    base_endpoint = API_ENDPOINTS.get('Log Upload Base Endpoint')

    # Ensure the base endpoint exists
    if not base_endpoint:
        log_to_general_log("Log Upload Base Endpoint is not configured.")
        return False

    for log_file in xlsx_log_files:
        try:
            # Check if the file exists
            if not os.path.exists(log_file):
                log_to_general_log(f"Excel log file not found: {log_file}, skipping upload.")
                continue

            # Dynamically create the endpoint URL with the log file name
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            file_name = os.path.basename(log_file)
            file_name_without_extension, file_extension = os.path.splitext(file_name)
            # Append the timestamp to the file name
            new_file_name = f"{file_name_without_extension}_{timestamp}{file_extension}"
            # Use the new file name in the endpoint URL
            endpoint = f"{base_endpoint}{new_file_name}?force=true"

            with open(log_file, 'rb') as f:
                # Upload the file as binary data
                response = requests.post(endpoint, headers=headers, data=f)

                # Log the response status
                if response.status_code in [200, 201]:
                    log_to_general_log(f"Successfully uploaded Excel log file: {log_file}")
                else:
                    log_to_general_log(f"Failed to upload Excel log file: {log_file}. Status code: {response.status_code} - Response: {response.text}")

        except Exception as e:
            log_to_general_log(f"Error uploading Excel log file: {log_file} - {str(e)}")
            continue

    return True