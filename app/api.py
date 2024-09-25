import requests
import json
from config import API_ENDPOINTS, OUTPUT_PATHS
from logger import log_to_general_log, handle_success, handle_failure, log_duplicate
from sdf_processing import construct_payload, check_uniqueness
import os
from datetime import datetime

# API functions
def post_to_api(molecule_data, fragment_data, file, callback, api_key, OUTPUT_PATHS):
    log_to_general_log(f"Checking uniqueness of molecule: {molecule_data.get('MolecularFormula', '')} in file: {file}")
    # Check uniqueness of the compound
    uniqueness_result = check_uniqueness(molecule_data, api_key)
    
    if uniqueness_result is None:
        log_to_general_log(f"Uniqueness result is None. Failed to verify uniqueness for {file}, aborting upload.")
        callback(f"Uniqueness result is None. Could not verify uniqueness for {file}. Aborting upload.")
        return False
    
    if uniqueness_result.get("data"):  # If the result is not empty, it's a duplicate
        orm_code = uniqueness_result["data"][0]["attributes"].get("name", "Unknown")
        log_duplicate(file, molecule_data, callback, orm_code)
        return True
    
    # Prepare api headers
    headers = {
        'Content-Type': 'application/vnd.api+json',
        'accept': 'application/vnd.api+json',
        'x-api-key': api_key
    }
    
    # Upload fragments (salts)
    salt_id = upload_fragments(fragment_data, callback, headers, api_key)
    
    # Construct the payload
    payload = construct_payload(molecule_data, salt_id, fragment_data)
    
    # Send the request
    success = send_request(payload, file, callback, API_ENDPOINTS['Compound Endpoint'], headers, OUTPUT_PATHS, molecule_data)
    print(f"Payload sent for file {file}: {payload}")
    # write payload to json
    with open('payload.json', 'w') as f:
        json.dump(payload, f, indent=4)
    return success

def send_request(data, file, callback, endpoint, headers, OUTPUT_PATHS, molecule_data):
    data_json = json.dumps(data)
    try:
        response = requests.post(endpoint, data=data_json, headers=headers)
        # write response to json file
        with open('response.json', 'w') as f:
            json.dump(response.json(), f, indent=4)
        
        # Extract orm_code from the response
        try:
            orm_code = response.json().get("data", {}).get("attributes", {}).get("name", "Unknown")
        except:
            callback(f"ORM code not found in request: {str(e)}")
        
        if response.status_code in [200, 201]:
            handle_success(file, data, orm_code, OUTPUT_PATHS, callback, response, molecule_data)
            return True
        handle_failure(file, data, response, orm_code, OUTPUT_PATHS, callback, response, molecule_data)
        return False
    except requests.exceptions.RequestException as e:
        callback(f"Network error during request: {str(e)}")
        return False

def upload_fragment(fragment_data, fragment_type, callback, headers):
    fragment_endpoint = f"{API_ENDPOINTS['Fragment Endpoint']}/{fragment_type}"
    log_to_general_log(f"Attempting to upload fragment to {fragment_endpoint} with fragment data: {fragment_data}")
    response = requests.post(fragment_endpoint, data=json.dumps(fragment_data), headers=headers)
    log_to_general_log(f"Received response for fragment upload: {response.status_code} - {response.text}")
    if response.status_code in [200, 201]:
        fragment_id = response.json()['data']['id']
        return fragment_id
    callback(f"Failed to upload {fragment_type}. Status code: {response.status_code}, response: {response.text}")
    return None

def upload_fragments(fragment_data, callback, headers, api_key):
    salt_details = None
    salt_mf = fragment_data.get('MolecularFormula', '').strip().upper()  # Use MF for detection

    # Ensure salt_mf is a valid string
    if not salt_mf:
        print("No valid salt molecular formula found, skipping fragment upload.")
        return None

    # Debug log for salt molecular formula
    print(f"Debug: Retrieved salt molecular formula from fragment_data: '{salt_mf}'")

    callback(f"Salt detected: Molecular Formula {salt_mf}")
    print(f"Salt detected in fragment data: Molecular Formula {salt_mf}")  # Debug

    salt_details = get_existing_fragment_details(salt_mf, "salts", api_key)
    print(f"Debug: Retrieved salt details from API: '{salt_details}' for salt MF: '{salt_mf}'")

    if salt_details and salt_details.get('id'):
        callback(f"Salt with Molecular Formula '{salt_mf}' already exists with ID: {salt_details['id']}, not uploading duplicate")
    else:
        salt_data = {
            "data": {
                "type": "salt",
                "attributes": {
                    "mf": salt_mf,
                    "mw": f"{fragment_data.get('MW_salt', '')} g/mol"
                }
            }
        }
        salt_id = upload_fragment(salt_data, "salts", callback, headers)
        if salt_id:
            callback(f"Successfully uploaded salt with Molecular Formula: {salt_mf}")
            salt_details = {'id': salt_id, 'mf': salt_mf, 'mw': salt_data['data']['attributes']['mw']}
        else:
            callback(f"Failed to upload salt with Molecular Formula: {salt_mf}")
    
    return salt_details

def get_existing_fragment_details(fragment_mf, fragment_type, api_key):
    headers = {
        'accept': 'application/vnd.api+json',
        'x-api-key': api_key
    }
    params = {
        'filter[mf]': fragment_mf  # Use 'mf' for filtering
    }
    fragment_endpoint = f"{API_ENDPOINTS['Fragment Endpoint']}/{fragment_type}"
    response = requests.get(fragment_endpoint, headers=headers, params=params)

    if response.status_code in [200, 201]:
        results = response.json().get('data', [])
        for fragment_data in results:
            # Ensure the molecular formula matches exactly
            if fragment_data['attributes'].get('mf', '').upper() == fragment_mf.upper():
                print(f"Exact match found: {fragment_data}")
                return {
                    'id': fragment_data['id'],
                    'mf': fragment_data['attributes'].get('mf', ''),
                    'mw': fragment_data['attributes'].get('mw', '')
                }
        print(f"No exact match found for salt molecular formula: {fragment_mf}")
    return None

def upload_txt_logs(api_key):
    """
    Uploads log files to the specified endpoint as plain text.
    """
    headers = {
        'accept': 'application/vnd.api+json',
        'Content-Type': 'text/plain',  # Updated to plain text
        'x-api-key': api_key
    }

    # List of log files to upload
    log_files = [
        OUTPUT_PATHS['success_log'],
        OUTPUT_PATHS['duplicate_log'],
        OUTPUT_PATHS['failed_log'],
        OUTPUT_PATHS['general_log']
    ]

    # Base endpoint URL for uploading logs from the configuration
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

import requests
import os
from config import API_ENDPOINTS, OUTPUT_PATHS
from logger import log_to_general_log

def upload_xlsx_logs(api_key):
    """
    Uploads xlsx log files to the specified endpoint.
    """
    headers = {
        'accept': 'application/vnd.api+json',
        'Content-Type': 'application/octet-stream',  # Correct content type for xlsx files
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