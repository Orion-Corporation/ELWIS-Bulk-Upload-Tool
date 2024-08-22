import requests
import json
from config import API_ENDPOINTS
from logger import log_to_general_log, handle_success, handle_failure, log_duplicate
from sdf_processing import construct_payload, check_uniqueness

# API functions
def post_to_api(molecule_data, fragment_data, file, callback, api_key, OUTPUT_PATHS):
    log_to_general_log(f"Checking uniqueness of molecule in file: {file}")
    # Check uniqueness of the compound
    uniqueness_result = check_uniqueness(molecule_data, api_key)
    
    if uniqueness_result is None:
        log_to_general_log(f"Failed to verify uniqueness for {file}, aborting upload.")
        callback(f"Could not verify uniqueness for {file}. Aborting upload.")
        return False
    
    if uniqueness_result.get("data"):  # If the result is not empty, it's a duplicate
        orm_code = uniqueness_result["data"][0]["attributes"].get("name", "Unknown")
        log_duplicate(file, molecule_data, callback, orm_code)
        return False
    
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
    success = send_request(payload, file, callback, API_ENDPOINTS['Compound Endpoint'], headers, OUTPUT_PATHS)
    print(f"Payload sent for file {file}: {payload}")
    # write payload to json
    with open('payload.json', 'w') as f:
        json.dump(payload, f, indent=4)
    return success

def upload_fragment(fragment_data, fragment_type, callback, headers):
    fragment_endpoint = f"{API_ENDPOINTS['Fragment Endpoint']}/{fragment_type}"
    log_to_general_log(f"Attempting to upload fragment to {fragment_endpoint} with data: {fragment_data}")
    response = requests.post(fragment_endpoint, data=json.dumps(fragment_data), headers=headers)
    log_to_general_log(f"Received response for fragment upload: {response.status_code} - {response.text}")
    if response.status_code in [200, 201]:
        fragment_id = response.json()['data']['id']
        return fragment_id
    callback(f"Failed to upload {fragment_type}. Status code: {response.status_code}, response: {response.text}")
    return None

def get_existing_fragment_details(fragment_name, fragment_type, api_key):
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
        for fragment_data in results:
            # Ensure the name matches exactly
            if fragment_data['attributes'].get('name', '').lower() == fragment_name.lower():
                print(f"Exact match found: {fragment_data}")
                return {
                    'id': fragment_data['id'],
                    'name': fragment_data['attributes'].get('name', ''),
                    'mf': fragment_data['attributes'].get('mf', ''),
                    'mw': fragment_data['attributes'].get('mw', '')
                }
        print(f"No exact match found for salt name: {fragment_name}")
    return None

def upload_fragments(fragment_data, callback, headers, api_key):
    salt_details = None
    salt_name = fragment_data.get('Salt_name', '') or fragment_data.get('Salt_Name', '')

    # Ensure salt_name is a string
    if isinstance(salt_name, (list, tuple)):
        salt_name = salt_name[0] if salt_name else ''
    else:
        salt_name = str(salt_name)

    # Debug log for salt name
    print(f"Debug: Retrieved salt name from fragment_data: '{salt_name}'")

    if salt_name:
        callback(f"Salt detected: {salt_name}")
        print(f"Salt detected in fragment data: {salt_name}")  # Debug

        salt_details = get_existing_fragment_details(salt_name, "salts", api_key)
        print(f"Debug: Retrieved salt details from API: '{salt_details}' for salt_name: '{salt_name}'")

        if salt_details and salt_details.get('id'):
            callback(f"Salt '{salt_name}' already exists with ID: {salt_details['id']}, not uploading duplicate")
        else:
            salt_data = {
                "data": {
                    "type": "salt",
                    "attributes": {
                        "name": salt_name,
                        "mf": fragment_data.get('MolecularFormula', ''),
                        "mw": f"{fragment_data.get('MW_salt', '')} g/mol"
                    }
                }
            }
            salt_id = upload_fragment(salt_data, "salts", callback, headers)
            if salt_id:
                callback(f"Successfully uploaded salt: {salt_name}")
                salt_details = {'id': salt_id, 'name': salt_name, 'mf': salt_data['data']['attributes']['mf'], 'mw': salt_data['data']['attributes']['mw']}
            else:
                callback(f"Failed to upload salt: {salt_name}")
    else:
        print("No valid salt name found, skipping fragment upload.")
    
    return salt_details

def send_request(data, file, callback, endpoint, headers, OUTPUT_PATHS):
    data_json = json.dumps(data)
    try:
        response = requests.post(endpoint, data=data_json, headers=headers)
        # write response to json file
        with open('response.json', 'w') as f:
            json.dump(response.json(), f, indent=4)
        
        # Extract orm_code from the response
        orm_code = response.json().get("data", {}).get("attributes", {}).get("name", "Unknown")
        
        if response.status_code in [200, 201]:
            handle_success(file, data, orm_code, OUTPUT_PATHS, callback)
            return True
        handle_failure(file, data, response, orm_code, OUTPUT_PATHS, callback)
        return False
    except requests.exceptions.RequestException as e:
        callback(f"Network error during request: {str(e)}")
        return False
    
    