import json
import requests
from logger import log_to_general_log
from config import API_ENDPOINTS, OUTPUT_PATHS, BATCH_FIELDS_CONFIG

def post_to_api(molecule_data, fragment_data, file, callback, api_key, OUTPUT_PATHS):
    log_to_general_log(f"Checking uniqueness of molecule in file: {file}")
    # Check uniqueness of the compound
    uniqueness_result = check_uniqueness(molecule_data, api_key)
    
    if uniqueness_result is None:
        log_to_general_log(f"Failed to verify uniqueness for {file}, aborting upload.")
        callback(f"Could not verify uniqueness for {file}. Aborting upload.")
        return False
    
    if uniqueness_result.get("data"):  # If the result is not empty, it's a duplicate. Tested.
        orm_code = uniqueness_result["data"][0]["attributes"].get("name", "Unknown")
        log_duplicate(file, molecule_data, callback, orm_code)
        return False
    
    # Prepare headers
    headers = {
        'Content-Type': 'application/vnd.api+json',
        'accept': 'application/vnd.api+json',
        'x-api-key': api_key
    }
    
    # Upload fragments (salts/solvates)
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

# Other functions such as upload_fragments, send_request, check_uniqueness, etc.
