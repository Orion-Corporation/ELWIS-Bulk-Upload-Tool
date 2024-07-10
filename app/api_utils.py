import requests
import json

# Load configuration files
def load_config(file_path, default):
    try:
        with open(file_path, 'r') as file:
            return json.load(file)
    except FileNotFoundError:
        return default
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from {file_path}: {e}")
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
            if salt_id:
                callback(f"Successfully uploaded salt: {molecule_data.get('Salt_name', '')}")
            else:
                callback(f"Failed to upload salt: {molecule_data.get('Salt_name', '')}")
    except Exception as e:
        callback(f"No Salt_name found: {str(e)}")
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
            if solvate_id:
                callback(f"Successfully uploaded solvate: {molecule_data.get('Solvate_name', '')}")
            else:
                callback(f"Failed to upload solvate: {molecule_data.get('Solvate_name', '')}")
    except Exception as e:
        callback(f"No Solvate_name found: {str(e)}")

    return salt_id, solvate_id

def upload_fragment(fragment_data, fragment_type, callback, headers):
    fragment_endpoint = f"{API_ENDPOINTS['Fragment Endpoint']}/{fragment_type}"
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
                        "id": "5d6e0287ee35880008c18d6d",
                        "value": molecule_data.get("cdxml", "")
                    },
                    {
                        "id": "62f9fe5b74770f14d1de43a8",
                        "value": "No stereochemistry"
                    }
                ]
            },
            "relationships": {
                "batch": {
                    "data": {
                        "type": "batch",
                        "attributes": {
                            "fields": [
                                {
                                "id": "62fcceeb19660304d1e5beee",
                                "value": "Internal"
                                },
                                {
                                "id": "63469c69ed8a726a31923537",
                                "value": "Unspecified"
                                },
                                {
                                "id": "62fcceeb19660304d1e5bef1",
                                "value": "2011-10-10T14:48:00Z"
                                },
                                {
                                "id": "6384a1270d28381d21deaca7",
                                "value": "TestUser MCChemist"
                                },
                                {
                                "id": "62fa096d19660304d1e5b2db",
                                "value": "Dummy compound"
                                },
                                {
                                "id": "62fa096d19660304d1e5b2da",
                                "value": "Discovery"
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
                "coefficient": 1 # Where does this value come from?
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
                "coefficient": 1 # Where does this value come from?
            }]
        data["data"]["relationships"]["batch"]["data"]["attributes"]["fragments"] = fragments

    return data

def send_request(data, file, callback, endpoint, headers, OUTPUT_PATHS):
    # Write data to json file
    # with open('data.json', 'w') as f:
    #     f.write(json.dumps(data, indent=4))
    data_json = json.dumps(data)
    # callback(f"Sending POST request for molecule: {data['data']['attributes']['synonyms'][0]} to endpoint: {endpoint}")
    # callback(f"POST payload: {data_json}")

    try:
        response = requests.post(endpoint, data=data_json, headers=headers)
        # callback(f"Response status code: {response.status_code}, response text: {response.text}")

        if response.status_code == 200:
            handle_success(file, data, OUTPUT_PATHS, callback)
            return True
        else:
            handle_failure(file, data, response, OUTPUT_PATHS, callback)
            return False
    except requests.exceptions.RequestException as e:
        callback(f"Network error during request: {str(e)}")
        return False

def handle_success(file, data, OUTPUT_PATHS, callback):
    if not os.path.exists(OUTPUT_PATHS['successful_files']):
        os.makedirs(OUTPUT_PATHS['successful_files'])
    shutil.copy(file, OUTPUT_PATHS['successful_files'])
    
    with open(OUTPUT_PATHS['success_log'], 'a') as success_log:
        success_log.write(f"File: {file}, Molecule: {data['data']['attributes']['synonyms'][0]}, API Response Status Code: 200\n")
    
    callback(f"Successfully uploaded compound from {file}. Copying file to '{OUTPUT_PATHS['successful_files']}' folder.")

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
