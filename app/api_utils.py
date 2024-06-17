import json
import requests
import os
import shutil
from app.key import api_key

def post_to_api(molecule_data, file, callback, endpoint, param_set):
    print(f"PARAMS: {param_set}")
    if not api_key:
        callback("Error: API key not found. Cannot proceed with the upload.")
        return False

    headers = {
        'accept': 'application/vnd.api+json',
        'x-api-key': api_key,
        'Content-Type': 'application/vnd.api+json',
    }

    # Construct the API request payload based on the parameter set
    if param_set == "Create Compound (Does not work)":
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
                            "id": "00011",  # Can this be assigned by the api? 
                                            # How do we know the batch id? 
                                            # Do we need to create it separately?
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

        fragments = {}
        if 'Salt_name' in molecule_data:
            fragments["salts"] = [
                {
                    "type": "SALT",
                    "name": molecule_data.get('Salt_name', ''),
                    "mf": molecule_data.get('Salt smiles', ''),
                    "mw": {
                        "rawValue": molecule_data.get('MW_salt', ''),
                        "displayValue": f"{molecule_data.get('MW_salt', '')} g/mol"
                    },
                    "id": "SALT:17",
                    "coefficient": 1
                }
            ]
        if 'Solvate_name' in molecule_data:
            fragments["solvates"] = [
                {
                    "type": "SOLVATE",
                    "name": molecule_data.get('Solvate_name', ''),
                    "mf": molecule_data.get('Solvate smiles', ''),
                    "mw": {
                        "rawValue": molecule_data.get('MW_solvate', ''),
                        "displayValue": f"{molecule_data.get('MW_solvate', '')} g/mol"
                    },
                    "id": "Solvate:17",
                    "coefficient": 1
                }
            ]
        
        if fragments:
            data["data"]["relationships"]["batch"]["data"]["attributes"]["fragments"] = fragments

    elif param_set == "Create Salt":
        data = {
            "data": {
                "type": "salt",
                "attributes": {
                    "name": molecule_data.get('Salt_name', ''),
                    "mf": molecule_data.get('Formula', ''),
                    "mw": f"{molecule_data.get('MW_salt', '')} g/mol"
                }
            }
        }
    elif param_set == "Create Solvate":
        data = {
            "data": {
                "type": "solvate",
                "attributes": {
                    "name": molecule_data.get('Solvate_name', ''),
                    "mf": molecule_data.get('Formula', ''),
                    "mw": f"{molecule_data.get('MW', '')} g/mol"
                }
            }
        }
    else:
        callback(f"Invalid parameter in molecule {molecule_data}. Cannot proceed with the upload.")

    # write data to json file
    with open('data.json', 'w') as f:
        f.write(json.dumps(data, indent=4))
    data_json = json.dumps(data)
    callback(f"Sending POST request for molecule: {molecule_data.get('Chemical name', 'Unknown')} to endpoint: {endpoint}")
    callback(f"POST payload: {data_json}")

    response = requests.post(endpoint, data=data_json, headers=headers)
    callback(f"Response status code: {response.status_code}, response text: {response.text}")

    if response.status_code == 200:
        if not os.path.exists(OUTPUT_PATHS['successful_files']):
            os.makedirs(OUTPUT_PATHS['successful_files'])
        shutil.copy(file, OUTPUT_PATHS['successful_files'])
        
        # Log success to success log specified in OUTPUT_PATHS
        with open(OUTPUT_PATHS['success_log'], 'a') as success_log:
            success_log.write(f"File: {file}, Molecule: {molecule_data.get('Chemical name', 'Unknown')}, {molecule_data.get('Formula', 'Unknown')}, API Response Status Code: {response.status_code}, response text: {response.text})\n")
        
        callback(f"Successfully uploaded {file} with molecule data {molecule_data}. Copying file to '{OUTPUT_PATHS['successful_files']}' folder.")
        return True
    else:
        if "duplicate" in response.text.lower():
            duplicate_folder = OUTPUT_PATHS['duplicate_files']
            if not os.path.exists(duplicate_folder):
                os.makedirs(duplicate_folder)
            shutil.copy(file, duplicate_folder)
            
            # Log duplicate to duplicates log specified in OUTPUT_PATHS
            with open(OUTPUT_PATHS['duplicate_log'], 'a') as duplicate_log:
                duplicate_log.write(f"File: {file}, Molecule: {molecule_data.get('Chemical name', 'Unknown')}, {molecule_data.get('Formula', 'Unknown')}, API Response Status Code: {response.status_code}, response text: {response.text})\n")
            
            callback(f"API request failed for {file} with status code {response.status_code}, response: {response.text}. Copying file to '{OUTPUT_PATHS['duplicate_files']}' folder.")
        else:
            if not os.path.exists(OUTPUT_PATHS['failed_files']):
                os.makedirs(OUTPUT_PATHS['failed_files'])
            shutil.copy(file, OUTPUT_PATHS['failed_files'])

        return False
