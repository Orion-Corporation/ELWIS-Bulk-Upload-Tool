import datetime
import json
from config import OUTPUT_PATHS

def handle_success(file, data, orm_code, OUTPUT_PATHS, callback):
    molecular_formula = data['data']['attributes']['synonyms'][1]
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')  # Get the current timestamp
    log_to_general_log(f"Successfully uploaded molecule: {molecular_formula} - ORM Code: {orm_code} - File: {file}")
    
    # Log the success information with timestamp
    with open(OUTPUT_PATHS['success_log'], 'a') as success_log:
        success_log.write(f"Timestamp: {timestamp} - File: {file} - Molecule: {molecular_formula} - ORM Code: {orm_code} - API Response Status Code: 200\n")
    
    callback(f"Success: {molecular_formula} - ORM Code: {orm_code} - File: {file} - logged successfully.")

def handle_failure(file, data, response, orm_code, OUTPUT_PATHS, callback):
    molecular_formula = data['data']['attributes']['synonyms'][1]
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')  # Get the current timestamp
    log_to_general_log(f"Failed to upload molecule: {molecular_formula} - ORM Code: {orm_code} - File: {file} - Status Code: {response.status_code} - Response: {response.text}")
    log_file = OUTPUT_PATHS['failed_log']
    
    # Log the failure information with timestamp
    with open(log_file, 'a') as log:
        log.write(f"Timestamp: {timestamp} - File: {file} - Molecule: {molecular_formula} - ORM Code: {orm_code} - API Response Status Code: {response.status_code} - response text: {response.text}\n")
    
    callback(f"Failed: {molecular_formula} - ORM Code: {orm_code} - with status code {response.status_code} - Logged as failure.")

    try:
        response_json = response.json()
        with open(f"{log_file}_response.json", 'w') as f:
            json.dump(response_json, f, indent=4)
    except ValueError:
        with open(f"{log_file}_response.json", 'w') as f:
            f.write(response.text)

def log_duplicate(file, molecule_data, callback, orm_code):
    molecular_formula = molecule_data.get('MolecularFormula', '')
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')  # Get the current timestamp
    log_to_general_log(f"Detected duplicate molecule: {molecular_formula} - ORM Code: {orm_code} - File: {file}")
    
    # Log the duplicate information with timestamp
    with open(OUTPUT_PATHS['duplicate_log'], 'a') as duplicate_log:
        duplicate_log.write(f"Timestamp: {timestamp} - File: {file} - Molecule: {molecular_formula} - ORM Code: {orm_code}\n")
    
    callback(f"Duplicate: {molecular_formula} - ORM Code: {orm_code} - Logged as duplicate.")

def log_to_general_log(message):
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    log_message = f"{timestamp} - {message}"
    with open(OUTPUT_PATHS['general_log'], 'a') as f:
        f.write(log_message + '\n')
        f.flush()