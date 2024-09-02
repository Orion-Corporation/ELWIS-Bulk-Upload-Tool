import datetime
import json
from config import OUTPUT_PATHS

# logger functions
def handle_success(file, data, orm_code, OUTPUT_PATHS, callback, response, molecule_data):
    molecular_formula = data['data']['attributes']['synonyms'][1]
    supplier_code = molecule_data.get('Supplier code', '')  # Extract the supplier code
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')  # Get the current timestamp
    log_to_general_log(f"Successfully uploaded: - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From file: {file} - API response code: {response.status_code} - API response: {response.text}")
    
    # Log the success information with timestamp
    with open(OUTPUT_PATHS['success_log'], 'a') as success_log:
        success_log.write(f"Timestamp: {timestamp} - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula}  - File: {file} - API response code: {response.status_code} - API Response: {response.text}\n")
    
    callback(f"Successfully uploaded: - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From file: {file} - API response code: {response.status_code}")

def handle_failure(file, data, orm_code, OUTPUT_PATHS, callback, response, molecule_data):
    molecular_formula = data['data']['attributes']['synonyms'][1]
    supplier_code = molecule_data.get('Supplier code', '')  # Extract the supplier code
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')  # Get the current timestamp
    log_to_general_log(f"Failed to upload: - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From File: {file} - API Response Code: {response.status_code} - Response: {response.text}")
    
    # Log the failure information with timestamp
    with open(OUTPUT_PATHS['failed_log'], 'a') as log:
        log.write(f"Timestamp: {timestamp} - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From file: {file} - API Response Status Code: {response.status_code} - API response text: {response.text}\n")

    try:
        response_json = response.json()
        with open(f"{OUTPUT_PATHS['failed_log']}_response.json", 'w') as f:
            json.dump(response_json, f, indent=4)
    except ValueError:
        with open(f"{OUTPUT_PATHS['failed_log']}_response.json", 'w') as f:
            f.write(response.text)
    
    callback(f"Failed to upload: - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From file: {file} - API Response Status Code: {response.status_code}\n")


def log_duplicate(file, molecule_data, callback, orm_code):
    molecular_formula = molecule_data.get('MolecularFormula', '')
    supplier_code = molecule_data.get('Supplier code', '')  # Extract the supplier code
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')  # Get the current timestamp
    log_to_general_log(f"Detected duplicate molecule: - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From File: {file}")
    
    # Log the duplicate information with timestamp
    with open(OUTPUT_PATHS['duplicate_log'], 'a') as duplicate_log:
        duplicate_log.write(f"Timestamp: {timestamp} - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From File: {file}\n")
    
    callback(f"Duplicate: - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From File: {file}")

def log_to_general_log(message):
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    log_message = f"{timestamp} - {message}"
    with open(OUTPUT_PATHS['general_log'], 'a') as f:
        f.write(log_message + '\n')
        f.flush()

def log_failed_upload(file, molecule_data):
    molecular_formula = molecule_data.get('MolecularFormula', '')
    supplier_code = molecule_data.get('Supplier code', '')  # Extract the supplier code
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')  # Get the current timestamp
    log_message = f"Failed to upload molecule: Timestamp: {timestamp} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From File: {file}"
    
    # Log the failed upload with timestamp
    with open(OUTPUT_PATHS['failed_log'], 'a') as failed_log:
        failed_log.write(log_message + '\n')
    
    log_to_general_log(log_message)
