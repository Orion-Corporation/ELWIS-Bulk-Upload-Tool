import datetime
import json
from config import OUTPUT_PATHS
import pandas as pd

# logger functions
def handle_success(file, data, orm_code, OUTPUT_PATHS, callback, response, molecule_data):
    molecular_formula = data['data']['attributes']['synonyms'][1]
    supplier_code = molecule_data.get('Supplier code', '')  # Extract the supplier code
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')  # Get the current timestamp
    log_to_general_log(f"Successfully uploaded: - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From file: {file} - API response code: {response.status_code} - API response: {response.text}\n")
    
    # Log the success information with timestamp
    with open(OUTPUT_PATHS['success_log'], 'a') as success_log:
        success_log.write(f"Timestamp: {timestamp} - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula}  - File: {file} - API response code: {response.status_code} - API Response: {response.text}\n")
    
    success_log_entry = {
        "Timestamp": timestamp,
        "ORM Code": orm_code,
        "Supplier Code": supplier_code,
        "Molecular Formula": molecular_formula,
        "File": file,
        "Reason": "Success",
        "API Response Code": response.status_code,
        "API Response": response.text
    }
    append_to_excel(OUTPUT_PATHS['success_log_excel'], success_log_entry)
    append_to_excel(OUTPUT_PATHS['general_log_excel'], success_log_entry)
    callback(f"Successfully uploaded: - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From file: {file} - API response code: {response.status_code}\n")

def handle_failure(file, data, orm_code, OUTPUT_PATHS, callback, response, molecule_data):
    molecular_formula = data['data']['attributes']['synonyms'][1]
    supplier_code = molecule_data.get('Supplier code', '')  # Extract the supplier code
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')  # Get the current timestamp
    log_to_general_log(f"Failed to upload: - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From File: {file} - API Response Code: {response.status_code} - Response: {response.text}\n")
    
    failure_log_entry = {
        "Timestamp": timestamp,
        "ORM Code": orm_code,
        "Supplier Code": supplier_code,
        "Molecular Formula": molecular_formula,
        "File": file,
        "Reason": "Failed",
        "API Response Status Code": response.status_code,
        "API Response": response.text
    }
    append_to_excel(OUTPUT_PATHS['failed_log_excel'], failure_log_entry)
    append_to_excel(OUTPUT_PATHS['general_log_excel'], failure_log_entry)
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
    log_to_general_log(f"Detected duplicate molecule: - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From File: {file}\n")

    duplicate_log_entry = {
        "Timestamp": timestamp,
        "ORM Code": orm_code,
        "Supplier Code": supplier_code,
        "Molecular Formula": molecular_formula,
        "Reason": "Duplicate",
        "File": file
    }
    append_to_excel(OUTPUT_PATHS['duplicate_log_excel'], duplicate_log_entry)
    append_to_excel(OUTPUT_PATHS['general_log_excel'], duplicate_log_entry)
    # Log the duplicate information with timestamp
    with open(OUTPUT_PATHS['duplicate_log'], 'a') as duplicate_log:
        duplicate_log.write(f"Timestamp: {timestamp} - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From File: {file}\n")
    
    callback(f"Duplicate: - ORM Code: {orm_code} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From File: {file}\n")

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
    log_message = f"Failed to upload molecule: Timestamp: {timestamp} - Supplier Code: {supplier_code} - Molecular Formula: {molecular_formula} - From File: {file}\n"
    
    # Log the failed upload with timestamp
    with open(OUTPUT_PATHS['failed_log'], 'a') as failed_log:
        failed_log.write(log_message + '\n')
    
    failure_log_entry = {
        "Timestamp": timestamp,
        "Supplier Code": supplier_code,
        "Molecular Formula": molecular_formula,
        "Reason": "Failed upload",
        "File": file
    }
    append_to_excel(OUTPUT_PATHS['failed_log_excel'], failure_log_entry)
    append_to_excel(OUTPUT_PATHS['general_log_excel'], failure_log_entry)
    log_to_general_log(log_message)

def append_to_excel(file_path, data_dict):
    # Append data to an Excel file
    try:
        # Read the existing file if it exists
        try:
            df = pd.read_excel(file_path)
        except FileNotFoundError:
            df = pd.DataFrame()

        # Convert data_dict to a DataFrame and append to the existing data
        new_data = pd.DataFrame([data_dict])
        df = pd.concat([df, new_data], ignore_index=True)

        # Save to the Excel file
        df.to_excel(file_path, index=False)
    except Exception as e:
        print(f"Error writing to Excel file {file_path}: {str(e)}")