import os
import json
import time
from dotenv import load_dotenv
from rdkit import Chem
from rdkit.Chem.MolStandardize.rdMolStandardize import FragmentParent
from rdkit.Chem.rdMolDescriptors import CalcNumAtomStereoCenters
import requests

# Load API key from .env file
load_dotenv()
api_key = os.getenv('API_KEY')

# API endpoint and headers
base_url = "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0"
bulk_import_url = f"{base_url}/materials/Compounds/bulkImport?rule=USE_MATCHES&importType=zip"
eid = "test"
update_compound_with_salt = f"{base_url}/materials/{eid}?force=true'"

headers = {
    'Content-Type': 'application/vnd.api+json',
    'accept': 'application/vnd.api+json',
    'x-api-key': api_key
}

# Load salt formulas from salts_dictionary.json
with open("config/salts_dictionary.json", "r") as f:
    salts_dict = json.load(f)

# Process molecules from SDF file
def process_sdf(file_path):    
    libraryID = "testlibrary1"
    project = "Unspecified"
    synthesis_date = "2011-11-10T13:37:00Z"
    supplier_name = "Enamine"
    source = "Acquired"
    chemist = "External chemist"
    batch_purpose = "Test compound"
    batch_type = "Discovery"

    file = Chem.SDMolSupplier(file_path)
    sdf_writer = Chem.SDWriter("summary.sdf")
    
    for mol in file:
        if mol is None:
            print("Error reading molecule.")
            continue

        # Apply FragmentParent to remove salts from the molecule
        mol = FragmentParent(mol)
        stereochemistry = "No stereochemistry" if CalcNumAtomStereoCenters(mol) == 0 else "Unresolved stereochemistry"

        # Set other molecule properties
        mol.SetProp("Supplier_Product_Code", mol.GetProp("ID"))
        mol.SetProp("Plate_ID", mol.GetProp("Plate_ID"))
        mol.SetProp("Well", mol.GetProp("Well"))
        mol.SetProp("Batch_Code", mol.GetProp("Barcode"))
        mol.SetProp("Amount", mol.GetProp("Amount_mg") + " mg")
        mol.SetProp("Library_ID", libraryID)
        mol.SetProp("Project", project)
        mol.SetProp("Synthesis_Date", synthesis_date)
        mol.SetProp("Stereochemistry", stereochemistry)
        mol.SetProp("Supplier_Name", supplier_name)
        mol.SetProp("Source", source)
        mol.SetProp("Chemist", chemist)
        mol.SetProp("Batch_Purpose", batch_purpose)
        mol.SetProp("Batch_Type", batch_type)
        
        # Collect salt properties
        if mol.HasProp("Salt_ratio") and mol.HasProp("Salt_Name") and mol.HasProp("MW_salt"):
            salt_ratio = float(mol.GetProp("Salt_ratio"))
            salt_name = mol.GetProp("Salt_Name")
            mw_salt = mol.GetProp("MW_salt")
            # print(f"Salt ratio: {salt_ratio}")
            # print(f"Salt name: {salt_name}")
            # print(f"MW of salt: {mw_salt}")
            
            # Save supplier code : salt data in hashmap "supplier_code_salt_data"
            supplier_code = mol.GetProp("Supplier_Product_Code")
            supplier_code_salt_data[supplier_code] = {
                "Salt_ratio": salt_ratio,
                "Salt_name": salt_name,
                "MW_salt": mw_salt
            }
        else: print("No salt properties found.")
        
        # print hashmap for debugging
        # print(supplier_code_salt_data)
        
        sdf_writer.write(mol)

    sdf_writer.close()
    os.system("rm -f summary.zip")
    os.system("zip summary.zip summary.sdf")
    os.remove("summary.sdf")
    return "summary.zip", libraryID

def send_to_api(zip_file_name):
    with open(zip_file_name, 'rb') as zip_file_contents:
        payload = zip_file_contents.read()
    response = requests.post(bulk_import_url, headers=headers, data=payload)
    if response.status_code in [200, 201]:
        print("Data successfully sent to the API.")
        job_id = response.json().get("data", {}).get("id")
        if job_id:
            check_job_status(job_id)
    else:
        print(f"Failed to send data. Status code: {response.status_code}")
        print(response.text)
        
def check_job_status(job_id):
    # Polls the job status until completion and saves failure report if available.
    status_url = f"{base_url}/materials/bulkImport/jobs/{job_id}"
    failures_url = f"{status_url}/failures"

    while True:
        response = requests.get(status_url, headers=headers)
        job_data = response.json()

        # Get job status from response
        status = job_data.get("data", {}).get("attributes", {}).get("status", "")
        print(f"Current job status: {status}")
        # print(f"Response: {job_data}")

        if status == "COMPLETED":
            # print("Import job COMPLETED successfully.")
            print(f"Response: {job_data}")
            report = job_data.get("data", {}).get("attributes", {}).get("report", {})
            duplicated = report.get("duplicated", 0)
            failed = report.get("failed", 0)

            # Only download the failure report if there are any duplicated or failed records
            if duplicated > 0 or failed > 0:
                print("Detected duplicated or failed items. Downloading failure report...")
                download_failure_report(failures_url)
            else:
                print("No duplicated or failed items. No failure report to download.")
            break

        elif status == "FAILED":
            print("Import job failed. No failure report seems to be available for failed jobs.")
            break

        # Wait before polling again
        time.sleep(5)

# Fetch all compounds and retrieve batch EIDs for those with salts
def fetch_batches(library_id):
    query_url = f"{base_url}/entities/search"
    query_payload = {
        "query": {
            "$and": [
                {
                    "$match": {
                        "field": "materials.Library ID",
                        "in": "tags",
                        "as": "text",
                        "value": library_id,
                        "mode": "keyword"
                    }
                }
            ]
        },
        "options": {
            "offset": 0,
            "limit": 1000,
            "sort": {
                "createdAt": "desc"
            }
        },
        "meta": {
            "reason": "Advanced Search"
        }
    }
    
    response = requests.post(query_url, headers=headers, data=json.dumps(query_payload))
    if response.status_code == 200:
        batches = {}
        data = response.json().get("data", [])
        
        # Iterate over batches, ensuring we only get the latest unique batches per supplier code
        for item in data:
            attributes = item.get("attributes", {})
            supplier_code = attributes.get("tags", {}).get("materials.Supplier Product Code")
            batch_eid = attributes.get("eid")
            
            # Only add batch_eid if it hasn't been processed before for this supplier code
            if supplier_code and batch_eid and supplier_code in supplier_code_salt_data:
                # Check if this supplier_code is already mapped; skip if already processed
                if supplier_code not in batches:
                    batches[supplier_code] = batch_eid
        return batches
    else:
        print("Failed to fetch batches.")
        print(response.text)
        return {}

# Upload salts to the correct batch based on supplier product code
def upload_salts(batches):
    salt_id_counter = 1  # Unique counter for salt IDs

    for supplier_code, batch_eid in batches.items():
        salt_data = supplier_code_salt_data.get(supplier_code)
        if salt_data:
            # Lookup molecular formula in salts_dictionary.json
            mf_formula = salts_dict.get(salt_data["Salt_name"].lower(), "")
            if not mf_formula:
                print(f"No molecular formula found for salt '{salt_data['Salt_name']}'.")

            salt_payload = {
                "data": {
                    "attributes": {
                        "fragments": {
                            "salts": [
                                {
                                    "id": f"SALT:{salt_id_counter}",
                                    "coefficient": int(salt_data["Salt_ratio"]),
                                    "mf": mf_formula,  # Use the mapped molecular formula
                                    "mw": {
                                        "displayValue": f"{salt_data['MW_salt']} g/mol",
                                        "rawValue": str(salt_data["MW_salt"])
                                    },
                                    "name": salt_data["Salt_name"],
                                    "type": "SALT"
                                }
                            ],
                            "solvates": []
                        }
                    }
                }
            }
            patch_url = f"{base_url}/materials/{batch_eid}?force=true"
            response = requests.patch(patch_url, headers=headers, data=json.dumps(salt_payload))
            if response.status_code in [200, 201]:
                print(f"Salt data for {supplier_code} uploaded successfully to batch {batch_eid}.")
            else:
                print(f"Failed to upload salt data for {supplier_code}. Status code: {response.status_code}")
                print(response.text)
            
            # Increment salt ID counter
            salt_id_counter += 1

def download_failure_report(failures_url, retries=3, delay=10):
    """Attempts to download the failure report, with retries if not immediately available."""
    for attempt in range(retries):
        print(f"Attempt {attempt + 1} to download failure report from: {failures_url}")
        response = requests.get(failures_url, headers=headers)

        # Check response status and log details
        print(f"Failure report download response code: {response.status_code}")
        print(f"Failure report download response text: {response.text}")

        if response.status_code == 200 or response.status_code == 201:
            with open("bulk_import_enamine/failure_report.json", "w") as f:
                f.write(response.text)
            print("Failure report saved to failure_report.json.")
            return
        elif response.status_code == 404:
            print("Failure report not found. Retrying...") if attempt < retries - 1 else print("Failure report is not available.")
            time.sleep(delay)  # Wait before retrying
        else:
            print("Failed to download failure report for an unknown reason.")
            break
        
# Main function
# Fetch all compounds and retrieve batch EIDs for those with salts
def fetch_batches(library_id):
    query_url = f"{base_url}/entities/search"
    query_payload = {
        "query": {
            "$and": [
                {
                    "$match": {
                        "field": "materials.Library ID",
                        "in": "tags",
                        "as": "text",
                        "value": library_id,
                        "mode": "keyword"
                    }
                }
            ]
        },
        "options": {
            "offset": 0,
            "limit": 1000,
            "sort": {
                "createdAt": "desc"
            }
        },
        "meta": {
            "reason": "Advanced Search"
        }
    }
    
    response = requests.post(query_url, headers=headers, data=json.dumps(query_payload))
    if response.status_code == 200:
        batches = {}
        data = response.json().get("data", [])
        
        # Iterate over batches, ensuring we only get the latest unique batches per supplier code
        for item in data:
            attributes = item.get("attributes", {})
            supplier_code = attributes.get("tags", {}).get("materials.Supplier Product Code")
            batch_eid = attributes.get("eid")
            
            # Only add batch_eid if it hasn't been processed before for this supplier code
            if supplier_code and batch_eid and supplier_code in supplier_code_salt_data:
                # Check if this supplier_code is already mapped; skip if already processed
                if supplier_code not in batches:
                    batches[supplier_code] = batch_eid
        return batches
    else:
        print("Failed to fetch batches.")
        print(response.text)
        return {}

# Main function
def main():
    # Clear any stateful data before running the script
    global supplier_code_salt_data
    supplier_code_salt_data = {}

    # Process and upload compounds
    sdf_file_path = 'examples/compound_test.sdf'
    zip_file_name, library_id = process_sdf(sdf_file_path)  # Get library ID from process_sdf
    send_to_api(zip_file_name)
    
    # Fetch batches and upload salts to matched compounds
    batches = fetch_batches(library_id)  # Pass library ID dynamically
    upload_salts(batches)

if __name__ == "__main__":
    main()