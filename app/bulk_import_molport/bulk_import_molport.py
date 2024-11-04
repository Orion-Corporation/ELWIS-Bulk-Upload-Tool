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

# Load salt formulas from salts_dictionary.json
with open("config/salts_dictionary.json", "r") as f:
    salts_dict = json.load(f)
    
# Load Molport supplier synonyms
with open("config/supplier_synonyms_molport.json", "r") as f:
    supplier_synonyms = json.load(f)

# API endpoint and headers
base_url = "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0"

headers = {
    'Content-Type': 'application/vnd.api+json',
    'accept': 'application/vnd.api+json',
    'x-api-key': api_key
}

# Process molecules from SDF file
def process_sdf(file_path):    
    libraryID = "testlibrary1"
    project = "Unspecified"
    synthesis_date = "2011-11-10T13:37:00Z"
    # supplier_name = "Enamine" # Set dynamically from Molport SDF
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
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True) 
        stereochemistry = "No stereochemistry" if CalcNumAtomStereoCenters(mol) == 0 else "Unresolved stereochemistry"

        # Map Molport Supplier Name using supplier synonyms
        original_supplier_name = mol.GetProp("SUPPLIER NAME")
        supplier_product_code = mol.GetProp("CATALOG NUMBER")
        
        # Check if supplier name exists in synonyms, otherwise skip molecule and report an error
        if original_supplier_name not in supplier_synonyms:
            print(f"Error: Supplier name '{original_supplier_name}' not found in supplier_synonyms_molport.json. Skipping molecule {supplier_product_code}")
            continue
        supplier_name = supplier_synonyms[original_supplier_name]
        
        # Generate Well value from BOX_ROW and BOX_COLUMN
        box_row = mol.GetProp("BOX_ROW")
        box_column = int(mol.GetProp("BOX_COLUMN"))
        if box_column < 10:
            well = box_row + "0" + str(box_column) # Add leading zero for single digit columns
            # print(well)
        else:
            well = box_row + str(box_column)
            # print(well)

    
        # Set other molecule properties
        mol.SetProp("Supplier_Product_Code", mol.GetProp("CATALOG NUMBER"))
        mol.SetProp("Plate_ID", mol.GetProp("BOX_NAME"))
        mol.SetProp("Well", well)
        mol.SetProp("Batch_Code", mol.GetProp("MOLPORTID"))
        mol.SetProp("Amount", mol.GetProp("AMOUNT_IN_VIAL") + " mg")
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
        if mol.HasProp("SALT RATIO") and mol.HasProp("SALT DATA") and mol.HasProp("MOL WEIGHT SALT"):
            print(f"Salt properties found for {mol.GetProp('Supplier_Product_Code')}, collecting salt data...")
            salt_ratio = float(mol.GetProp("SALT RATIO"))
            salt_name = mol.GetProp("SALT DATA")
            mw_salt = mol.GetProp("MOL WEIGHT SALT")
            
            # Save supplier code : salt data in hashmap "supplier_code_salt_data"
            supplier_code = mol.GetProp("Supplier_Product_Code")
            supplier_code_salt_data[supplier_code] = {
                "Salt_ratio": salt_ratio,
                "Salt_name": salt_name,
                "MW_salt": mw_salt
            }
        else: print(f"No salt properties found for {mol.GetProp('Supplier_Product_Code')}.")
        
        sdf_writer.write(mol)

    sdf_writer.close()
    os.system("rm -f summary.zip")
    os.system("zip summary.zip summary.sdf")
    os.remove("summary.sdf")
    return "summary.zip", libraryID

def send_to_api(zip_file_name):
    bulk_import_url = f"{base_url}/materials/Compounds/bulkImport?rule=USE_MATCHES&importType=zip"
    with open(zip_file_name, 'rb') as zip_file_contents:
        payload = zip_file_contents.read()
    response = requests.post(bulk_import_url, headers=headers, data=payload)
    if response.status_code in [200, 201]:
        print("Bulk import job submitted successfully... Checking job status.")
        job_id = response.json().get("data", {}).get("id")
        if job_id:
            check_job_status(job_id)
            # Clean up the job on server
            delete_url = f"{base_url}/materials/bulkImport/jobs/{job_id}"
            response = requests.delete(delete_url, headers=headers)
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

            # Complete with duplicates or failed items
            if duplicated > 0 or failed > 0:
                print("Detected duplicated or failed items. Downloading failure report...")
                download_failure_report(failures_url)
            else: break

        elif status == "FAILED":
            print("Import job failed. Downloading failure report...")
            # download_failure_report(failures_url)
            break

        # Wait before polling again
        time.sleep(5)

# Fetch all compounds and retrieve batch EIDs for those with salts
def fetch_batches(library_id):
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
    
    query_url = f"{base_url}/entities/search"
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
                # Check if this supplier_code is already mapped, skip if already processed
                if supplier_code not in batches:
                    batches[supplier_code] = batch_eid
        return batches
    else:
        print("Failed to fetch batches.")
        print(response.text)
        return {}

# Upload salts to the correct batch based on supplier product code
def upload_salts(batches):
    for supplier_code, batch_eid in batches.items():
        salt_data = supplier_code_salt_data.get(supplier_code)
        if salt_data:
            # Lookup molecular formula in salts_dictionary.json
            mf_formula = salts_dict.get(salt_data["Salt_name"].lower(), "")
            if not mf_formula:
                print(f"Incomplete salt data in salts_dictionary.json for '{salt_data['Salt_name']}'.")

            salt_payload = {
                "data": {
                    "attributes": {
                        "fragments": {
                            "salts": [
                                {
                                    "id": f"SALT:{supplier_code}", # OK?? Will create duplicate salts?
                                    "coefficient": int(salt_data["Salt_ratio"]),
                                    "mf": mf_formula, # Mapped MF based on salt name --> MF in salts_dictionary.json
                                    "mw": {
                                        "displayValue": f"{salt_data['MW_salt']} g/mol",
                                        "rawValue": str(salt_data["MW_salt"])
                                    },
                                    "name": salt_data["Salt_name"],
                                    "type": "SALT"
                                }
                            ],
                            "solvates": [] # Needs to be here in payloads structure
                        }
                    }
                }
            }
            patch_url = f"{base_url}/materials/{batch_eid}?force=true"
            response = requests.patch(patch_url, headers=headers, data=json.dumps(salt_payload))
            if response.status_code in [200, 201]:
                print(f"Salt data for {supplier_code} uploaded successfully to batch {batch_eid}")
            else:
                print(f"Failed to upload salt data for {supplier_code}. Status code: {response.status_code}")
                print(response.text)

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
def main():
    global supplier_code_salt_data
    supplier_code_salt_data = {}

    # Process and upload compounds
    sdf_file_path = 'MOLPORT/molport_test.sdf'
    zip_file_name, library_id = process_sdf(sdf_file_path)  # Get library ID from process_sdf
    send_to_api(zip_file_name)
    
    # Fetch batches for specific library ID (in this case "testlibrary1")
    batches = fetch_batches(library_id)
    
    # Upload salts to the correct ORM-code-batch based on supplier product code
    upload_salts(batches)
    os.system("rm -f summary.zip")

if __name__ == "__main__":
    main()
