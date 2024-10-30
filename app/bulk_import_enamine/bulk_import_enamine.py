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

headers = {
    'Content-Type': 'application/vnd.api+json',
    'accept': 'application/vnd.api+json',
    'x-api-key': api_key
}

# supplier product code - salt data hashmap for later upload of salts to the right cmpound and batch
supplier_code_salt_data = {}

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
            print(f"Salt ratio: {salt_ratio}")
            print(f"Salt name: {salt_name}")
            print(f"MW of salt: {mw_salt}")
            
            # Save supplier code : salt data in hashmap "supplier_code_salt_data"
            supplier_code = mol.GetProp("Supplier_Product_Code")
            supplier_code_salt_data[supplier_code] = {
                "Salt_ratio": salt_ratio,
                "Salt_name": salt_name,
                "MW_salt": mw_salt
            }
        else: continue
        
        # print hashmap for debugging
        print(supplier_code_salt_data)
        
        sdf_writer.write(mol)

    sdf_writer.close()
    os.system("rm -f summary.zip")
    os.system("zip summary.zip summary.sdf")
    os.remove("summary.sdf")
    
    return "summary.zip"

def send_to_api(zip_file_name):
    with open(zip_file_name, 'rb') as zip_file_contents:
        payload = zip_file_contents.read()
    response = requests.post(bulk_import_url, headers=headers, data=payload)
    response_data = response.json()
    
    if response.status_code in [200, 201]:
        print("Data successfully sent to the API.")
        job_id = response_data.get("data", {}).get("id")
        if job_id:
            print(f"Job ID: {job_id}")
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
    sdf_file_path = 'examples/compound_test.sdf'
    zip_file_name = process_sdf(sdf_file_path)
    # send_to_api(zip_file_name)

if __name__ == "__main__":
    main()
