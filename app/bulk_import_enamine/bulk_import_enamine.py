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
with open("config/salts_dictionary_enamine.json", "r") as f:
    salts_dict = json.load(f)

# API endpoint and headers
base_url = "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0"

headers = {
    'Content-Type': 'application/vnd.api+json',
    'accept': 'application/vnd.api+json',
    'x-api-key': api_key
}

def process_sdf(molecules):
    libraryID = "testlibrary1"
    project = "Unspecified"
    synthesis_date = "2011-11-10T13:37:00Z"
    supplier_name = "Enamine"
    source = "Acquired"
    chemist = "External chemist"
    batch_purpose = "Test compound"
    batch_type = "Discovery"

    sdf_writer = Chem.SDWriter("summary.sdf")
    for mol in molecules:
        mol = FragmentParent(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        stereochemistry = "No stereochemistry" if CalcNumAtomStereoCenters(mol) == 0 else "Unresolved stereochemistry"

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

        if mol.HasProp("Salt_ratio") and mol.HasProp("Salt_Name") and mol.HasProp("MW_salt"):
            salt_ratio = float(mol.GetProp("Salt_ratio"))
            salt_name = mol.GetProp("Salt_Name")
            mw_salt = mol.GetProp("MW_salt")
            supplier_code = mol.GetProp("Supplier_Product_Code")
            supplier_code_salt_data[supplier_code] = {
                "Salt_ratio": salt_ratio,
                "Salt_name": salt_name,
                "MW_salt": mw_salt
            }
        sdf_writer.write(mol)

    sdf_writer.close()
    os.system("zip summary.zip summary.sdf")
    os.remove("summary.sdf")
    return "summary.zip", libraryID

def send_to_api(processed_zip_file):
    bulk_import_url = f"{base_url}/materials/Compounds/bulkImport?rule=USE_MATCHES&importType=zip"
    with open(processed_zip_file, 'rb') as zip_file_contents:
        payload = zip_file_contents.read()
    response = requests.post(bulk_import_url, headers=headers, data=payload)
    if response.status_code in [200, 201]:
        job_id = response.json().get("data", {}).get("id")
        if job_id:
            check_job_status(job_id)
            delete_url = f"{base_url}/materials/bulkImport/jobs/{job_id}"
            requests.delete(delete_url, headers=headers)
    else:
        print(f"Failed to send data. Status code: {response.status_code}")
        print(response.text)
    os.remove(processed_zip_file)

def check_job_status(job_id):
    status_url = f"{base_url}/materials/bulkImport/jobs/{job_id}"
    failures_url = f"{status_url}/failures"
    while True:
        response = requests.get(status_url, headers=headers)
        job_data = response.json()
        status = job_data.get("data", {}).get("attributes", {}).get("status", "")
        print(f"Current job status: {status}")

        if status == "COMPLETED":
            report = job_data.get("data", {}).get("attributes", {}).get("report", {})
            print(f"Job Data: {job_data}")
            if report.get("duplicated", 0) > 0 or report.get("failed", 0) > 0:
                download_failure_report(failures_url)
            break
        elif status == "FAILED":
            download_failure_report(failures_url)
            break
        time.sleep(5)

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
    batches = {}
    
    if response.status_code == 200:
        data = response.json().get("data", [])
        for item in data:
            attributes = item.get("attributes", {})
            supplier_code = attributes.get("tags", {}).get("materials.Supplier Product Code")
            batch_eid = attributes.get("eid")
            orm_name = attributes.get("name")
            if supplier_code and batch_eid and supplier_code in supplier_code_salt_data and supplier_code not in batches:
                batches[supplier_code] = {"batch_eid": batch_eid, "name": orm_name}
    else:
        print("Failed to fetch batches.")
        print(response.text)
    return batches

def fetch_salts():
    salt_fetch_url = f"{base_url}/fragments/salts"
    response = requests.get(salt_fetch_url, headers=headers)
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Failed to fetch salts. Status code: {response.status_code}")
        return None

def upload_salts(batches, elwis_salts_list):
    for supplier_code, batch_info in batches.items(): # supplier_code is the key, batch_info is the value, supplier_code mapped to salt data in fetch_batches()
        batch_eid = batch_info["batch_eid"]
        orm_name = batch_info["name"]
        
        salt_data = supplier_code_salt_data.get(supplier_code) # Get salt data based on supplier product code
        
        if salt_data: # If False, skip salt upload
            salt_code = salts_dict.get(salt_data["Salt_name"].lower(), "") # Get salt name from supplier_code_salt_data and map to salt code in salts dictionary
            
            elwis_salt = next((item for item in elwis_salts_list["data"] if item["attributes"]["name"] == salt_code), None) # Match sdf salt code to ELWIS salt code
            
            if elwis_salt is None:
                print(f"Salt '{salt_code}' not found in ELWIS for '{supplier_code}', name '{salt_data['Salt_name']}' missing in salt dictionary")
                continue
            
            # Construct payload, coefficient from sdf, other data from ELWIS
            elwis_salt_id = elwis_salt["id"]
            elwis_salt_type = elwis_salt["attributes"]["type"]
            elwis_salt_name = elwis_salt["attributes"]["name"]
            elwis_mf = elwis_salt["attributes"]["mf"]
            elwis_mw_salt = elwis_salt["attributes"]["mw"]

            salt_payload = {
                "data": {
                    "attributes": {
                        "fragments": {
                            "salts": [
                                {
                                    "id": f"SALT:{elwis_salt_id}",
                                    "coefficient": int(salt_data["Salt_ratio"]),
                                    "mf": elwis_mf,
                                    "mw": {
                                        "displayValue": elwis_mw_salt,
                                        "rawValue": elwis_mw_salt.split()[0]
                                    },
                                    "name": elwis_salt_name,
                                    "type": elwis_salt_type
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
                # Print success message with supplier code and ORM name
                print(f"Salt data for {supplier_code} ({orm_name}) uploaded successfully to batch {batch_eid}")
            else:
                print(f"Failed to upload salt data for {supplier_code} ({orm_name}). Status code: {response.status_code}")
                print(response.text)

def download_failure_report(failures_url, retries=3, delay=10):
    for attempt in range(retries):
        response = requests.get(failures_url, headers=headers)
        if response.status_code in [200, 201]:
            with open("bulk_import_enamine/failure_report.json", "w") as f:
                f.write(response.text)
            print("Failure report saved to failure_report.json.")
            return
        elif response.status_code == 404:
            if attempt < retries - 1:
                time.sleep(delay)
        else:
            print("Failed to download failure report.")
            break

def main():
    global supplier_code_salt_data
    elwis_salts_list = fetch_salts()

    file_path = 'examples/compound_test.sdf' # One SDF file per run enough? Otherwise loop over files
    sdf_supplier = Chem.SDMolSupplier(file_path)
    molecules = [mol for mol in sdf_supplier if mol is not None]
    chunk_size = 2

    for i in range(0, len(molecules), chunk_size): # from 0 to len(molecules), incrementing by chunk_size
        mol_chunk = molecules[i:i + chunk_size] # slices to sub-lists of size chunk_size, ending at but excluding i + chunk_size
        
        ''' Example:
        If molecules = [m0, m1, m2] and chunk_size = 2:
        1st iteration (i=0): mol_chunk = molecules[0:2] → [m0, m1]
        2nd iteration (i=2): mol_chunk = molecules[2:4] → [m2] (last smaller chunk, handled by python automatically in "molecules[i:i + chunk_size]")
        '''
        
        # Reset supplier_code_salt_data for each chunk
        supplier_code_salt_data = {}
        
        # Process the chunk
        processed_zip_file, library_id = process_sdf(mol_chunk)
        
        # Send the processed chunk to the API
        send_to_api(processed_zip_file)

        # Fetch batches for the processed chunk
        batches = fetch_batches(library_id)
        
        # Upload salts for the fetched batches
        upload_salts(batches, elwis_salts_list)

if __name__ == "__main__":
    main()
