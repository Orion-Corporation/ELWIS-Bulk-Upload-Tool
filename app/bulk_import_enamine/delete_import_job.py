import os
import requests
from dotenv import load_dotenv

# Load API key from .env file
load_dotenv()
api_key = os.getenv('API_KEY')

# API base URL and headers
base_url = "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0"
headers = {
    'Content-Type': 'application/vnd.api+json',
    'accept': 'application/vnd.api+json',
    'x-api-key': api_key
}

def get_job_status(job_id):
    """Fetch the status of a specific import job."""
    status_url = f"{base_url}/materials/bulkImport/jobs/{job_id}"
    response = requests.get(status_url, headers=headers)
    
    if response.status_code == 200:
        job_data = response.json()
        status = job_data.get("data", {}).get("attributes", {}).get("status", "")
        print(f"Job ID {job_id} status: {status}")
        return status
    else:
        print(f"Failed to retrieve job status. Status code: {response.status_code}")
        return None

def delete_import_job(job_id):
    """Delete a specific import job if it is completed or failed."""
    # Check job status
    status = get_job_status(job_id)
    
    if status in ["COMPLETED", "FAILED", "COMPLETED_WITH_ERRORS"]:
        # Proceed to delete the job
        delete_url = f"{base_url}/materials/bulkImport/jobs/{job_id}"
        response = requests.delete(delete_url, headers=headers)
        
        if response.status_code == 204:
            print(f"Job ID {job_id} deleted successfully.")
        else:
            print(f"Failed to delete job ID {job_id}. Status code: {response.status_code}")
            print(response.text)
    else:
        print(f"Job ID {job_id} is not finished (status: {status}). Cannot delete unfinished jobs.")

# Main function to delete a job by job ID
def main():
    # Input job ID
    job_id = input("Enter the job ID to delete: ")
    delete_import_job(job_id)

if __name__ == "__main__":
    main()
