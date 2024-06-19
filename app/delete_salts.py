import requests
from dotenv import load_dotenv
import os

# Load environment variables from .env file
load_dotenv()

# Get the API key from environment variable
API_KEY = os.getenv("API_KEY")

# Define the list of salts to delete
salts_to_delete = [

]

# Base URL for the API
base_url = "https://orionsandbox.signalsresearch.revvitycloud.eu/api/rest/v1.0/fragments/salts/"

# Function to delete a salt by ID
def delete_salt(salt_id):
    url = f"{base_url}{salt_id}"
    headers = {
        'accept': 'application/vnd.api+json',
        'x-api-key': API_KEY,
        'Content-Type': 'application/vnd.api+json',
    }
    
    response = requests.delete(url, headers=headers)
    return response.status_code, response.text

# Delete the salts and print the results
for salt_name, salt_id in salts_to_delete:
    status_code, response_text = delete_salt(salt_id)
    print(f"Deleting {salt_name} with ID {salt_id}: Status Code: {status_code}, Response: {response_text}")

