import os
import requests
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Get API key from environment variable
api_key = os.getenv('API_KEY')

# Define the endpoint URL
url = "https://orion.signalsresearch.revvitycloud.eu/api/rest/v1.0/fragments/salts"

# Define headers with the API key
headers = {
    "accept": "application/vnd.api+json",
    "x-api-key": api_key
}

# Make the GET request
response = requests.get(url, headers=headers)

# Check for successful response
if response.status_code == 200:
    # Print the JSON response
    print(response.json())
else:
    # Print error message if request fails
    print(f"Error: {response.status_code}")
    print(response.text)
