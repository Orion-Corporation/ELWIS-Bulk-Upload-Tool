import json

# Sample JSON data
json_data = {
}

# Initialize an empty dictionary to store the results
processed_data = {}

# Loop through each item in the "data" list
for item in json_data["data"]:
    attributes = item["attributes"]
    name = attributes["name"]
    formula = attributes["mf"]

    # Generate a generic name by processing the original name (e.g., lowercase, simple words)
    generic_name = name.lower()
    
    # Add to dictionary with generic name as key and original name as value
    processed_data[generic_name] = name

# Output the processed data
print(processed_data)
