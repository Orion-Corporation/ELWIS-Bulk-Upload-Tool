## Enamine Bulk Import Script 
https://orionfi-my.sharepoint.com/:u:/g/personal/robert_kottelin_orionpharma_com/EfdN7WO_oNpOiYUMHWypGcIBkk0s_8Uia4qQM1HSi9-8Hg?e=DkKIJS


### 1. Process Compounds (`process_sdf`)
   - Reads molecules from an SDF file, removes salts using `FragmentParent`, and sets metadata for each compound (e.g., `Supplier_Product_Code`, `Amount`).
   - Extracts salt properties (ratio, name, molecular weight) if available and stores them in `supplier_code_salt_data`, keyed by `Supplier_Product_Code`.

### 2. Upload Compounds (`send_to_api`)
   - Compresses processed compounds into `summary.zip` and uploads to the API.
   - Polls the upload status, and if there are any errors, downloads a failure report.

### 3. Fetch Batches (`fetch_batches`)
   - Retrieves batch IDs associated with the specified `library_id`.
   - Filters the results to ensure only the latest unique batch ID is captured for each `Supplier_Product_Code` present in `supplier_code_salt_data`.

### 4. Upload Salts (`upload_salts`)
   - For each batch ID retrieved, finds corresponding salt data, looks up the molecular formula from `salts_dict`, and builds a structured payload.
   - Sends this payload to the API to attach the salts to the respective batch.
