## Enamine Bulk Import Script 
https://orionfi-my.sharepoint.com/:u:/g/personal/robert_kottelin_orionpharma_com/EfdN7WO_oNpOiYUMHWypGcIBkk0s_8Uia4qQM1HSi9-8Hg?e=DkKIJS


# ELWIS Bulk Upload Tool

This tool automates the process of reading chemical compounds from an SDF file, removing salts, setting compound metadata, and uploading the processed data to the ELWIS API. It also associates salt data with batches if available, handling each compound in manageable chunks for efficiency.

## Steps Overview

### 1. Process Compounds (`process_sdf`)

- Reads molecules from an SDF file.
- Removes salts using `FragmentParent` from RDKit, then assigns metadata properties to each molecule, including `Supplier_Product_Code`, `Amount`, `Plate_ID`, etc.
- If available, extracts salt properties (`Salt_ratio`, `Salt_Name`, `MW_salt`) and stores these in a global dictionary, `supplier_code_salt_data`, using `Supplier_Product_Code` as the key.
- Writes the processed molecules to an output file, `summary.sdf`, which is compressed into `summary.zip`.

### 2. Upload Compounds (`send_to_api`)

- Compresses the processed SDF file (`summary.sdf`) into `summary.zip` and sends it to the ELWIS API.
- Polls the upload job status with `check_job_status` until completion or failure.
- If any errors are detected, downloads a failure report for diagnostic purposes.

### 3. Check Job Status (`check_job_status`)

- Continuously checks the status of an ongoing upload job.
- If the job completes but has duplicated or failed items, a failure report is downloaded for further inspection.
- Manages retries for status polling until the job succeeds or fails.

### 4. Fetch Batches (`fetch_batches`)

- Retrieves batches from the ELWIS API based on `library_id`.
- For each batch, extracts `batch_eid` (batch ID) and `name` (ORM name).
- Ensures that only the latest batch ID for each `Supplier_Product_Code` in `supplier_code_salt_data` is retained in `batches`.
- Returns a dictionary where each `supplier_code` maps to a dictionary with `batch_eid` and `name`.

### 5. Fetch Salts (`fetch_salts`)

- Retrieves all available salts from the ELWIS API, storing them in `elwis_salts_list`.
- This list is used to validate and map salts extracted from the SDF file against ELWIS’s records during the salt upload process.

### 6. Upload Salts (`upload_salts`)

- Iterates over each `supplier_code` in `batches` and checks if corresponding salt data exists in `supplier_code_salt_data`.
- For each compound with salt data:
  - Looks up the salt’s formula using `salts_dict` to find the ELWIS salt code and retrieve the salt's molecular data (`mf`, `type`, `name`).
  - Constructs a payload with both SDF and ELWIS salt details, then sends this data to the ELWIS API to attach the salts to the corresponding batch.
- If the ELWIS salt code is missing, it logs an error for the respective compound and skips the upload.

### 7. Download Failure Report (`download_failure_report`)

- Attempts to download a failure report from the ELWIS API if there are duplicated or failed items during the upload.
- Retries up to 3 times, with a delay between each attempt, in case of temporary issues.

## Chunk Processing in `main`

The main function handles processing in manageable chunks to optimize memory usage and API calls:
- The list of molecules is divided into sub-lists (`mol_chunk`) of size `chunk_size`.
- Each chunk goes through:
  - **Processing**: Compounds are processed with `process_sdf`.
  - **Upload**: Compounds are uploaded via `send_to_api`.
  - **Batch Retrieval**: Batch information is retrieved for the uploaded compounds using `fetch_batches`.
  - **Salt Upload**: Salt data for each batch is uploaded with `upload_salts`.
  
Example chunking behavior:
- If `molecules = [m0, m1, m2]` and `chunk_size = 2`:
  - **1st iteration (i=0)**: `mol_chunk = [m0, m1]`
  - **2nd iteration (i=2)**: `mol_chunk = [m2]` (the last smaller chunk is handled automatically by Python’s slicing)

## Usage

1. Place the SDF file for the compounds under `examples/compound_test.sdf`.
2. Ensure `config/salts_dictionary_enamine.json` contains the salt mappings.
3. Set up the `.env` file with the API key.
4. Run the script: python bulk_import_enamine.py
