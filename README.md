# ERAT: ELWIS Registration Tool

The ELWIS Registration Tool (ERAT) helps users process SDF files and sort compounds based on ELWIS API upload responses (successfull/failed/duplicate).

API docs: https://orionsandbox.signalsresearch.revvitycloud.eu/docs/extapi/swagger/index.html

## Todo:
- [X] Add logic for salt duplicate detection based on SMILES instead of salt names - Not possible since API "get_salts" endpoint only return name, mf and mw. Added duplicate detection logic based on salt molecular formula instead.
- [X] Update logs to include: timestamp, ORM code, molecule, file, and API response (more important variables first)
    - [X] Add SDF supplier code (e.g., > <ID> Z2754556176, Query Mcule ID, and MOLPORTID) to extraction of properties and to logs
- [X] Write local logs as Excel files
- [X] Retain general logs when deleting logs after upload
- [X] Save logs to an ELWIS notebook via API
    - [X] .txt
    - [X] .xlsx
- [ ] Add validation schema step for SDF file property mappings and report any errors in properties
- [ ] Add config json file for SDF properties - refactor current hard coding
- [ ] Implement rollback function based on database IDs recorded in logs
- [ ] Authentication? mitigated by installation on the user's laptop?
- [ ] For one bulk upload only fetch salts once and store them for next compound check - don't fetch per compound to spare API overload
- [x] Support for supplier SDF-files
    - [x] ENAMINE
    - [ ] MOLPORT
    - [X] MCULE
        Query Mcule ID -> batch field
        Ordered Mcule ID_single container ID
    - [ ] Salts & Solvates
- [ ] Library ID input from user 
- [ ] User input for batch properties
    - [ ] Use Materials libraries for lists
    - [ ] Update list from API "/materials/libraries"
        Map input data to pre-defined limited lists expected by the API
- [ ] Remove salt name from Compound name
- [ ] Upload compound without chemical name -> elwis adds this automatically?
- [ ] >  <Supplier name> ChemBridge, Can vary, add vendor dictionary to map different names - based on Materials.json. Should be Easy to change as a user
- [ ] enamine ID = Mcule Product ID, when expanding to Mcule
- [ ] Document test cases: 
    - [ ] Upload new compounds
    - [ ] Handle duplicates
    - [ ] Compound contains fragment / compound does not contain fragment / contains multiple fragments
    - [ ] Repeat test cases for ENAMINE, MOLPORT, MCULE

## Discussion points with Revvity:
- [X] SDF file ingestion, how to handle varying formats and data?
    Lift out to config file and hardcode based on supplier ID
- [X] General API usage
- [X] New bulk import function, would it make sense to use that instead?  
- [X] More robust Mol --> CDXML conversion
    Tool exists but proprietary
    RDKit is fine

- [X] API rate limitations? Batch size limitations? Bulk Upload Optimization: Is there a way to optimize bulk uploads for large data sets to reduce processing times on your end?
    keep under 200 api calls per minute, error code 429 check for rate limitation error, if sleep else not
    no need for bulk size limitation
- [X] Does the API support additional custom fields or metadata within payloads (specifically fragments) that may not be documented? POST schema in docs only describe name, mf, mw.
    Materials support custom fields (through admin panel configuration) but not fragments


