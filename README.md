# ELWIS-Bulk-Upload-Tool
Tool to bulk upload compounds and associated salts from SDF files to Revvitys Signals Inventory

API docs: https://orionsandbox.signalsresearch.revvitycloud.eu/docs/extapi/swagger/index.html
Visio: https://orionfi-my.sharepoint.com/:u:/g/personal/robert_kottelin_orionpharma_com/EfdN7WO_oNpOiYUMHWypGcIBkk0s_8Uia4qQM1HSi9-8Hg?e=DkKIJS

## Todo:
- [X] Add logic for salt duplicate detection based on SMILES instead of salt names - Not possible since API "get_salts" endpoint only return name, mf and mw. Added duplicate detection logic based on salt molecular formula instead.
- [X] Update logs to include: timestamp, ORM code, molecule, file, and API response (more important variables first)
    - [X] Add SDF supplier code (e.g., > <ID> Z2754556176, Query Mcule ID, and MOLPORTID) to extraction of properties and to logs
- [X] Write local logs as Excel files
- [X] Retain general logs when deleting logs after upload
- [X] Save logs to an ELWIS notebook via API
    - [X] .txt
    - [X] .xlsx
- [X] Add config json file for SDF property names - refactor current hard coding, possible to add more fields without updating code
- [ ] Add validation schema step for SDF file property mappings and report any errors in properties
- [ ] Implement rollback function based on database IDs recorded in logs
- [x] Support for supplier SDF-files
    - [x] ENAMINE
    - [X] MOLPORT
        - [ ] Problems with reading all molecules in file? SDF files sometimes have odd delimiters?
    - [X] MCULE
        Query Mcule ID -> batch field
        Ordered Mcule ID_single container ID
        enamine ID = Mcule Product ID
    - [ ] Salts & Solvates
    - [ ] Refactor process_sdf function to vendor specific?
- [ ] User input for batch properties
    - [X] Library ID input from user (user interactive prompt, copy paste)
    - [X] Use Materials table for projects
    - [X] Update list from API "/materials/libraries" automatically
        Map input data to pre-defined limited lists expected by the API
    - [ ] UI input for the rest
- [X] Upload compound without chemical name -> elwis adds this automatically?
- [X] >  <Supplier name> ChemBridge, Can vary, add vendor dictionary to map different names - based on Materials.json. Should be Easy to change as a user
    - [ ] Fill out / complete manual mapping
    - [X] Automatic normalization and mapping
        - [ ] Still issues with "Inc, Corp, Ltd..."
- [ ] Authentication? mitigated by installation on the user's laptop?
    - [ ] Encrypted build of python executable with hidden api key
- [ ] Optimize:
    - [ ] For one bulk upload only fetch salts once and store them for next compound check - don't fetch per compound speed up
- [ ] Document test cases: 
    - [ ] Upload new compounds
    - [ ] Handle duplicates
    - [ ] Compound contains fragment / compound does not contain fragment / contains multiple fragments
    - [ ] Repeat test cases for ENAMINE, MOLPORT, MCULE

## Tuomo discussion
- [ ] Investigate bulk import function, delete Openbabel'
    - [ ] need to create summary file which is problematic?
- [X] List projects (API call Materials.json from Elwis) and map selected project to batch payload
- [X] Input library ID (user pastes text string (no need to fetch with API)) and map to batch payload
- [ ] Dictionary
    - [X] supplier_synonyms.json created and started, implemented in logic
        - [ ] fill out / complete with problematic mappings (Inc, Corp, Ltd...)
    - [X] Automatic normalization and mapping of supplier names to fetched materials_table
- [ ] Build app to encrypted python executable with hidden api key

## Discussion points with Revvity:
- [X] SDF file ingestion, how to handle varying formats and data?
    Lift out to config file and hardcode based on supplier ID, needs to be hardcoded...
- [X] General API usage
- [X] New bulk import function, would it make sense to use that instead?
    No.. maybe in the future
- [X] More robust Mol --> CDXML conversion
    Tool exists but proprietary
    RDKit is fine
    Openbabel or https://github.com/kienerj/pycdxml only option ...
- [X] API rate limitations? Batch size limitations? Bulk Upload Optimization: Is there a way to optimize bulk uploads for large data sets to reduce processing times on your end?
    keep under 200 api calls per minute, error code 429 check for rate limitation error, if sleep else not
    no need for bulk size limitation
- [X] Does the API support additional custom fields or metadata within payloads (specifically fragments) that may not be documented? POST schema in docs only describe name, mf, mw.
    Materials support custom fields (through admin panel configuration) but not fragments



# Git Cheatsheet
git status
git pull origin <branch-name> (Fetches and merges changes from the remote branch to your current local branch)
git fetch origin (Downloads changes from the remote branch but doesn’t merge them. Use git status to see the changes ready to be merged.)
git add <file-name> (or . to add all)
git commit -m ""
git push origin <branch-name>
git push origin <branch-name> --force
git merge <branch-name> (Merges the specified branch into your current branch.)
git checkout <branch-name> (Switches to a new branch)
git checkout -b <new-branch-name> (Creates and switches to a new branch.)
git checkout -- <file-name> (Discards local changes to a specific file.)
git reset --hard (Resets all local changes to the last commit. Warning: This will erase all changes that haven’t been committed.)
git stash (Temporarily saves all modified tracked files, allowing you to switch branches without committing
git stash apply - Reapplies the last stashed changes)

