# ERAT: ELWIS Registration Tool

The ELWIS Registration Tool (ERAT) helps users process SDF files and sort compounds based on ELWIS API upload responses (successfull/failed/duplicate).

## Todo:
- [X] Add logic for salt duplicate detection based on SMILES instead of salt names - Not possible since API get salts only return name, mf and mw. Added duplicate logic based on salt molecular formula instead. 
- [X] Update logs to include: timestamp, ORM code, molecule, file, and API response (more important variables first)
- [X] Add SDF supplier code (e.g., > <ID> Z2754556176, Query Mcule ID, and MOLPORTID) to extraction of properties and to logs
- [X] Write local logs as Excel files
- [X] Retain general logs when deleting logs after upload
- [ ] Save logs to an ELWIS notebook via API
    - [X] .txt
    - [ ] .xlsx
- [ ] Add validation schema step for SDF file property mappings and report any errors in properties
- [ ] Add config json file for SDF properties - refactor current hard coding
- [ ] Implement rollback function based on database IDs recorded in logs
- [ ] Document test cases: 
    - [ ] Upload new compounds
    - [ ] Handle duplicates
    - [ ] Compound contains fragment / compound does not contain fragment / contains multiple fragments
    - [ ] Repeat test cases for ENAMINE, MOLPORT, MCULE
- [ ] Authentication? mitigated by installation on the user's laptop?
- [ ] For one bulk upload only fetch salts once and store them for next compound check - don't fetch per compound to spare API overload


- [x] Support for SDF-files
    - [x] ENAMINE
    - [ ] MOLPORT
    - [X] MCULE
    - [ ] Salts & Solvates

https://orionsandbox.signalsresearch.revvitycloud.eu/docs/extapi/swagger/index.html

