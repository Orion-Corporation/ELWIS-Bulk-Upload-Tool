# ERAT: ELWIS Registration Tool

The ELWIS Registration Tool (ERAT) is designed to help users process SDF files and sort them based on ELWIS API responses. This tool simplifies the process of uploading SDF files to an API and categorizing the responses.

## Features
- **Select Multiple SDF Files**: Choose multiple SDF files for processing.
- **Upload to API**: Upload selected SDF files to a specified API endpoint (salt, solvate, compound).
- **File Sorting**: Automatically sort input files based on the API response into "Successfully Uploaded" and "Failed" categories. Will also log molecules based on successful/failed API responses.

Todo:
- [X] Add upload functions for salts
- [X] Add upload functions for solvates
- [ ] Add upload functions for molecules
    - What properties does ELWIS want?
- [X] Separate salts from moleculs
- [X] Let user choose output path of successfull/failed files 


https://orionsandbox.signalsresearch.revvitycloud.eu/docs/extapi/swagger/index.html
