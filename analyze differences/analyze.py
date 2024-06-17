import pandas as pd
import numpy as np
from rdkit import Chem

# Function to load SDF file into DataFrame
def sdf_to_df(file_path):
    supplier = Chem.SDMolSupplier(file_path)
    data = []
    for mol in supplier:
        if mol is not None:
            props = mol.GetPropsAsDict()
            props['SMILES'] = Chem.MolToSmiles(mol)
            data.append(props)
    return pd.DataFrame(data)

# Load each SDF file
df_enamine = sdf_to_df('SDF_enamine.SDF')
df_mcule = sdf_to_df('SDF_MCULE.sdf')
df_molport = sdf_to_df('SDF_MOLPORT.sdf')

# Get the columns in each DataFrame
columns_enamine = set(df_enamine.columns)
columns_mcule = set(df_mcule.columns)
columns_molport = set(df_molport.columns)

# Find unique columns in each DataFrame
unique_enamine = columns_enamine - columns_mcule - columns_molport
unique_mcule = columns_mcule - columns_enamine - columns_molport
unique_molport = columns_molport - columns_enamine - columns_mcule

# Find common columns
common_columns = columns_enamine & columns_mcule & columns_molport

# Create a summary DataFrame to show the differences
summary_df = pd.DataFrame({
    'Enamine Unique Columns': list(unique_enamine) + [np.nan] * (max(len(unique_enamine), len(unique_mcule), len(unique_molport)) - len(unique_enamine)),
    'Mcule Unique Columns': list(unique_mcule) + [np.nan] * (max(len(unique_enamine), len(unique_mcule), len(unique_molport)) - len(unique_mcule)),
    'Molport Unique Columns': list(unique_molport) + [np.nan] * (max(len(unique_enamine), len(unique_mcule), len(unique_molport)) - len(unique_molport)),
    'Common Columns': list(common_columns) + [np.nan] * (max(len(unique_enamine), len(unique_mcule), len(unique_molport)) - len(common_columns))
})

# Save summary to a CSV file
summary_df.to_csv('sdf_files_structure_differences.csv', index=False)

# Print the summary DataFrame
print(summary_df)
