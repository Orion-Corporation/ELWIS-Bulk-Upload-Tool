from openbabel import openbabel

def convert_sdf_to_cdxml(input_sdf, output_cdxml):
    obConversion = openbabel.OBConversion()
    if not obConversion.SetInAndOutFormats("sdf", "cdxml"):
        print("Error: Could not set input/output formats to sdf/cdxml")
        return

    mol = openbabel.OBMol()
    builder = openbabel.OBBuilder()

    # Read the SDF file
    if not obConversion.ReadFile(mol, input_sdf):
        print(f"Error: Could not read the input file '{input_sdf}'")
        return

    # Check if the molecule has 3D coordinates
    if not mol.Has3D():
        print("No 3D coordinates found. Generating 3D coordinates...")
        builder.Build(mol)
        mol.AddHydrogens()

    # Write the output CDXML file
    if not obConversion.WriteFile(mol, output_cdxml):
        print(f"Error: Could not write the output file '{output_cdxml}'")
        return

    print(f"Successfully converted '{input_sdf}' to '{output_cdxml}'")

# Replace 'input.sdf' and 'output.cdxml' with your actual file names
convert_sdf_to_cdxml('/home/robert/ERAT/ENAMINE/45841679_FL10231460_Orion_2.SDF', 'output.cdxml')
