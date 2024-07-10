def convert_mol_to_cdxml(molecule_data):
    # Template for the CDXML format
    cdxml_template = """<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd" >
<CDXML
 CreationProgram="ChemDraw JS 23.2.0.0"
 BoundingBox="213.06 68.93 368.28 126.88"
 WindowPosition="0 0"
 WindowSize="0 0"
 FractionalWidths="yes"
 InterpretChemically="yes"
 ShowAtomQuery="yes"
 ShowAtomStereo="no"
 ShowAtomEnhancedStereo="yes"
 ShowAtomNumber="no"
 ShowResidueID="no"
 ShowBondQuery="yes"
 ShowBondRxn="yes"
 ShowBondStereo="no"
 ShowTerminalCarbonLabels="no"
 ShowNonTerminalCarbonLabels="no"
 HideImplicitHydrogens="no"
 Magnification="666"
 LabelFont="24"
 LabelSize="10"
 LabelFace="96"
 CaptionFont="24"
 CaptionSize="10"
 HashSpacing="2.50"
 MarginWidth="1.60"
 LineWidth="0.60"
 BoldWidth="2"
 BondLength="14.40"
 BondSpacing="18"
 ChainAngle="120"
 LabelJustification="Auto"
 CaptionJustification="Left"
 AminoAcidTermini="HOH"
 ShowSequenceTermini="yes"
 ShowSequenceBonds="yes"
 ShowSequenceUnlinkedBranches="no"
 ResidueWrapCount="40"
 ResidueBlockCount="10"
 PrintMargins="36 36 36 36"
 MacPrintInfo="0003000000480048000000000318026400000000031802640367052803FC00020000004800480000000003180264000100000064000000010001010100000001270F000100010000000000000000000000000002001901900000000000400000000000000000000100000000000000000000000000000000"
 ChemPropName=""
 ChemPropFormula="Chemical Formula: "
 ChemPropExactMass="Exact Mass: "
 ChemPropMolWt="Molecular Weight: "
 ChemPropMOverZ="m/z: "
 ChemPropAnalysis="Elemental Analysis: "
 ChemPropBoilingPt="Boiling Point: "
 ChemPropMeltingPt="Melting Point: "
 ChemPropCritTemp="Critical Temp: "
 ChemPropCritPres="Critical Pres: "
 ChemPropCritVol="Critical Vol: "
 ChemPropGibbs="Gibbs Energy: "
 ChemPropLogP="Log P: "
 ChemPropMR="MR: "
 ChemPropHenry="Henry&apos;s Law: "
 ChemPropEForm="Heat of Form: "
 ChemProptPSA="tPSA: "
 ChemPropID=""
 ChemPropFragmentLabel=""
 color="0"
 bgcolor="1"
 RxnAutonumberStart="1"
 RxnAutonumberConditions="no"
 RxnAutonumberStyle="Roman"
 RxnAutonumberFormat="(#)"
 MonomerRenderingStyle="graphic"
><colortable>
<color r="1" g="1" b="1"/>
<color r="0" g="0" b="0"/>
<color r="1" g="0" b="0"/>
<color r="1" g="1" b="0"/>
<color r="0" g="1" b="0"/>
<color r="0" g="1" b="1"/>
<color r="0" g="0" b="1"/>
<color r="1" g="0" b="1"/>
</colortable><fonttable>
<font id="24" charset="utf-8" name="Arial"/>
</fonttable><page
 id="51"
 BoundingBox="0 0 581.33 196"
 Width="581.33"
 Height="196"
 HeaderPosition="36"
 FooterPosition="36"
 PageOverlap="0"
 PrintTrimMarks="yes"
 HeightPages="1"
 WidthPages="2"
 DrawingSpace="poster"
><fragment
 id="1"
 BoundingBox="213.06 69.12 368.28 126.88"
 Z="1"
>{atoms}
{bonds}
</fragment></page></CDXML>
"""
    
    atoms = []
    bonds = []
    
    # Add atoms
    atom_template = '<n\n id="{id}"\n p="{x} {y}"\n Z="{z}"\n AS="N"\n AtomID="{atom_id}"\n/>'
    for i, atom in enumerate(molecule_data["atom_block"]):
        atoms.append(atom_template.format(id=i+1, x=atom["x"], y=atom["y"], z=1, atom_id=i+1))
    
    # Add bonds
    bond_template = '<b\n id="{id}"\n Z="{z}"\n B="{begin}"\n E="{end}"\n BS="N"\n Order="{order}"\n/>'
    for i, bond in enumerate(molecule_data["bond_block"]):
        bonds.append(bond_template.format(id=i+1, z=i+1, begin=bond["begin_atom_idx"]+1, end=bond["end_atom_idx"]+1, order=int(bond["bond_type"])))

    # Combine atoms and bonds into the final CDXML content
    cdxml_content = cdxml_template.format(atoms="\n".join(atoms), bonds="\n".join(bonds))
    
    return cdxml_content
