import xml.etree.ElementTree as ET
import xml.dom.minidom
from rdkit import Chem
from rdkit.Chem import AllChem

def mol_to_cdxml(mol):
    # Create the root CDXML element with all attributes
    cdxml = ET.Element("CDXML")
    cdxml.set("CreationProgram", "ChemDraw JS 23.2.0.0")
    cdxml.set("BoundingBox", "237.08 69.10 392.30 126.86")
    cdxml.set("WindowPosition", "0 0")
    cdxml.set("WindowSize", "0 0")
    cdxml.set("FractionalWidths", "yes")
    cdxml.set("InterpretChemically", "yes")
    cdxml.set("ShowAtomQuery", "yes")
    cdxml.set("ShowAtomStereo", "no")
    cdxml.set("ShowAtomEnhancedStereo", "yes")
    cdxml.set("ShowAtomNumber", "no")
    cdxml.set("ShowResidueID", "no")
    cdxml.set("ShowBondQuery", "yes")
    cdxml.set("ShowBondRxn", "yes")
    cdxml.set("ShowBondStereo", "no")
    cdxml.set("ShowTerminalCarbonLabels", "no")
    cdxml.set("ShowNonTerminalCarbonLabels", "no")
    cdxml.set("HideImplicitHydrogens", "no")
    cdxml.set("Magnification", "666")
    cdxml.set("LabelFont", "24")
    cdxml.set("LabelSize", "10")
    cdxml.set("LabelFace", "96")
    cdxml.set("CaptionFont", "24")
    cdxml.set("CaptionSize", "10")
    cdxml.set("HashSpacing", "2.50")
    cdxml.set("MarginWidth", "1.60")
    cdxml.set("LineWidth", "0.60")
    cdxml.set("BoldWidth", "2")
    cdxml.set("BondLength", "14.40")
    cdxml.set("BondSpacing", "18")
    cdxml.set("ChainAngle", "120")
    cdxml.set("LabelJustification", "Auto")
    cdxml.set("CaptionJustification", "Left")
    cdxml.set("AminoAcidTermini", "HOH")
    cdxml.set("ShowSequenceTermini", "yes")
    cdxml.set("ShowSequenceBonds", "yes")
    cdxml.set("ShowSequenceUnlinkedBranches", "no")
    cdxml.set("ResidueWrapCount", "40")
    cdxml.set("ResidueBlockCount", "10")
    cdxml.set("PrintMargins", "36 36 36 36")
    cdxml.set("MacPrintInfo", "0003000000480048000000000318026400000000031802640367052803FC00020000004800480000000003180264000100000064000000010001010100000001270F000100010000000000000000000000000002001901900000000000400000000000000000000100000000000000000000000000000000")
    cdxml.set("ChemPropName", "")
    cdxml.set("ChemPropFormula", "Chemical Formula: ")
    cdxml.set("ChemPropExactMass", "Exact Mass: ")
    cdxml.set("ChemPropMolWt", "Molecular Weight: ")
    cdxml.set("ChemPropMOverZ", "m/z: ")
    cdxml.set("ChemPropAnalysis", "Elemental Analysis: ")
    cdxml.set("ChemPropBoilingPt", "Boiling Point: ")
    cdxml.set("ChemPropMeltingPt", "Melting Point: ")
    cdxml.set("ChemPropCritTemp", "Critical Temp: ")
    cdxml.set("ChemPropCritPres", "Critical Pres: ")
    cdxml.set("ChemPropCritVol", "Critical Vol: ")
    cdxml.set("ChemPropGibbs", "Gibbs Energy: ")
    cdxml.set("ChemPropLogP", "Log P: ")
    cdxml.set("ChemPropMR", "MR: ")
    cdxml.set("ChemPropHenry", "Henry's Law: ")
    cdxml.set("ChemPropEForm", "Heat of Form: ")
    cdxml.set("ChemProptPSA", "tPSA: ")
    cdxml.set("ChemPropID", "")
    cdxml.set("ChemPropFragmentLabel", "")
    cdxml.set("color", "0")
    cdxml.set("bgcolor", "1")
    cdxml.set("RxnAutonumberStart", "1")
    cdxml.set("RxnAutonumberConditions", "no")
    cdxml.set("RxnAutonumberStyle", "Roman")
    cdxml.set("RxnAutonumberFormat", "(#)")
    cdxml.set("MonomerRenderingStyle", "graphic")

    # Add colortable
    colortable = ET.SubElement(cdxml, "colortable")
    colors = [
        (1, 1, 1), (0, 0, 0), (1, 0, 0), (1, 1, 0),
        (0, 1, 0), (0, 1, 1), (0, 0, 1), (1, 0, 1)
    ]
    for r, g, b in colors:
        color = ET.SubElement(colortable, "color")
        color.set("r", str(r))
        color.set("g", str(g))
        color.set("b", str(b))

    # Add fonttable
    fonttable = ET.SubElement(cdxml, "fonttable")
    font = ET.SubElement(fonttable, "font")
    font.set("id", "24")
    font.set("charset", "utf-8")
    font.set("name", "Arial")

    # Add page
    page = ET.SubElement(cdxml, "page")
    page.set("id", "481")
    page.set("BoundingBox", "0 0 629.33 196")
    page.set("Width", "629.33")
    page.set("Height", "196")
    # ... (add other page attributes)

    # Add fragment
    fragment = ET.SubElement(page, "fragment")
    fragment.set("id", "1")
    fragment.set("BoundingBox", "237.08 69.10 392.30 126.86")
    fragment.set("Z", "1")

    # Generate 2D coordinates if not present
    try:
        conf = mol.GetConformer()
    except ValueError:
        # No conformer exists, so we need to generate 2D coordinates
        mol = Chem.AddHs(mol)
        AllChem.Compute2DCoords(mol)
        mol = Chem.RemoveHs(mol)
    else:
        if conf.Is3D():
            # If it's a 3D conformer, we still want 2D coordinates
            mol = Chem.RemoveHs(mol)
            AllChem.Compute2DCoords(mol)

    # Add atoms
    for i, atom in enumerate(mol.GetAtoms()):
        n = ET.SubElement(fragment, "n")
        n.set("id", str(i+2))  # Start from 2 to match example
        pos = mol.GetConformer().GetAtomPosition(i)
        n.set("p", f"{pos.x:.2f} {pos.y:.2f}")
        n.set("Z", str(i+2))
        n.set("AS", "N")
        n.set("AtomID", str(i+1))

        if atom.GetSymbol() != "C":
            n.set("Element", atom.GetSymbol())
            n.set("NumHydrogens", "0")
            n.set("NeedsClean", "yes")
            
            t = ET.SubElement(n, "t")
            t.set("p", f"{pos.x-3.61:.2f} {pos.y+3.62:.2f}")
            t.set("BoundingBox", f"{pos.x-3.61:.2f} {pos.y-5.62:.2f} {pos.x+3.61:.2f} {pos.y+3.62:.2f}")
            t.set("LabelJustification", "Left")
            
            s = ET.SubElement(t, "s")
            s.set("font", "24")
            s.set("size", "10")
            s.set("color", "0")
            s.set("face", "96")
            s.text = atom.GetSymbol()

    # Add bonds
    for i, bond in enumerate(mol.GetBonds()):
        b = ET.SubElement(fragment, "b")
        b.set("id", str(i+25))  # Start from 25 to match example
        b.set("Z", str(i+25))
        b.set("B", str(bond.GetBeginAtomIdx()+2))
        b.set("E", str(bond.GetEndAtomIdx()+2))
        b.set("BS", "N")
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            b.set("Order", "2")
        elif bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            b.set("Order", "3")

    # Convert to string and prettify
    xml_str = ET.tostring(cdxml, encoding='unicode')
    dom = xml.dom.minidom.parseString(xml_str)
    pretty_xml_str = dom.toprettyxml(indent="  ")

    # Add DOCTYPE and XML declaration
    doctype = '<!DOCTYPE CDXML SYSTEM "https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd" >'
    xml_declaration = '<?xml version="1.0" encoding="UTF-8" ?>'
    final_xml = f"{xml_declaration}\n{doctype}\n{pretty_xml_str}"

    return final_xml

# Example usage (you would replace this with your actual mol object)
mol = Chem.MolFromSmiles("CN1C(SCC(=O)N2CCOCC2)=NN=C1C=3C=CC(F)=CC3")
cdxml_string = mol_to_cdxml(mol)
print(cdxml_string)
# save file
with open('output.cdxml', 'w') as f:
    f.write(cdxml_string)