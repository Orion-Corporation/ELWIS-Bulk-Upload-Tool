import xml.etree.ElementTree as ET
import re

def molfile_to_cdxml(molfile_content):
    # Parse the molfile content
    lines = molfile_content.split('\n')
    atom_count = int(lines[3].split()[0])
    bond_count = int(lines[3].split()[1])

    # Create the root CDXML element
    root = ET.Element('CDXML')
    root.set('CreationProgram', 'ChemDraw JS 23.2.0.0')
    # Add other attributes as needed

    # Create the page element
    page = ET.SubElement(root, 'page')
    page.set('id', '481')

    # Create the fragment element
    fragment = ET.SubElement(page, 'fragment')
    fragment.set('id', '1')

    # Process atoms
    for i in range(atom_count):
        line = lines[4 + i]
        x, y, z, symbol = line.split()[:4]
        
        n = ET.SubElement(fragment, 'n')
        n.set('id', str(i + 2))
        n.set('p', f"{float(x) * 40:.2f} {-float(y) * 40:.2f}")
        n.set('Z', str(i + 2))
        n.set('AS', 'N')
        n.set('AtomID', str(i + 1))
        
        if symbol != 'C':
            n.set('Element', symbol)
            n.set('NumHydrogens', '0')
            n.set('NeedsClean', 'yes')
            
            t = ET.SubElement(n, 't')
            t.set('p', f"{float(x) * 40 - 3.61:.2f} {-float(y) * 40 + 3.62:.2f}")
            s = ET.SubElement(t, 's')
            s.set('font', '24')
            s.set('size', '10')
            s.set('color', '0')
            s.set('face', '96')
            s.text = symbol

    # Process bonds
    for i in range(bond_count):
        line = lines[4 + atom_count + i]
        atom1, atom2, bond_type = map(int, line.split()[:3])
        
        b = ET.SubElement(fragment, 'b')
        b.set('id', str(atom_count + i + 2))
        b.set('Z', str(atom_count + i + 2))
        b.set('B', str(atom1))
        b.set('E', str(atom2))
        b.set('BS', 'N')
        
        if bond_type == 2:
            b.set('Order', '2')

    # Convert to string and return
    return ET.tostring(root, encoding='unicode')

# Example usage
molfile_content = """
  Mrv1925 04252409172D          

 23 25  0  0  0  0            999 V2000
    0.7145   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145   -1.2375    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.3819   -1.7224    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1665   -1.4675    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    2.7796   -2.0195    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5642   -1.7646    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7358   -0.9576    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1773   -2.3166    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.9620   -2.0617    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5751   -2.6137    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.4035   -3.4207    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6189   -3.6756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0058   -3.1236    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1270   -2.5070    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.3020   -2.5070    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0470   -1.7224    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7376   -1.4675    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3507   -2.0195    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1353   -1.7646    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3068   -0.9576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0915   -0.7027    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6937   -0.4056    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9091   -0.6605    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6  7  2  0  0  0  0
  6  8  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  1  0  0  0  0
 10 11  1  0  0  0  0
 11 12  1  0  0  0  0
 12 13  1  0  0  0  0
  8 13  1  0  0  0  0
  3 14  2  0  0  0  0
 14 15  1  0  0  0  0
 15 16  2  0  0  0  0
  2 16  1  0  0  0  0
 16 17  1  0  0  0  0
 17 18  1  0  0  0  0
 18 19  2  0  0  0  0
 19 20  1  0  0  0  0
 20 21  1  0  0  0  0
 20 22  2  0  0  0  0
 22 23  1  0  0  0  0
 17 23  2  0  0  0  0
M  END
"""

cdxml_content = molfile_to_cdxml(molfile_content)
print(cdxml_content)
with open('output.cdxml', 'w') as f:
    f.write(cdxml_content)