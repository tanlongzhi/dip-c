import sys
from pdbx.reader.PdbxReader  import PdbxReader
from pdbx.writer.PdbxWriter  import PdbxWriter
from pdbx.reader.PdbxContainers import *

# open mmCIF file to write
myDataList = []
curContainer = DataContainer("myblock")
aCat = DataCategory("atom_site")
aCat.appendAttribute("group_PDB")
aCat.appendAttribute("type_symbol")
aCat.appendAttribute("id")
aCat.appendAttribute("label_asym_id")
aCat.appendAttribute("label_comp_id")
aCat.appendAttribute("label_seq_id")  
aCat.appendAttribute("label_atom_id")  
aCat.appendAttribute("Cartn_x")  
aCat.appendAttribute("Cartn_y")  
aCat.appendAttribute("Cartn_z")
aCat.appendAttribute("B_iso_or_equiv")

sCat = DataCategory("struct_conn")
sCat.appendAttribute("id")
sCat.appendAttribute("conn_type_id")
sCat.appendAttribute("ptnr1_label_asym_id")
sCat.appendAttribute("ptnr1_label_comp_id")
sCat.appendAttribute("ptnr1_label_seq_id")
sCat.appendAttribute("ptnr1_label_atom_id")
sCat.appendAttribute("ptnr2_label_asym_id")
sCat.appendAttribute("ptnr2_label_comp_id")
sCat.appendAttribute("ptnr2_label_seq_id")
sCat.appendAttribute("ptnr2_label_atom_id")

# write atoms
atom_id = 0
for cen_tel_file_line in open(sys.argv[1], "rb"):
    atom_id += 1
    if cen_tel_file_line.endswith("None\n"):
        continue
    cen_tel_file_line_data = cen_tel_file_line.strip().split("\t")
    chain_id = cen_tel_file_line_data[0] # column 1: name -> chain_id
    b_factor = float(cen_tel_file_line_data[1]) # column 2: color -> b_factor
    x, y, z = map(float, cen_tel_file_line_data[2:5]) # columns 3-5: x, y, z
    aCat.append(("HETATM", ".", atom_id, chain_id, 1, 1, 1, x, y, z, b_factor))

# write output
curContainer.append(sCat)
curContainer.append(aCat)
myDataList.append(curContainer)
pdbxW = PdbxWriter(sys.stdout)
pdbxW.write(myDataList)    