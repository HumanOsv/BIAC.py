


import re
import numpy as np
import os
import sys
import itertools
import pandas as pd

##########################################################
#
#       FUNCTIONS
#
##########################################################

def coordsXYZ(fileName):
    filexyz = str(fileName) + "_opt.xyz"
    tmp_xyz = []
    axis_x  = []
    axis_y  = []
    axis_z  = []
    FI = open(filexyz, 'r')
    for line in FI:
        line = line.replace('\n', '').replace('\r', '')
        tmp_xyz.append(line)
    FI.close()
    del tmp_xyz[:2]
    length = len(tmp_xyz)
    for i in range(0,length):
        parameter_mol = tmp_xyz[i].split()
        coord_x   = float(parameter_mol[1])
        coord_y   = float(parameter_mol[2])
        coord_z   = float(parameter_mol[3])
        axis_x.append(coord_x)
        axis_y.append(coord_y)
        axis_z.append(coord_z)
    return axis_x, axis_y, axis_z

def Gaussian_coords(g09file,fileName,NAtoms):
    ## Rangsiman Ketkaew, MSc student in physical chemistry
    ## Computational Chemistry Research Unit, Thammasat University, Thailand
    ## https://github.com/rangsimanketkaew
    start = 0
    end = 0
    ## nutt: updated 03 March
    newfile = str(fileName) + "_opt.xyz"
    openold = open(g09file,"r")
    opennew = open(newfile,"w")
    rline = openold.readlines()
    for i in range (len(rline)):
        if "Standard orientation:" in rline[i]:
            start = i
    for m in range (start + 5, len(rline)):
        if "---" in rline[m]:
            end = m
            break
    ## Convert to Cartesian coordinates format
    ## convert atomic number to atomic symbol
    opennew.write(str(NAtoms))
    opennew.write("\n")
    opennew.write("opt molecule")
    opennew.write("\n")
    for line in rline[start+5 : end] :
        words = line.split()
        word1 = int(words[1])
        word3 = str(words[3])
        ## Periodic table supported
        ## Your molecule should comprise only atoms between 1-90
        if   word1 ==   1 : word1 = "H"
        elif word1 ==   2 : word1 = "He"
        elif word1 ==   3 : word1 = "Li"
        elif word1 ==   4 : word1 = "Be"
        elif word1 ==   5 : word1 = "B"
        elif word1 ==   6 : word1 = "C"
        elif word1 ==   7 : word1 = "N"
        elif word1 ==   8 : word1 = "O"
        elif word1 ==   9 : word1 = "F"
        elif word1 ==  10 : word1 = "Ne"
        elif word1 ==  11 : word1 = "Na"
        elif word1 ==  12 : word1 = "Mg"
        elif word1 ==  13 : word1 = "Al"
        elif word1 ==  14 : word1 = "Si"
        elif word1 ==  15 : word1 = "P"
        elif word1 ==  16 : word1 = "S"
        elif word1 ==  17 : word1 = "Cl"
        elif word1 ==  18 : word1 = "Ar"
        elif word1 ==  19 : word1 = "K"
        elif word1 ==  20 : word1 = "Ca"
        elif word1 ==  21 : word1 = "Sc"
        elif word1 ==  22 : word1 = "Ti"
        elif word1 ==  23 : word1 = "V"
        elif word1 ==  24 : word1 = "Cr"
        elif word1 ==  25 : word1 = "Mn"
        elif word1 ==  26 : word1 = "Fe"
        elif word1 ==  27 : word1 = "Co"
        elif word1 ==  28 : word1 = "Ni"
        elif word1 ==  29 : word1 = "Cu"
        elif word1 ==  30 : word1 = "Zn"
        elif word1 ==  31 : word1 = "Ga"
        elif word1 ==  32 : word1 = "Ge"
        elif word1 ==  33 : word1 = "As"
        elif word1 ==  34 : word1 = "Se"
        elif word1 ==  35 : word1 = "Br"
        elif word1 ==  36 : word1 = "Kr"
        elif word1 ==  37 : word1 = "Rb"
        elif word1 ==  38 : word1 = "Sr"
        elif word1 ==  39 : word1 = "Y"
        elif word1 ==  40 : word1 = "Zr"
        elif word1 ==  41 : word1 = "Nb"
        elif word1 ==  42 : word1 = "Mo"
        elif word1 ==  43 : word1 = "Tc"
        elif word1 ==  44 : word1 = "Ru"
        elif word1 ==  45 : word1 = "Rh"
        elif word1 ==  46 : word1 = "Pd"
        elif word1 ==  47 : word1 = "Ag"
        elif word1 ==  48 : word1 = "Cd"
        elif word1 ==  49 : word1 = "In"
        elif word1 ==  50 : word1 = "Sn"
        elif word1 ==  51 : word1 = "Sb"
        elif word1 ==  52 : word1 = "Te"
        elif word1 ==  53 : word1 = "I"
        elif word1 ==  54 : word1 = "Xe"
        elif word1 ==  55 : word1 = "Cs"
        elif word1 ==  56 : word1 = "Ba"
        elif word1 ==  57 : word1 = "La"
        elif word1 ==  58 : word1 = "Ce"
        elif word1 ==  59 : word1 = "Pr"
        elif word1 ==  60 : word1 = "Nd"
        elif word1 ==  61 : word1 = "Pm"
        elif word1 ==  62 : word1 = "Sm"
        elif word1 ==  63 : word1 = "Eu"
        elif word1 ==  64 : word1 = "Gd"
        elif word1 ==  65 : word1 = "Tb"
        elif word1 ==  66 : word1 = "Dy"
        elif word1 ==  67 : word1 = "Ho"
        elif word1 ==  68 : word1 = "Er"
        elif word1 ==  69 : word1 = "Tm"
        elif word1 ==  70 : word1 = "Yb"
        elif word1 ==  71 : word1 = "Lu"
        elif word1 ==  72 : word1 = "Hf"
        elif word1 ==  73 : word1 = "Ta"
        elif word1 ==  74 : word1 = "W"
        elif word1 ==  75 : word1 = "Re"
        elif word1 ==  76 : word1 = "Os"
        elif word1 ==  77 : word1 = "Ir"
        elif word1 ==  78 : word1 = "Pt"
        elif word1 ==  79 : word1 = "Au"
        elif word1 ==  80 : word1 = "Hg"
        elif word1 ==  81 : word1 = "Tl"
        elif word1 ==  82 : word1 = "Pb"
        elif word1 ==  83 : word1 = "Bi"
        elif word1 ==  84 : word1 = "Po"
        elif word1 ==  85 : word1 = "At"
        elif word1 ==  86 : word1 = "Rn"
        elif word1 ==  87 : word1 = "Fe"
        elif word1 ==  88 : word1 = "Ra"
        elif word1 ==  89 : word1 = "Ac"
        elif word1 ==  90 : word1 = "Th"
        ## copy from atom list.
        s = ("%s%s" % (word1,line[30:-1]))
        opennew.write(s)
        opennew.write("\n")
    openold.close()
    opennew.close()

def vector_points(end_x,end_y,end_z,start_x,start_y,start_z,numb_div):
    #/// Returns an array of (`n` + 1) equidistant points from `start` to `end`.
    x = []
    y = []
    z = []
    delta_x = (end_x - start_x) / float(numb_div)
    delta_y = (end_y - start_y) / float(numb_div)
    delta_z = (end_z - start_z) / float(numb_div)
    for k in range(0,numb_div):
        new_x =  start_x + (delta_x * float(k))
        new_y =  start_y + (delta_y * float(k))
        new_z =  start_z + (delta_z * float(k))
        x.append(new_x)
        y.append(new_y)
        z.append(new_z)
    return x,y,z

def PDBfile (fileName,index_elements, index_elements_num, axis_x, axis_y, axis_z, charge_NPA):
    newfile = str(fileName) + "_wiberg.pdb"
    opennew = open(newfile,"w")
    length  = len(axis_x)
    count   = 0
    opennew.write ("COMPND    " + str(fileName) + "\n")
    opennew.write ("AUTHOR    GENERATED BY Osvaldo Yanez-Osses\n")
    for i in range(0,length):
        count = count + 1
        arr = ['ATOM',str(count), index_elements_num[i], 'MOL', 'A', '1', axis_x[i], axis_y[i], axis_z[i], charge_NPA[i], 0.00, index_elements[i] ]
        arr[0] = arr[0].ljust(6)  #atom#6s
        arr[1] = arr[1].rjust(5)  #aomnum#5d
        arr[2] = arr[2].center(4) #atomname$#4s
        arr[3] = arr[3].ljust(3)  #resname#1s
        arr[4] = arr[4].rjust(1)  #Astring
        arr[5] = arr[5].rjust(4)  #resnum
        arr[6] = str('%8.3f' % (float(arr[6]))).rjust(8)  #x
        arr[7] = str('%8.3f' % (float(arr[7]))).rjust(8)  #y
        arr[8] = str('%8.3f' % (float(arr[8]))).rjust(8)  #z
        arr[9] = str('%6.2f' % (float(arr[9]))).rjust(6)  #occ
        arr[10]= str('%6.2f' % (float(arr[10]))).ljust(6) #temp
        arr[11]= arr[11].rjust(12)  #elname
        opennew.write ("{0}{1} {2} {3} {4}{5}    {6}{7}{8}{9}{10}{11}\n" .format(arr[0],arr[1],arr[2],arr[3],arr[4],arr[5],arr[6],arr[7],arr[8],arr[9],arr[10],arr[11]) )
    opennew.close()

def PDBfile_wiberg (fileName, axis_x, axis_y, axis_z, values, NAtoms ):
    newfile = str(fileName) + "_wiberg.pdb"
    opennew = open(newfile,"a")
    length  = len(axis_x)
    count   = NAtoms
    for i in range(0,length):
        count = count + 1
        arr = ['ATOM',str(count), 'XX', 'MOL', 'B', '2', axis_x[i], axis_y[i], axis_z[i], 1.00, values[i], 'Xx' ]
        arr[0] = arr[0].ljust(6)  #atom#6s
        arr[1] = arr[1].rjust(5)  #aomnum#5d
        arr[2] = arr[2].center(4) #atomname$#4s
        arr[3] = arr[3].ljust(3)  #resname#1s
        arr[4] = arr[4].rjust(1)  #Astring
        arr[5] = arr[5].rjust(4)  #resnum
        arr[6] = str('%8.3f' % (float(arr[6]))).rjust(8)  #x
        arr[7] = str('%8.3f' % (float(arr[7]))).rjust(8)  #y
        arr[8] = str('%8.3f' % (float(arr[8]))).rjust(8)  #z
        arr[9] = str('%6.2f' % (float(arr[9]))).rjust(6)  #occ
        arr[10]= str('%6.2f' % (float(arr[10]))).ljust(6) #temp
        arr[11]= arr[11].rjust(12)  #elname
        opennew.write ("{0}{1} {2} {3} {4}{5}    {6}{7}{8}{9}{10}{11}\n" .format(arr[0],arr[1],arr[2],arr[3],arr[4],arr[5],arr[6],arr[7],arr[8],arr[9],arr[10],arr[11]) )
    opennew.write ("END\n")
    opennew.close()

def VMD_representations (fileName,index_elements,unique_wiberg):
    fileName = fileName.replace(".", "")
    fileName = fileName.replace("\\", "")
    newfile  = str(fileName) + "_wiberg.pdb"
    newfilevmd = str(fileName) + "_VMD.vmd"
    opennew    = open(newfilevmd, "w")
    unique_elements = pd.unique(index_elements)
    length     = len(unique_elements)
    #
    opennew.write ("display resetview\n")
    opennew.write ("display projection   Orthographic\n")
    opennew.write ("display depthcue   off\n")
    opennew.write ("display shadows off\n")
    opennew.write ("display ambientocclusion off\n")
    opennew.write ("animate goto 0\n")
    opennew.write ("axes location off\n")
    opennew.write ("color Display Background white\n")
    opennew.write ("scale to 0.18\n")
    opennew.write ("mol load pdb " + str(newfile) + " \n")
    opennew.write ("mol delrep 0 top\n")
    for i in range(0,length):
        opennew.write ("mol representation CPK 0.800000 0.000000 600.000000 600.000000\n")
        opennew.write ("mol color Name\n")
        opennew.write ("mol selection {element " + str(unique_elements[i]) + " }\n")
        opennew.write ("mol material Glossy\n")
        opennew.write ("mol addrep top\n")
        opennew.write ("mol selupdate " + str(i) + " top 0\n")
        opennew.write ("mol colupdate " + str(i) + " top 0\n")
        opennew.write ("mol scaleminmax top " + str(i) + " 0.000000 0.000000\n")
    opennew.write ("mol representation CPK 0.200000 0.000000 600.000000 600.000000\n")
    opennew.write ("mol color Beta\n")
    opennew.write ("mol selection {element X}\n")
    opennew.write ("mol material AOChalky\n")
    opennew.write ("mol addrep top\n")
    opennew.write ("mol selupdate " + str(length) + " top 0\n")
    opennew.write ("mol colupdate " + str(length) + " top 0\n")
    opennew.write ("mol scaleminmax top " + str(length) + " 0.000000 " + str(unique_wiberg) + " \n")


'''
    Main
'''
# parameters wiberg
value_wiberg_threshold = 0.05
numb_div = 500

# Determine logfile
if len(sys.argv) == 1:
    print ("\nUsage:\n")
    print ("\tEnter name output gaussian 16 .log or .out\n")
    print ("\tpython Wiberg_VMD.py outputfileGaussian\n")
    sys.exit(1)
g09file = sys.argv[1]
# Nombre del archivo y extension
fileName, fileExtension = os.path.splitext(g09file)

np.set_printoptions(precision=2)

# Extract Wiberg matrix and Index from G16 .log file
get_Wiberg_matrix = False
get_Wiberg_index  = False
NAtoms            = 0

logfile = open(g09file,'r')
for text in logfile:
	words = text.split()
	if all(x in words for x in ['Wiberg','matrix','NAO','basis:']):
		get_Wiberg_matrix = True
	if all(x in words for x in ['Wiberg','Totals', 'atom:']):
		get_Wiberg_index = True
	# We check the number of atoms, which will indicate the size of the wiberg matrix.
	if all(x in words for x in ['NAtoms=','NQM=']):
		# Atom number
		NAtoms = int(words[1])
logfile.close()

'''
	 Build matrix
'''
logfile = open(g09file,'r')
data = logfile.read()
# Guarda solo lo que este entre las palabras claves "Overlap ***" y "*** Kinetic"
raw_string_matrix = re.findall(r'Wiberg bond index matrix in the NAO basis(.*?)Wiberg bond index, Totals by atom',data,re.DOTALL)
raw_string_matrix = itertools.chain(*[x.split(':\n\n') for x in raw_string_matrix])
tmp_raw_string_matrix = []
for value in raw_string_matrix:
    tmp_raw_string_matrix.append(value)
length = len (tmp_raw_string_matrix)
raw_string_matrix = tmp_raw_string_matrix[length-1]
raw_string_matrix = raw_string_matrix.split('\n')
#
index_segment       = []
raw_matrix_elements = []
count_segment = 0
#
for matrix_value in raw_string_matrix:
    if all(x in matrix_value for x in ['Atom']):
        index_segment.append(count_segment)
        continue
    count_segment = count_segment + 1
    raw_matrix_elements.append(matrix_value)
wiberg_matrix_total = []
for i in index_segment:
    segment_matrix = []
    for j in range (i+1,i+NAtoms+1):
        raw_wiberg_elements = raw_matrix_elements[j].split()
        del raw_wiberg_elements[:2]
        tmp_matrix = "\t".join(raw_wiberg_elements)
        segment_matrix.append(tmp_matrix)
    wiberg_matrix_total.append(segment_matrix)
length = len(wiberg_matrix_total)
#
# Funcion para juntar segementos de matrices
wiberg_matrix = []
for k in range(0,NAtoms):
    values_w = ''
    for i in range(0,length):
        for j in range(0,NAtoms):
            values_w+=(" " + wiberg_matrix_total[i][k])
            break
    list_values = values_w.split()
    wiberg_matrix.append(list_values)
matrix_elements = []
# Tamaño de la matrix segun numero de atomos
for i in range(0,NAtoms):
    for j in range(0,NAtoms):
        matrix_elements.append(wiberg_matrix[i][j])
unique_wiberg_values = pd.unique(matrix_elements)
unique_wiberg_values.sort()
#
raw_string_index   = re.findall(r'Wiberg bond index, Totals by atom(.*?)NBI: Natural Binding',data,re.DOTALL)
raw_string_index = itertools.chain(*[x.split(':\n\n') for x in raw_string_index])
tmp_raw_string_index = []
for value in raw_string_index:
    tmp_raw_string_index.append(value)
length = len (tmp_raw_string_index)
raw_string_index = tmp_raw_string_index[length-1]
raw_elements_index = raw_string_index.split('\n')
index_elements     = []
index_elements_num = []
count_element      = 1
length = len(raw_elements_index)
for i in range(2,length-3):
    tmp_index = raw_elements_index[i].split()
    numb = str(tmp_index[1]) + str(count_element)
    #print (numb)
    index_elements_num.append(numb)
    index_elements.append(tmp_index[1])
    count_element = count_element + 1
logfile.close()
#
# Guarda solo lo que este entre las palabras claves "Summary of Natural Population Analysis:" y "* Total"
raw_NPA_string = re.findall(r'Summary of Natural Population Analysis:(.*?)\* Total \*',data,re.DOTALL)
raw_NPA_string = itertools.chain(*[x.split(':\n\n') for x in raw_NPA_string])
tmp_raw_NPA_string = []
for value in raw_NPA_string:
    tmp_raw_NPA_string.append(value)
length = len (tmp_raw_NPA_string)
raw_NPA_string = tmp_raw_NPA_string[length-1]
raw_NPA_string = raw_NPA_string.split('\n')


Charge_NPA = []
length = len(raw_NPA_string)
for i in range (6,length-2):
    #print (raw_NPA_string[i])
    parameter_mol = raw_NPA_string[i].split()
    Charge_NPA.append(parameter_mol[2])
# Coordinates
Gaussian_coords(g09file,fileName,NAtoms)
axis_x, axis_y, axis_z = coordsXYZ(fileName)
# Generated PDB file
PDBfile (fileName,index_elements, index_elements_num, axis_x, axis_y, axis_z, Charge_NPA)
# Recorrer matrix
x_arr = []
y_arr = []
z_arr = []
values_arr = []
length = len(wiberg_matrix)
for i in range(0,length):
    for j in range(0,length):
        if i < j:
            value = float(wiberg_matrix[i][j])
            #print (" " + str(value) + " " + index_elements_num[i] + " " + index_elements_num[j] )
            if (value > value_wiberg_threshold):
                x,y,z = vector_points(axis_x[j],axis_y[j],axis_z[j],axis_x[i],axis_y[i],axis_z[i],numb_div)
                for k in range(0,len(x)):
                    x_arr.append(x[k])
                    y_arr.append(y[k])
                    z_arr.append(z[k])
                    values_arr.append(value)

PDBfile_wiberg (fileName, x_arr, y_arr, z_arr, values_arr, NAtoms )
VMD_representations (fileName,index_elements, unique_wiberg_values[-1])


print ("Scale Wiberg index bond:\n")
print ("\t" + str(unique_wiberg_values[1]) + " to " + str (unique_wiberg_values[-1]) + "\n")
print ("\tWiberg threshold: " + str(value_wiberg_threshold) + "\n")

quit()
