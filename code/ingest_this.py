import os
import pandas as pd
import fileinput
from pypif import pif
from pypif.obj import *
            
def remove_empty_lines(filename):
    """removes empty lines and \n in a file"""
    if not os.path.isfile(filename):
        print("{} does not exist ".format(filename))
        return
    with open(filename) as filehandle:
        lines = filehandle.readlines()

    with open(filename, 'w') as filehandle:
        lines = filter(lambda x: x.strip() and x.strip('\n'), lines)
        filehandle.writelines(lines)
        
def check_int(obj,fdir,filename):
    """checks if the object is an integer"""
    try:
        int(obj)
    except ValueError:
        print('Warning: "%s" in line %d in "%s" is not an integer' %(obj,fdir.filelineno(),filename))
        
def check_float(obj,fdir,filename):
    """checks if the object is a float"""
    try:
        float(obj)
    except ValueError:
        print('Warning: "%d"in line %d in "%s" is not a float' %(obj,fdir.filelineno(),filename))

def check_notafloat(obj):
    """checks if the object is a str"""
    try: float(obj)
    except ValueError:
        return 1
    else:
        return 0


"""loads files in a directory into a list"""
file_list = os.listdir('.') #define target directory
new_file_list = [] 
fdir = []

for names in file_list:
    if names.endswith(".xyz"):
        remove_empty_lines(names)
        new_file_list.append(names)
        fdir = fileinput.FileInput(new_file_list)
        
#print(fdir.filelineno())

def read_xyz(filename):
    """Read filename in XYZ format and return series of identifiers, properties, frequencies,
    atoms,coordinates and partial charges.
    
    Warnings are raised if:
    If the first line of a file doesn't have a single field.
    If the second line doesn't have 17 fields or does not contain the tag "gdb"
    If the frequency value in the first field of the corresponding line is not a float.
    If the SMILES line doesn't have 2 fields
    If the InChI line doesn't include the word InChI
    If no. of coordinates does not agree with the no. of atoms stated in the file."""

    atoms = []
    coordinates = []
    dcharge = []
    prop = []
    freq = []  
    
    n_a = fdir.readline()   #define 1st line as no of atoms
#    print(fdir.filelineno())
    check_int(n_a,fdir,filename)    #check if value is an integer
    if len(n_a.split()) != 1:
        print ('Warning: Expected line %d in "%s" to have 1 field "No. of atoms"' %(fdir.filelineno(),filename))
        return
    else:
        n_atoms = n_a.split()
                
    properties = fdir.readline()    #define 2nd line as properties
#    print(fdir.filelineno())
    if len(properties.split()) != 17:
        print ('Warning: Expected line %d in "%s" to have 17 fields "Properties"' %(fdir.filelineno(),filename))
        return
    elif 'gdb' not in properties.split()[0]:
        print('Warning: Expected line %d field %d in "%s" to be "gdb"' %(fdir.filelineno(),int(prop[0],filename)))
        return
    else:
        prop = properties.split()
        
    #lists next lines with 5 fields in atoms, with the 1st field not a float, coordinates and partial charge
    for line in fdir:
        if len(line.split()) == 5 and check_notafloat(line[0]) == 1:
            atom,x,y,z,de = line.split()
            atoms.append(atom)
            coordinates.append([float(x), float(y), float(z)])
            dcharge.append(float(de))
#            print(fdir.filelineno())
        else: break
#    print(fdir.filelineno())
    if int(n_atoms[0]) != len(coordinates):
        raise ValueError('File says %d atoms but read %d points' % (int(n_atoms[0]), len(coordinates)))
    
    freq = line.split()     #lists next line as frequencies
    check_float(freq[0],fdir,filename)  #Checks if the first index is a float

    smiles = fdir.readline() #defines next line as SMILES
#    print(fdir.filelineno())
    if len(smiles.split()) !=2:
        print ('Warning: Expected line %d in "%s" to have 2 fields "SMILES"' %(fdir.filelineno(),filename))
        return
    elif check_notafloat(smiles.split()[0]) != 1:
        print ('Warning: Expected line %d in "%s" to be a string' %(fdir.filelineno(),filename))
        return
    else:
        sm = smiles.split()[0]
    
    idchem = fdir.readline() #defines last line as InChI
    if 'InChI' not in idchem.split()[0]:
        print ('Warning: line',fdir.filelineno(),'in',filename,'should have include "InChI"')
        return
    else:   
        inchi = idchem.split()[0]

    fileinput.nextfile()    #close line-by-line reading iterations, begins for next file
    
    return n_atoms, prop, atoms, coordinates, dcharge, freq, sm, inchi
	
"""Converts .xyz lists into series"""
def Series_me(filename):

    n_atoms, prop, atoms, coordinates, dcharge, freq, sm, inchi = read_xyz(filename)
    s1 = pd.Series(n_atoms)
    s1 = s1.replace('\n','', regex=True)    #removes \n for clarity  
    s2 = pd.Series(sm)
    s3 = pd.Series(inchi)
    s4 = pd.Series(prop)
    s5 = pd.Series(freq)
    s6 = pd.Series(atoms)
    s7 = pd.Series(coordinates)
    s8 = pd.Series(dcharge)

    return s1,s2,s3,s4,s5,s6,s7,s8

"""Combines series into a single dataframe, one for each file"""
def DataFrame_me(filename):
    s1,s2,s3,s4,s5,s6,s7,s8 = Series_me(filename)
    df = pd.concat([s1,s2,s3,s4,s5,s6,s7,s8],axis=1,ignore_index=True)
    df.columns = ['no. of atoms', 'SMILES','InChI','properties','frequencies','atoms','XYZ','partial charge']
    return df

"""combines dataframes from multiple files in a directory"""
dfs=[]
for j in new_file_list:
    dfs.append(DataFrame_me(j))
alldf = pd.concat(dfs,axis=1)

"""Ingester function converting each series in the dataframe into a PIF object,
    returns JSON """
def ingester(dataframe):
    out=[]
    rot=[]
    fr=[]
    ap=[]
    
    for i,j in zip(range(len(new_file_list)),range(0,dataframe.shape[1],8)):    #simultaneous loop for filename tagging
        syst = System()
        chem_syst = ChemicalSystem()
        syst.tags = new_file_list[i]    #tags filename
            
        chem_syst.tags = dataframe.iloc[1,j+3]  #"gdb9"; string constant to ease extraction via grep
        chem_syst.uid = dataframe.iloc[0,j+3]   #Consecutive, 1-based integer identifier of molecule
        chem_syst.ids = dataframe.iloc[0,j+2]   #InChI for GDB9 and for relaxed geometry
        chem_syst.names = dataframe.iloc[0,j+1] #SMILES from GDB9 and for relaxed geometry
        
        #Chemical System properties 
        no_atoms = Property()                   
        no_atoms.name = 'Number of atoms'
        no_atoms.scalars = dataframe.iloc[0,j]

        dip_mom = Property()
        dip_mom.name = 'Dipole moment'
        dip_mom.scalars = dataframe.iloc[0,j+3]
        dip_mom.units = 'Debye'
        
        homo = Property()
        homo.name = 'HOMO'
        homo.scalars = dataframe.iloc[7,j+3]
        homo.units = 'Hartree'
        
        lumo = Property()
        lumo.name = 'LUMO'
        lumo.scalars = dataframe.iloc[8,j+3]
        lumo.units = 'Hartree'
        
        band_gap = Property()
        band_gap.name = 'Band gap'
        band_gap.scalars = dataframe.iloc[9,j+3]
        band_gap.units = 'eV'
        
        zpve = Property()
        zpve.name = 'Zero point vibrational energy'
        zpve.scalars = dataframe.iloc[11,j+3]
        zpve.units = 'Hartree'
        
        u_zero = Property()
        u_zero.name = 'Internal energy at 0 K'
        u_zero.scalars = dataframe.iloc[12,j+3]
        u_zero.units = 'Hartree'
        
        u = Property()
        u.name = 'Internal energy at 298.15 K'
        u.scalars = dataframe.iloc[13,j+3]
        u.units = 'Hartree'
        
        H = Property()
        H.name = 'Enthalpy at 298.15 K'
        H.scalars = dataframe.iloc[14,j+3]
        H.units = 'Hartree'
        
        G = Property()
        G.name = 'Free energy at 298.15 K'
        G.scalars = dataframe.iloc[15,j+3]
        G.units = 'Hartree'
        
        Cv = Property()
        Cv.name = 'Heat capacity at 298.15 K'
        Cv.scalars = dataframe.iloc[16,j+3]
        Cv.units = 'cal/(mol K)'
        
        atom_pos = Property()
        atom_pos.name  = 'Element, XYZ-coordinates, partial charge'
        for k in range(5,8):
            ap.append(dataframe.iloc[:dataframe.count()[j+5],j+k].tolist()) #combines element, XYZ and partial charge
        atom_pos.matrices = ap
        ap =[]
        
        #System properties
        rotconst = Property()
        rotconst.name  = 'Rotational constants'
        rot.append(dataframe.iloc[2:5,j+3].tolist())
        rotconst.vectors = [['a','b','c'],rot[int(j/8)]]    #label the rotation vectors in a, b and c direction
        
        isopol = Property()
        isopol.name  = 'Isotropic polarizability'
        isopol.scalars = dataframe.iloc[6,j+3]
        isopol.units = 'Bohr^3'
        
        espat = Property()
        espat.name = 'Electronic spatial extent'
        espat.scalars = dataframe.iloc[10,j+3]
        espat.units = 'Bohr^3'
        
        freqs = Property()
        freqs.name = 'Molecular vibrational frequencies'
        fr.append(dataframe.iloc[:dataframe.count()[j+4],j+4].tolist())
        freqs.vectors = fr[int(j/8)]
        freqs.units = 'cm^-1'
         
        syst.properties = [rotconst,isopol,espat,freqs]
        chem_syst.properties = [no_atoms,dip_mom,homo,lumo,band_gap,u_zero,u,H,G,Cv,atom_pos]
        
        #converts list of system and chemical system PIF objects into JSON string
        out.append(pif.dumps([syst,chem_syst], indent=4))   
        str = ''.join(out)  #outputs a string

    return str
        
""" run Dataframe as a parameter into ingester function and outputs a json text file """
jsonData = ingester(alldf)

with open('test.json', 'w') as f:
    f.write(jsonData)