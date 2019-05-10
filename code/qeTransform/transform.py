import os
import re
from math import sqrt
import numpy as np
from numpy.linalg import inv
import pandas as pd
from collections import Counter


class PwInput():
    _pw = 'pw.x'

    def __init__(self, filename=None):
        """initial defs"""
        # kpoints
        self.ktype = "automatic"
        self.kpoints = [1, 1, 1]
        self.shiftk = [0, 0, 0]
        self.klist = []
        # dictionaries
        self.control = dict()
        self.system = dict()
        self.electrons = dict()
        self.ions = dict()
        self.cell = dict()
        self.atypes = dict()
        self.atoms = []
        self.atyp = []
        self.cell_parameters = []
        self.cell_units = 'angstrom'
        self.atomic_pos_type = 'angstrom'

        # from a reference file
        if filename:
            f = open(filename, "r")
            self.file_name = filename
            self.file_lines = f.readlines()
            self.singlespaces()
            self.remove_empty_lines()
            self.store(self.control, "control")  # read &control as dicts
            self.store(self.system, "system")  # read &system ..
            self.store(self.electrons, "electrons")  # read &electrons ..
            self.store(self.ions, "ions")  # read &ions ..
            self.store(self.cell, "cell")  # read &cells ..
            # read ATOMIC_SPECIES
            self.read_atomicspecies()
            # read ATOMIC_POSITIONS
            self.read_atoms()
            # read K_POINTS
            self.read_kpoints()
            # read CELL_PARAMETERS
            self.read_cell_parameters()

    def get_atoms(self):
        """ Get the lattice parameters, postions of the atoms and chemical symbols
        """
        self.read_cell_parameters()
        cell = self.cell_parameters
        sym = [atom[0] for atom in self.atoms]
        pos = [atom[1] for atom in self.atoms]
        if self.atomic_pos_type == 'bohr':
            pos = car_red(pos, cell)
        return cell, pos, sym

    def get_masses(self):
        """ Get an array with the masses of all the atoms
        """
        masses = []
        for atom in self.atoms:
            atype = self.atypes[atom[0]]
            mass = float(atype[0])
            masses.append(mass)
        return masses

    def get_symmetry(self):
        """get the symmetry group of this system"""
        import spglib

        lat, positions, atypes = self.get_atoms()
        lat = np.array(lat)

        at = np.unique(atypes)
        an = dict(list(zip(at, list(range(len(at))))))
        atypes = [an[a] for a in atypes]

        cell = (lat, positions, atypes)

        spacegroup = spglib.get_spacegroup(cell, symprec=1e-5)
        return spacegroup

    def read_atomicspecies(self):
        """find ATOMIC_SPECIES keyword in file line by line"""
        lines = iter(self.file_lines)
        for line in lines:
            if "ATOMIC_SPECIES" in line:
                for i in range(int(self.system["ntyp"])):
                    atype, mass, psp = next(lines).split()
                    self.atypes[atype] = [mass, psp]  # e.g. atypes['C']=['12.017',C_1star-blabla.UPF]

    def read_atoms(self):
        """find ATOMIC_POSITIONS keyword in file and read next lines"""
        lines = iter(self.file_lines)
        for line in lines:
            if "ATOMIC_POSITIONS" in line:
                atomic_pos_type = line
                self.atomic_pos_type = re.findall('([A-Za-z]+)', line)[-1]  # last word in line: angstrom or crystal
                for i in range(int(self.system["nat"])):
                    atype, x, y, z = next(lines).split()
                    self.atyp.append(np.array(atype))
                    self.atoms.append(np.array([float(i) for i in (x, y, z)]))
        self.atomic_pos_type = atomic_pos_type.replace('{', '').replace('}', '').strip().split()[1]

    def read_cell_parameters(self):
        ibrav = int(self.system['ibrav'])

        def rmchar(string, symbols):
            return ''.join([c for c in string if c not in symbols])

        if ibrav == 0:
            if 'celldm(1)' in list(self.system.keys()):
                a = float(self.system['celldm(1)'])
            else:
                a = 1
            lines = iter(self.file_lines)
            for line in lines:
                if "CELL_PARAMETERS" in line:
                    units = rmchar(line.strip(), '{}()').split()
                    self.cell_parameters = [[], [], []]
                    if len(units) > 1:
                        self.cell_units = units[1]
                    else:
                        self.cell_units = 'bohr'
                    for i in range(3):
                        self.cell_parameters[i] = np.array([float(x) * a for x in next(lines).split()])
            if self.cell_units == 'angstrom' or self.cell_units == 'bohr':
                if 'celldm(1)' in self.system: del self.system['celldm(1)']
            if 'celldm(1)' not in list(self.system.keys()):
                a = np.linalg.norm(self.cell_parameters[0])
        elif ibrav == 1:
            a = float(self.system['celldm(1)'])
            self.cell_parameters = [[a, 0, 0],
                                    [0, a, 0],
                                    [0, 0, a]]
        elif ibrav == 2:
            a = float(self.system['celldm(1)'])
            self.cell_parameters = [[-a / 2, 0, a / 2],
                                    [0, a / 2, a / 2],
                                    [-a / 2, a / 2, 0]]
        elif ibrav == 3:
            a = float(self.system['celldm(1)'])
            self.cell_parameters = [[a / 2, a / 2, a / 2],
                                    [-a / 2, a / 2, a / 2],
                                    [-a / 2, -a / 2, a / 2]]
        elif ibrav == 4:
            a = float(self.system['celldm(1)'])
            c = float(self.system['celldm(3)'])
            self.cell_parameters = [[a, 0, 0],
                                    [-a / 2, sqrt(3) / 2 * a, 0],
                                    [0, 0, c * a]]
        elif ibrav == 6:
            a = float(self.system['celldm(1)'])
            c = float(self.system['celldm(3)'])
            self.cell_parameters = [[a, 0, 0],
                                    [0, a, 0],
                                    [0, 0, c * a]]
        else:
            raise NotImplementedError('ibrav = %d not implemented' % ibrav)
        self.alat = a

    def read_kpoints(self):
        """find K_POINTS keyword in file and read next lines"""
        lines = iter(self.file_lines)
        for line in lines:
            if "K_POINTS" in line:
                if "automatic" in line:
                    self.ktype = "automatic"
                    vals = list(map(float, next(lines).split()))
                    self.kpoints, self.shiftk = vals[0:3], vals[3:6]
                # otherwise read a list
                else:
                    # read number of kpoints
                    nkpoints = int(next(lines).split()[0])
                    self.klist = []
                    self.ktype = ""
                    try:
                        lines_list = list(lines)
                        for n in range(nkpoints):
                            vals = lines_list[n].split()[:4]
                            self.klist.append(list(map(float, vals)))
                    except IndexError:
                        print("wrong k-points list format")
                        exit()

    def remove_key(self, group, key):
        """ if a certain key exists in the group, remove it"""
        if key in list(group.items()):
            del group[key]

    def remove_empty_lines(self):
        """removes empty lines and \n in a file"""
        if not os.path.isfile(self.file_name):
            print("{} does not exist ".format(self.file_name))
            return
        with open(self.file_name) as filehandle:
            lines = filehandle.readlines()
        with open(self.file_name, 'w') as filehandle:
            lines = filter(lambda x: x.strip() and x.strip('\n'), lines)
            filehandle.writelines(lines)

    def run(self, filename, procs=1, folder='.'):
        """ this function is used to run this job locally
        """
        os.system('mkdir -p %s' % folder)
        self.write("%s/%s" % (folder, filename))
        if procs == 1:
            os.system('cd %s; OMP_NUM_THREADS=1 %s -inp %s > %s.log' % (folder, self._pw, filename, filename))
        else:
            os.system('cd %s; OMP_NUM_THREADS=1 mpirun -np %d %s -inp %s > %s.log' % (
                folder, procs, self._pw, filename, filename))

    def set_atoms_string(self, string):
        """
        set the atomic postions using string of the form
        Si 0.0 0.0 0.0
        Si 0.5 0.5 0.5
        """
        atoms_str = [line.strip().split() for line in string.strip().split('\n')]
        self.atoms = []
        for atype, x, y, z in atoms_str:
            self.atoms.append([atype, list(map(float, [x, y, z]))])

    def set_atoms(self, atoms):
        """ set the atomic postions using a Atoms datastructure from ase
        """
        # we will write down the cell parameters explicitly
        self.system['ibrav'] = 0
        if 'celldm(1)' in self.system: del self.system['celldm(1)']
        self.cell_parameters = atoms.get_cell()
        self.atoms = list(zip(atoms.get_chemical_symbols(), atoms.get_scaled_positions()))
        self.system['nat'] = len(self.atoms)

    def singlespaces(self):
        for line in self.file_lines:
            # ' '.join(line.split())
            re.sub(' +', ' ', line)

    def slicefile(self, keyword):
        """find lines after keyword between & and /"""
        lines = re.findall('&%s(?:.?)+\n((?:.+\n)+?)(?:\s+)?[\/&]' % keyword, "".join(self.file_lines), re.MULTILINE)
        return lines

    def store(self, group, name):
        """
        Save the variables specified in each of the groups on the structure, separated by =
        """
        for file_slice in self.slicefile(name):
            for keyword, value in re.findall('([a-zA-Z_0-9_\(\)]+)(?:\s+)?=(?:\s+)?([a-zA-Z\'"0-9_.-]+)', file_slice):
                group[keyword.strip()] = value.strip()

    def stringify_group(self, keyword, group):
        if group != {}:
            string = '&%s\n' % keyword
            for keyword in sorted(group):  # Py2/3 discrepancy in keyword order
                string += "%20s = %s\n" % (keyword, group[keyword])
            string += "/&end\n"
            return string
        else:
            return ''

    def write(self, filename):
        f = open(filename, 'w')
        f.write(str(self))
        f.close()

    def __str__(self):
        """Output the file in the form of a string"""
        string = ''
        string += self.stringify_group("control", self.control)  # print control
        string += self.stringify_group("system", self.system)  # print system
        string += self.stringify_group("electrons", self.electrons)  # print electrons
        string += self.stringify_group("ions", self.ions)  # print ions
        string += self.stringify_group("cell", self.cell)  # print ions

        # print atomic species
        string += "ATOMIC_SPECIES\n"
        for atype in self.atypes:
            string += " %3s %8s %20s\n" % (atype, self.atypes[atype][0], self.atypes[atype][1])
        # print atomic positions
        string += "ATOMIC_POSITIONS { %s }\n" % self.atomic_pos_type
        for atom in self.atoms:
            string += "%3s %14.10lf %14.10lf %14.10lf\n" % (atom[0], atom[1][0], atom[1][1], atom[1][2])
        # print kpoints
        if self.ktype == "automatic":
            string += "K_POINTS { %s }\n" % self.ktype
            string += ("%3d" * 6 + "\n") % tuple(self.kpoints + self.shiftk)
        elif self.ktype == "crystal":
            string += "K_POINTS { %s }\n" % self.ktype
            string += "%d\n" % len(self.klist)
            for i in self.klist:
                string += ('%12.8lf ' * 4 + '\n') % tuple(i)
        else:
            string += "K_POINTS { }\n"
            string += "%d\n" % len(self.klist)
            for i in self.klist:
                string += (("%12.8lf " * 4) + "\n") % tuple(i)
        if self.system['ibrav'] == 0 or self.system['ibrav'] == '0':
            string += "CELL_PARAMETERS %s\n" % self.cell_units
            for i in range(3):
                string += ("%14.10lf " * 3 + "\n") % tuple(self.cell_parameters[i])
        return string

    def df_file(self):
        """combines all input parameters from one file to one dataframe """
        ctrl_df = pd.DataFrame.from_dict(self.control, orient='index')
        syst_df = pd.DataFrame.from_dict(self.system, orient='index')
        electron_df = pd.DataFrame.from_dict(self.electrons, orient='index')
        atom_species = pd.DataFrame.from_dict(self.atypes, orient='index', columns=['mass', 'psp'])
        atom_pos = pd.DataFrame(self.atoms, columns=['atom', 'xyz'])
        k_pts = pd.DataFrame(self.kpoints, index=['ka', 'kb', 'kc'], columns=['k'])
        cell_param = pd.DataFrame(self.cell_parameters, index=['a', 'b', 'c'], columns=['va', 'vb', 'vc'])
        groups = [ctrl_df, syst_df, electron_df, atom_species, atom_pos, k_pts, cell_param]
        return pd.concat(groups,
                         keys=['ctrl_df', 'syst_df', 'electron_df', 'atom_species', 'atom_pos', 'k_pts', 'cell_param'])


def rename_duplicates(mylist):
    counts = Counter(mylist)  # so we have: {'name':3, 'state':1, 'city':1, 'zip':2}
    for s, num in counts.items():
        if num > 1:  # ignore strings that only appear once
            for suffix in range(1, num + 1):  # suffix starts at 1 and increases by 1 each time
                mylist[mylist.index(s)] = s + '(' + str(suffix) + ')'  # replace each appearance of s


def cart_to_cryst(atname, atoms, lattice):
    tfnm = pd.DataFrame(np.dot(atoms, lattice))
    at_n = pd.DataFrame(atname)
    df = pd.concat([at_n, tfnm], axis=1)
    return df


# Mymatrix = np.array([["a", "b"], ["c", "d"]])
# Mycol = np.array([0.4, 0.6])
#
# dt = np.dtype([('col0', 'U1'), ('col1', 'U1'), ('col2', float)])
# new_recarr = np.empty((2,), dtype=dt)
# new_recarr['col0'] = Mymatrix[:, 0]
# new_recarr['col1'] = Mymatrix[:, 1]
# new_recarr['col2'] = Mycol[:]
# print(new_recarr)


def cryst_to_cart(atname, atoms, lattice): return np.dot(atoms, inv(lattice))


if __name__ == '__main__':
    path = r'C:\Users\hwahab\Documents\qe-stuff\data\kgec\g_copy\pw.in'
    data = PwInput(path)

    new_atpos = cart_to_cryst(data.atyp, data.atoms, data.cell_parameters)
    # print(new_atpos.values)
# new_coords = cart_to_cryst(data.atyp, data.atoms, data.cell_parameters)

for item in new_atpos.values:
    print(item[0], ' '.join(map(str, item[1:])))

# r'C:\Users\hwahab\Documents\qe-stuff\data\kgec\g_copy\pw.in'

# f1_n = os.path.basename(os.path.dirname(path))
# data.control['prefix'] = "'gc'"
# data.system['nat'] = 2
# data.system['ntyp'] = 2
# data.system['tot_charge'] = 1
# data.kpoints = [1,1,1]


#
