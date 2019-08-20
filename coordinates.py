import os
import numpy as np
import database as db
from database import Atom

class XYZfile():
    """
    Class containing the XYZ trajectory file.
    """
    def __init__(self, filename, periodic_table, boxside = 0):
        self._periodic_table = periodic_table
        self.atoms, self.frames, self.energies, self.nat = self.loadall(filename)
        self.nframes = self.nframes()
        self.atomtypes = self.atomtypes()
        self.atcount = self.atcount()
        self.pbc = boxside

    def loadall(self, filename):
        """
        Loads the XYZ trajectory file.
        
        Parameters
        ----------
        filename: str
                Name of the XYZ file.

        Returns
        -------
        atoms: list
                List of int of the atomic numbers in order of appearance.
        frames: list
                List of np.arrays of the coordinates.
        energies: list
                List of str containing the comment line of the XYZ file.
        nat: int
                Number of atoms.
        """
        with open(filename, 'r') as xyzfile:
            spl_file = xyzfile.readlines()
        atoms, frames, energies = [], [], []
        nat = int(spl_file[0])
        splitted_traj = [spl_file[line_n*(nat + 2):line_n*(nat + 2) + (nat + 2)] for line_n in
                range(int(len(spl_file) / (nat + 2)))]
        for at in splitted_traj[0][2:]:
            atom = at.split()
            atoms.append(self._periodic_table[db.get_name(atom[0], self._periodic_table)].atnum)
        for frame in splitted_traj:
            conf = Frame(frame, en = True)
            frames.append(conf.coordset)
            energies.append(conf.energy)
        return atoms, frames, energies, nat

    def nframes(self):
        """
        Counts the frames in the trajectory.

        Returns
        -------
        int
                Number of frames in the XYZ file.
        """
        return len(self.frames)

    def atomtypes(self):
        return set(self.atoms)

    def atcount(self):
        """ 
        Counts the atom frequency.

        Returns
        -------
        atcount: dict
                Dict with the atom labels as keys and the frequencies as values.
        """
        at_cnt = {}
        for at_typ in self.atomtypes:
            cnt = 0
            for atom in self.atoms:
                if atom == at_typ:
                    cnt += 1
            at_cnt[db.get_name(at_typ, self._periodic_table)] = cnt
            del cnt
        return at_cnt

    def coordprint(self, frame, save = False):
        if save:
            with open('output.xyz', 'a+') as output:
                output.write(F" {self.nat}\n {' '.join(self.energies[frame].split())}\n")
                for at, coord in enumerate(self.frames[frame]):
                    label = self._periodic_table[db.get_name(self.atoms[at], self._periodic_table)].label
                    output.write(F'{label} {coord[0]} {coord[1]} {coord[2]}\n')
        else:
            print(F" {self.nat}\n {' '.join(self.energies[frame].split())}")
            for at, coord in enumerate(self.frames[frame]):
                label = self._periodic_table[db.get_name(self.atoms[at], self._periodic_table)].label
                print(F'{label} {coord[0]} {coord[1]} {coord[2]}')

class Frame():
    def __init__(self, conf, en = False):
        self.coordset = self.coordset(conf[2:])
        if en == True:
            self.energy = conf[1]
        else:
            self.energy = 0

    def coordset(self, coords):
        coordset = []
        for at in coords:
            atom = at.split()
            coordset.append(np.array([float(coord) for coord in atom[1:]]))
        return coordset
