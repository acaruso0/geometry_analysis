import os
import numpy as np
import database as db
from database import periodic_table
import math

class XYZfile():
    """
    Class containing the XYZ trajectory file.
    """
    def __init__(self, filename, boxside = (0, 0, 0)):
        self.pbc = boxside
        self.atoms, self.frames, self.energies, self.nat = self.loadall(filename)
        self.nframes = self.nframes()
        self.atomtypes = self.atomtypes()
        self.atcount = self.atcount()

    class XYZFrame():
        def __init__(self, conf, boxside, en = False):
            self.pbc = boxside
            self.coordset = self.coordset(conf[2:])
            if en == True:
                self.energy = conf[1]
            else:
                self.energy = 0

        def coordset(self, coords):
            coordset = []
            for at in coords:
                atom = at.split()
                coord = np.array([float(coord) for coord in atom[1:]])

                ################################
                ### MINIMUM IMAGE CONVENTION ###
                ################################
                
                if self.pbc != (0, 0, 0):
                    for q in range(3):
                        coord[q] -= math.floor(coord[q] / self.pbc[q]) * self.pbc[q]

                ################################
                ################################

                #coordset.append(np.array([float(coord) for coord in atom[1:]]))
                coordset.append(coord)
            return coordset

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
            atoms.append(periodic_table[db.get_name(atom[0])].atnum)
        for frame in splitted_traj:
            conf = self.XYZFrame(frame, self.pbc, en = True)
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
            at_cnt[db.get_name(at_typ)] = cnt
            del cnt
        return at_cnt

    def wrap_h2o(self, frame):

        def pbc(oxygen, hydrogen):
            xx = np.linalg.norm(oxygen[0] - hydrogen[0])
            yy = np.linalg.norm(oxygen[1] - hydrogen[1])
            zz = np.linalg.norm(oxygen[2] - hydrogen[2])
            
            delta = [xx, yy, zz]
            delta = np.array([delta[i] - round(delta[i] / self.pbc[i]) * self.pbc[i] for i in range(3)])
            print(delta)

            hydrogen = oxygen + delta
            print(hydrogen)
            return hydrogen

        for at in range(self.nat):
            if self.atoms[at] == 8:
                self.frames[frame][at + 1] = pbc(self.frames[frame][at], self.frames[frame][at + 1])
                self.frames[frame][at + 2] = pbc(self.frames[frame][at], self.frames[frame][at + 2])

        with open('wrapped.xyz', 'a+') as output:
            output.write(F" {self.nat}\n {' '.join(self.energies[frame].split())}\n")
            for at, coord in enumerate(self.frames[frame]):
                label = periodic_table[db.get_name(self.atoms[at])].label
                output.write(F'{label} {coord[0]} {coord[1]} {coord[2]}\n')

    def coordprint(self, frame, save = False):
        if save:
            with open('output.xyz', 'a+') as output:
                output.write(F" {self.nat}\n {' '.join(self.energies[frame].split())}\n")
                for at, coord in enumerate(self.frames[frame]):
                    label = periodic_table[db.get_name(self.atoms[at])].label
                    output.write(F'{label:<2} {coord[0]:13.8f} {coord[1]:13.8f} {coord[2]:13.8f}\n')
        else:
            print(F" {self.nat}\n {' '.join(self.energies[frame].split())}")
            for at, coord in enumerate(self.frames[frame]):
                label = periodic_table[db.get_name(self.atoms[at])].label
                print(F'{label:<2} {coord[0]:13.8f} {coord[1]:13.8f} {coord[2]:13.8f}')

class History():
    def __init__(self, filename):
        self.atoms, self.frames, self.nat, self.pbc = self.loadall(filename)

    def loadall(self, filename):
        with open(filename, 'r') as hisfile:
            spl_file = hisfile.readlines()
        atoms, frames = [], []
        pbc = (float(spl_file[1].split()[0]), float(spl_file[2].split()[1]), float(spl_file[3].split()[2]))
        nat = int(spl_file[0].split()[2])
        splitted_traj = [spl_file[line_n*(2*nat + 4):line_n*(2*nat + 4) + (2*nat + 4)] for line_n in range(int(len(spl_file) / (2*nat + 4)))]
        for at_n in range(nat):
            atom = splitted_traj[0][4:][2*at_n].split()[0]

            if atom == 'OW' or atom == 'HW':
                atom = atom[0]

            atoms.append(periodic_table[db.get_name(atom)].atnum)
        for frame in splitted_traj:
            coordset = []
            for at_n in range(nat):
                coords = frame[4:][2*at_n + 1].split()
                coordset.append(np.array([float(coord) for coord in coords]))
            frames.append(coordset)
        return atoms, frames, nat, pbc
