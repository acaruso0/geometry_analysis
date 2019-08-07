import numpy as np

class XYZfile():
    def __init__(self, filename):
        self._frames, self.energies, self.nat = self.loadall(filename)
        self.nframes = self.nframes()
        self.atomtypes = self.atomtypes()

    def loadall(self, filename):
        with open(filename, 'r') as xyzfile:
            spl_file = xyzfile.readlines()
        frames, energies = [], []
        nat = int(spl_file[0])
        splitted_traj = [spl_file[line_n*(nat + 2):line_n*(nat + 2) + (nat + 2)] for line_n in
                range(int(len(spl_file) / (nat + 2)))]
        for frame in splitted_traj:
            conf = Frame(frame, en = True)
            frames.append(conf.coordset)
            energies.append(conf.energy)
        return frames, energies, nat

    def nframes(self):
        return len(self._frames)

    def atomtypes(self):
        if not self._frames:
            print("Make sure to load a trajectory first!")
        else:
            atomtypes = set()
            for at in self._frames[0]:
                atomtypes.add(at.split('_')[0])
            return atomtypes

    def coordinates(self, frame):
        print(F" {self.nat}\n {' '.join(self.energies[frame].split())}")
        for at in self._frames[frame]:
            print(F"{at.split('_')[0].ljust(3)} {self._frames[frame][at][0]} {self._frames[frame][at][1]} {self._frames[frame][at][2]}")#   {:.3f}")

#    def __del__(self):
#        self.handle.close()

class Frame():
    def __init__(self, conf, en = False):
        self.coordset = self.coordset(conf[2:])
        if en == True:
            self.energy = conf[1]
        else:
            self.energy = 0

    def coordset(self, coords):
        coordset = {}
        for at in coords:
            count = 1
            atom = at.split()
            for label in coordset:
                if atom[0] + "_" in label:
                    count += 1
            coordset[atom[0] + "_" + str(count)] = np.array([float(coord) for coord in atom[1:]])

        return coordset
