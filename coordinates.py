import numpy as np

class XYZfile():
    def __init__(self, filename):
        self.frames, self.energies, self.nat = self.loadall(filename)
        self.nframes = self.nframes()
        self.atomtypes = self.atomtypes()
        self.atcount = self.atcount()

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
        return len(self.frames)

    def atcount(self):
        atcount = {}
        for atomtype in self.atomtypes:
            it = 0
            for at in self.frames[0]:
                if atomtype == at.split('_')[0]:
                    it += 1
            atcount[atomtype] = it
        return atcount

    def atomtypes(self):
        if not self.frames:
            print("Make sure to load a trajectory first!")
        else:
            atomtypes = set()
            for at in self.frames[0]:
                atomtypes.add(at.split('_')[0])
            return atomtypes

    def coordprint(self, frame, label = False):
        print(F" {self.nat}\n {' '.join(self.energies[frame].split())}")
        for at in self.frames[frame]:
            if label == False:
                print(F"{at.split('_')[0].ljust(3)} {self.frames[frame][at][0]} {self.frames[frame][at][1]} {self.frames[frame][at][2]}")#   {:.3f}")
            else:
                print(F"{at.split('_')[0]}.{at.split('_')[1]} {self.frames[frame][at][0]} {self.frames[frame][at][1]} {self.frames[frame][at][2]}")

    def getatom(self, in_frame, end_frame, *args):
        arglist, sliced = [], []
        if args:
            for arg in args:
                arg = arg.split('.')
                if ':' in arg[1]:
                    argrange = range(int(arg[1].split(':')[0]), int(arg[1].split(':')[1]) + 1)
                arglist += [arg[0] + "_" + str(n) for n in list(argrange)]
        else:
            arglist = list(self.frames[0].keys())#[in_frame:end_frame+1][0].keys())

        for frame in range(in_frame, end_frame + 1):
            sl_atoms = {}
            for atlabel in arglist:
                sl_atoms[atlabel] = self.frames[frame][atlabel]
            sliced.append(sl_atoms)
            
        return sliced

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
