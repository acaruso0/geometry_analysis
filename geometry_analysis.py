import numpy as np
import coordinates as crd

if __name__ == "__main__":
    testXYZ = crd.XYZfile("input2.xyz")

#    print(testXYZ.frames)
#    print(testXYZ.atomtypes)
#    print(testXYZ.nframes)
#    for n in range(testXYZ.nframes):
#        testXYZ.energies[n] = F'{testXYZ.energies[n]} {testXYZ.energies[n]} 0.0 0.0'
#        testXYZ.coordprint(n)
#        testXYZ.coordprint(n, label = True)
#    print(testXYZ.atcount)
#    print(testXYZ.getatom([1], 'O.2'))
    print(testXYZ.getatom(2, 2, 'H.1:2'))
    print(testXYZ.getatom(1, 2))
