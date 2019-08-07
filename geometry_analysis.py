import numpy as np
import coordinates as crd

if __name__ == "__main__":
    testXYZ = crd.XYZfile("input2.xyz")

#    print(testXYZ._frames)
#    print(testXYZ.atomtypes)
#    print(testXYZ.nframes)
    for n in range(testXYZ.nframes):
        testXYZ.energies[n] = F'{testXYZ.energies[n]} {testXYZ.energies[n]} 0.0 0.0'
        testXYZ.coordinates(n)
