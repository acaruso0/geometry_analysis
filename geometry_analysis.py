import numpy as np
import coordinates as crd

if __name__ == "__main__":
    testXYZ = crd.XYZfile("input.xyz")

#    print(testXYZ._frames)
#    print(testXYZ.atomtypes)
    print(testXYZ.nframes)
    for n in range(testXYZ.nframes):
        testXYZ.coordinates(n)
