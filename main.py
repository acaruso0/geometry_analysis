import numpy as np
import geometry_analysis as geom
from geometry_analysis import Atom
import math

if __name__ == "__main__":
    periodic_table = geom.loaddb()

    testXYZ = geom.History("history", periodic_table)
#    testXYZ = geom.XYZfile("p_mbtr_5k.xyz", periodic_table, (19.37, 19.37, 22.34))
#    for n in range(834):
#        testXYZ.coordprint(n, save=True)
#    geom.rdf(testXYZ, 'Cs', 'O', 0, 10, 300)
