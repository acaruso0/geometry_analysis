import numpy as np
import geometry_analysis as geom
from geometry_analysis import Atom

if __name__ == "__main__":
    periodic_table = geom.loaddb()

#    print(periodic_table["oxygen"].label)
#    print(db.label_to_name("Cl", periodic_table))

#    testXYZ = geom.XYZfile("input.xyz", periodic_table, 1.5)
    testXYZ = geom.XYZfile("simulation.pos_0.xyz", periodic_table, 12.42855)
#    print(testXYZ.atoms)
#    print(np.asarray(testXYZ.frames[0]))
#    print(testXYZ.atcount)
#    testXYZ.coordprint(1)#, save=True)
    geom.rdf(testXYZ, 'O', 'O', 0, 15, 1000)
