import geometry_analysis as geom

if __name__ == "__main__":

#    testXYZ = geom.History("history")
    testXYZ = geom.XYZfile("./data/simulation.pos_0.xyz", (12.42855, 12.42855, 12.42855))
#    testXYZ = geom.XYZfile("test.xyz")
#    print(testXYZ.frames)
#    for n in range(834):
#        testXYZ.coordprint(n, save=True)
    geom.rdf(testXYZ, 'O', 'O', 0, 10, 300)
