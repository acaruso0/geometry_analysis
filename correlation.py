import math
import numpy as np
from database import Atom, get_name

def rdf(xyzfile, atom1, atom2, start, stop, nbins):
    assert start >= 0, "The starting value of the rdf has to be positive."
    assert stop > start, "The final value of the rdf has to be greater than the starting one."

    #! ENTER PBC CONDITION, THROW ERROR

    atom1 = get_name(atom1, xyzfile._periodic_table)
    atom2 = get_name(atom2, xyzfile._periodic_table)

    at1_cnt = xyzfile.atcount[atom1]
    at2_cnt = xyzfile.atcount[atom2]
    
    atom1 = xyzfile._periodic_table[atom1].atnum
    atom2 = xyzfile._periodic_table[atom2].atnum

    dr = (stop - start) / nbins

    couples = []
    for at1_n, at1 in enumerate(xyzfile.atoms):
        if at1 == atom1:
            for at2_n, at2 in enumerate(xyzfile.atoms):
                if at2 == atom2 and at2_n != at1_n: #this last bool checks for rdf with itself
                    couples.append((at1_n, at2_n))

    dists = []
    for frame in xyzfile.frames:
        for couple in couples:
            dx = frame[couple[0]][0] - frame[couple[1]][0]
            dy = frame[couple[0]][1] - frame[couple[1]][1]
            dz = frame[couple[0]][2] - frame[couple[1]][2]

            ############################
            # Minimum Image Convention #
            ############################

            dx = dx - round(dx / xyzfile.pbc) * xyzfile.pbc
            dy = dy - round(dy / xyzfile.pbc) * xyzfile.pbc
            dz = dz - round(dz / xyzfile.pbc) * xyzfile.pbc

            ############################

            delta = dx**2 + dy**2 + dz**2
            d_ij = math.sqrt(delta)

            dists.append(d_ij)

    rdf = list(np.histogram(dists, nbins, range=(start, stop)))

    #################
    # Normalization #
    #################

    rho = at2_cnt / pow(xyzfile.pbc, 3)
    norm = 4 * math.pi * rho * dr * xyzfile.nframes * at1_cnt
    
    if start == 0:
        rdf[0] = list(map(lambda x, y: x / (norm * y**2), rdf[0][1:], rdf[1][1:-1]))
    else:
        rdf[0] = list(map(lambda x, y: x / (norm * y**2), rdf[0], rdf[1][:-1]))


    #################

    with open('output.dat', 'w+') as outfile:
        for n in range(len(rdf[0])):
            outfile.write(F'{rdf[1][n]} {rdf[0][n]}\n')
