import math
import numpy as np
from database import Atom, get_name
import geometry_analysis as geom

import matplotlib.pyplot as plt

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

            dx = dx - round(dx / xyzfile.pbc[0]) * xyzfile.pbc[0]
            dy = dy - round(dy / xyzfile.pbc[1]) * xyzfile.pbc[1]
            dz = dz - round(dz / xyzfile.pbc[2]) * xyzfile.pbc[2]

            ############################

            delta = dx**2 + dy**2 + dz**2
            d_ij = math.sqrt(delta)

            dists.append(d_ij)

    rdf = list(np.histogram(dists, nbins, range=(start, stop)))

    #################
    # Normalization #
    #################

    rho = at2_cnt / (xyzfile.pbc[0] * xyzfile.pbc[1] * xyzfile.pbc[2])
    norm = 4 * math.pi * rho * dr * xyzfile.nframes * at1_cnt
    
    if start == 0:
        rdf[0] = list(map(lambda x, y: x / (norm * y**2), rdf[0][1:], rdf[1][1:-1]))
    else:
        rdf[0] = list(map(lambda x, y: x / (norm * y**2), rdf[0], rdf[1][:-1]))


    #################

    with open('out_rdf.dat', 'w+') as outfile:
        for n in range(len(rdf[0])):
            outfile.write(F'{rdf[1][n]} {rdf[0][n]}\n')

#def cdf(xyzfile, atom1, atom2, r_range, th_range):
def cdf(xyzfile, atom1, atom2, atom3, r_range, th_range):
    atom1 = get_name(atom1, xyzfile._periodic_table)
    atom2 = get_name(atom2, xyzfile._periodic_table)
    atom3 = get_name(atom3, xyzfile._periodic_table)

    at1_cnt = xyzfile.atcount[atom1]
    at2_cnt = xyzfile.atcount[atom2]
    at3_cnt = xyzfile.atcount[atom3]
    
    atom1 = xyzfile._periodic_table[atom1].atnum
    atom2 = xyzfile._periodic_table[atom2].atnum
    atom3 = xyzfile._periodic_table[atom3].atnum
    
    triples = []
    for at1_n, at1 in enumerate(xyzfile.atoms):
        if at1 == atom1:
            for at2_n, at2 in enumerate(xyzfile.atoms):
                if at2 == atom2 and at2_n != at1_n:
#                    for at3_n, at3 in enumerate(xyzfile.atoms):
#                        if at3 == atom2 and at3_n != at1_n and at3_n != at2_n:
                    for at3_n, at3 in enumerate(xyzfile.atoms):
                        if at3 == atom3 and at3_n != at1_n and at3_n != at2_n:
                            triples.append((at1_n, at2_n, at3_n))

    print("Triples generated.")
    print(triples)
    dists = []
    angles = []
    for frame in xyzfile.frames:
        for triple in triples:
            dx12 = frame[triple[0]][0] - frame[triple[1]][0]
            dy12 = frame[triple[0]][1] - frame[triple[1]][1]
            dz12 = frame[triple[0]][2] - frame[triple[1]][2]

            dx13 = frame[triple[0]][0] - frame[triple[2]][0]
            dy13 = frame[triple[0]][1] - frame[triple[2]][1]
            dz13 = frame[triple[0]][2] - frame[triple[2]][2]

            dx23 = frame[triple[1]][0] - frame[triple[2]][0]
            dy23 = frame[triple[1]][1] - frame[triple[2]][1]
            dz23 = frame[triple[1]][2] - frame[triple[2]][2]

            ############################
            # Minimum Image Convention #
            ############################

            dx13 = dx13 - round(dx12 / xyzfile.pbc[0]) * xyzfile.pbc[0]
            dy13 = dy13 - round(dy12 / xyzfile.pbc[1]) * xyzfile.pbc[1]
            dz13 = dz13 - round(dz12 / xyzfile.pbc[2]) * xyzfile.pbc[2]

            dx12 = dx12 - round(dx12 / xyzfile.pbc[0]) * xyzfile.pbc[0]
            dy12 = dy12 - round(dy12 / xyzfile.pbc[1]) * xyzfile.pbc[1]
            dz12 = dz12 - round(dz12 / xyzfile.pbc[2]) * xyzfile.pbc[2]

            ############################

            delta12 = dx12**2 + dy12**2 + dz12**2
            delta13 = dx13**2 + dy13**2 + dz13**2
            delta23 = dx23**2 + dy23**2 + dz23**2

            d_12 = math.sqrt(delta12)
            d_13 = math.sqrt(delta13)
            d_23 = math.sqrt(delta23)

            dists.append(d_12)
            print(d_13, d_12, d_23)
            angles.append(geom.cosrule(d_23, d_12, d_13, deg = True))

    print("Distances and angles calculated.")

    print(len(angles))
    print(angles)

    cdf, dist, angle = np.histogram2d(dists, angles, bins=(r_range[2], th_range[2]), range=(r_range[:2], th_range[:2]))

#    with open('before_norm.dat', 'w+') as outfile:
#        for row in range(len(cdf)):
#            for col in range(len(cdf[0])):
#                outfile.write(F'{dist[row]} {angle[col]} {cdf[row][col]}\n')

    #################
    # Normalization #
    #################

    def volume(r_range, th_range):
        th_range = list(map(lambda x: x * (math.pi / 180), th_range))
        volume = (2 / 3) * math.pi * (pow(r_range[1], 3) - pow(r_range[0], 3)) * (math.cos(th_range[0]) - math.cos(th_range[1]))
        return volume
    
    rho = at2_cnt / (xyzfile.pbc[0] * xyzfile.pbc[1] * xyzfile.pbc[2])

    norm = 2 * xyzfile.nframes * at1_cnt * rho

#    print(norm, xyzfile.nframes, at1_cnt, at2_cnt, rho, xyzfile.pbc)
    plt_dist = []
    plt_angle = []
    plt_cdf = []
    for row in range(len(cdf)):
        for col in range(len(cdf[0])):
#            cdf[row][col] /= (norm*volume((dist[row], dist[row + 1]), (angle[col], angle[col + 1])))
            plt_dist.append(dist[row])
            plt_angle.append(angle[col])
            plt_cdf.append(cdf[row][col])
    print("Normalized.")
    #################

    print(sum(plt_cdf))

#    cdf = cdf.T
    fig, ax = plt.subplots()
    ax.scatter(plt_dist, plt_angle, c=plt_cdf)#, s=30)
#    ax.imshow(cdf, interpolation='nearest', origin='low', extent=[dist[0], dist[-1], angle[0], angle[-1]])
    plt.show()

#    with open('out_cdf.dat', 'w+') as outfile:
#        for row in range(len(cdf)):
#            for col in range(len(cdf[0])):
#                outfile.write(F'{dist[row]} {angle[col]} {cdf[row][col]}\n')
