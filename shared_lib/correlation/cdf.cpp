#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

std::vector<std::vector<unsigned int>> triples(std::vector<int> atoms, int atom1, int atom2, int atom3) {
    std::vector<std::vector<unsigned int>> triples;
    for(unsigned int i = 0; i < atoms.size(); i++) {
        if(atoms.at(i) == atom1) {
            for(unsigned int j = 0; j < atoms.size(); j++) {
                if(atoms.at(j) == atom2 && j != i) {
                    for(unsigned int k = 0; k < atoms.size(); k++) {
                        if(atoms.at(k) == atom3 && k != i && k != j) {
                            std::vector<unsigned int> triple{i, j, k};
                            triples.push_back(triple);
                        }
                    }
                }
            }
        }
    }
    return triples;
}

std::vector<std::vector<double>> cdf_distang(std::vector<std::vector<Eigen::VectorXd>> frames,
                std::vector<double> pbc,
                std::vector<std::vector<unsigned int>> triples) {

    std::vector<std::vector<double>> distang;
    std::vector<double> dists, angs;
    const double PI = 3.141592653589793;
    
    for(std::vector<Eigen::VectorXd> frame : frames) {
        for(std::vector<unsigned int> triple : triples) {
            double dx12 = frame.at(triple.at(0))(0) - frame.at(triple.at(1))(0);
            double dy12 = frame.at(triple.at(0))(1) - frame.at(triple.at(1))(1);
            double dz12 = frame.at(triple.at(0))(2) - frame.at(triple.at(1))(2);

            double dx13 = frame.at(triple.at(0))(0) - frame.at(triple.at(2))(0);
            double dy13 = frame.at(triple.at(0))(1) - frame.at(triple.at(2))(1);
            double dz13 = frame.at(triple.at(0))(2) - frame.at(triple.at(2))(2);
            
            double dx23 = frame.at(triple.at(1))(0) - frame.at(triple.at(2))(0);
            double dy23 = frame.at(triple.at(1))(1) - frame.at(triple.at(2))(1);
            double dz23 = frame.at(triple.at(1))(2) - frame.at(triple.at(2))(2);

            // ######################## //
            // Minimum Image Convention //
            // ######################## //
            
            dx13 = dx13 - std::round(dx13 / pbc.at(0)) * pbc.at(0);
            dy13 = dy13 - std::round(dy13 / pbc.at(1)) * pbc.at(1);
            dz13 = dz13 - std::round(dz13 / pbc.at(2)) * pbc.at(2);

            dx12 = dx12 - std::round(dx12 / pbc.at(0)) * pbc.at(0);
            dy12 = dy12 - std::round(dy12 / pbc.at(1)) * pbc.at(1);
            dz12 = dz12 - std::round(dz12 / pbc.at(2)) * pbc.at(2);
            
            dx23 = dx23 - std::round(dx23 / pbc.at(0)) * pbc.at(0);
            dy23 = dy23 - std::round(dy23 / pbc.at(1)) * pbc.at(1);
            dz23 = dz23 - std::round(dz23 / pbc.at(2)) * pbc.at(2);
            
            // ######################## //
            
            double delta12 = std::pow(dx12, 2) + pow(dy12, 2) + pow(dz12, 2);
            double delta13 = std::pow(dx13, 2) + pow(dy13, 2) + pow(dz13, 2);
            double delta23 = std::pow(dx23, 2) + pow(dy23, 2) + pow(dz23, 2);
            
            double d_12 = std::sqrt(delta12);
            double d_13 = std::sqrt(delta13);
            double d_23 = std::sqrt(delta23);

            dists.push_back(d_12);
            angs.push_back(std::acos((std::pow(d_12, 2) + std::pow(d_13, 2) - std::pow(d_23, 2)) / (2 * d_12 * d_13)) * (180 / PI));
        }
    }
    distang.push_back(dists);
    distang.push_back(angs);
    
    return distang;
}
