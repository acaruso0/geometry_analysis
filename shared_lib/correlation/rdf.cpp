#include <vector>
#include <Eigen/Dense>
#include <cmath>

//using namespace std;

//typedef Eigen::VectorXd eigvec;

std::vector<std::vector<unsigned int>> couples(std::vector<int> atoms, int atom1, int atom2) {
    std::vector<std::vector<unsigned int>> couples;
    for(unsigned int i = 0; i < atoms.size(); i++) {
        if(atoms.at(i) == atom1) {
            for(unsigned int j = 0; j < atoms.size(); j++) {
                if(atoms.at(j) == atom2 && j != i) {
                    std::vector<unsigned int> couple{i, j};
                    couples.push_back(couple);
                }
            }
        }
    }
    return couples;
}

std::vector<double> rdf_dists(std::vector<std::vector<Eigen::VectorXd>> frames,
                std::vector<double> pbc,
                std::vector<std::vector<unsigned int>> couples) {
    std::vector<double> dists;
    for(std::vector<Eigen::VectorXd> frame : frames) {
        for(std::vector<unsigned int> couple : couples) {
            double dx = frame.at(couple.at(0))(0) - frame.at(couple.at(1))(0);
            double dy = frame.at(couple.at(0))(1) - frame.at(couple.at(1))(1);
            double dz = frame.at(couple.at(0))(2) - frame.at(couple.at(1))(2);

            dx = dx - std::round(dx / pbc.at(0)) * pbc.at(0);
            dy = dy - std::round(dy / pbc.at(1)) * pbc.at(1);
            dz = dz - std::round(dz / pbc.at(2)) * pbc.at(2);

            double delta = std::pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
            double d_ij = std::sqrt(delta);

            dists.push_back(d_ij);
        }
    }
    return dists;
}
