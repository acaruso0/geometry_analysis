#include <vector>
#include <Eigen/Dense>

std::vector<std::vector<unsigned int>> couples(std::vector<int>, int, int);
std::vector<double> rdf_dists(std::vector<std::vector<Eigen::VectorXd>>, std::vector<double>, std::vector<std::vector<unsigned int>>);
