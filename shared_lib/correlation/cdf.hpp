#include <vector>
#include <Eigen/Dense>

std::vector<std::vector<unsigned int>> triples(std::vector<int>, int, int, int);
std::vector<std::vector<double>> cdf_distang(std::vector<std::vector<Eigen::VectorXd>>, std::vector<double>, std::vector<std::vector<unsigned int>>);
