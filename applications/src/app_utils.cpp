#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "app_utils.h"

namespace po = boost::program_options;

std::tuple<po::variables_map, po::options_description> parse_command_line(int ac, char **av)
{
    po::options_description desc(
        "\nNAME\n\tMultitensor factorization algorithm\n\nUSAGE:\n\t" + std::string(av[0]) + " <option> where <option> is one or more of");
    desc.add_options()("help", "");
    desc.add_options()("a", po::value<std::string>()->default_value("adjacency.dat"), "filename for adjacency matrix");
    desc.add_options()("l", po::value<int>()->default_value(4), "number of layers");
    desc.add_options()("k", po::value<int>()->default_value(5), "number of groups");
    desc.add_options()("r", po::value<int>()->default_value(1), "number of realizations");
    desc.add_options()("maxit", po::value<int>()->default_value(500), "maximum number of iterations");
    desc.add_options()("y", po::value<int>()->default_value(10), "number of positive checks required for reaching convergence");
    desc.add_options()("with-w-init", po::value<bool>()->default_value(false), "whether to initialize the affinity matrix");
    desc.add_options()("w", po::value<std::string>()->default_value("w.dat"), "filename for initiliazing the affinity matrix");
    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);
    return std::make_tuple(vm, desc);
}