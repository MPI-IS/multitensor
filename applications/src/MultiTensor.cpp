#include <iostream>
#include <string>
#include <boost/program_options.hpp>

#include "app_utils.h"
#include "multitensor.hpp"

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
    using namespace multitensor;

    auto parser = parse_command_line(argc, argv);
    po::variables_map vm = std::get<0>(parser);
    po::options_description desc = std::get<1>(parser);

    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 0;
    }

    // Read-in data
    std::vector<size_t> edges_in, edges_out;
    std::vector<unsigned int> edges_weight;
    size_t nof_nodes, nof_layers;

    read_adjacency_data(edges_in, edges_out, edges_weight, nof_nodes, nof_layers,
                        vm["a"].as<std::string>());

    // Call algorithm
    multitensor_algo(edges_in, edges_out, edges_weight,
                     nof_nodes, nof_layers,
                     static_cast<size_t>(vm["k"].as<int>()),
                     static_cast<unsigned int>(vm["r"].as<int>()),
                     static_cast<unsigned int>(vm["maxit"].as<int>()),
                     static_cast<unsigned int>(vm["y"].as<int>()));
}