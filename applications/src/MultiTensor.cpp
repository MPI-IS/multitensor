#include <iostream>
#include <string>

#include "app_utils.h"
#include "multitensor/main.hpp"

int main(int argc, char *argv[])
{
    using namespace multitensor;

    // Input adjacency file
    std::string adjacency_filename = "adjacency.dat";

    // Output U file
    std::string u_output_filename = "u_out.dat";

    // Output V file
    std::string v_output_filename = "v_out.dat";

    // Output affinity file
    std::string affinity_output_filename = "w_out.dat";

    // Number of groups
    size_t nof_groups;

    // Number of realizations
    size_t nof_realizations = 1;

    // Maximum number of iterations
    size_t max_nof_iterations = 500;

    // Number of positive checks required for reaching convergence
    size_t nof_convergences = 10;

    // Help option
    if (cmd_option_exists(argv, argv + argc, "--help"))
    {
        std::cout
            << "\nNAME\n"
            << "\tMultitensor factorization algorithm\n\n"
            << "DESCRIPTION\n"
            << "\tCommand line interface to run the Multitensor factorization algorithm\n\n"
            << "SYNOPSIS\n"
            << "\t./Multitensor --k <nof_groups> --a <input-adjacency-file> \n\n"
            << "OPTIONS\n"
            << "\t--k <nof_groups>\n"
            << "\t\t Number of groups (required)\n\n"
            << "\t--a <input-adjacency-file>\n"
            << "\t\t Adjacency file (default: " << adjacency_filename << ")\n\n"
            << "\t--r <nof_realizations>\n"
            << "\t\t Number of realizations (default : " << nof_realizations << ")\n\n"
            << "\t--maxit <max_nof_iterations>\n"
            << "\t\t Maximum number of iterations (default : " << max_nof_iterations << ")\n\n"
            << "\t--y <nof_convergences>\n"
            << "\t\t Number of positive checks required for reaching convergence (default : "
            << nof_convergences << ")\n\n"
            << "\t--wo <affinity-output-file>\n"
            << "\t\t Affinity output file (default: " << affinity_output_filename << ")\n\n"
            << "\t--uo <u-output-file>\n"
            << "\t\t Output file for outgoing vertices (default: " << u_output_filename << ")\n\n"
            << "\t--vo <v-output-file>\n"
            << "\t\t Output file for incoming vertices (default: " << v_output_filename << ")\n\n"
            << std::endl;
        return 0;
    }

    // Version option
    if (cmd_option_exists(argv, argv + argc, "--version"))
    {
        std::cout << MULTITENSOR_VERSION << std::endl;
        return 0;
    }

    // Read in command line arguments
    if (cmd_option_exists(argv, argv + argc, "--a"))
    {
        adjacency_filename = get_cmd_option(argv, argv + argc, "--a");
    }
    if (cmd_option_exists(argv, argv + argc, "--o"))
    {
        affinity_output_filename = get_cmd_option(argv, argv + argc, "--o");
    }
    if (cmd_option_exists(argv, argv + argc, "--k"))
    {
        nof_groups = std::stoi(get_cmd_option(argv, argv + argc, "--k"));
    }
    else
    {
        throw std::runtime_error("Please specify the number of groups (--k option).");
    }

    // Read-in data
    std::vector<size_t> edges_in, edges_out;
    std::vector<unsigned int> edges_weight;
    size_t nof_nodes, nof_layers;

    read_adjacency_data(adjacency_filename, edges_in, edges_out, edges_weight, nof_nodes, nof_layers);

    // Call algorithm
    multitensor_factorization(edges_in, edges_out, edges_weight,
                              affinity_output_filename, u_output_filename, v_output_filename,
                              nof_nodes, nof_layers, nof_groups,
                              nof_realizations, max_nof_iterations, nof_convergences);
}