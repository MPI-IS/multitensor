#include <iostream>
#include <string>
#include <ctime>
#include <boost/graph/adjacency_list.hpp>

#include "app_utils.hpp"
#include "multitensor/main.hpp"
#include "multitensor/utils.hpp"

int main(int argc, char *argv[])
{
    using namespace multitensor;

    // Input adjacency file
    std::string adjacency_filename = "adjacency.dat";

    // Output directory
    std::string output_directory = "results";

    // Graph type used
    std::string graph_type = "directed";

    // Number of groups
    size_t nof_groups = 0;

    // Number of realizations
    size_t nof_realizations = 1;

    // Maximum number of iterations
    size_t max_nof_iterations = 500;

    // Number of positive checks required for reaching convergence
    size_t nof_convergences = 10;

    // Seed for the random generator
    std::string seed = "random";

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
            << "\t--undirected\n"
            << "\t\t Use undirected graphs (instead of the default directed graphs)\n\n"
            << "\t--r <nof_realizations>\n"
            << "\t\t Number of realizations (default : " << nof_realizations << ")\n\n"
            << "\t--maxit <max_nof_iterations>\n"
            << "\t\t Maximum number of iterations (default : " << max_nof_iterations << ")\n\n"
            << "\t--y <nof_convergences>\n"
            << "\t\t Number of positive checks required for reaching convergence (default : "
            << nof_convergences << ")\n\n"
            << "\t--s <seed>\n"
            << "\t\t Seed for the random generator. Options are 'random' (default) or an integer\n\n"
            << "\t--o <output-dir>\n"
            << "\t\t Output directory (default: '" << output_directory << "')\n\n"
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
    if (cmd_option_exists(argv, argv + argc, "--k"))
    {
        nof_groups = std::stoi(get_cmd_option(argv, argv + argc, "--k"));
    }
    else
    {
        throw std::runtime_error("Please specify the number of groups (--k option).");
    }
    if (cmd_option_exists(argv, argv + argc, "--a"))
    {
        adjacency_filename = get_cmd_option(argv, argv + argc, "--a");
    }
    if (cmd_option_exists(argv, argv + argc, "--undirected"))
    {
        graph_type = "undirected";
    }
    if (cmd_option_exists(argv, argv + argc, "--o"))
    {
        output_directory = get_cmd_option(argv, argv + argc, "--o");
    }
    if (cmd_option_exists(argv, argv + argc, "--r"))
    {
        nof_realizations = std::stoi(get_cmd_option(argv, argv + argc, "--r"));
    }
    if (cmd_option_exists(argv, argv + argc, "--s"))
    {
        seed = get_cmd_option(argv, argv + argc, "--s");
    }

    // Read-in data
    std::vector<size_t> edges_start, edges_end, edges_weight;
    read_adjacency_data(adjacency_filename, edges_start, edges_end, edges_weight);
    const size_t nof_vertices = utils::get_num_vertices(edges_start, edges_end);
    const size_t nof_layers = edges_weight.size() / edges_start.size();

    // Prepare output data
    tensor::Tensor<double> w(nof_groups, nof_groups, nof_layers), u(nof_vertices, nof_groups), v;
    std::vector<size_t> labels;

    // Create random generator
    std::time_t seed_value;
    if (seed == "random")
    {
        seed_value = std::time(nullptr);
    }
    else
    {
        seed_value = std::stoi(seed);
    }
    utils::RandomGenerator random_generator{seed_value};

    // Call algorithm
    utils::Report results;
    if (graph_type == "directed")
    {
        // We need v
        v.resize(nof_vertices, nof_groups);
        results = multitensor_factorization(
            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, w, random_generator);
    }
    else
    {
        results = multitensor_factorization<boost::undirectedS>(
            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, w, random_generator);
    }

    // Write output files
    boost::filesystem::create_directory(output_directory);
    auto dpath = boost::filesystem::canonical(output_directory);
    std::cout << "Writing output files in directory " << dpath << std::endl;
    write_affinity_file(dpath / WOUT_FILENAME, w, results);
    write_membership_file(dpath / UOUT_FILENAME, labels, u, results);
    write_membership_file(dpath / VOUT_FILENAME, labels, v, results);
    write_info_file(dpath / INFO_FILENAME, results);

    return 0;
}