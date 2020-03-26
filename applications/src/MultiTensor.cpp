#include <iostream>
#include <string>
#include <ctime>
#include <boost/graph/adjacency_list.hpp>

#include "app_utils.hpp"
#include "multitensor/main.hpp"
#include "multitensor/utils.hpp"

int main(int argc, char *argv[])
{
    using namespace boost;
    using namespace multitensor::initialization;
    using namespace multitensor::tensor;

    // Input adjacency file
    std::string adjacency_filename = "adjacency.dat";

    // Input affinity file
    std::string affinity_filename = "";

    // Output directory
    std::string output_directory = "results";

    // Graph type used
    bool directed_graph = true;

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

    // Option for assortative groups
    bool assortative = false;

    // Help option
    if (cmd_option_exists(argv, argv + argc, "--help"))
    {
        std::cout
            << "\nNAME\n"
            << "\tMultitensor factorization algorithm\n\n"
            << "DESCRIPTION\n"
            << "\tCommand line interface to run the Multitensor factorization algorithm\n\n"
            << "SYNOPSIS\n"
            << "\t./Multitensor --k <nof_groups> [options]\n\n"
            << "OPTIONS\n"
            << "\t--k <nof_groups>\n"
            << "\t\t Number of groups (required)\n\n"
            << "\t--a <input-adjacency-file>\n"
            << "\t\t Adjacency file (default: " << adjacency_filename << ")\n\n"
            << "\t--w <input-affinity-file>\n"
            << "\t\t Affinity file for initialization (optional)\n\n"
            << "\t--assortative\n"
            << "\t\t Use assortative model\n\n"
            << "\t--undirected\n"
            << "\t\t Use undirected graphs (instead of the default directed graphs)\n\n"
            << "\t--r <nof_realizations>\n"
            << "\t\t Number of realizations (default: " << nof_realizations << ")\n\n"
            << "\t--maxit <max_nof_iterations>\n"
            << "\t\t Maximum number of iterations (default: " << max_nof_iterations << ")\n\n"
            << "\t--y <nof_convergences>\n"
            << "\t\t Number of positive checks required for reaching convergence (default: "
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
    if (cmd_option_exists(argv, argv + argc, "--w"))
    {
        affinity_filename = get_cmd_option(argv, argv + argc, "--w");
    }
    if (cmd_option_exists(argv, argv + argc, "--undirected"))
    {
        directed_graph = false;
    }
    if (cmd_option_exists(argv, argv + argc, "--assortative"))
    {
        assortative = true;
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
    if (cmd_option_exists(argv, argv + argc, "--maxit"))
    {
        max_nof_iterations = std::stoi(get_cmd_option(argv, argv + argc, "--maxit"));
    }

    // Read-in data
    // Adjacency data
    std::vector<size_t> edges_start, edges_end, edges_weight;
    read_adjacency_data(adjacency_filename, edges_start, edges_end, edges_weight);
    const size_t nof_vertices = utils::get_num_vertices(edges_start, edges_end);
    const size_t nof_layers = edges_weight.size() / edges_start.size();

    // Affinity tensor
    std::vector<double> affinity;
    // Different size expected if assortative or not
    const size_t affinity_size = assortative ? nof_groups * nof_layers
                                             : nof_groups * nof_groups * nof_layers;
    affinity.resize(affinity_size);

    // Read-in file if requested
    const bool w_init_defined = (affinity_filename != "");
    if (w_init_defined)
    {
        read_affinity_data(affinity_filename, assortative, affinity);
    }

    // Prepare output data
    // We do not set the size of v right away - we might not need it if we use an undirected network
    Matrix<double> u(nof_vertices, nof_groups), v;
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

    // Merge the booleans in a single integer
    const size_t selection = directed_graph + 2 * assortative + 4 * w_init_defined;

    // Call algorithm
    utils::Report results;

    switch (selection)
    {
    case 0:
        // Undirected + non-assortative + w random
        results = multitensor_factorization<undirectedS>(
            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity, random_generator);
        break;
    case 1:
        // Directed + non-assortative + w random
        v.resize(nof_vertices, nof_groups); // we need v
        results = multitensor_factorization(
            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity, random_generator);
        break;
    case 2:
        // Undirected + assortative + w random
        results = multitensor_factorization<undirectedS, DiagonalTensor<double>>(
            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity, random_generator);
        break;
    case 3:
        // Directed + assortative + w random
        v.resize(nof_vertices, nof_groups); // we need v
        results = multitensor_factorization<bidirectionalS, DiagonalTensor<double>>(
            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity, random_generator);
        break;
    case 4:
        // Undirected + non-assortative + w from file
        results = multitensor_factorization<undirectedS, SymmetricTensor<double>,
                                            init_symmetric_tensor_from_initial<SymmetricTensor<double>>>(
            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity, random_generator);
        break;
    case 5:
        // Directed + non-assortative + w from file
        v.resize(nof_vertices, nof_groups); // we need v
        results = multitensor_factorization<bidirectionalS, SymmetricTensor<double>,
                                            init_symmetric_tensor_from_initial<SymmetricTensor<double>>>(
            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity, random_generator);
        break;
    case 6:
        // Undirected + assortative + w from file
        results = multitensor_factorization<undirectedS, DiagonalTensor<double>,
                                            init_symmetric_tensor_from_initial<DiagonalTensor<double>>>(
            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity, random_generator);
        break;
    case 7:
        // Directed + assortative + w from file
        v.resize(nof_vertices, nof_groups); // we need v
        results = multitensor_factorization<bidirectionalS, DiagonalTensor<double>,
                                            init_symmetric_tensor_from_initial<DiagonalTensor<double>>>(
            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity, random_generator);
        break;
    default:
        throw std::runtime_error(
            "[multitensor] Something went wrong with the algorithm type. It should be >= 0 and < 8, intead got " +
            std::to_string(selection) + "\n");
    }

    // Write output files
    filesystem::create_directory(output_directory);
    const auto dpath = filesystem::canonical(output_directory);
    std::cout << "Writing output files in directory " << dpath << std::endl;
    write_info_file(dpath / INFO_FILENAME, results);
    write_affinity_file(dpath / WOUT_FILENAME, affinity, results, nof_groups, nof_layers);
    write_membership_file(dpath / UOUT_FILENAME, labels, u, results);
    if (directed_graph)
    {
        write_membership_file(dpath / VOUT_FILENAME, labels, v, results);
    }

    return 0;
}