// Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
// Distributed under the GNU GPL license version 3
// See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md

/*!
 * @file
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 * @author Caterina De Bacco (caterina.debacco@tuebingen.mpg.de)
 */

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "app_utils.hpp"

char *
get_cmd_option(char **begin, char **end, const std::string &option)
{
    char **itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return nullptr;
}

bool cmd_option_exists(char **begin, char **end, const std::string &option)
{
    return std::find(begin, end, option) != end;
}

void read_affinity_data(const boost::filesystem::path &filename,
                        const bool &assortative,
                        std::vector<double> &w)
{
    std::ifstream in(filename.string());
    if (in.fail())
    {
        throw std::runtime_error(
            std::string("In read_affinity_data, failed to open ") + filename.string());
    }

    std::cout << "Reading affinity file " << filename << std::endl;

    std::string line;

    // First parse the file to get dimensions
    size_t nof_groups(0), nof_layers(0);
    std::string tok;
    double value;
    while (!in.eof())
    {
        std::getline(in, line);
        // skip over empty lines
        if (line.size() == 0)
        {
            continue;
        }

        // Remove trailing whitespaces
        line.erase(line.find_last_not_of(" ") + 1);

        // Count number of groups and layers
        std::istringstream is(line);
        size_t current_nof_groups(0);

        // First character - could be # or layer id
        is >> tok;
        if (tok == "#")
        {
            continue;
        }

        // Groups values
        while (is >> value)
        {
            current_nof_groups++;
        }
        if (nof_groups == 0)
        {
            nof_groups = current_nof_groups;
        }
        assert(current_nof_groups = nof_groups);
        nof_layers++;
    }

    // Now build vector...
    in.clear();
    in.seekg(0);

    while (!in.eof())
    {
        std::getline(in, line);
        if (line.size() == 0)
            continue; // skip over empty lines and comments

        // Remove trailing whitespaces
        line.erase(line.find_last_not_of(" ") + 1);

        std::istringstream is(line);

        // Layer ID
        size_t layer;
        is >> layer;

        // Groups values - only diagnoal terms
        size_t group(1), index(0);
        while (is >> value)
        {
            index = assortative ? group
                                : group * group - 1 + layer * nof_groups * nof_groups;
            w[index] = value;
            group++;
        }
    }
}

void write_info_file(const boost::filesystem::path &output_filename,
                     const utils::Report &results)
{
    std::ofstream stream_out(output_filename.string());
    if (stream_out.fail())
    {
        throw(std::runtime_error(
            std::string("In write_info_file, failed to open ") + output_filename.string()));
    }

    stream_out << "# Number of realization = " << results.nof_realizations << std::endl;
    stream_out << "# Maximum Likelihood = " << results.max_L2() << std::endl;
    stream_out << "# Duration (s) = " << results.duration << std::endl;
    stream_out << "# Seed = " << results.seed << std::endl;
    stream_out << "# real  num_iters  term_reason  L2" << std::endl;
    for (size_t i = 0; i < results.vec_iter.size(); i++)
    {
        stream_out << "  " << i << "  " << results.vec_iter[i] << "  " << results.vec_term_reason[i]
                   << "  " << results.vec_L2[i] << std::endl;
    }
}
