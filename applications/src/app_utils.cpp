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
