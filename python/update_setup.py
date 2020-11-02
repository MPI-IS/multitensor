# Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
# Distributed under the GNU GPL license version 3
# See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md

"""Script for including the relevant extensions to the setup.py"""

import sys

if __name__ == '__main__':
    file_setup = sys.argv[1]
    version = sys.argv[2]
    libname = sys.argv[3]

    # Read
    with open(file_setup, 'r') as fin:
        content = fin.read()
    # Replace content
    content = content.replace('{{version}}', version)
    content = content.replace('{{libname}}', libname)

    # Write over
    with open(file_setup, 'w') as fout:
        fout.write(content)
