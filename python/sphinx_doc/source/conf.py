# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

from pathlib import Path
import os
import sys
import sphinx_bootstrap_theme


# -- Project information -----------------------------------------------------
project = "@CMAKE_PROJECT_NAME@"
copyright = (
    "2019, Max Planck Society / Software Workshop - "
    "Max Planck Institute for Intelligent Systems"
)
author = (
    "Caterina De Bacco, "
    "Jean-Claude Passy, "
    "Ivan Oreshnikov, "
)

# The full version, including alpha/beta/rc tags
version = "@CMAKE_PROJECT_VERSION@"
release = "@CMAKE_PROJECT_VERSION@"


# -- General configuration ---------------------------------------------------

# those informations are passed by the cmake
if os.getenv('PYTHON_ADDITIONAL_FOLDERS'):
    sys.path.insert(0, os.getenv('PYTHON_ADDITIONAL_FOLDERS'))
else:
    sys.path.insert(0, os.path.abspath('.'))

# Make sure __init__ is always documented
# from https://stackoverflow.com/a/5599712/1617295


def skip(app, what, name, obj, skip, options):
    if name == "__init__":
        return not (hasattr(obj, "__doc__") and obj.__doc__)
    return skip


def setup(app):
    app.connect("autodoc-skip-member", skip)


# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'bootstrap'
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []
