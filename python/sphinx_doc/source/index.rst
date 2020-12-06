.. MultiTensor documentation master file, created by
   sphinx-quickstart on Fri Oct 30 19:44:15 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MultiTensor's documentation!
=======================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Definitions<definitions>

Introduction
============

This is the documentation for the Python extension of the MultiTensor,
a library for multilayer network tensor factorization that can be used
for community detection, link prediction and measure layer interdependence.

Requirements
============

* Python 3.6 or higher
* `NumPy <https://numpy.org/>`_ installed

Installation
============

The package can be easily installed using `pip <https://pypi.python.org/pypi/pip>`_ on a distribution (e.g. wheel, source tarball):

.. code::

	$ pip install multitensor_dist

Alternatively, you can also build the python extension yourself
in your build directory:

.. code::

    $ pip install -U numpy cython # cython is optional
    $ make multitensor_py


Getting started
===============

To use the algorithm, import the module inside you Python session
and use the :func:`multitensor.run` function:

.. code-block:: python

    import multitensor

    adjacency_file = '$MULTITENSOR_SRC_DIR/data/main/adjacency.dat'
    num_groups = 2

    u, v, affinity, report = multitensor.run(adjacency_file, num_groups)


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
