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

The package can be easily installed using `pip <https://pypi.python.org/pypi/pip>`_ on a distribution:

.. code::

	$ pip install multitensor_dist

Alternatively, you can also build the python extension yourself
in your build directory:

.. code::

    $ make multitensor_py


Getting started
===============

To use the algorithm, import the module inside you *Python* session:

.. code:: python

	import multitensor


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
