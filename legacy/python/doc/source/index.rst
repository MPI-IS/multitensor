.. MultiTensor documentation master file, created by
   sphinx-quickstart on Mon Feb 18 13:27:35 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MultiTensor's documentation!
=======================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
Introduction
============

This is the documentation for the Python extension of the multi-tensor factorization algorithm 
developed by `Caterina de Bacco <https://is.tuebingen.mpg.de/person/cdebacco>`__ and
the `Software Workshop <https://is.tuebingen.mpg.de/en/software-workshop>`_
at the `Max Planck Institute for Intelligent Systems <https://is.mpg.de/>`_.

Installation
============

The package can be easily installed using `pip <https://pypi.python.org/pypi/pip>`_ directly on a distribution
(source, binary, wheel):

.. code::
	
	$ pip install multi_tensor_dist

or use the setup.py

.. code::
	
	$ cd /path/to/multi/tensor/setup.py
	$ pip install .
	

Getting started
===============

To use the algorithm, import the module inside you *Python* session:

.. code:: python

	import MultiTensor

It comes with a demo which can be run using the following command:

.. code:: python

	MultiTensor.demo()

Reference
=========

De Bacco, C., Power, E. A., Larremore, D. B., & Moore, C. (2017). 
*Community detection, link prediction, and layer interdependence in multilayer networks.* 
Physical Review E, 95(4), 042317.

Authors
=======

`Caterina De Bacco <mailto:caterina.debacco@tuebingen.mpg.de>`_

`Jean-Claude Passy <mailto:jean-claude.passy@tuebignen.mpg.de>`_

License
=======

GNU GPL version 3 (see LICENSE.md)

Copyright
=========

\(c\) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
