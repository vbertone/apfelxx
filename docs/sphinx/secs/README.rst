|CircleCI| |CodeFactor|

.. image:: https://github.com/vbertone/apfelxx/raw/master/resources/logo.png

**A PDF evolution library in C++**

Introduction
============

APFEL++ is a C++ rewriting of the Fortran 77 evolution code
`APFEL <https://github.com/vbertone/apfel>`__. However, APFEL++ is based
on a completely new code design and guarantees a better performance
along with a more optimal memory management. The new modular structure
allows for better maintainability and easier extensibility. This makes
APFEL++ suitable for a wide range of tasks: from the solution of the
DGLAP evolution equations to the computation of
deep-inelastic-scattering (DIS) and single-inclusive-annihilation cross
sections. Also more complex computations, such as differential
semi-inclusive DIS and Drell-Yan cross sections, are easily
implementable in APFEL++.

APFEL++ is used as a prediction engine in
`NangaParbat <https://github.com/vbertone/NangaParbat>`__, a code
devoted to the extraction of Trasverse-Momentum-Dependent (TMD)
distributions, and is also currently interfaced to
`PARTONS <http://partons.cea.fr/partons/doc/html/index.html>`__, a
software dedicated to the phenomenology of Generalised Parton
Distributions (GPDs) and TMDs, and to
`xFitter <https://www.xfitter.org/xFitter/>`__, an open source fit
framework devoted to the extraction of collinear distributions and to
the assessment of the impact of new experimental data.

Homebrew (macOS)
----------------

APFEL++ is available through the
`Homebrew-hep <https://github.com/davidchall/homebrew-hep>`__ tap of
`Homebrew <https://brew.sh>`__. After having installed Homebrew and
tapped Homebrew-hep, the latest stable release of APFEL++ can be
downloaded and installed with:

.. code:: shell

   brew install apfelxx

Download
--------

You can obtain APFEL++ directly from the github repository:

https://github.com/vbertone/apfelxx/releases

For the latest development branch you can clone the master code:

.. code:: shell

   git clone https://github.com/vbertone/apfelxx.git

If you instead want to download a specific tag:

.. code:: shell

   git tag -l
   git checkout tags/tag_name

Installation
------------

The code can be compiled using the following procedure:

.. code:: shell

   cd apfelxx
   mkdir build && cd build
   cmake -DCMAKE_INSTALL_PREFIX=/your/installation/path/ ..
   make && make install

By the default, if no prefix is given, the program will be installed in
the /usr/local folder. If you want (or need) to use a different path,
remember to export the APFEL++ ``/lib`` and ``/bin`` folders into the
``$LD_LIBRARY_PATH`` and ``$PATH``, respectively. More configuration
options can be accessed by typing:

.. code:: shell

   ccmake .

Code documentation
------------------

Code documentation automatically generated with
`Doxygen <https://www.doxygen.nl/index.html>`__ can be found at
`this <https://vbertone.github.io/apfelxx/html/index.html>`__ link.

Relevant references
-------------------

-  V. Bertone, *APFEL++: A new PDF evolution library in C++*,
   `arXiv:1708.00911 <https://arxiv.org/pdf/1708.00911.pdf>`__.
-  V. Bertone, S. Carrazza, J. Rojo, *APFEL: A PDF Evolution Library
   with QED corrections*,
   `arXiv:1310.1394 <http://arxiv.org/abs/arXiv:1310.1394>`__.

Contact information
-------------------

-  Valerio Bertone: valerio.bertone@cern.ch

.. |CircleCI| image:: https://circleci.com/gh/vbertone/apfelxx.svg?style=svg&circle-token=079509ce5a2d4da15596d7812b33dca48eec8dc7
   :target: https://circleci.com/gh/vbertone/apfelxx
.. |CodeFactor| image:: https://www.codefactor.io/repository/github/vbertone/apfelxx/badge
   :target: https://www.codefactor.io/repository/github/vbertone/apfelxx
