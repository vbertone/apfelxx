[![CircleCI](https://circleci.com/gh/vbertone/apfelxx.svg?style=svg&circle-token=079509ce5a2d4da15596d7812b33dca48eec8dc7)](https://circleci.com/gh/vbertone/apfelxx)
[![CodeFactor](https://www.codefactor.io/repository/github/vbertone/apfelxx/badge)](https://www.codefactor.io/repository/github/vbertone/apfelxx)

![alt text](https://github.com/vbertone/apfelxx/raw/master/resources/logo.png "Logo APFEL")

# APFEL++: A PDF evolution library in C++
 
APFEL++ is a C++ rewriting of the Fortran 77 evolution code APFEL (see
http://apfel.hepforge.org and
https://github.com/vbertone/apfel). APFEL++ is based on a completely
new code design and guarantees a better performance along with a more
optimal memory management. The new modular structure allows for a
better maintainability and an easier extensibility. This makes APFEL++
suitable for a wide range of tasks: from the solution of the DGLAP
evolution equations to the computation of deep-inelastic-scattering
(DIS) and single-inclusive-annihilation cross sections. Also more
complex computations, such as semi-inclusive DIS and Drell-Yan cross
sections, are easily implementable in APFEL++.

## Homebrew (macOS)

APFEL++ is available through the [Homebrew-hep](https://github.com/davidchall/homebrew-hep) tap of [Homebrew](https://brew.sh). After having installed Homebrew and tapped Homebrew-hep, the latest stable release of APFEL++ can be downloaded and installed with:
```Shell
brew install apfelxx
```

## Download

You can obtain APFEL++ directly from the github repository:

https://github.com/vbertone/apfelxx/releases

For the last development branch you can clone the master code:

```Shell
git clone https://github.com/vbertone/apfelxx.git
```

If you instead want to download a specific tag:

```Shell
git tag -l
git checkout tags/tag_name
```
## Installation 

The code can be compiled using the following procedure:

```Shell
cd apfelxx
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/your/installation/path/ ..
make && make install
```
By the default, if no prefix specification is given, the program will
be installed in the /usr/local folder. If you want (or need) to use a
different path, remember to export the APFEL++ /lib folder into the
LD_LIBRARY_PATH. More configuration options can be accessed by typing:

```Shell
ccmake .
```

## Documetation

Code documentation generated with Doxygen can be found here:
https://vbertone.github.io/apfelxx/html/index.html.

## Relevant references

- V. Bertone, *APFEL++: A new PDF evolution library in C++*, [arXiv:1708.00911](https://arxiv.org/pdf/1708.00911.pdf).
- V. Bertone, S. Carrazza, J. Rojo, *APFEL: A PDF Evolution Library with QED corrections*, [arXiv:1310.1394](http://arxiv.org/abs/arXiv:1310.1394).

## Contact Information

- Valerio Bertone: valerio.bertone@cern.ch
