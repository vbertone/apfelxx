[![CircleCI](https://circleci.com/gh/vbertone/apfelxx.svg?style=svg&circle-token=079509ce5a2d4da15596d7812b33dca48eec8dc7)](https://circleci.com/gh/vbertone/apfelxx)

![alt text](https://github.com/vbertone/apfelxx/raw/master/resources/logo.png "Logo APFEL")

# APFEL++: A new PDF evolution library in C++
 
APFEL++ is a C++ rewriting of the Fortran 77 evolution code APFEL (see
http://apfel.hepforge.org and
https://github.com/scarrazza/apfel). APFEL++ is based on a completely
new code design guarantees a better performance and a more optimal
memory management. The new modular structure allows for a better
maintainability and an easier extensibility. This makes APFEL++
suitable for a wide range of tasks: for the solution of the DGLAP
evolution equations to the computation of deep-inelastic-scattering
(DIS) and single-inclusive-annihilation cross sections. Also more
complex computations, like semi-inclusive DIS and Drell-Yan cross
sections, are easily implementable in APFEL++.

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
mkdir build;
cd build;
cmake ..
make && make install
```
By the default, if no further specification is given, the program will
be installed in the /usr/local folder. If you want (or need) to use a
different path, remember to export the APFEL++ /lib folder into the
LD_LIBRARY_PATH.

## Documetation

Code documentation generated with Doxygen can be found here:
https://vbertone.github.io/apfelxx/html/index.html.

Visit also the HepForge webpage http://apfel.hepforge.org and the
APFEL Web interface http://apfel.mi.infn.it/.

## Relevant references

- V. Bertone, *APFEL++: A new PDF evolution library in C++*, [arXiv:1708.00911](https://arxiv.org/pdf/1708.00911.pdf).
- V. Bertone, S. Carrazza, J. Rojo, *Doped Parton Distributions*, [arXiv:1509.04022](http://arxiv.org/abs/1509.04022).
- V. Bertone, S. Carrazza, D. Pagani, M. Zaro, *On the Impact of Lepton PDFs*, [arXiv:1508.07002](http://arxiv.org/abs/1508.07002).
- V. Bertone, S. Carrazza, E. Nocera, *Reference results for time-like evolution up to O(α_s^3)*, [arXiv:1501.00494](http://arxiv.org/abs/1501.00494).
- S. Carrazza, A. Ferrara, D. Palazzo, J. Rojo, *APFEL Web: a web-based application for the graphical visualization of parton distribution functions*, [arXiv:1410.5456](http://arxiv.org/abs/1410.5456).
- V. Bertone, S. Carrazza, J. Rojo, *APFEL: A PDF Evolution Library with QED corrections*, [arXiv:1310.1394](http://arxiv.org/abs/arXiv:1310.1394).

## Contact Information

- Valerio Bertone: valerio.bertone@cern.ch
- Stefano Carrazza: stefano.carrazza@cern.ch
