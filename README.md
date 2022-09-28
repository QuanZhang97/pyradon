**Pyradon**
======

## Description

**Pyradon* is a python package for denoising and interpolation of multi-channel seismic data. This package has a variety of applications in both exploration and earthquake seismology.

## Reference
Zhang et al., 2023, Pyradon: a python package of Radon transform for multi-channel seismic data processing, in preparation. 

BibTeX:

	@article{pyradon,
	  title={Pyradon: a python package of Radon transform for multi-channel seismic data processing},
	  author={Quan Zhang et al},
	  journal={TBD},
	  volume={1},
	  number={1},
	  pages={1-10},
	  year={2022}
	}

-----------
## Copyright
    The pyradon developing team, 2021-present
-----------

## License
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)   

-----------

## Install
Using the latest version

    git clone https://github.com/QuanZhang97/pyradon
    cd pyradon
    pip install -v -e .
or using Pypi

    pip install pyradon

-----------
## Examples
    The "demo" directory contains all runable scripts to demonstrate different applications of pyradon. 

-----------
## Gallery
The gallery figures of the pyradon package can be found at
    https://github.com/QuanZhang97/gallery/tree/main/pyradon
Each figure in the gallery directory corresponds to a DEMO script in the "demo" directory with the exactly the same file name.

-----------
## Dependence Packages
* scipy 
* numpy 
* matplotlib

-----------
## Modules
    xxx.py  -> description
    
-----------
## Development
    The development team welcomes voluntary contributions from any open-source enthusiast. 
    If you want to make contribution to this project, feel free to contact the development team. 

-----------
## Contact
    Regarding any questions, bugs, developments, collaborations, please contact  
    Quan Zhang
    quanzhang1997@gmail.com

-----------
## Examples

A simple example of the adjoint Radon transform

<img src='./test.png' alt='Adj' width=960/>
	Figure 1. Original Z component data with its Radon spectrum.

<img src='./test_i.png' alt='Adj' width=960/>
	Figure 2. HRT denoised Z component data with its Radon spectrum.
