# [Birkhoff.py](https://github.com/yonesuke/Birkhoff.py)

Implementation of [Birkhoff normal form](https://encyclopediaofmath.org/wiki/Birkhoff_normal_form) in 
Python

## Example

For Hamiltonian

<img src="https://latex.codecogs.com/gif.latex?H=\frac{1}{4}p_{r}^{2}+\frac{\omega^{2}}{4r^{2}}+\frac{\alpha+2}{4\alpha}p_{z}^{2}-\frac{1}{2\alpha r}-\frac{2}{\sqrt{r^{2}+z^{2}}},"/>

the point <img src="https://latex.codecogs.com/gif.latex?(p_r,p_z,r,z)=(0,0,\omega^{2}\alpha/(4\alpha+1),0)"/> is an equilibrium point.
By using [Birkhoff.py](https://github.com/yonesuke/Birkhoff.py), we can calculate the Birkhoff normal form as

![Birkhoff normal form](figure/shibayama2009.png)

See [shibayama2009.py](shibayama2009.py) for a detailed code.

## Installation
```bash
git clone https://github.com/yoensuke/Birkhoff.py
```

## Requirements

You will need
 `Python` 3
 `SymPy`

If you do not have `SymPy`, install it via pip.
```bash
pip3 install sympy
```
For another installation, check [here](https://docs.sympy.org/latest/install.html).

## Usage

Place [`Birkhoff.py`](Birkhoff.py) and [`Hamilton.py`](Hamilton.py) to a directory you work for.
For the rest, just bring Hamiltonian and its equilibrium point and mimic [shibayama2009.py](shibayama2009.py).



## [License](LICENSE)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
