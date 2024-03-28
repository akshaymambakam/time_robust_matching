# The robustTRE package
This repository is for sharing the source code of the robustTRE package. The package is for doing time-robust pattern matching using timed regular expressions in Python. We make use of **Parma Polyhedra Library (PPL)** for implementing the algorithm for getting the time-robust match-set from a given match-set and a given robustness value. We use **pybind11** to create **Python3** bindings.

## Installation
We give installation instructions for Ubuntu.
1. Install Parma Polyhedra Library (PPL).
```
sudo apt install ppl-dev
```
2. Install C++ compiler.
```
sudo apt install g++
```
4. Install the GNU Multiple Precision (GMP) Arithmetic Library.
```
sudo apt install libgmp-dev
```
5. Install
```
./build.sh
```
6. Run
```
python3 robustness_map.py
```

## Authors and acknowledgment
This package is built on top of the work done at the VERIMAG laboratory located in the Grenoble city of France. This work is based on the theory of timed pattern matching developed by [Dogan Ulus](https://www.cmpe.boun.edu.tr/tr/people/dogan.ulus). Check out his [github page](https://github.com/doganulus). Timed pattern matching has been implemented in [montre](https://github.com/doganulus/montre) and [timedrel](https://github.com/doganulus/timedrel).

## License
For open source projects, say how it is licensed.