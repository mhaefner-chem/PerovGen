# PerovGen
![GitHub Release](https://img.shields.io/github/v/release/mhaefner-chem/PerovGen?include_prereleases) ![GitHub License](https://img.shields.io/github/license/mhaefner-chem/PerovGen)

A python program with a GUI that turns a Jagodzinski sequence into a symmetrized CIF-file of the corresponding perovskite.

## Table of contents

- [How to use PerovGen?](#how-to-use-perovgen)
- [Requirements & Installation](#requirements-and-installation)
- [How to cite?](#how-to-cite)
- [Further license information](#further-license-information)

## How to use PerovGen?
PerovGen performs a symmetry analysis on a Jagodzinski sequence and turns it into a symmetrized CIF structure file of the corresponding perovskite. As such, only the sequence needs to be provided into the field "Jagodzinski sequence" as a string of "h" and "c" denoting a hexagonal and cubic ordering of the stacking sequence of the associated perovskite. The program performs all necessary reductions and expansions of the sequence to turn it into a valid structure, and in addition to generating the structure file, it also outputs the stacking sequence of layers A, B, and C and interlayers α, β, and γ as well as the location of symmetry centers along the layers.

Additionally, the elements comprising the perovskite can be specified in the fields A, B, and X, with A being the in-layer cation, B the interstitial cation, and X the anion according to the generalized sum formula ABX<sub>3</sub>.

## Requirements and installation

So far, the program has been successfully tested with `python 3.10` and `3.12` and `spglib 2.4.0`. The compiled binaries have been tested on Windows 10-64bit.
The source code for the program is obtained with the command

```console
$ git clone https://github.com/mhaefner-chem/PerovGen
```

Running the program with python requires the python package `spglib` and all its dependencies for standardizing the generated unit cell. It can be installed via pip

```console
$ pip install spglib
```
or via conda as

```console
$ conda install -c conda-forge spglib
```
For further details on installing spglib, please consult the [github page](https://github.com/spglib/spglib) or the [documentation](https://spglib.readthedocs.io) of spglib.

The compiled binaries are attached to the assets of the released version of the program and should work out of the box on Windows systems.

## How to cite?
Not yet published.

## Further license information

**Icon**

The icon for this software was created with VESTA (version 3.5.8)

K. Momma and F. Izumi, 
"VESTA 3 for three-dimensional visualization of crystal, volumetric and morphology data," 
J. Appl. Crystallogr., 44, 1272-1276 (2011).

**spglib**

The binaries were compiled with spglib for symmetry analysis.

Copyright (c) 2024, Spglib team
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
* Neither the name of the <organization> nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT
HOLDER> BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

<!--stackedit_data:
eyJoaXN0b3J5IjpbMTc3MzE0MTE1MiwtNTE5NDYwNTk4LC05OD
I4MjI4NzNdfQ==
-->
