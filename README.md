seis-stack-td
=============

Andy Nowacki, University of Leeds <a.nowacki@leeds.ac.uk>

A set of Fortran programs to compute vespagrams and f-k spectra for broadband
seismic data in SAC format.

The computations are all done in the time domain, making these routines very
slow.  They are recommended only for educational purposes (and were written for
that reason).


Installation
============

Required are:

- A Fortran compiler
- The Generic Mapping Tools for plotting
- The NetCDF library for writing files as NetCDF grid files
- FFTW

Optional dependency:
- SAC's XAPIIR filtering library (libxapiir)

To build, adjust the variables in `Makefile`, and the `make`.

If using the XAPIIR library for filtering, do `make filt`.  This is recommended
if you have SAC installed as the library is distributed with SAC.


Usage
=====

The programs `plot_fk.sh` and `plot_vesp.sh` do stacking and create plots.
Run these programs with the `-h` command line option to print out instructions.

They call the programs `bin/stack_fk` and `bin/stack_vespa`, which should be
available in your `$PATH`.  If you are current;y in the root directory of the
repo, this can be done by

```sh
$ export PATH="$PATH:$PWD/bin:$PWD/script"
```


Licence
=======

This software is distributed under the terms of the three-clause BSD licence:

Copyright (c) 2016-17, Andy Nowacki
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the University of Leeds nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
