# Variable Projection for DMD

Code for creating the examples in the paper "Variable
projection methods for debiasing and generalizing
the dynamic mode decomposition"

## Set-up

These codes make use of Timothy Toolan's MatLab
scripts for calling LAPACK functions. These require
that you compile a mex binary, but the process is
largely automated. Open MatLab and run setup.m
to install. Enter "y" when prompted.

## How to use

If you'd like to see how to use the OPTDMD wrappers
(for computing the optimized dynamic mode decomposition)
the best place to start is to check out simple_example.m
(in "examples" folder)

## To run the figure generation (may take a while)

Be sure to add the "src" directory to MatLab's path.
The scripts for running the examples are found in the
"examples" folder. To run all of them at once (just
generating dat), see "run_all_examples.m". Then, to
create figures, see the individual *figs.m files.

## Updates

This directory is meant to be something of a snapshot
of the codes as they were used for generating the
examples. For the most up-to-date versions, see 
[the optdmd repository](github.com/askhamwhat/optdmd).

## Citing this software

To cite this figure generation package specifically
(not recommended) follow this Zenodo link: 
[![DOI](https://zenodo.org/badge/82845075.svg)](https://zenodo.org/badge/latestdoi/82845075)

## License 

The files in the "src" and "examples" directories are available under the MIT license unless noted otherwise (see license* files in src directory).

The MIT License (MIT)

Copyright (c) 2017 Travis Askham

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.