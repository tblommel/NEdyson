NEdyson is a program which solves the Kadanoff-Baym equations on the three legged Keldysh 
contour.  These coupled integro-differential equations are solved using a timestepping
scheme described in the NESSi paper: https://arxiv.org/abs/1911.01211.  One difference from
the NESSi implementation is that we use a Legendre grid instead of a uniform grid to store 
data along the imaginary axis.

NEdyson depends on the gfmol package in order to solve the equilibrium system, this package
can be installed from https://github.com/CQMP/gfmol.
NEdyson also depends on cpsitop, which handles the Legendre grid convolutions and 
integrations and can be installed from https://github.com/xinyangd/cpsitop.

The build script provided will build the program molNEdyson in the build directory.
To run the program, one must provide an input file which contains all of the parameters
necessary for the calculation.  An example input file is given in data/NEdyson-input.ini.
For information on what the input parameters are, one can call the program with the -h
flag for a list of parameters and their meaning.

This package also forms the basis for the compression-based integration package: hodlr.  As well as NEadapt,
an adaptive implementation which interfaces with the SUNDIALS package.  For further information on either of these
implementations, please contact the author.
