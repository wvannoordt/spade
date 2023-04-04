# Getting Started

This tutorial will demonstrate setting up a basic executable with SPADE.
SPADE requires MPI to operate, so it's recommended that you install MPI
on your machine as executables should be built with the `mpicxx` compiler
wrapper.

SPADE is written in C++20, meaning you will need at least G++ version 10.3+
to compile any executables that use SPADE. For this tutorial, it is assumed
that you have set the environment variable SPADE to the root directory of
the main SPADE repository, i.e. `/my/path/to/spade`. It is also assumed that
`mpicxx` is on your path.

`main.cc` contains some simple code to check if SPADE is working.