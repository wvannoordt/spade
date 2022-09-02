# SPADE (Static Polymorphic Algorithms for Differential Equations)

SPADE is a framework for the solution of partial differential equations in space and
time. It provides a series of high-order algorithms that are intended to make the
the implementation of solvers as smooth and efficient as possible.

## Building and Installation

SPADE is currently a header-only, single-include implementation that utilized
features available in the c++20 standard. To build a project using SPADE, simply
include the header `src/spade.h` and ensure that the c++20 standard is used.

## Examples

Examples are currently found in the `runs` directory. `runs/channel` contains an
example solver for a compressible turbulent channel flow, solved using the
Navier-stokes equations.