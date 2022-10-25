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

## Namespaces

Everything that SPADE can do is organized into its various namespaces, all of which
fall under the `spade` namespace, e.g. the `io` namespace would be used as
`spade::io`. The list of namespaces within SPADE are as follows:

- `timing`: contains various utilities for measuring the performance of implemented
applications.

- `algs`: i.e. "algorithms". Contains a series of general-purpose algorithms that
are useful for managine the state of arrays or performing loops compactly.

- `aliases`: an intermediate namespace that forwards STL types for use in CPU
implementations and providing polymorphism for GPU implementations. For example,
`aliases::vector<double>` will be `std::vector<double>` when compiling as a CPU
application, but will be something else when compiling with CUDA.

- `array_container`: A specialized namespace that handles underlying storage for
grid arrays.

- `detail` or `[namespace]::detail`: a quarantined space for highly specialized
implementations that should not be called or manipulated outside very specific
situations. Don't come in here!

- `block_config`: provides various implementations of block configurations for
grid arrays.

- `cli_args`: contains utilities for handling command-line arguments

- `consts`: contains mathematical constants

- `coords`: contains everything to do with coordinate systems and their jacobians
for use with generalized-coordinates formulations

- `ctrs`: i.e. "containers". Contains various flavors of basic arrays with arithmetic
operations enabled.

- `dims`: soon to be deprecated!

- `finite_diff`: algoriths for the calculation of finite-difference coefficients

- `grid`: contains everything to do with the semantics of handling data stored on a grid

- `io`: contains functions for outputting data to files

- `linear_algebra`: contains implementations for numerical linear algebra algorithms
for small-size problems

- `parallel`: contains types and functions dealing with the use of a parallel computing
context

- `utils`: contains general-purpose utilities

- `partition`: contains implementations of partitioning algorithms for load-balancing

- `pde_algs`: contains algorithms specifically intended for the numerical solution of
partial differential equations

- `reduce_ops`: contains general definitions of reduction operations 

- `static_math`: contains implementations for compile-time mathematical calculations

- `time_integration`: contains algorithms for causal integration in a single dimension,
e.g. through time

- `fetch`: contains specialized algorithms for buffering grid array data for the 
calculation of PDE terms. Careful in here!

- `convective`: contains implementations of inviscid fluxes for the solution of
the compressible Navier-Stokes equations (CNSE)

- `viscous`: contains implementations of viscous fluxes for the solution of
CNSE

- `fluid_state`: contains implementations of and conversions between fluid states
of a gas

- `navier_stokes_mms`: soon to be deprecated!

- `viscous_laws`: contains state-dependent viscous laws for CNSE

- `amr`: contains algorithms for adaptive mesh refinement