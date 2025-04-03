# MoT-Voellmy –– a Numerical Model for Avalanche Motion

## Introduction

"MoT-Voellmy" stands for "Method of Transport" and "Voellmy-type avalanche
model". This is a code for simulating the propagation of a (granular) mass down
a slope:
- The slope can be an arbitrary topography described by a Digital Terrain Model
  (DTM), provided as a raster file with quadratic cells spanning a rectangle.
- The flow is described by the depth-averaged balance equations for mass and
  momentum. This means that the 3D conservation equations for mass and momentum
  are integrated across the depth of the flow so that the resulting equations
  are effectively 2D (in the local tangential plane to the topography), with an
  extra equation for the spatial distribution and temporal evolution of the
  flow depth.
- The avalanche moves under the action of gravity and resistive forces along
  the bottom of the flow.
- The Voellmy avalanche model from 1955 assumes that the bed shar stress,
  __*τ*__<sub>b</sub>, has two components, one proportional to the normal
  stress at the bed, *σ<sub>n</sub>*, and independent of the velocity __*u*__,
  and another independent of the normal stress and proportional to the square
  of the velocity:<br>
  __*τ*__<sub>*b*</sub> = &minus;__*u*__/||__*u*__||
                         [*μσ*<sub>*n*</sub> + *kρ*__*u*__<sup>2</sup>].<br>
  *ρ* is the density of the avalanche flow and is assumed constant in this
  model. The first term can be considered as a Coulomb failure criterion, the
  second as some kind of "turbulent" drag or the collisional shear stress in a
  confined and rapidly sheared granular material. The corresponding
  dimensionless coefficients *μ* and *k* are to be considered empirical
  paramaters.

The *Method of Transport* is a particular method for discretizing hyperbolic
(advection-dominated) conservation equations. It is implemented in its very
simplest variant in this code, using only one wave component (corresponding to
the flow velocity _*u*_) and constant reconstruction of field values. In this
way, it is very similar to a first-order upwind finite-difference scheme. As a
consequence, the code is rather diffusive but also very fast.

The code has the tendency to develop so-called chequerboard oscillations, the
reason for which is not entirely clear yet. They may lead to instability. To
counter this, there are a few ad-hoc measures the user can take:
- The so-called CFL number can be chosen smaller than its default value 0.8.
- The minimum allowed time step (typically set to 1 ms for avalanches on a
  5 m grid) can be reduced.
- The pressure may be multiplied with a coefficient < 1 to reduce its impact on
  the simulation. The next stage of development should aim at eliminating this
  weakness.

## Getting Started

Getting MoT-Voellmy to run on your machine is as simple as copying the
executable to some directory. Sorry, no fancy installer for Windows... There
are no extra libraries (DLLs on Windows) to install. However, MoT-Voellmy will
not be known to the Windows registry.

On Linux and macOS, you need to make sure that the file permissions allow
executing the code. On Windows, you can simply run the command from the command
shell (cmd.exe) if the executable is in a directory in the user's path.

__Command synopsis:__<br><br>
`MoT-Voellmy.<version date> <(path)name of simulation control file>`<br><br>
`<version date>` is to be replaced by the date characterizing the version you
have downloaded and wish to use. At the time of writing, the most recent
version is 2025-02-10. Windows users may have to write<br><br>
`MoT-Voellmy.<version date>.exe <(path)name of simulation control file>`

The simulation control file (SCF) is a simple text file in a specific format.
In the course of development of MoT-Voellmy, this format has evolved somewhat.
A template with the most recent version (which may be older than the
executable!) is available in the main branch of the repository.

At NGI, MoT-Voellmy can be run from within ArcGIS Pro, which makes it easier to
prepare the necessary input files. The GUI also takes care of writing the SCF.
The form for specifying the input data also offers a help facility, which
explains the meaning of the parameters that the user can set.

## Build and Test

The code is written in ISO-C without C99 extensions, etc. It should therefore
compile out of the box with any standard-compliant C compiler. However, this
has not been tested; in particular, there might be adjustments needed if one
wishes to use Microsoft's C compiler instead of GCC. The Windows executable in
this repository is compiled under Linux using the MinGW tool chain.

The command lines for compiling with gcc are:<br><br>
*Linux*:<br>
    `gcc -Wall -pedantic -o MoT-Voellmy-l64 MoT-Voellmy.c -lm`<br><br>
*macOS*:<br>
    `gcc -Wall -pedantic -o MoT-Voellmy-m64 MoT-Voellmy.c -lm`<br><br>
*Cross-compilation on Linux for MS Windows*:<br>
    `x86_64-w64-mingw32-gcc -Wall -pedantic -o MoT-Voellmy-w64.exe
 MoT-Voellmy.c -lm`

## Further development

At this point, DI's prioritized list of further developments is the following:
- Improve the stability and eliminate chequerboard oscillations by adding a
  mass diffusion term that mimicks the neglected gravity waves.
- Implement the deposition model built into MoT-PSA by Hervé Vicari.
- Replace the input and output functions by more robust and general ones.
- Merge MoT-Voellmy with Callum Tregaskis' fork MoT-muI, which implements the
  *μ*(*I*) rheology.
- Implement impermeable cells, i.e., obstacles to the flow.

If you want to contribute to the development of the code, please create a new
branch and explain the goals of the development.
