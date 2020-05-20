# ChocoboJ
Chocobo CFD Code - Julia Version

## Building ChocoboJ
ChocoboJ is intended to be built into an executable using the PackageCompiler Julia package. Compilation is performed using the `compile.jl` file, but generally the provided makefile should be used. The standard build process. Generally, the following command should be used:

    make; make install

This will build ChocoboJ, and install it. By default, the installation path is `/prod/ChocoboJ/bin`, with a simlink placed in `/prod/bin`, although these can be changed in the makefile.

## Running ChocoboJ
To run ChocoboJ, assuming its installation location is in your `$PATH`, simply type:

    chocoboj

This will assume that the following input files are located in the current directory:

| File        | Purpose                          |
|-------------|----------------------------------|
| input.in    | Main input file                  |
| mesh.key    | Mesh file (in LSDYNA KEY format) |
| eos_data.jl | EoS data file                    |

