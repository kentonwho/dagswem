# DaGSWEM

**D**iscontinuous **A**daptive **G**alerkin (**D**irected **A**cyclic **G**raph) **S**hallow **W**ater **E**quation **M**odel

This project is a fork of the main DGSWEM project that implements both an adaptive local timestepping scheme and a task-parallel (or directed acyclic task graph) execution mode, hence the double acronym. The project was chosen to be forked completely since the end product has a completely different execution architecture, and the upstream is from a decade+ old untracked branch (DGSWEM-LTS), so merging the two now is an impossibility. 

## Prerequisites

Before you begin, ensure you have the following installed:
*   A Fortran compiler (e.g., `gfortran`)
*   `make` build automation tool
*   An MPI implementation (e.g., OpenMPI, MPICH) for parallel execution

## Building

The project uses `make` for compilation. Compiler flags and library paths can be configured in [`cmplrflags.mk`](cmplrflags.mk).

1.  **Build METIS library**: The METIS partitioning library is included as a dependency.
    ```bash
    cd metis/
    make
    cd ..
    ```

2.  **Build the main project**: Compile the preprocessing utilities and the main simulation code.
    ```bash
    make
    ```
    This will use the main [`makefile`](makefile) to build the executables from the source files in [`prep/`](prep/) and [`src/`](src/).

## Usage

The typical workflow involves a preprocessing step followed by the main simulation run.

1.  **Preprocessing**: Use the tools built from the [`prep/`](prep/) directory (e.g., `adcprep`, `decomp`) to prepare your input mesh and data.

2.  **Running the simulation**: Execute the main program, which appears to be `adcirc` from the [`src/`](src/) directory. For parallel runs, use `mpirun`.
    ```bash
    mpirun -np <number_of_processes> ./adcirc
    ```

## Directory Structure

*   `src/`: Contains the main source code for the simulation (e.g., [`adcirc.F`](src/adcirc.F)).
*   `prep/`: Contains source code for preprocessing and utility programs (e.g., [`adcprep.F`](prep/adcprep.F), [`decomp.F`](prep/decomp.F)).
*   `metis/`: Source code and compiled library for the METIS graph partitioning software.
*   `makefile`: The main makefile for building the project.
*   `cmplrflags.mk`: Makefile include for compiler flags and settings.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

[MIT](https://choosealicense.com/licenses/mit/)