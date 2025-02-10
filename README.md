# Trimmed Mesher
![C](https://img.shields.io/badge/language-C-00599C.svg) ![CMake](https://img.shields.io/badge/build-CMake-brightgreen.svg)

This repository contains a C implementation of a **trimmed mesher** designed for **Computational Fluid Dynamics (CFD)** simulations. The mesher generates a structured **Cartesian 2D grid** composed of square cells and modifies it to accommodate an **embedded closed curve**, which represents an obstacle in the flow domain.

## Requirements

- C compiler (e.g., `gcc`, `clang`, etc.)
- Standard C libraries (`stdio.h`, `stdlib.h`, `math.h`)
- [CMake](https://cmake.org/download/) (version 3.10 or higher)
- [git](https://git-scm.com/downloads) (optional, for version control)

## Compilation

To compile the code, navigate to the project directory and use the following commands:

```bash
mkdir build
cd build
cmake ..
make
```
The compiled executable will be generated inside the `build` folder.

## Usage
Run the compiled executable using the following command:
```bash
./mesher <curve.csv> <cell_size> <rows> [columns] [X0] [Y0]
```
### Positional Arguments
- `<curve.csv>`: Path to the CSV file containing the **closed curve** (list of X, Y coordinates).
- `<cell_size>`: Size of each square cell in the Cartesian grid.
- `<rows>`: Number of rows in the initial structured grid.

### Optional Arguments
- `[columns]` (default = same as `rows`): Number of columns in the initial structured grid.
- `[X0]` (default = `0.0`): X-coordinate of the center of the grid.
- `[Y0]` (default = `0.0`): Y-coordinate of the center of the grid.

### Example usage
1) **Basic Usage**: Generate a mesh with a **100x100 grid** and a cell size of `0.05`, using `curve.csv` as input:
    >`./mesher curve.csv 0.05 100`

2) **Custom Grid Size**: Generate a **150x200 grid** with cell size `0.02`:
    >`./mesher curve.csv 0.02 150 200`

3) **Offset Grid Position**: Generate a **200x200 grid** and move the center of the grid to `(-5.0, -5.0)`:
    >`./mesher curve.csv 0.01 200 200 -5.0 -5.0`

## Input Format
The input curve must be provided as a CSV file with **X, Y coordinates**, representing a **closed polygonal curve**.

## Output
The output `mesh.msh` file follows the **GMSH ASCII format (version 2.2)** and includes both **node** and **elements**.

## License
This project is licensed under the Apache 2.0 License. See the [LICENSE](LICENSE) file for details.