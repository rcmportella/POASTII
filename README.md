# BOAST II - Python Implementation

Python conversion of BOAST II (Black Oil Applied Simulation Tool), a three-dimensional, three-phase (oil, gas, and water) black oil reservoir simulator originally developed in Fortran 77.

## Overview

BOAST II uses the IMPES (Implicit Pressure, Explicit Saturation) method to simulate oil reservoir behavior:
- **Implicit Pressure**: Solves pressure equations using numerical methods
- **Explicit Saturation**: Updates saturations based on calculated flow rates

### Original Source
- **Language**: Fortran 77
- **Release**: Version 1.2
- **Author**: J.R. Fanchi (April 1987)

## Project Structure

### Core Simulation Files

#### `main.py`
Main entry point and simulator orchestration.
- `BOASTSimulator` class: Central coordinator for all simulation components
- `run_simulation()`: Main simulation loop implementing IMPES method
- Manages time stepping and dataset processing
- Coordinates all block modules
- Handles input/output file operations

### Block Modules (Fortran → Python)

#### `block1.py`
Aquifer models, material balance, and well management.
- `AquiferModel`: Carter-Tracy analytical aquifer
- `MaterialBalance`: Track oil, water, and gas volumes
- `WellManager`: Read and manage well data from input files

#### `block2.py`
Grid setup, solution control, and utility functions.
- `GridSetup`: Read and initialize grid dimensions and properties
- `SolutionControl`: Read simulation control parameters (time steps, solution method, etc.)
- `Interpolation`: Linear interpolation utilities for PVT and relative permeability tables
- `LSORSolver`: Line Successive Over-Relaxation iterative solver

#### `block3.py`
Output generation and visualization.
- `OutputVisualization`: Format and write simulation output
- `RockProperties`: Manage rock property distributions
- `PostProcessing`: Generate summary reports and pressure/saturation maps

#### `block4.py`
Gas properties and tridiagonal solver.
- `GasProperties`: Calculate gas viscosity, Z-factor, and FVF using correlations
- `Repressurization`: Handle repressurization calculations
- `TridiagonalSolver`: Direct solution for 1D problems (GAUS1D)

#### `block5.py`
Well rate calculations and boundary conditions.
- `WellRates`: Calculate production/injection rates for all well types
  - Rate-controlled wells (oil, water, gas)
  - Pressure-controlled wells (BHP)
  - GOR/WOR constraints
  - Gas wells with back-pressure equation
- `qrate()`: Main well rate calculation routine
- `pratei()`: Modify pressure matrix for implicit wells
- `prateo()`: Calculate implicit well rates after pressure solution

#### `block6.py`
PVT tables, transmissibilities, and initialization.
- `PVTTables`: Read and store fluid property tables (Bo, Bg, Bw, Rs, viscosities)
- `Transmissibility`: Calculate geometric transmissibilities between grid blocks
- `Initialization`: Initialize pressure and saturation distributions
- Helper functions: `zandc()`, `viscy()`, `trikro()` (Stone's 3-phase relative permeability)

#### `block7.py`
Flow equation assembly (one-point upstream weighting).
- `FlowEquation`: Build 7-diagonal coefficient matrix for pressure equation
- `solone()`: Single-point upstream weighting scheme
- Calculates flow terms between neighboring blocks
- Handles gravity, capillary pressure, and mobility

#### `block8.py`
Flow equation assembly (two-point weighting).
- `FlowEquationTwoPoint`: Alternative flow equation formulation
- `soltwo()`: Two-point weighted scheme
- More accurate than single-point for heterogeneous reservoirs

## Input/Output Files

### Input Files
- `EX1.DAT`: Example reservoir data file containing:
  - Grid dimensions and properties
  - Rock properties (porosity, permeability)
  - PVT tables (oil, water, gas properties)
  - Relative permeability curves
  - Well locations and rates
  - Time stepping parameters

### Output Files
- `OUT1_*.DAT`: Simulation results including:
  - Initialization summary
  - Time step summaries (rates, pressures, material balance)
  - Detailed reports at specified times
  - Pressure and saturation distributions
  - Production/injection summaries

## Running the Simulator

### Prerequisites
```bash
python 3.10+
numpy
```

### Basic Execution

#### Interactive Mode
```bash
python main.py
```
Then enter:
- Input file name: `EX1.DAT`
- Output file name: `OUT1.DAT`

#### Automated Mode
```bash
python run_test.py
```
Uses predefined input/output filenames.

#### Command Line
```bash
echo -e "EX1.DAT\nOUT1.DAT" | python main.py
```

## Current Implementation Status

### ✅ Fully Implemented and Tested
- Complete input reading (all blocks)
- Grid initialization and property assignment
- PVT table interpolation with bubble point handling
- Rock property distribution (relative permeability, capillary pressure)
- Total compressibility calculation (CT = CR + CO*SO + CW*SW + CG*SG)
- Well rate calculations (QRATE) for all well types
- Pressure equation assembly (SOLONE/SOLTWO)
- Implicit well control (PRATEI/PRATEO)
- Pressure solution methods:
  - Direct solver (GAUS1D) for 1D problems
  - LSOR iterative solvers (X, Y, Z directions)
- Explicit saturation updates using IMPES method
- Material balance calculations (production/injection rates)
- Time stepping with proper state array updates
- Output formatting:
  - Time step summary tables
  - Detailed pressure and saturation maps (table and contour formats)
  - Production/injection rate summaries
  - Material balance error tracking
- Debug output control (matrix printing, convergence details)

### ✅ Recent Fixes (Feb 2026)
- Fixed array indexing conversions (Fortran 1-based → Python 0-based)
  - PVT table access in saturation calculations
  - Formation volume factor interpolation in well calculations
- Corrected debug output control (ksn1, ksm1, kco1, kcoff parameters)
- Added missing matbal() call for production/injection rate calculations
- Implemented old-time array updates (pn, son, swn, sgn) after each time step
- Fixed saturation map output formatting consistency
- Corrected time step summary printing:
  - Added one-line summary for each time step
  - Header reprints after each detailed report
  - Actual calculated values (not placeholders)

### 🔄 In Progress
- Material balance error calculations (currently showing -100%)
- Initial fluid in place calculations (TOOIP, TOWIP, TOGIP)

### 📋 To Do
- Aquifer influx calculations (AQUI module)
- Restart capability (read/write restart files)
- Additional solution methods (D4 Gaussian elimination)
- Performance optimization and vectorization
- Validation against original Fortran results for multiple test cases

## Simulation Method: IMPES

The simulator uses IMPES (Implicit Pressure, Explicit Saturation):

1. **Calculate Well Rates** (`qrate`)
   - Determine production/injection rates based on constraints

2. **Build Pressure Matrix** (`solone` or `soltwo`)
   - Assemble 7-diagonal coefficient matrix
   - Include transmissibilities, mobilities, gravity, and capillary pressure

3. **Apply Well Boundary Conditions** (`pratei`)
   - Modify matrix for implicit well treatment

4. **Solve Pressure Equation** (`gaus1d` or `lsor`)
   - Compute new pressure distribution

5. **Calculate Final Well Rates** (`prateo`)
   - Update rates based on new pressures

6. **Update Saturations** (explicit)
   - Calculate saturation changes based on flow rates
   - Update oil, water, and gas saturations

7. **Advance Time**
   - Move to next time step
   - Repeat cycle

## Key Features

- **3D Grid**: Fully three-dimensional Cartesian grid
- **Three-Phase**: Oil, water, and gas with solution gas
- **Black Oil**: Uses standard black oil PVT correlations
- **Multiple Wells**: Supports producers and injectors
- **Well Controls**: Rate or pressure constraints with GOR/WOR limits
- **Flexible Grid**: Variable block sizes and properties
- **Rock Types**: Multiple rock regions with different properties
- **PVT Regions**: Multiple PVT regions for fluid property variations

## Code Organization

```
source/
├── main.py           # Main simulator (IMPES loop)
├── block1.py         # Aquifer, material balance, well management
├── block2.py         # Grid, solution control, interpolation
├── block3.py         # Output, visualization, post-processing
├── block4.py         # Gas properties, tridiagonal solver
├── block5.py         # Well rates and boundary conditions
├── block6.py         # PVT tables, transmissibilities, initialization
├── block7.py         # Flow equation (single-point)
├── block8.py         # Flow equation (two-point)
├── run_test.py       # Test runner with automatic inputs
├── EX1.DAT           # Example input data
└── OUT*.DAT          # Output files

Fortran originals:
├── MAIN.FOR          # Original Fortran main program
├── BLOCK1.FOR - BLOCK8.FOR
└── PARAMS.FOR        # Parameter definitions
```

## Development Notes

### Array Indexing
- **Fortran**: 1-based indexing
- **Python**: 0-based indexing
- Conversion requires careful index management

### Common Blocks
Fortran COMMON blocks converted to:
- Class attributes in `BOASTSimulator`
- Shared references passed between modules

### Data Structures
- Multi-dimensional arrays: `numpy` arrays
- Dynamic sizing based on input parameters
- Pre-allocated for performance

## Troubleshooting

### Common Issues

**Issue**: Uniform pressure/saturation (no variation between blocks)  
**Cause**: Flow arrays not properly initialized  
**Solution**: Verify transmissibility and mobility calculations

**Issue**: Simulation hangs or runs very slowly  
**Cause**: Iterative solver not converging  
**Solution**: Check convergence parameters (omega, tolerance)

**Issue**: Array index errors  
**Cause**: 0-based vs 1-based indexing mismatch  
**Solution**: Review index conversions in specific routine

## References

1. Fanchi, J.R. (1987). "BOAST II: A Three-Dimensional, Three-Phase Black Oil Applied Simulation Tool", DOE Report
2. Aziz, K. and Settari, A. (1979). "Petroleum Reservoir Simulation", Applied Science Publishers
3. Peaceman, D.W. (1977). "Fundamentals of Numerical Reservoir Simulation", Elsevier

## License

Original BOAST II was released as public domain software by the U.S. Department of Energy.

## Contributing

This is a research/educational project. Contributions should maintain compatibility with the original Fortran version for validation purposes.
