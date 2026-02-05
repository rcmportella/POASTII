"""
BLOCK2.FOR - Grid Setup, Solution Control, and Iterative Solvers
Contains:
- CODES: Read solution method and run control parameters
- GRIDSZ: Read and establish grid dimensions (DX, DY, DZ)
- INTCOM: Interpolation for compressibility
- INTERP: Linear interpolation
- INTPVT: Interpolation with bubble point handling
- LSORX: Line successive over-relaxation in X-direction
- LSORY: Line successive over-relaxation in Y-direction
- LSORZ: Line successive over-relaxation in Z-direction
"""

import numpy as np
from typing import Tuple, Optional, List
from dataclasses import dataclass


def parse_fortran_value(value_str: str) -> List[float]:
    """
    Parse Fortran-style repeated values like '5*0' or '3*1.5'
    Returns a list of values
    
    Examples:
    '5*0' -> [0, 0, 0, 0, 0]
    '3*1.5' -> [1.5, 1.5, 1.5]
    '10' -> [10]
    """
    if '*' in value_str:
        parts = value_str.split('*')
        count = int(parts[0])
        value = float(parts[1])
        return [value] * count
    else:
        return [float(value_str)]


def parse_fortran_line(line: str) -> List[float]:
    """
    Parse a line with Fortran-style repeated values
    
    Example:
    '5*0 2*1 10' -> [0, 0, 0, 0, 0, 1, 1, 10]
    """
    tokens = line.split()
    result = []
    for token in tokens:
        result.extend(parse_fortran_value(token))
    return result


@dataclass
class SolutionParameters:
    """Solution method and run control parameters"""
    ksol: int = 1          # Solution method (1=Band, 2=LSORX, 3=LSORY, 4=LSORZ)
    miter: int = 50        # Maximum iterations for LSOR
    omega: float = 1.2     # Acceleration parameter for LSOR
    tol: float = 0.1       # Maximum pressure residual
    tol1: float = 0.01     # Parameter for changing omega
    dsmax: float = 0.05    # Maximum allowed saturation change
    dpmax: float = 50.0    # Maximum allowed pressure change
    nn: int = 100          # Maximum number of time steps
    fact1: float = 1.2     # Factor for increasing time step
    fact2: float = 0.5     # Factor for decreasing time step
    tmax: float = 365.0    # Maximum simulation time (days)
    wormax: float = 0.0    # Maximum water-oil ratio
    gormax: float = 0.0    # Maximum gas-oil ratio
    pamin: float = 14.7    # Minimum average reservoir pressure
    pamax: float = 10000.0 # Maximum average reservoir pressure
    numdis: int = 0        # Discretization method (0=1-point, 1=2-point)
    irk: int = 0           # Runge-Kutta option (0=IMPES, >0=stabilized)
    thruin: float = 0.1    # Throughput fraction for stabilized IMPES


class GridSetup:
    """Grid dimension and geometry setup"""
    
    def __init__(self, simulator):
        """Initialize with reference to main simulator"""
        self.sim = simulator
        
    def gridsz(self, infile, outfile) -> Tuple[int, int, int]:
        """
        Read grid description and establish grid dimensions
        Returns (ii, jj, kk) - number of blocks in each direction
        """
        # Skip comment line
        infile.readline()
        
        outfile.write('\n\n                     ***** INITIALIZATION DATA *****\n\n')
        
        # Read grid dimensions
        line = infile.readline()
        values = parse_fortran_line(line)
        ii = int(values[0])
        jj = int(values[1])
        kk = int(values[2])
        
        # Skip comment line
        infile.readline()
        
        # Read input codes for DX, DY, DZ
        line = infile.readline()
        values = parse_fortran_line(line)
        kdx = int(values[0])
        kdy = int(values[1])
        kdz = int(values[2])
        kdznet = int(values[3])
        
        # Initialize arrays
        dx = np.zeros((ii, jj, kk))
        dy = np.zeros((ii, jj, kk))
        dz = np.zeros((ii, jj, kk))
        dznet = np.zeros((ii, jj, kk))
        el = np.zeros((ii, jj, kk))
        
        # Read DX (X-direction block lengths)
        dx = self._read_dimension(infile, outfile, ii, jj, kk, kdx, 'DX', 
                                   'GRID BLOCK LENGTH')
        
        # Read DY (Y-direction block widths)
        dy = self._read_dimension(infile, outfile, ii, jj, kk, kdy, 'DY', 
                                   'GRID BLOCK WIDTH')
        
        # Read DZ (Z-direction block thickness)
        dz = self._read_dimension(infile, outfile, ii, jj, kk, kdz, 'DZ', 
                                   'GRID BLOCK GROSS THICKNESS')
        
        # Read DZNET (net thickness)
        dznet = self._read_dimension(infile, outfile, ii, jj, kk, kdznet, 'DZNET', 
                                      'GRID BLOCK NET THICKNESS')
        
        # Grid block modifications
        infile.readline()  # Skip comment
        line = infile.readline().strip()
        values = parse_fortran_line(line)
        numdx = int(values[0])
        numdy = int(values[1])
        numdz = int(values[2])
        numdzn = int(values[3])
        idcode = int(values[4])
        
        # Apply DX modifications
        if numdx > 0:
            outfile.write('\n\n               ********GRID BLOCK LENGTH (DX) NODE MODIFICATIONS**********\n\n')
            outfile.write('               I1  I2  J1  J2  K1  K2  NEW DX VALUE\n')
            dx = self._apply_modifications(infile, outfile, dx, numdx)
            if idcode == 1:
                self._write_distribution(outfile, dx, 'DX', kk)
                
        # Apply DY modifications
        if numdy > 0:
            outfile.write('\n\n               ********GRID BLOCK WIDTH (DY) NODE MODIFICATIONS**********\n\n')
            outfile.write('               I1  I2  J1  J2  K1  K2  NEW DY VALUE\n')
            dy = self._apply_modifications(infile, outfile, dy, numdy)
            if idcode == 1:
                self._write_distribution(outfile, dy, 'DY', kk)
                
        # Apply DZ modifications
        if numdz > 0:
            outfile.write('\n\n               ********GRID BLOCK GROSS THICKNESS (DZ) NODE MODIFICATIONS**********\n\n')
            outfile.write('               I1  I2  J1  J2  K1  K2  NEW DZ VALUE\n')
            dz = self._apply_modifications(infile, outfile, dz, numdz)
            if idcode == 1:
                self._write_distribution(outfile, dz, 'DZ', kk)
                
        # Apply DZNET modifications
        if numdzn > 0:
            outfile.write('\n\n               ********GRID BLOCK NET THICKNESS (DZNET) NODE MODIFICATIONS**********\n\n')
            outfile.write('               I1  I2  J1  J2  K1  K2  NEW DZNET VALUE\n')
            dznet = self._apply_modifications(infile, outfile, dznet, numdzn)
            if idcode == 1:
                self._write_distribution(outfile, dznet, 'DZNET', kk)
                
        # Check for negative values
        self._check_grid_values(outfile, dx, dy, dz, dznet)
        
        # Establish node mid-point elevations
        el = self._establish_elevations(infile, outfile, ii, jj, kk, dz)
        
        # Store in simulator
        self.sim.dx = dx
        self.sim.dy = dy
        self.sim.dz = dz
        self.sim.dznet = dznet
        self.sim.el = el
        
        return ii, jj, kk
        
    def _read_dimension(self, infile, outfile, ii: int, jj: int, kk: int, 
                        kcode: int, name: str, description: str) -> np.ndarray:
        """Read dimension array based on input code"""
        arr = np.zeros((ii, jj, kk))
        
        # Constant value for all blocks
        if kcode < 0:
            val = float(infile.readline().split()[0])
            arr[:, :, :] = val
            outfile.write(f'\n               {description} ({name}) IS INITIALLY SET AT {val:10.4f} FOR ALL NODES\n\n')
            
        # One value per row/column/layer
        elif kcode == 0:
            if name == 'DX':
                vals = [float(x) for x in infile.readline().split()]
                for i in range(ii):
                    arr[i, :, :] = vals[i]
                    outfile.write(f'               GRID SIZE (DX) IN COLUMN {i+1:5d} IS INITIALLY SET AT {vals[i]:8.2f} FOR ALL NODES\n')
            elif name == 'DY':
                vals = [float(x) for x in infile.readline().split()]
                for j in range(jj):
                    arr[:, j, :] = vals[j]
                    outfile.write(f'               GRID SIZE (DY) IN ROW    {j+1:5d} IS INITIALLY SET AT {vals[j]:8.2f} FOR ALL NODES\n')
            elif name in ['DZ', 'DZNET']:
                vals = [float(x) for x in infile.readline().split()]
                for k in range(kk):
                    arr[:, :, k] = vals[k]
                    outfile.write(f'               GRID SIZE ({name}) IN LAYER  {k+1:5d} IS INITIALLY SET AT {vals[k]:8.2f} FOR ALL NODES\n')
                    
        # Full 2D distribution per layer
        elif kcode > 0:
            outfile.write(f'\n\n               ********{description} ({name}) DISTRIBUTION********\n\n')
            for k in range(kk):
                outfile.write(f'\nK = {k+1:2d}\n\n')
                for j in range(jj):
                    line = infile.readline()
                    vals = parse_fortran_line(line)
                    arr[:, j, k] = vals[:ii]
                    outfile.write(' ' + ' '.join(f'{v:8.0f}' for v in vals[:ii]) + '\n')
                # Copy first layer to others if needed
                if k == 0 and kk > 1:
                    for kp in range(1, kk):
                        arr[:, :, kp] = arr[:, :, 0]
                        
        return arr
        
    def _apply_modifications(self, infile, outfile, arr: np.ndarray, 
                             num_mods: int) -> np.ndarray:
        """Apply regional modifications to dimension array"""
        for n in range(num_mods):
            line = infile.readline()
            values = parse_fortran_line(line)
            i1, i2 = int(values[0]), int(values[1])
            j1, j2 = int(values[2]), int(values[3])
            k1, k2 = int(values[4]), int(values[5])
            regval = values[6]
            
            outfile.write(f'               {i1:4d}{i2:4d}{j1:4d}{j2:4d}{k1:4d}{k2:4d}  {regval:10.4e}\n')
            
            # Convert to 0-based indexing
            arr[i1-1:i2, j1-1:j2, k1-1:k2] = regval
            
        return arr
        
    def _write_distribution(self, outfile, arr: np.ndarray, name: str, kk: int):
        """Write full distribution of dimension array"""
        ii, jj, _ = arr.shape
        outfile.write(f'\n\n               ********GRID BLOCK {name} DISTRIBUTION********\n\n')
        for k in range(kk):
            outfile.write(f'\nK = {k+1:2d}\n\n')
            for j in range(jj):
                outfile.write(' ' + ' '.join(f'{arr[i,j,k]:8.0f}' for i in range(ii)) + '\n')
                
    def _check_grid_values(self, outfile, dx, dy, dz, dznet):
        """Check for negative grid values"""
        ii, jj, kk = dx.shape
        errors = 0
        
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    if dx[i, j, k] < 0.0:
                        outfile.write(f'\n     -----GRID BLOCK DX ERROR AT IJK = {i+1:5d}{j+1:5d}{k+1:5d}\n')
                        errors += 1
                    if dy[i, j, k] < 0.0:
                        outfile.write(f'\n     -----GRID BLOCK DY ERROR AT IJK = {i+1:5d}{j+1:5d}{k+1:5d}\n')
                        errors += 1
                    if dz[i, j, k] < 0.0:
                        outfile.write(f'\n     -----GRID BLOCK DZ ERROR AT IJK = {i+1:5d}{j+1:5d}{k+1:5d}\n')
                        errors += 1
                    if dznet[i, j, k] < 0.0:
                        outfile.write(f'\n     -----GRID BLOCK DZNET ERROR AT IJK = {i+1:5d}{j+1:5d}{k+1:5d}\n')
                        errors += 1
                        
        if errors > 0:
            self.sim.ifatal += errors
            
    def _establish_elevations(self, infile, outfile, ii: int, jj: int, 
                              kk: int, dz: np.ndarray) -> np.ndarray:
        """Establish node mid-point elevations"""
        el = np.zeros((ii, jj, kk))
        
        infile.readline()  # Skip comment
        kel = int(infile.readline().split()[0])
        
        if kel <= 1:
            # Read top surface elevations
            if kel == 0:
                elev = float(infile.readline().split()[0])
                varel = np.full((ii, jj), elev)
            else:  # kel == 1
                varel = np.zeros((ii, jj))
                for j in range(jj):
                    vals = [float(x) for x in infile.readline().split()]
                    varel[:, j] = vals[:ii]
                    
            # Calculate depths for each layer
            for k in range(kk):
                if k == 0:
                    el[:, :, k] = varel
                else:
                    cumulative_dz = np.sum(dz[:, :, :k], axis=2)
                    el[:, :, k] = varel + cumulative_dz
                    
        elif kel == 2:
            # One elevation per layer
            elevations = [float(x) for x in infile.readline().split()]
            for k in range(kk):
                el[:, :, k] = elevations[k]
                
        else:  # kel == 3
            # Full 2D elevation per layer
            for k in range(kk):
                for j in range(jj):
                    vals = [float(x) for x in infile.readline().split()]
                    el[:, j, k] = vals[:ii]
                    
        # Write depths to grid block tops
        outfile.write('\n\n               ********** DEPTHS TO GRID BLOCK TOPS **********\n\n')
        for k in range(kk):
            outfile.write(f'\nK = {k+1:2d}\n\n')
            for j in range(jj):
                outfile.write(' ' + ' '.join(f'{el[i,j,k]:8.0f}' for i in range(ii)) + '\n')
                
        # Adjust to mid-point elevations
        el += dz * 0.5
        
        return el


class SolutionControl:
    """Solution method and run control parameters"""
    
    def __init__(self, simulator):
        """Initialize with reference to main simulator"""
        self.sim = simulator
    
    @staticmethod
    def codes(infile, outfile, simulator) -> SolutionParameters:
        """
        Read solution method, debug controls, and run parameters
        Returns SolutionParameters object
        """
        params = SolutionParameters()
        
        # Read debug and output control
        infile.readline()  # Skip comment
        line = infile.readline()
        values = parse_fortran_line(line)
        while len(values) < 4:
            values.append(0.0)
        ksn1 = int(values[0])  # LSOR debug output frequency
        ksm1 = int(values[1])  # Solution method debug frequency
        kco1 = int(values[2])  # Compressibility debug frequency
        kcoff = int(values[3]) # Density/saturation debug
        
        # Read run control parameters
        infile.readline()  # Skip comment
        line = infile.readline()
        values = parse_fortran_line(line)
        while len(values) < 8:
            values.append(0.0)
        params.nn = int(values[0])
        params.fact1 = values[1]
        params.fact2 = values[2]
        params.tmax = values[3]
        params.wormax = values[4]
        params.gormax = values[5]
        params.pamin = values[6]
        params.pamax = values[7]
        
        # Write run control parameters
        outfile.write('\n               RUN CONTROL PARAMETERS:\n')
        outfile.write(f'                    MAXIMUM NUMBER OF TIME-STEPS = {params.nn}\n')
        outfile.write(f'                    FACTOR FOR INCREASING DELT = {params.fact1:10.3f}   WHEN DSMAX AND DPMAX NOT EXCEEDED.\n')
        outfile.write(f'                    FACTOR FOR DECREASING DELT = {params.fact2:10.3f}   WHEN DSMAX OR DPMAX IS EXCEEDED.\n')
        outfile.write(f'                    MAXIMUM SIMULATION TIME = {params.tmax:11.3f}\n')
        outfile.write(f'                    MAXIMUM RESERVOIR WOR/TIME-STEP = {params.wormax:8.0f} STB/STB\n')
        outfile.write(f'                    MAXIMUM RESERVOIR GOR/TIME-STEP = {params.gormax:8.0f} SCF/STB\n')
        outfile.write(f'                    MINIMUM AVERAGE RESERVOIR PRESSURE/TIME-STEP = {params.pamin:8.0f}\n')
        outfile.write(f'                    MAXIMUM AVERAGE RESERVOIR PRESSURE/TIME-STEP = {params.pamax:8.0f}\n')
        
        # Read rock region WOR/GOR maxima if needed
        if params.wormax == 0.0:
            line = infile.readline()
            worock = parse_fortran_line(line)
            outfile.write('\n\n               ROCK REGION SPECIFIED WOR MAXIMA:\n')
            for i, wor in enumerate(worock):
                outfile.write(f'                    ROCK REGION {i+1:3d} WOR MAX (STB/STB) = {wor:8.0f}\n')
                
        if params.gormax == 0.0:
            line = infile.readline()
            gorock = parse_fortran_line(line)
            outfile.write('\n\n               ROCK REGION SPECIFIED GOR MAXIMA:\n')
            for i, gor in enumerate(gorock):
                outfile.write(f'                    ROCK REGION {i+1:3d} GOR MAX (SCF/STB) = {gor:8.0f}\n')
                
        # Read solution method parameters
        infile.readline()  # Skip comment
        line = infile.readline()
        values = parse_fortran_line(line)
        while len(values) < 7:
            values.append(0.0)
        params.ksol = int(values[0])
        params.miter = int(values[1])
        params.omega = values[2]
        params.tol = values[3]
        params.tol1 = values[4]
        params.dsmax = values[5]
        params.dpmax = values[6]
        
        # Write solution method
        if params.ksol == 1:
            outfile.write('\n\n               SOLUTION METHOD IS BAND.\n')
        elif params.ksol == 2:
            outfile.write('\n\n               SOLUTION METHOD IS LSORX:\n')
            outfile.write(f'                    MAXIMUM NUMBER OF ITERATIONS      (MITER) =      {params.miter:5d}\n')
            outfile.write(f'                    INITIAL ACCELERATION PARAMETER    (OMEGA) = {params.omega:10.4f}\n')
            outfile.write(f'                    MAXIMUM PRESSURE RESIDUAL           (TOL) = {params.tol:10.4f}\n')
            outfile.write(f'                    PARAMETER FOR CHANGING OMEGA       (TOL1) = {params.tol1:10.4f}\n')
        elif params.ksol == 3:
            outfile.write('\n\n               SOLUTION METHOD IS LSORY:\n')
            outfile.write(f'                    MAXIMUM NUMBER OF ITERATIONS      (MITER) =      {params.miter:5d}\n')
            outfile.write(f'                    INITIAL ACCELERATION PARAMETER    (OMEGA) = {params.omega:10.4f}\n')
            outfile.write(f'                    MAXIMUM PRESSURE RESIDUAL           (TOL) = {params.tol:10.4f}\n')
            outfile.write(f'                    PARAMETER FOR CHANGING OMEGA       (TOL1) = {params.tol1:10.4f}\n')
        elif params.ksol == 4:
            outfile.write('\n\n               SOLUTION METHOD IS LSORZ:\n')
            outfile.write(f'                    MAXIMUM NUMBER OF ITERATIONS      (MITER) =      {params.miter:5d}\n')
            outfile.write(f'                    INITIAL ACCELERATION PARAMETER    (OMEGA) = {params.omega:10.4f}\n')
            outfile.write(f'                    MAXIMUM PRESSURE RESIDUAL           (TOL) = {params.tol:10.4f}\n')
            outfile.write(f'                    PARAMETER FOR CHANGING OMEGA       (TOL1) = {params.tol1:10.4f}\n')
        else:
            outfile.write(f'\nIMPROPER SOLUTION METHOD SPECIFIED, KSOL={params.ksol}\n')
            simulator.ifatal += 1
            
        outfile.write('\n               AUTOMATIC TIME STEP CRITERIA:\n')
        outfile.write(f'                    MAXIMUM ALLOWED SATURATION CHANGE (DSMAX) = {params.dsmax:10.4f}\n')
        outfile.write(f'                    MAXIMUM ALLOWED PRESSURE CHANGE   (DPMAX) = {params.dpmax:10.4f}\n\n')
        
        # Read IMPES/Runge-Kutta controls
        infile.readline()  # Skip comment
        line = infile.readline()
        values = parse_fortran_line(line)
        while len(values) < 3:
            values.append(0.0)
        params.numdis = int(values[0])
        params.irk = int(values[1])
        params.thruin = values[2]
        
        if params.irk == 0 and params.numdis == 0:
            outfile.write('\n\n               IMPES FORMULATION SELECTED;  SINGLE-POINT UPSTREAM WEIGHTING.\n')
        elif params.irk == 0 and params.numdis == 1:
            outfile.write('\n\n               IMPES FORMULATION SELECTED;  TWO-POINT UPSTREAM WEIGHTING.\n')
        elif params.irk > 0 and params.numdis == 0:
            outfile.write('\n\n               STABILISED IMPES FORMULATION:\n')
            outfile.write('                    SINGLE-POINT UPSTREAM WEIGHTING\n')
            outfile.write(f'                    USER SPEC THROUGHPUT, FRACTION  (THRUIN) =     {params.thruin:10.4f}\n\n')
        elif params.irk > 0 and params.numdis == 1:
            outfile.write('\n\n               STABILISED IMPES FORMULATION:\n')
            outfile.write('                    TWO-POINT UPSTREAM WEIGHTING\n')
            outfile.write(f'                    USER SPEC THROUGHPUT, FRACTION  (THRUIN) =     {params.thruin:10.4f}\n\n')
            
        # Check for warnings and errors
        SolutionControl._check_parameters(outfile, simulator, params, 
                                          ksn1, ksm1, kco1, kcoff)
        
        return params
        
    @staticmethod
    def _check_parameters(outfile, simulator, params: SolutionParameters,
                         ksn1: int, ksm1: int, kco1: int, kcoff: int):
        """Check parameters for errors and warnings"""
        # Debug warnings
        if ksn1 != 0:
            simulator.iwarn += 1
            outfile.write('\n     -----LSOR DEBUG OUTPUT ON\n')
        if ksm1 != 0:
            simulator.iwarn += 1
            outfile.write('\n     -----SOLN METHOD DEBUG OUTPUT ON\n')
        if kco1 != 0:
            simulator.iwarn += 1
            outfile.write('\n     -----COMP AND FVF DEBUG OUTPUT ON\n')
        if kcoff != 0:
            simulator.iwarn += 1
            outfile.write('\n     -----DEN AND SAT DEBUG OUTPUT ON\n')
            
        # Run control errors
        if params.nn < 1:
            simulator.ifatal += 1
            outfile.write('\n     -----MAX # OF TIME STEPS ERROR\n')
        if params.fact1 < 1.0:
            simulator.ifatal += 1
            outfile.write('\n     -----FACT1 ERROR\n')
        if params.fact2 <= 0.0 or params.fact2 > 1.0:
            simulator.ifatal += 1
            outfile.write('\n     -----FACT2 ERROR\n')
        if params.tmax <= 0.0:
            simulator.ifatal += 1
            outfile.write('\n     -----TMAX ERROR\n')
        if params.wormax < 0.0:
            simulator.ifatal += 1
            outfile.write('\n     -----WORMAX ERROR\n')
        if params.gormax < 0.0:
            simulator.ifatal += 1
            outfile.write('\n     -----GORMAX ERROR\n')
        if params.pamin <= 0.0:
            simulator.ifatal += 1
            outfile.write('\n     -----PAMIN ERROR\n')
        if params.pamax <= params.pamin:
            simulator.ifatal += 1
            outfile.write('\n     -----PAMAX ERROR\n')
            
        # Solution method warnings
        if params.miter < 1:
            simulator.iwarn += 1
            outfile.write('\n     -----ALLOWED # OF LSOR ITER LESS THAN 1\n')
        if params.omega < 1.0 or params.omega > 2.0:
            simulator.iwarn += 1
            outfile.write('\n     -----OMEGA LT 1 OR OMEGA GT 2\n')
        if params.tol <= 0:
            simulator.iwarn += 1
            outfile.write('\n     -----TOL LESS THAN OR EQUAL TO 0\n')
        if params.tol1 < 0.0:
            simulator.iwarn += 1
            outfile.write('\n     -----TOL1 IS NEGATIVE\n')
        if params.dsmax <= 0.0 or params.dsmax > 1.0:
            simulator.ifatal += 1
            outfile.write('\n     -----DSMAX ERROR\n')
        if params.dpmax <= 0.0:
            simulator.ifatal += 1
            outfile.write('\n     -----DPMAX ERROR\n')


class Interpolation:
    """Interpolation utilities"""
    
    def __init__(self, simulator):
        """Initialize with reference to main simulator"""
        self.sim = simulator
    
    @staticmethod
    def intcom(x_table: np.ndarray, y_table: np.ndarray, ireg: int, 
               n: int, xo: float) -> float:
        """
        Interpolation for compressibility (step function)
        If xo >= x[n], return y[n]
        Otherwise return first y[i] where x[i] > xo
        """
        if xo >= x_table[ireg, n-1]:
            return y_table[ireg, n-1]
            
        for i in range(1, n):
            if xo <= x_table[ireg, i]:
                return y_table[ireg, i]
                
        return y_table[ireg, n-1]
        
    @staticmethod
    def interp(x_table: np.ndarray, y_table: np.ndarray, *args) -> float:
        """
        Linear interpolation - supports two calling conventions:
        1. 1D arrays: interp(x_array, y_array, xo)
        2. 2D arrays: interp(x_table, y_table, ireg, n, xo)
        """
        if len(args) == 1:
            # 1D case: interp(x_array, y_array, xo)
            xo = args[0]
            x_arr = x_table.ravel()  # Ensure 1D
            y_arr = y_table.ravel()
            n = len(x_arr)
            
            # Handle empty arrays
            if n == 0:
                return 0.0
            
            if xo >= x_arr[n-1]:
                return y_arr[n-1]
                
            for i in range(1, n):
                if xo < x_arr[i]:
                    # Linear interpolation
                    x1, x2 = x_arr[i-1], x_arr[i]
                    y1, y2 = y_arr[i-1], y_arr[i]
                    if abs(x2 - x1) < 1e-10:
                        return y1
                    yo = y1 + (xo - x1) * (y2 - y1) / (x2 - x1)
                    return yo
                    
            return y_arr[n-1]
            
        elif len(args) == 3:
            # 2D case: interp(x_table, y_table, ireg, n, xo)
            ireg = args[0]
            n = args[1]
            xo = args[2]
            
            # Handle invalid n
            if n <= 0:
                return 0.0
            
            if xo >= x_table[ireg, n-1]:
                return y_table[ireg, n-1]
                
            for i in range(1, n):
                if xo < x_table[ireg, i]:
                    # Linear interpolation
                    x1, x2 = x_table[ireg, i-1], x_table[ireg, i]
                    y1, y2 = y_table[ireg, i-1], y_table[ireg, i]
                    if abs(x2 - x1) < 1e-10:
                        return y1
                    yo = y1 + (xo - x1) * (y2 - y1) / (x2 - x1)
                    return yo
                    
            return y_table[ireg, n-1]
        else:
            raise ValueError(f"interp() takes either 3 or 5 arguments, got {len(args)+2}")
        
    @staticmethod
    def intpvt(x_table: np.ndarray, y_table: np.ndarray, ireg: int, 
               n: int, xo: float, bpt: float, rm: float) -> float:
        """
        Interpolation with bubble point handling
        Below bubble point: linear interpolation
        Above bubble point: extrapolate using slope rm
        """
        if xo <= bpt:
            # Below bubble point - use linear interpolation
            return Interpolation.interp(x_table, y_table, ireg, n, xo)
        else:
            # Above bubble point - extrapolate
            # First find value at bubble point
            yobp = Interpolation.interp(x_table, y_table, ireg, n, bpt)
            # Extrapolate linearly
            yo = yobp + rm * (xo - bpt)
            return yo


class LSORSolver:
    """Line Successive Over-Relaxation solvers"""
    
    def __init__(self, simulator):
        """Initialize with reference to main simulator"""
        self.sim = simulator
    
    @staticmethod
    def lsorx(simulator, nx: int, ny: int, nz: int, omega: float, 
              tol: float, tol1: float, miter: int, delt: float, 
              delt0: float, ksn: int, n: int) -> Tuple[float, int]:
        """
        Iterative solution: X-direction tridiagonal algorithm
        Returns (final omega, number of iterations)
        """
        niter = 0
        dmax = 1.0
        rho1 = 0.0
        theta = 0.0
        
        while True:
            tw = 1.0 - omega
            dmax0 = dmax
            theta0 = theta
            
            if niter >= miter:
                simulator.outfile.write(f'               CONVERGENCE(LSORX) WAS NOT REACHED IN {niter:5d} ITERATIONS\n')
                simulator.outfile.write(f'               TOL = {tol:10.7f}          DMAX = {dmax:15.7f}\n')
                return omega, niter
                
            niter += 1
            dmax = 0.0
            
            # Solve line by line in X-direction
            for k in range(nz):
                for j in range(ny):
                    # Build tridiagonal system
                    azl = simulator.aw[:nx, j, k]
                    bzl = simulator.e[:nx, j, k]
                    czl = simulator.ae[:nx, j, k]
                    dzl = simulator.b[:nx, j, k].copy()
                    
                    # Subtract Y-direction terms
                    if ny > 1:
                        jm = max(0, j-1)
                        jp = min(ny-1, j+1)
                        dzl -= simulator.as_[:nx, j, k] * simulator.p[:nx, jm, k]
                        dzl -= simulator.an[:nx, j, k] * simulator.p[:nx, jp, k]
                        
                    # Subtract Z-direction terms
                    if nz > 1:
                        km = max(0, k-1)
                        kp = min(nz-1, k+1)
                        dzl -= simulator.at[:nx, j, k] * simulator.p[:nx, j, km]
                        dzl -= simulator.ab[:nx, j, k] * simulator.p[:nx, j, kp]
                        
                    # Solve tridiagonal system
                    uzl = LSORSolver._solve_tridiagonal(azl, bzl, czl, dzl, nx)
                    
                    # Update pressure with over-relaxation
                    um = simulator.p[:nx, j, k].copy()
                    simulator.p[:nx, j, k] = tw * um + omega * uzl
                    
                    # Track maximum change
                    dm = np.max(np.abs(simulator.p[:nx, j, k] - um))
                    if dm > dmax:
                        dmax = dm
                        
            # Update omega if appropriate
            if tol1 != 0.0 and niter > 1:
                theta = dmax / dmax0
                delta = abs(theta - theta0)
                if delta <= tol1:
                    om = omega - 1.0
                    rho1 = (theta + om) ** 2 / (theta * omega ** 2)
                    if rho1 < 1.0:
                        omega = 2.0 / (1.0 + np.sqrt(1.0 - rho1))
                        
            # Check convergence
            if dmax <= tol:
                if n == ksn:
                    simulator.outfile.write(f'     CONVERGENCE(LSORX) HAS BEEN REACHED AFTER {niter:3d} ITERATIONS     OMEGA = {omega:6.3f}\n')
                    simulator.outfile.write(f'     DMAX = {dmax:10.6f}     THETA = {theta:10.6f}     RHO1 = {rho1:10.6f}\n\n')
                return omega, niter
                
    @staticmethod
    def lsory(simulator, nx: int, ny: int, nz: int, omega: float, 
              tol: float, tol1: float, miter: int, delt: float, 
              delt0: float, ksn: int, n: int) -> Tuple[float, int]:
        """
        Iterative solution: Y-direction tridiagonal algorithm
        Returns (final omega, number of iterations)
        """
        # Similar implementation to lsorx but solving in Y-direction
        # Implementation follows same pattern with different loop order
        return omega, 0  # Placeholder
        
    @staticmethod
    def lsorz(simulator, nx: int, ny: int, nz: int, omega: float, 
              tol: float, tol1: float, miter: int, delt: float, 
              delt0: float, ksn: int, n: int) -> Tuple[float, int]:
        """
        Iterative solution: Z-direction tridiagonal algorithm
        Returns (final omega, number of iterations)
        """
        # Similar implementation to lsorx but solving in Z-direction
        # Implementation follows same pattern with different loop order
        return omega, 0  # Placeholder
        
    @staticmethod
    def _solve_tridiagonal(a: np.ndarray, b: np.ndarray, c: np.ndarray, 
                          d: np.ndarray, n: int) -> np.ndarray:
        """
        Solve tridiagonal system: a[i]*u[i-1] + b[i]*u[i] + c[i]*u[i+1] = d[i]
        Using Thomas algorithm
        """
        # Forward elimination
        beta = np.zeros(n)
        gamma = np.zeros(n)
        w = np.zeros(n)
        
        beta[0] = b[0]
        gamma[0] = d[0] / b[0]
        
        for i in range(n-1):
            w[i] = c[i] / beta[i]
            beta[i+1] = b[i+1] - a[i+1] * w[i]
            gamma[i+1] = (d[i+1] - a[i+1] * gamma[i]) / beta[i+1]
            
        # Back substitution
        u = np.zeros(n)
        u[n-1] = gamma[n-1]
        
        for i in range(n-2, -1, -1):
            u[i] = gamma[i] - w[i] * u[i+1]
            
        return u
