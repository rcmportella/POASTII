"""
BLOCK3.FOR - Output Visualization, Rock Properties, and Post-Processing

This module contains:
1. OutputVisualization class: Digital contour plotting and line printer plots
2. RockProperties class: Porosity and permeability input/validation
3. PostProcessing class: Summary reports and time series plots

Converted from BOAST II (Release 1.2) Fortran code
"""

import numpy as np
from typing import TextIO, List, Tuple
from dataclasses import dataclass


class OutputVisualization:
    """Handles visualization output including contour plots and line printer plots"""
    
    def __init__(self, simulator):
        """Initialize with reference to main simulator"""
        self.sim = simulator
    
    @staticmethod
    def plot(arr: np.ndarray, ii: int, jj: int, kk: int, iocode: TextIO, 
             iplot: int) -> None:
        """
        Digital contour plotting routine for 3D arrays
        
        Parameters:
        -----------
        arr : numpy array
            Array to plot (P, SO, SW, SG, PBOT, or CUMAQW)
        ii, jj, kk : int
            Grid dimensions
        iplot : int
            Plot type (1=P, 2=SO, 3=SW, 4=SG, 5=PBOT, 6=CUMAQW)
        iocode : file object
            Output file handle
        """
        # Determine plot title and type
        titles = {
            1: '***** RESERVOIR PRESSURE (PSIA) *****',
            2: '***** OIL SATURATION *****',
            3: '***** WATER SATURATION *****',
            4: '***** GAS SATURATION *****',
            5: '***** BUBBLE POINT PRESSURE (PSIA) *****',
            6: '***** CUMULATIVE AQUIFER INFLUX (STB) *****'
        }
        
        iocode.write(f"\n\n{'':>14}{titles.get(iplot, 'UNKNOWN')}\n")
        
        # Process each layer
        for k in range(kk):
            iocode.write(f"\n   LAYER = {k+1}\n\n")
            
            # Find min and max values for this layer
            vmin = float('inf')
            vmax = float('-inf')
            
            for j in range(jj):
                for i in range(ii):
                    val = arr[i, j, k]
                    if val < vmin:
                        vmin = val
                    if val > vmax:
                        vmax = val
            
            # Calculate contour intervals (11 levels: -, 1-9, T)
            if vmax > vmin:
                dv = (vmax - vmin) / 10.0
            else:
                dv = 1.0
            
            # Create contour levels
            levels = [vmin + i * dv for i in range(11)]
            
            # Print scale
            iocode.write("   CONTOUR SCALE:\n")
            iocode.write(f"   -  : {levels[0]:10.3f} to {levels[1]:10.3f}\n")
            for i in range(1, 10):
                iocode.write(f"   {i}  : {levels[i]:10.3f} to {levels[i+1]:10.3f}\n")
            iocode.write(f"   T  : {levels[10]:10.3f} (max)\n\n")
            
            # Print contour map (Y is vertical, X is horizontal)
            for j in range(jj-1, -1, -1):  # Print from top to bottom
                line = f"   J={j+1:2d}  "
                for i in range(ii):
                    val = arr[i, j, k]
                    # Determine contour character
                    if val <= levels[0]:
                        char = '-'
                    elif val >= levels[10]:
                        char = 'T'
                    else:
                        for lev in range(1, 10):
                            if levels[lev] <= val < levels[lev+1]:
                                char = str(lev)
                                break
                        else:
                            char = '9'
                    line += char
                iocode.write(line + "\n")
            
            # Print X axis labels
            iocode.write("        ")
            for i in range(ii):
                iocode.write(str((i+1) % 10))
            iocode.write("\n")
            iocode.write(f"        I = 1 to {ii}\n\n")
    
    @staticmethod
    def ploti(title: str, xdata: List[float], ydata: List[float], 
              npts: int, iocode: TextIO, xlabel: str = "TIME",
              ylabel: str = "VALUE") -> None:
        """
        Line printer plot for time series data
        
        Parameters:
        -----------
        title : str
            Plot title
        xdata, ydata : list of float
            Data arrays to plot
        npts : int
            Number of data points
        iocode : file object
            Output file handle
        xlabel, ylabel : str
            Axis labels
        """
        if npts < 1:
            return
        
        # Find data ranges
        xmin = min(xdata[:npts])
        xmax = max(xdata[:npts])
        ymin = min(ydata[:npts])
        ymax = max(ydata[:npts])
        
        # Handle constant data
        if abs(ymax - ymin) < 1.0e-10:
            ymin = ymax - 1.0
            ymax = ymax + 1.0
        
        # Set up plot dimensions (101 characters wide, 51 lines high)
        nwidth = 101
        nheight = 51
        
        # Create plot grid
        grid = [[' ' for _ in range(nwidth)] for _ in range(nheight)]
        
        # Draw axes
        for i in range(nheight):
            grid[i][0] = '|'
        for i in range(nwidth):
            grid[nheight-1][i] = '-'
        grid[nheight-1][0] = '+'
        
        # Scale and plot data
        dy = (ymax - ymin) / (nheight - 2)
        dx = (xmax - xmin) / (nwidth - 2) if xmax > xmin else 1.0
        
        for n in range(npts):
            ix = int((xdata[n] - xmin) / dx) + 1
            iy = int((ydata[n] - ymin) / dy)
            
            if 1 <= ix < nwidth and 0 <= iy < nheight-1:
                grid[nheight-2-iy][ix] = '*'
        
        # Print title
        iocode.write(f"\n\n{title:^101s}\n")
        iocode.write(f"{ylabel:^10s}\n")
        
        # Print grid with Y-axis labels
        for i in range(nheight):
            y_val = ymin + (nheight - 1 - i) * dy
            if i % 5 == 0:
                iocode.write(f"{y_val:10.2e} ")
            else:
                iocode.write(f"{'':11s}")
            iocode.write(''.join(grid[i]) + "\n")
        
        # Print X-axis label
        iocode.write(f"{'':>11s}{xlabel:^101s}\n")
        iocode.write(f"{'':>11s}{xmin:10.2e}{'':>80s}{xmax:10.2e}\n\n")


class RockProperties:
    """Handles porosity and permeability input and validation"""
    
    def __init__(self, simulator):
        """
        Initialize with reference to main simulator
        
        Parameters:
        -----------
        simulator : BOASTSimulator
            Reference to main simulator object
        """
        self.sim = simulator
    
    def porprm(self) -> None:
        """
        Read and validate porosity and permeability distributions
        
        Reads PHI (porosity) and KX, KY, KZ (permeability) with multiple input modes:
        - Code -1: Constant value (value on next line)
        - Code 0: Constant value (value on same line)
        - Code 1: Values by layer
        - Code 2: Full 3D distribution
        
        Validates: 0 <= PHI <= 1, K >= 0
        """
        ii, jj, kk = self.sim.ii, self.sim.jj, self.sim.kk
        iread = self.sim.iread
        iocode = self.sim.iocode
        
        # Read codes for PHI, KX, KY, KZ
        iread.readline()  # Skip comment line "POROSITY AND PERMEABILITY"
        line = iread.readline()
        codes = [int(x) for x in line.split()]
        iphi = codes[0]
        ikx = codes[1]
        iky = codes[2]
        ikz = codes[3]
        
        # Read porosity
        iocode.write("\n   READING POROSITY DISTRIBUTION...\n")
        
        if iphi == -1 or iphi == 0:
            # Constant porosity
            if iphi == -1:
                line = iread.readline()
                phi_const = float(line.strip())
            else:
                phi_const = float(line.split()[1])
            self.sim.phi.fill(phi_const)
            iocode.write(f"   CONSTANT POROSITY = {phi_const:.4f}\n")
        
        elif iphi == 1:
            # Porosity by layer
            iocode.write("   POROSITY BY LAYER:\n")
            for k in range(kk):
                line = iread.readline()
                phi_k = float(line.strip())
                self.sim.phi[:, :, k] = phi_k
                iocode.write(f"   LAYER {k+1}: PHI = {phi_k:.4f}\n")
        
        elif iphi == 2:
            # Full 3D porosity distribution
            iocode.write("   READING FULL 3D POROSITY DISTRIBUTION\n")
            for k in range(kk):
                for j in range(jj):
                    line = iread.readline()
                    values = [float(x) for x in line.split()]
                    for i in range(ii):
                        self.sim.phi[i, j, k] = values[i]
        
        # Validate porosity
        ierr = 0
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    phi = self.sim.phi[i, j, k]
                    if phi < 0.0 or phi > 1.0:
                        iocode.write(f"   ERROR: PHI({i+1},{j+1},{k+1}) = {phi:.4f} OUT OF RANGE [0,1]\n")
                        ierr += 1
        
        if ierr > 0:
            raise ValueError(f"POROSITY VALIDATION FAILED: {ierr} errors found")
        
        # Read X-direction permeability
        iocode.write("\n   READING X-DIRECTION PERMEABILITY...\n")
        
        if ikx == -1 or ikx == 0:
            # Constant KX
            if ikx == -1:
                line = iread.readline()
                kx_const = float(line.strip())
            else:
                line = iread.readline()
                kx_const = float(line.split()[1])
            self.sim.kx.fill(kx_const)
            iocode.write(f"   CONSTANT KX = {kx_const:.2f} md\n")
        
        elif ikx == 1:
            iocode.write("   KX BY LAYER:\n")
            for k in range(kk):
                line = iread.readline()
                kx_k = float(line.strip())
                self.sim.kx[:, :, k] = kx_k
                iocode.write(f"   LAYER {k+1}: KX = {kx_k:.2f} md\n")
        
        elif ikx == 2:
            iocode.write("   READING FULL 3D KX DISTRIBUTION\n")
            for k in range(kk):
                for j in range(jj):
                    line = iread.readline()
                    values = [float(x) for x in line.split()]
                    for i in range(ii):
                        self.sim.kx[i, j, k] = values[i]
        
        # Read Y-direction permeability
        iocode.write("\n   READING Y-DIRECTION PERMEABILITY...\n")
        
        if iky == -1 or iky == 0:
            # Constant KY
            if iky == -1:
                line = iread.readline()
                ky_const = float(line.strip())
            else:
                line = iread.readline()
                ky_const = float(line.split()[1])
            self.sim.ky.fill(ky_const)
            iocode.write(f"   CONSTANT KY = {ky_const:.2f} md\n")
        
        elif iky == 1:
            iocode.write("   KY BY LAYER:\n")
            for k in range(kk):
                line = iread.readline()
                ky_k = float(line.strip())
                self.sim.ky[:, :, k] = ky_k
                iocode.write(f"   LAYER {k+1}: KY = {ky_k:.2f} md\n")
        
        elif iky == 2:
            iocode.write("   READING FULL 3D KY DISTRIBUTION\n")
            for k in range(kk):
                for j in range(jj):
                    line = iread.readline()
                    values = [float(x) for x in line.split()]
                    for i in range(ii):
                        self.sim.ky[i, j, k] = values[i]
        
        # Read Z-direction permeability
        iocode.write("\n   READING Z-DIRECTION PERMEABILITY...\n")
        
        if ikz == -1 or ikz == 0:
            # Constant KZ
            if ikz == -1:
                line = iread.readline()
                kz_const = float(line.strip())
            else:
                line = iread.readline()
                kz_const = float(line.split()[1])
            self.sim.kz.fill(kz_const)
            iocode.write(f"   CONSTANT KZ = {kz_const:.2f} md\n")
        
        elif ikz == 1:
            iocode.write("   KZ BY LAYER:\n")
            for k in range(kk):
                line = iread.readline()
                kz_k = float(line.strip())
                self.sim.kz[:, :, k] = kz_k
                iocode.write(f"   LAYER {k+1}: KZ = {kz_k:.2f} md\n")
        
        elif ikz == 2:
            iocode.write("   READING FULL 3D KZ DISTRIBUTION\n")
            for k in range(kk):
                for j in range(jj):
                    line = iread.readline()
                    values = [float(x) for x in line.split()]
                    for i in range(ii):
                        self.sim.kz[i, j, k] = values[i]
        
        # Validate permeability
        ierr = 0
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    if self.sim.kx[i, j, k] < 0.0:
                        iocode.write(f"   ERROR: KX({i+1},{j+1},{k+1}) < 0\n")
                        ierr += 1
                    if self.sim.ky[i, j, k] < 0.0:
                        iocode.write(f"   ERROR: KY({i+1},{j+1},{k+1}) < 0\n")
                        ierr += 1
                    if self.sim.kz[i, j, k] < 0.0:
                        iocode.write(f"   ERROR: KZ({i+1},{j+1},{k+1}) < 0\n")
                        ierr += 1
        
        if ierr > 0:
            raise ValueError(f"PERMEABILITY VALIDATION FAILED: {ierr} errors found")
        
        # Read porosity and permeability modifications
        iread.readline()  # Skip comment "POROSITY AND PERMEABILITY MODIFICATIONS"
        line = iread.readline()
        from block2 import parse_fortran_line
        values = parse_fortran_line(line)
        while len(values) < 5:
            values.append(0.0)
        numphi = int(values[0])
        numkx = int(values[1])
        numky = int(values[2])
        numkz = int(values[3])
        idcode = int(values[4])
        
        # Apply PHI modifications
        if numphi > 0:
            iocode.write("\n\n               ********POROSITY NODE MODIFICATIONS**********\n\n")
            iocode.write("               I1  I2  J1  J2  K1  K2  NEW PHI VALUE\n")
            
            for n in range(numphi):
                line = iread.readline()
                values = parse_fortran_line(line)
                while len(values) < 7:
                    values.append(0.0)
                i1, i2, j1, j2, k1, k2 = [int(x) for x in values[:6]]
                phi_new = values[6]
                
                iocode.write(f"               {i1:4d}{i2:4d}{j1:4d}{j2:4d}{k1:4d}{k2:4d}  {phi_new:10.4f}\n")
                
                for k in range(k1-1, k2):
                    for j in range(j1-1, j2):
                        for i in range(i1-1, i2):
                            self.sim.phi[i, j, k] = phi_new
        
        # Apply KX modifications
        if numkx > 0:
            iocode.write("\n\n               ********KX NODE MODIFICATIONS**********\n\n")
            iocode.write("               I1  I2  J1  J2  K1  K2  NEW KX VALUE\n")
            
            for n in range(numkx):
                line = iread.readline()
                values = parse_fortran_line(line)
                while len(values) < 7:
                    values.append(0.0)
                i1, i2, j1, j2, k1, k2 = [int(x) for x in values[:6]]
                kx_new = values[6]
                
                iocode.write(f"               {i1:4d}{i2:4d}{j1:4d}{j2:4d}{k1:4d}{k2:4d}  {kx_new:10.2f}\n")
                
                for k in range(k1-1, k2):
                    for j in range(j1-1, j2):
                        for i in range(i1-1, i2):
                            self.sim.kx[i, j, k] = kx_new
        
        # Apply KY modifications
        if numky > 0:
            iocode.write("\n\n               ********KY NODE MODIFICATIONS**********\n\n")
            iocode.write("               I1  I2  J1  J2  K1  K2  NEW KY VALUE\n")
            
            for n in range(numky):
                line = iread.readline()
                values = parse_fortran_line(line)
                while len(values) < 7:
                    values.append(0.0)
                i1, i2, j1, j2, k1, k2 = [int(x) for x in values[:6]]
                ky_new = values[6]
                
                iocode.write(f"               {i1:4d}{i2:4d}{j1:4d}{j2:4d}{k1:4d}{k2:4d}  {ky_new:10.2f}\n")
                
                for k in range(k1-1, k2):
                    for j in range(j1-1, j2):
                        for i in range(i1-1, i2):
                            self.sim.ky[i, j, k] = ky_new
        
        # Apply KZ modifications
        if numkz > 0:
            iocode.write("\n\n               ********KZ NODE MODIFICATIONS**********\n\n")
            iocode.write("               I1  I2  J1  J2  K1  K2  NEW KZ VALUE\n")
            
            for n in range(numkz):
                line = iread.readline()
                values = parse_fortran_line(line)
                while len(values) < 7:
                    values.append(0.0)
                i1, i2, j1, j2, k1, k2 = [int(x) for x in values[:6]]
                kz_new = values[6]
                
                iocode.write(f"               {i1:4d}{i2:4d}{j1:4d}{j2:4d}{k1:4d}{k2:4d}  {kz_new:10.2f}\n")
                
                for k in range(k1-1, k2):
                    for j in range(j1-1, j2):
                        for i in range(i1-1, i2):
                            self.sim.kz[i, j, k] = kz_new
        
        # CRITICAL: After all porosity modifications are complete, copy porosity to vp array
        # The Fortran code stores porosity in VP, then TRANS multiplies it by block volume
        # to compute pore volume. We must do the same.
        np.copyto(self.sim.vp, self.sim.phi)
        
        # DEBUG: Verify porosity was copied
        iocode.write(f"\n   DEBUG: Porosity copied to VP array\n")
        iocode.write(f"   DEBUG: VP min={self.sim.vp.min():.6f}, max={self.sim.vp.max():.6f}, mean={self.sim.vp.mean():.6f}\n")
        iocode.write(f"   DEBUG: PHI min={self.sim.phi.min():.6f}, max={self.sim.phi.max():.6f}, mean={self.sim.phi.mean():.6f}\n")
        
        iocode.write("\n   ROCK PROPERTIES SUCCESSFULLY READ AND VALIDATED\n")


class PostProcessing:
    """Handles post-processing output and summary reports"""
    
    def __init__(self, simulator):
        """
        Initialize with reference to main simulator
        
        Parameters:
        -----------
        simulator : BOASTSimulator
            Reference to main simulator object
        """
        self.sim = simulator
    
    def postp(self) -> None:
        """
        Post-processing package - generates time series plots
        
        Creates plots for:
        - Oil, gas, water production rates
        - Gas and water injection rates
        - Cumulative production and injection
        - GOR and WOR ratios
        - Average reservoir pressure
        - Aquifer influx
        """
        iocode = self.sim.iocode
        nsteps = self.sim.nsteps  # Current time step
        
        if nsteps < 1:
            return
        
        iocode.write("\n\n")
        iocode.write("*" * 80 + "\n")
        iocode.write("*" + " " * 78 + "*\n")
        iocode.write("*" + " " * 20 + "POST-PLOT PACKAGE: TIME SERIES PLOTS" + " " * 22 + "*\n")
        iocode.write("*" + " " * 78 + "*\n")
        iocode.write("*" * 80 + "\n\n")
        
        # Create visualization object
        viz = OutputVisualization()
        
        # Get time array
        time = self.sim.time_array[:nsteps]
        
        # Plot 1: Oil production rate
        if hasattr(self.sim, 'opr_array'):
            viz.ploti("OIL PRODUCTION RATE (STB/D)", 
                     time, self.sim.opr_array[:nsteps], nsteps, iocode,
                     "TIME (DAYS)", "OPR (STB/D)")
        
        # Plot 2: Gas production rate
        if hasattr(self.sim, 'gpr_array'):
            viz.ploti("GAS PRODUCTION RATE (MSCF/D)",
                     time, self.sim.gpr_array[:nsteps], nsteps, iocode,
                     "TIME (DAYS)", "GPR (MSCF/D)")
        
        # Plot 3: Water production rate
        if hasattr(self.sim, 'wpr_array'):
            viz.ploti("WATER PRODUCTION RATE (STB/D)",
                     time, self.sim.wpr_array[:nsteps], nsteps, iocode,
                     "TIME (DAYS)", "WPR (STB/D)")
        
        # Plot 4: Cumulative oil production
        if hasattr(self.sim, 'cop_array'):
            viz.ploti("CUMULATIVE OIL PRODUCTION (STB)",
                     time, self.sim.cop_array[:nsteps], nsteps, iocode,
                     "TIME (DAYS)", "CUM OIL (STB)")
        
        # Plot 5: Cumulative gas production
        if hasattr(self.sim, 'cgp_array'):
            viz.ploti("CUMULATIVE GAS PRODUCTION (MSCF)",
                     time, self.sim.cgp_array[:nsteps], nsteps, iocode,
                     "TIME (DAYS)", "CUM GAS (MSCF)")
        
        # Plot 6: Cumulative water production
        if hasattr(self.sim, 'cwp_array'):
            viz.ploti("CUMULATIVE WATER PRODUCTION (STB)",
                     time, self.sim.cwp_array[:nsteps], nsteps, iocode,
                     "TIME (DAYS)", "CUM WATER (STB)")
        
        # Plot 7: Gas injection rate
        if hasattr(self.sim, 'gir_array'):
            viz.ploti("GAS INJECTION RATE (MSCF/D)",
                     time, self.sim.gir_array[:nsteps], nsteps, iocode,
                     "TIME (DAYS)", "GIR (MSCF/D)")
        
        # Plot 8: Water injection rate
        if hasattr(self.sim, 'wir_array'):
            viz.ploti("WATER INJECTION RATE (STB/D)",
                     time, self.sim.wir_array[:nsteps], nsteps, iocode,
                     "TIME (DAYS)", "WIR (STB/D)")
        
        # Plot 9: Cumulative gas injection
        if hasattr(self.sim, 'cgi_array'):
            viz.ploti("CUMULATIVE GAS INJECTION (MSCF)",
                     time, self.sim.cgi_array[:nsteps], nsteps, iocode,
                     "TIME (DAYS)", "CUM GAS INJ (MSCF)")
        
        # Plot 10: Cumulative water injection
        if hasattr(self.sim, 'cwi_array'):
            viz.ploti("CUMULATIVE WATER INJECTION (STB)",
                     time, self.sim.cwi_array[:nsteps], nsteps, iocode,
                     "TIME (DAYS)", "CUM WATER INJ (STB)")
        
        # Plot 11: GOR
        if hasattr(self.sim, 'gor_array'):
            viz.ploti("GAS-OIL RATIO (SCF/STB)",
                     time, self.sim.gor_array[:nsteps], nsteps, iocode,
                     "TIME (DAYS)", "GOR (SCF/STB)")
        
        # Plot 12: WOR
        if hasattr(self.sim, 'wor_array'):
            viz.ploti("WATER-OIL RATIO (STB/STB)",
                     time, self.sim.wor_array[:nsteps], nsteps, iocode,
                     "TIME (DAYS)", "WOR (STB/STB)")
        
        # Plot 13: Average reservoir pressure
        if hasattr(self.sim, 'pavg_array'):
            viz.ploti("AVERAGE RESERVOIR PRESSURE (PSIA)",
                     time, self.sim.pavg_array[:nsteps], nsteps, iocode,
                     "TIME (DAYS)", "PRESSURE (PSIA)")
        
        # Plot 14: Aquifer influx rate
        if hasattr(self.sim, 'aqrate_array'):
            viz.ploti("AQUIFER INFLUX RATE (STB/D)",
                     time, self.sim.aqrate_array[:nsteps], nsteps, iocode,
                     "TIME (DAYS)", "AQ RATE (STB/D)")
        
        # Plot 15: Cumulative aquifer influx
        if hasattr(self.sim, 'cumaq_array'):
            viz.ploti("CUMULATIVE AQUIFER INFLUX (STB)",
                     time, self.sim.cumaq_array[:nsteps], nsteps, iocode,
                     "TIME (DAYS)", "CUM AQ (STB)")
    
    def prtps(self, nloop: int, time: float, delt: float, 
              oerror: float, gerror: float, werror: float,
              coerr: float, cgerr: float, cwerr: float) -> None:
        """
        Print summary report with material balance and production data
        
        Parameters:
        -----------
        nloop : int
            Time step loop counter
        time : float
            Current simulation time (days)
        delt : float
            Time step size (days)
        oerror, gerror, werror : float
            Oil, gas, water material balance errors (%)
        coerr, cgerr, cwerr : float
            Cumulative material balance errors (%)
        """
        ii, jj, kk = self.sim.ii, self.sim.jj, self.sim.kk
        iocode = self.sim.iocode
        
        # Print header
        iocode.write("\f")  # Form feed
        iocode.write(" " * 30 + "*" * 69 + "\n")
        iocode.write(" " * 30 + "*" + " " * 67 + "*\n")
        iocode.write(" " * 30 + "*" + " " * 67 + "*\n")
        iocode.write(" " * 30 + "* " + " " * 13 + 
                    "SUMMARY REPORT: BOAST II (RELEASE 1.2)  " + " " * 13 + "*\n")
        iocode.write(" " * 30 + "*" + " " * 67 + "*\n")
        iocode.write(" " * 30 + "*" + " " * 67 + "*\n")
        iocode.write(" " * 30 + "*" * 69 + "\n\n\n")
        
        # Print time step info
        iocode.write(f"\n TIME STEP NUMBER             = {nloop:8d}        ")
        iocode.write(f"CUMULATIVE TIME (DAYS)        = {time:10.2f}\n")
        iocode.write(f" TIME STEP SIZE (DAYS)        = {delt:10.4f}\n\n")
        
        # Print material balance errors
        iocode.write(f" OIL MATERIAL BALANCE (%)    = {oerror:9.6f}   ")
        iocode.write(f"GAS MATERIAL BALANCE (%)     = {gerror:9.6f}  ")
        iocode.write(f"WATER MATERIAL BALANCE (%)    = {werror:9.6f}\n")
        
        iocode.write(f" CUM. OIL MATERIAL BALANCE(%)= {coerr:9.6f}   ")
        iocode.write(f"CUM. GAS MATERIAL BALANCE(%) = {cgerr:9.6f}  ")
        iocode.write(f"CUM. WATER MATERIAL BALANCE(%)= {cwerr:9.6f}\n\n")
        
        # Print production/injection rates and cumulatives
        opr = self.sim.opr if hasattr(self.sim, 'opr') else 0.0
        cop = self.sim.cop if hasattr(self.sim, 'cop') else 0.0
        gpr = self.sim.gpr if hasattr(self.sim, 'gpr') else 0.0
        cgp = self.sim.cgp if hasattr(self.sim, 'cgp') else 0.0
        wpr = self.sim.wpr if hasattr(self.sim, 'wpr') else 0.0
        cwp = self.sim.cwp if hasattr(self.sim, 'cwp') else 0.0
        gir = self.sim.gir if hasattr(self.sim, 'gir') else 0.0
        cgi = self.sim.cgi if hasattr(self.sim, 'cgi') else 0.0
        wir = self.sim.wir if hasattr(self.sim, 'wir') else 0.0
        cwi = self.sim.cwi if hasattr(self.sim, 'cwi') else 0.0
        
        iocode.write(f" OIL PRODUCTION RATE (STB/D) = {opr:9.1f}   ")
        iocode.write(f"CUM. OIL PRODUCTION (STB)    = {cop:10.4e}\n")
        
        iocode.write(f" GAS PRODUCTION RATE (MSCF/D)= {gpr:9.1f}   ")
        iocode.write(f"CUM. GAS PRODUCTION (MSCF)   = {cgp:10.4e}\n")
        
        iocode.write(f" WATER PRODUCTION RATE(STB/D)= {wpr:9.1f}\n")
        iocode.write(f"{' ':>42}CUM. WATER PRODUCTION (STB)  = {cwp:10.4e}\n\n")
        
        iocode.write(f" GAS INJECTION RATE (MSCF/D) = {gir:9.1f}   ")
        iocode.write(f"CUM. GAS INJECTION (MSCF)    = {cgi:10.4e}\n")
        
        iocode.write(f" WATER INJECTION RATE (STB/D)= {wir:9.1f}   ")
        iocode.write(f"CUM. WATER INJECTION (STB)   = {cwi:10.4e}\n\n")
        
        # Calculate and print WOR and GOR
        wor = wpr / opr if opr > 0.001 else 0.0
        cwor = cwp / cop if cop > 0.001 else 0.0
        gor = gpr / opr if opr > 0.001 else 0.0
        cgor = cgp / cop if cop > 0.001 else 0.0
        
        iocode.write(f" PRODUCING WOR (STB/STB)     = {wor:9.3f}   ")
        iocode.write(f"CUM. WOR (STB/STB)           = {cwor:9.3f}\n")
        
        iocode.write(f" PRODUCING GOR (SCF/STB)     = {gor:9.1f}   ")
        iocode.write(f"CUM. GOR (SCF/STB)           = {cgor:9.1f}\n\n")
        
        # Call aquifer print routine if needed
        if hasattr(self.sim, 'iaqopt') and self.sim.iaqopt > 0:
            self.aqprnt()
        
        # Print pressure and saturation arrays based on map options
        iocode.write("\f")  # Form feed
        
        if nloop == 1:
            iocode.write(f"\n{'':>14}{'*' * 7} INITIAL ARRAYS {'*' * 7}\n\n")
        
        # Print arrays based on map options
        if self.sim.ipmap == 1:
            self.print_pressure_array()
        elif self.sim.ipmap == 2:
            OutputVisualization.plot(self.sim.p, ii, jj, kk, iocode, 1)
        
        if self.sim.isomap == 1:
            self.print_saturation_array('OIL', self.sim.so)
        elif self.sim.isomap == 2:
            OutputVisualization.plot(self.sim.so, ii, jj, kk, iocode, 2)
        
        if self.sim.iswmap == 1:
            self.print_saturation_array('WATER', self.sim.sw)
        elif self.sim.iswmap == 2:
            OutputVisualization.plot(self.sim.sw, ii, jj, kk, iocode, 3)
        
        if self.sim.isgmap == 1:
            self.print_saturation_array('GAS', self.sim.sg)
        elif self.sim.isgmap == 2:
            OutputVisualization.plot(self.sim.sg, ii, jj, kk, iocode, 4)
        
        if self.sim.ipbmap == 1 or nloop == 1:
            self.print_bubble_point_array()
        elif self.sim.ipbmap == 2:
            OutputVisualization.plot(self.sim.pbot, ii, jj, kk, iocode, 5)
        
        if hasattr(self.sim, 'iaqopt') and self.sim.iaqopt > 0:
            if self.sim.iaqmap == 1 or nloop == 1:
                self.print_aquifer_array()
            elif self.sim.iaqmap == 2:
                OutputVisualization.plot(self.sim.cumaqw, ii, jj, kk, iocode, 6)
        
        # Print footer
        if nloop != 1:
            iocode.write(f"\n\n\n{'':>3}{'*' * 49}  END OF SUMMARY REPORT  {'*' * 49}\n\n\n\n\n\n")
        else:
            iocode.write(f"\n\n\n{'':>3}{'*' * 49}  END OF INITIALIZATION  {'*' * 49}\n\n\n\n\n\n")
    
    def print_pressure_array(self) -> None:
        """Print reservoir pressure distribution in tabular format"""
        ii, jj, kk = self.sim.ii, self.sim.jj, self.sim.kk
        iocode = self.sim.iocode
        ihed = self.sim.ihed
        
        iocode.write(f"\n\n\n{'':>14}***** RESERVOIR PRESSURE DISTRIBUTION *****\n\n")
        
        # Print in blocks of 15 columns
        intrvl = 15
        ichop = ii // intrvl
        ir2 = 0
        
        for icnt in range(ichop + 1):
            ir1 = ir2
            ir2 = min(ir1 + intrvl, ii)
            if ir1 >= ir2:
                break
            
            for k in range(kk):
                iocode.write(f"\n K = {k+1:2d}\n")
                iocode.write("     ")
                for i in range(ir1, ir2):
                    iocode.write(f"  {ihed[i]:4d}  ")
                iocode.write("\n\n")
                
                for j in range(jj):
                    iocode.write(f" {j+1:3d} ")
                    for i in range(ir1, ir2):
                        iocode.write(f"{self.sim.p[i, j, k]:8.0f}")
                    iocode.write("\n")
    
    def print_saturation_array(self, phase: str, sat_array: np.ndarray) -> None:
        """Print saturation distribution in tabular format"""
        ii, jj, kk = self.sim.ii, self.sim.jj, self.sim.kk
        iocode = self.sim.iocode
        ihed = self.sim.ihed
        
        iocode.write(f"\n\n\n{'':>14}*********  {phase} SATURATION  *********\n\n")
        
        # Print in blocks of 15 columns
        intrvl = 15
        ichop = ii // intrvl
        ir2 = 0
        
        for icnt in range(ichop + 1):
            ir1 = ir2
            ir2 = min(ir1 + intrvl, ii)
            if ir1 >= ir2:
                break
            
            for k in range(kk):
                iocode.write(f"\n K = {k+1:2d}\n")
                iocode.write("     ")
                for i in range(ir1, ir2):
                    iocode.write(f"  {ihed[i]:4d}  ")
                iocode.write("\n\n")
                
                for j in range(jj):
                    iocode.write(f" {j+1:3d} ")
                    for i in range(ir1, ir2):
                        iocode.write(f"{sat_array[i, j, k]:8.3f}")
                    iocode.write("\n")
    
    def print_bubble_point_array(self) -> None:
        """Print bubble point pressure distribution"""
        ii, jj, kk = self.sim.ii, self.sim.jj, self.sim.kk
        iocode = self.sim.iocode
        ihed = self.sim.ihed
        
        iocode.write(f"\n\n\n{'':>14}***** BUBBLE POINT PRESSURE DISTRIBUTION *****\n\n")
        
        # Print in blocks of 15 columns
        intrvl = 15
        ichop = ii // intrvl
        ir2 = 0
        
        for icnt in range(ichop + 1):
            ir1 = ir2
            ir2 = min(ir1 + intrvl, ii)
            if ir1 >= ir2:
                break
            
            for k in range(kk):
                iocode.write(f"\n K = {k+1:2d}\n")
                iocode.write("     ")
                for i in range(ir1, ir2):
                    iocode.write(f"  {ihed[i]:4d}  ")
                iocode.write("\n\n")
                
                for j in range(jj):
                    iocode.write(f" {j+1:3d} ")
                    for i in range(ir1, ir2):
                        iocode.write(f"{self.sim.pbot[i, j, k]:8.0f}")
                    iocode.write("\n")
    
    def print_aquifer_array(self) -> None:
        """Print cumulative aquifer influx distribution"""
        ii, jj, kk = self.sim.ii, self.sim.jj, self.sim.kk
        iocode = self.sim.iocode
        ihed = self.sim.ihed
        
        iocode.write(f"\n\n\n{'':>14}***** CUM. AQ INFLUX (STB) DISTRIBUTION *****\n\n")
        
        # Print in blocks of 10 columns (wider format for scientific notation)
        intrvl = 10
        ichop = ii // intrvl
        ir2 = 0
        
        for icnt in range(ichop + 1):
            ir1 = ir2
            ir2 = min(ir1 + intrvl, ii)
            if ir1 >= ir2:
                break
            
            for k in range(kk):
                iocode.write(f"\n K = {k+1:2d}\n")
                iocode.write("     ")
                for i in range(ir1, ir2):
                    iocode.write(f"    {ihed[i]:4d}  ")
                iocode.write("\n\n")
                
                for j in range(jj):
                    iocode.write(f" {j+1:3d} ")
                    for i in range(ir1, ir2):
                        iocode.write(f"{self.sim.cumaqw[i, j, k]:10.4e}")
                    iocode.write("\n")
    
    def aqprnt(self) -> None:
        """Print aquifer information (placeholder for aquifer-specific output)"""
        iocode = self.sim.iocode
        
        if not hasattr(self.sim, 'iaqopt') or self.sim.iaqopt <= 0:
            return
        
        iocode.write("\n   AQUIFER INFORMATION:\n")
        iocode.write(f"   Aquifer option = {self.sim.iaqopt}\n")
        
        # Additional aquifer-specific output would go here
        # This would be implemented based on aquifer model details
