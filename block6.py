"""
BLOCK6.FOR - PVT Tables, Transmissibility, and Initialization

This module contains:
1. PVTTables class: Read and validate rock and PVT data
2. Transmissibility class: Calculate grid block transmissibilities
3. Initialization class: Initialize pressure and saturation distributions
4. Utility functions: VISCY, ZANDC, XLGR4, TRIKRO

Converted from BOAST II (Release 1.2) Fortran code
"""

import numpy as np
from typing import TextIO, Tuple
import math
from block2 import parse_fortran_line


class PVTTables:
    """
    Read and validate PVT and rock property tables
    
    Handles:
    - Rock region and PVT region distributions
    - Relative permeability and capillary pressure tables
    - Oil, water, gas PVT properties
    - Real gas properties (via correlations or tables)
    - Bubble point pressure initialization
    - Stock tank densities
    """
    
    def __init__(self, simulator):
        """
        Initialize with reference to main simulator
        
        Parameters:
        -----------
        simulator : BOASTSimulator
            Reference to main simulator object
        """
        self.sim = simulator
    
    def table(self, iocode: TextIO, ii: int, jj: int, kk: int) -> None:
        """
        Read rock and PVT data tables
        
        Parameters:
        -----------
        iocode : file object
            Output file handle
        ii, jj, kk : int
            Grid dimensions
        """
        iread = self.sim.iread
        
        iocode.write("\f\n")
        iocode.write(f"{'':>14}***** EMPIRICAL DATA TABLES *****\n")
        
        # Read number of regions
        iread.readline()  # Skip comment line
        line = iread.readline()
        values = parse_fortran_line(line)
        nrock, npvt = int(values[0]), int(values[1])
        
        self.sim.nrock = nrock
        self.sim.npvt = npvt
        
        iocode.write(f"\n\n{'':>14}*** NUMBER OF REGIONS ***\n")
        iocode.write(f"{'':>19}ROCK{'':>10}{nrock:5d}\n")
        iocode.write(f"{'':>19}PVT{'':>10}{npvt:5d}\n\n")
        
        # Initialize region arrays
        self.sim.irock.fill(1)
        self.sim.ipvt.fill(1)
        
        # Read rock regions
        if nrock > 1:
            self._read_rock_regions(iocode, ii, jj, kk)
        
        # Read PVT regions
        if npvt > 1:
            self._read_pvt_regions(iocode, ii, jj, kk)
        
        # Read relative permeability tables
        self._read_rel_perm_tables(iocode, nrock)
        
        # Read PVT tables for each region
        for np in range(npvt):
            self._read_pvt_region(iocode, np, ii, jj, kk)
    
    def _read_rock_regions(self, iocode: TextIO, ii: int, jj: int, kk: int) -> None:
        """Read rock region modifications and distribution"""
        iread = self.sim.iread
        
        iread.readline()  # Skip comment
        line = iread.readline()
        values = parse_fortran_line(line)
        numrok = int(values[0])
        
        if numrok > 0:
            iocode.write(f"\n\n{'':>14}***** ROCK REGION MODIFICATIONS *****\n\n")
            iocode.write(f"{'':>14}   I1  I2  J1  J2  K1  K2  IROCK\n")
            
            for l in range(numrok):
                line = iread.readline()
                values = parse_fortran_line(line)
                i1, i2, j1, j2, k1, k2, ival = [int(x) for x in values[:7]]
                iocode.write(f"{'':>14}{i1:4d}{i2:4d}{j1:4d}{j2:4d}{k1:4d}{k2:4d}   {ival:3d}\n")
                
                for k in range(k1-1, k2):
                    for j in range(j1-1, j2):
                        for i in range(i1-1, i2):
                            self.sim.irock[i, j, k] = ival
        
        # Print rock region distribution
        iocode.write(f"\n\n{'':>14}***** ROCK REGION DISTRIBUTION *****\n\n")
        for k in range(kk):
            iocode.write(f"\n K ={k+1:4d}\n\n")
            for j in range(jj):
                if numrok == 0:
                    line = iread.readline()
                    values = parse_fortran_line(line)
                    values = [int(x) for x in values]
                    for i in range(ii):
                        self.sim.irock[i, j, k] = values[i]
                
                iocode.write(" ")
                for i in range(min(20, ii)):
                    iocode.write(f"{self.sim.irock[i, j, k]:6d}")
                iocode.write("\n")
    
    def _read_pvt_regions(self, iocode: TextIO, ii: int, jj: int, kk: int) -> None:
        """Read PVT region modifications and distribution"""
        iread = self.sim.iread
        
        iread.readline()  # Skip comment
        line = iread.readline()
        values = parse_fortran_line(line)
        numpvt = int(values[0])
        
        if numpvt > 0:
            iocode.write(f"\n\n{'':>14}***** PVT REGION MODIFICATIONS *****\n\n")
            iocode.write(f"{'':>14}   I1  I2  J1  J2  K1  K2  IPVT\n")
            
            for l in range(numpvt):
                line = iread.readline()
                values = parse_fortran_line(line)
                i1, i2, j1, j2, k1, k2, ival = [int(x) for x in values[:7]]
                iocode.write(f"{'':>14}{i1:4d}{i2:4d}{j1:4d}{j2:4d}{k1:4d}{k2:4d}   {ival:3d}\n")
                
                for k in range(k1-1, k2):
                    for j in range(j1-1, j2):
                        for i in range(i1-1, i2):
                            self.sim.ipvt[i, j, k] = ival
        
        # Print PVT region distribution
        iocode.write(f"\n\n{'':>14}***** PVT REGION DISTRIBUTION *****\n\n")
        for k in range(kk):
            iocode.write(f"\n K ={k+1:4d}\n\n")
            for j in range(jj):
                if numpvt == 0:
                    line = iread.readline()
                    values = parse_fortran_line(line)
                    values = [int(x) for x in values]
                    for i in range(ii):
                        self.sim.ipvt[i, j, k] = values[i]
                
                iocode.write(" ")
                for i in range(min(20, ii)):
                    iocode.write(f"{self.sim.ipvt[i, j, k]:6d}")
                iocode.write("\n")
    
    def _read_rel_perm_tables(self, iocode: TextIO, nrock: int) -> None:
        """Read relative permeability and capillary pressure tables"""
        iread = self.sim.iread
        
        for nr in range(nrock):
            iread.readline()  # Skip comment
            iocode.write(f"\f\n{'':>14}ROCK REGION # {nr+1:5d}\n\n")
            iocode.write(f"  SATURATION    KROT      KRW       KRG        KROG     PCOW      PCGO\n")
            iocode.write(f"{'':>56}(PSI){'':>6}(PSI)\n\n")
            
            for i in range(25):
                line = iread.readline()
                values = parse_fortran_line(line)
                
                # Ensure we have at least 7 values (pad with zeros if needed)
                while len(values) < 7:
                    values.append(0.0)
                
                self.sim.sat[nr, i] = values[0]
                self.sim.krot[nr, i] = values[1]
                self.sim.krwt[nr, i] = values[2]
                self.sim.krgt[nr, i] = values[3]
                self.sim.krogt[nr, i] = values[4]
                self.sim.pcowt[nr, i] = values[5]
                self.sim.pcgot[nr, i] = values[6]
                
                iocode.write(f" {values[0]:10.4f}{values[1]:10.4f}{values[2]:10.4f}{values[3]:10.4f}" +
                           f"{values[4]:10.4f}{values[5]:10.4f}{values[6]:10.4f}\n")
                
                if values[0] >= 1.0001:
                    break
            
            self.sim.msat[nr] = i + 1
            
            # Read three-phase parameters
            iread.readline()  # Skip comment
            line = iread.readline()
            values = parse_fortran_line(line)
            ithree, swr = values[0], values[1]
            self.sim.ithree[nr] = int(ithree)
            self.sim.swr[nr] = swr
            
            if ithree == 0:
                iocode.write(f"\n\n{'':>9}THE THREE-PHASE REL. PERM. CALC. IS NOT USED.\n")
            else:
                iocode.write(f"\n\n{'':>9}THE THREE-PHASE REL. PERM. CALC. IS USED.\n")
            iocode.write(f"{'':>9}IRREDUCIBLE WATER SATURATION IS {swr:6.4f}\n")
            
            # Validate saturation data
            self._validate_sat_data(iocode, nr)
    
    def _validate_sat_data(self, iocode: TextIO, nr: int) -> None:
        """Validate saturation-dependent data"""
        ifatal = 0
        
        for i in range(self.sim.msat[nr]):
            sat = self.sim.sat[nr, i]
            if not (-0.1 <= sat <= 1.1):
                self.sim.ifatal += 1
                iocode.write(f"\n     -----SAT ERROR FOR ROCK REGION{nr+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if not (0.0 <= self.sim.krot[nr, i] <= 1.0):
                self.sim.ifatal += 1
                iocode.write(f"\n     -----KROW ERROR FOR ROCK REGION{nr+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if not (0.0 <= self.sim.krwt[nr, i] <= 1.0):
                self.sim.ifatal += 1
                iocode.write(f"\n     -----KRW ERROR FOR ROCK REGION{nr+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if not (0.0 <= self.sim.krgt[nr, i] <= 1.0):
                self.sim.ifatal += 1
                iocode.write(f"\n     -----KRG ERROR FOR ROCK REGION{nr+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if not (0.0 <= self.sim.krogt[nr, i] <= 1.0):
                self.sim.ifatal += 1
                iocode.write(f"\n     -----KROG ERROR FOR ROCK REGION{nr+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if self.sim.pcowt[nr, i] < 0.0:
                self.sim.ifatal += 1
                iocode.write(f"\n     -----PCOW ERROR FOR ROCK REGION{nr+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if self.sim.pcgot[nr, i] < 0.0:
                self.sim.ifatal += 1
                iocode.write(f"\n     -----PCGO ERROR FOR ROCK REGION{nr+1:5d} AND ENTRY # {i+1:5d}\n")
        
        # Check endpoints
        if self.sim.sat[nr, 0] >= 0.0:
            self.sim.ifatal += 1
            iocode.write(f"\n     -----FIRST SAT ENTRY FOR ROCK REGION{nr+1:5d} SHOULD EQUAL -0.10.\n")
            iocode.write(f"{'':>9}ADJUST REL PERM AND CAP PRESS CURVES ALSO.\n\n")
        
        if self.sim.sat[nr, self.sim.msat[nr]-1] < 1.05:
            self.sim.ifatal += 1
            iocode.write(f"\n     -----LAST SAT ENTRY FOR ROCK REGION{nr+1:5d} SHOULD EQUAL 1.10.\n")
            iocode.write(f"{'':>9}ADJUST REL PERM AND CAP PRESS CURVES ALSO.\n\n")
        
        if not (0.0 <= self.sim.swr[nr] <= 1.0):
            self.sim.ifatal += 1
            iocode.write(f"\n     -----3-PHASE SWR ERROR FOR ROCK REGION{nr+1:5d}\n")
    
    def _read_pvt_region(self, iocode: TextIO, np: int, ii: int, jj: int, kk: int) -> None:
        """Read PVT tables for one region"""
        iread = self.sim.iread
        
        # Read bubble point parameters
        iread.readline()  # Skip comment
        line = iread.readline()
        values = parse_fortran_line(line)
        pbo, pbodat, pbgrad = values[0], values[1], values[2]
        
        iocode.write(f"\f\n{'':>14}PVT REGION # {np+1:5d}\n\n")
        iocode.write(f"\n\n{'':>9}BUBBLE POINT PARAMETERS:\n")
        iocode.write(f"{'':>4}BUBBLE POINT PRESSURE (PSI){'':>14}{pbo:10.3f}\n")
        iocode.write(f"{'':>4}BUBBLE POINT DEPTH (FT){'':>18}{pbodat:10.3f}\n")
        iocode.write(f"{'':>4}BUBBLE POINT GRADIENT (PSI/FT){'':>10}{pbgrad:10.3f}\n")
        
        # Initialize bubble point pressures
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    if self.sim.ipvt[i, j, k] == np + 1:
                        pbot_ijk = pbo + (pbodat - self.sim.el[i, j, k]) * pbgrad
                        self.sim.pbot[i, j, k] = pbot_ijk
                        self.sim.pbotn[i, j, k] = pbot_ijk
                        
                        if pbot_ijk < 0.0:
                            self.sim.ifatal += 1
                            iocode.write(f"\n     -----BUBBLE POINT PRESSURE ERROR IN GRID BLOCK IJK = {i+1:5d}{j+1:5d}{k+1:5d}\n")
        
        # Read oil properties
        self._read_oil_properties(iocode, np)
        
        # Read water properties
        self._read_water_properties(iocode, np)
        
        # Read gas and rock properties
        self._read_gas_properties(iocode, np)
        
        # Calculate slopes for compressibility
        self._calculate_slopes(iocode, np)
    
    def _read_oil_properties(self, iocode: TextIO, np: int) -> None:
        """Read oil PVT properties"""
        iread = self.sim.iread
        
        iread.readline()  # Skip comment
        iocode.write("\n\n   VSLOPE     BSLOPE       RSLOPE     PMAX  REPRS\n")
        iocode.write("  (CP/PSI) (RB/STB/PSI)    (1/PSI)   (PSI)\n\n")
        
        line = iread.readline()
        values = parse_fortran_line(line)
        vslope = values[0]
        bslope = values[1]
        rslope = values[2]
        pmaxt = values[3]
        ireprs = int(values[4]) if len(values) > 4 else 0
        
        self.sim.vslope[np] = vslope
        self.sim.bslope[np] = bslope
        self.sim.rslope[np] = rslope
        self.sim.pmaxt = pmaxt
        self.sim.ireprs = ireprs
        
        iocode.write(f" {vslope:10.3e} {bslope:10.3e} {rslope:10.5f} {pmaxt:10.2f} {ireprs:5d}\n")
        
        # Validate slopes
        if vslope < 0.0:
            self.sim.iwarn += 1
            iocode.write(f"\n     -----VSLOPE FOR PVT REGION{np+1:5d} IS NEGATIVE\n")
        if bslope > 0.0:
            self.sim.iwarn += 1
            iocode.write(f"\n     -----BSLOPE FOR PVT REGION{np+1:5d} IS POSITIVE\n")
        if rslope != 0.0:
            self.sim.iwarn += 1
            iocode.write(f"\n     -----RSLOPE FOR PVT REGION{np+1:5d} IS NOT ZERO\n")
        
        # Read oil table
        iread.readline()  # Skip comment
        iocode.write(f"\n\n{'':>4} **** SATURATED OIL PVT PROPERTIES ****\n")
        iocode.write(f"{'':>4} PRESSURE  VISCOSITY    FVF    SOLN. GAS\n")
        iocode.write(f"{'':>4}  (PSI)       (CP)    (RB/STB)  (SCF/STB)\n\n")
        
        for i in range(25):
            line = iread.readline()
            values = parse_fortran_line(line)
            while len(values) < 4:
                values.append(0.0)
            
            self.sim.pot[np, i] = values[0]
            self.sim.muot[np, i] = values[1]
            self.sim.bot[np, i] = values[2]
            rso_input = values[3]  # Rs in SCF/STB from input
            self.sim.rsot[np, i] = rso_input * 0.17809  # Convert to CF/CF for internal use
            
            iocode.write(f"   {values[0]:10.1f} {values[1]:8.4f}  {values[2]:8.4f}  {rso_input:8.2f}\n")
            
            if values[0] >= pmaxt:
                break
        
        self.sim.mpot[np] = i + 1
        
        # Validate oil data
        self._validate_oil_data(iocode, np)
    
    def _read_water_properties(self, iocode: TextIO, np: int) -> None:
        """Read water PVT properties"""
        iread = self.sim.iread
        
        iread.readline()  # Skip comment
        iocode.write(f"\n\n{'':>4}     **** WATER PVT PROPERTIES ***\n")
        iocode.write(f"{'':>4} PRESSURE  VISCOSITY    FVF    SOLN. GAS\n")
        iocode.write(f"{'':>4}  (PSI)       (CP)    (RB/STB)  (SCF/STB)\n\n")
        
        for i in range(25):
            line = iread.readline()
            values = parse_fortran_line(line)
            while len(values) < 4:
                values.append(0.0)
            
            self.sim.pwt[np, i] = values[0]
            self.sim.muwt[np, i] = values[1]
            self.sim.bwt[np, i] = values[2]
            self.sim.rswt[np, i] = values[3] * 0.17809  # Convert to CF/CF
            
            iocode.write(f"   {values[0]:10.1f} {values[1]:8.4f}  {values[2]:8.4f}  {values[3]:8.2f}\n")
            
            if values[0] >= self.sim.pmaxt:
                break
        
        self.sim.mpwt[np] = i + 1
        
        # Validate water data
        self._validate_water_data(iocode, np)
    
    def _read_gas_properties(self, iocode: TextIO, np: int) -> None:
        """Read gas PVT and rock compressibility properties"""
        iread = self.sim.iread
        from block4 import GasProperties
        
        iread.readline()  # Skip comment
        line = iread.readline()
        values = parse_fortran_line(line)
        kgcor = int(values[0])
        
        iread.readline()  # Skip comment
        
        if kgcor != 1:
            # Read gas table
            iocode.write(f"\n\n{'':>5}  ****  GAS  PVT AND ROCK COMP. ****\n")
            iocode.write(f"{'':>4} PRESSURE  VISCOSITY    FVF     PSEUDO-PRS  ROCK COMP.\n")
            iocode.write(f"{'':>4}  (PSI)       (CP)   (RCF/SCF) (PSIA**2/CP)  (1/PSI)\n\n")
            
            for i in range(999):
                line = iread.readline()
                values = parse_fortran_line(line)
                while len(values) < 5:
                    values.append(0.0)
                
                self.sim.pgt[np, i] = values[0]
                self.sim.mugt[np, i] = values[1]
                self.sim.bgt[np, i] = values[2]
                self.sim.psit[np, i] = values[3]
                self.sim.crt[np, i] = values[4]
                self.sim.prt[np, i] = values[0]
                
                iocode.write(f"   {values[0]:10.1f} {values[1]:8.4f}  {values[2]:10.4e}  {values[3]:10.4e}  {values[4]:10.3e}\n")
                
                if values[0] >= self.sim.pmaxt:
                    break
            
            self.sim.mpgt[np] = i + 1
        else:
            # Calculate gas properties using correlations
            gas_props = GasProperties(self.sim)
            self.sim.ifatal = gas_props.pseudo(np, self.sim.ifatal, iocode)
            
            # Read rock compressibility
            iocode.write(f"\n\n{'':>5}  **** ROCK COMPRESSIBILITY ****\n")
            iocode.write(f"{'':>4} PRESSURE    ROCK COMP.\n")
            iocode.write(f"{'':>4}  (PSI)      (1/PSI)\n\n")
            
            iread.readline()  # Skip comment
            
            for ig in range(self.sim.mpgt[np]):
                if ig > 0 and self.sim.prt[np, 0] == self.sim.pmaxt:
                    self.sim.prt[np, ig] = self.sim.pgt[np, ig]
                    self.sim.crt[np, ig] = self.sim.crt[np, 0]
                    if ig == self.sim.mpgt[np] - 1:
                        self.sim.prt[np, 0] = self.sim.pgt[np, 0]
                else:
                    line = iread.readline()
                    values = parse_fortran_line(line)
                    while len(values) < 2:
                        values.append(0.0)
                    self.sim.prt[np, ig] = values[0]
                    self.sim.crt[np, ig] = values[1]
                
                iocode.write(f"   {self.sim.prt[np, ig]:10.1f}  {self.sim.crt[np, ig]:10.3e}\n")
        
        # Validate gas data
        self._validate_gas_data(iocode, np)
        
        # Read stock tank densities
        iread.readline()  # Skip comment
        iocode.write(f"\n\n{'':>7}*** DENSITIES AT STD. CONDITIONS ***\n")
        iocode.write(f"{'':>13}OIL      WATER      GAS\n")
        iocode.write(f"{'':>9} (LBM/SCF) (LBM/SCF) (LBM/SCF)\n\n")
        
        line = iread.readline()
        values = parse_fortran_line(line)
        while len(values) < 3:
            values.append(0.0)
        self.sim.rhosco[np] = values[0]
        self.sim.rhoscw[np] = values[1]
        self.sim.rhoscg[np] = values[2]
        
        iocode.write(f"{'':>8}{values[0]:10.4f}{values[1]:10.4f}{values[2]:10.4f}\n")
        
        # Validate densities
        if self.sim.rhosco[np] <= 0.0:
            self.sim.ifatal += 1
            iocode.write(f"\n     -----OIL DENSITY FOR PVT REGION{np+1:5d} IN ERROR\n")
        if self.sim.rhoscw[np] <= 0.0:
            self.sim.ifatal += 1
            iocode.write(f"\n     -----WATER DENSITY FOR PVT REGION{np+1:5d} IN ERROR\n")
        if self.sim.rhoscg[np] <= 0.0:
            self.sim.ifatal += 1
            iocode.write(f"\n     -----GAS DENSITY  FOR PVT REGION{np+1:5d} IN ERROR\n")
    
    def _calculate_slopes(self, iocode: TextIO, np: int) -> None:
        """Calculate slopes for compressibility calculations"""
        iocode.write(f"\n\n{'':>14}***** SLOPES FOR COMPRESSIBILITY CALCULATIONS ****\n")
        
        # Oil slopes
        iocode.write(f"\n\n{'':>11}PRESSURE{'':>13}BO{'':>9}DBO/DP{'':>13}RSO{'':>8}DRSO/DP\n")
        iocode.write(f"{'':>12}(PSI){'':>9}(RB/STB){'':>9}(RB/STB/PSI){'':>9}(CF/CF){'':>10}(1/PSI)\n\n")
        
        for i in range(1, self.sim.mpot[np]):
            div = 1.0 / (self.sim.pot[np, i] - self.sim.pot[np, i-1])
            self.sim.bopt[np, i] = (self.sim.bot[np, i] - self.sim.bot[np, i-1]) * div
            self.sim.rsopt[np, i] = (self.sim.rsot[np, i] - self.sim.rsot[np, i-1]) * div
            
            iocode.write(f"{'':>11}{self.sim.pot[np, i]:7.1f}{self.sim.bot[np, i]:10.4f}  " +
                        f"{self.sim.bopt[np, i]:11.4e}{self.sim.rsot[np, i]:9.1f}  {self.sim.rsopt[np, i]:11.4e}\n")
        
        # Water slopes
        iocode.write(f"\n\n{'':>11}PRESSURE{'':>13}BW{'':>9}DBW/DP{'':>13}RSW{'':>8}DRSW/DP\n")
        iocode.write(f"{'':>12}(PSI){'':>9}(RB/STB){'':>9}(RB/STB/PSI){'':>9}(CF/CF){'':>10}(1/PSI)\n\n")
        
        for i in range(1, self.sim.mpwt[np]):
            div = 1.0 / (self.sim.pwt[np, i] - self.sim.pwt[np, i-1])
            self.sim.bwpt[np, i] = (self.sim.bwt[np, i] - self.sim.bwt[np, i-1]) * div
            self.sim.rswpt[np, i] = (self.sim.rswt[np, i] - self.sim.rswt[np, i-1]) * div
            
            iocode.write(f"{'':>11}{self.sim.pwt[np, i]:7.1f}{self.sim.bwt[np, i]:10.4f}  " +
                        f"{self.sim.bwpt[np, i]:11.4e}{self.sim.rswt[np, i]:9.1f}  {self.sim.rswpt[np, i]:11.4e}\n")
        
        # Gas slopes
        iocode.write(f"\n\n{'':>11}PRESSURE{'':>9}BG{'':>12}DBG/DP\n")
        iocode.write(f"{'':>13}(PSI){'':>13}(RCF/SCF){'':>8}(RCF/SCF/PSI)\n\n")
        
        for i in range(1, self.sim.mpgt[np]):
            self.sim.bgpt[np, i] = (self.sim.bgt[np, i] - self.sim.bgt[np, i-1]) / \
                                   (self.sim.pgt[np, i] - self.sim.pgt[np, i-1])
            
            iocode.write(f"{'':>9}{self.sim.pgt[np, i]:9.1f} {self.sim.bgt[np, i]:15.4e}{self.sim.bgpt[np, i]:15.4e}\n")
    
    def _validate_oil_data(self, iocode: TextIO, np: int) -> None:
        """Validate oil PVT data"""
        for i in range(self.sim.mpot[np]):
            if i == 0:
                if not (0.0 <= self.sim.pot[np, 0] <= self.sim.pmaxt):
                    self.sim.ifatal += 1
                    iocode.write(f"\n     -----POT ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
            else:
                if not (self.sim.pot[np, i-1] <= self.sim.pot[np, i] <= self.sim.pmaxt):
                    self.sim.ifatal += 1
                    iocode.write(f"\n     -----POT ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if self.sim.muot[np, i] < 0.0:
                self.sim.ifatal += 1
                iocode.write(f"\n     -----MUO ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if self.sim.bot[np, i] < 0.0:
                self.sim.ifatal += 1
                iocode.write(f"\n     -----BO ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if self.sim.rsot[np, i] < 0.0:
                self.sim.ifatal += 1
                iocode.write(f"\n     -----RSO ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
    
    def _validate_water_data(self, iocode: TextIO, np: int) -> None:
        """Validate water PVT data"""
        for i in range(self.sim.mpwt[np]):
            if i == 0:
                if not (0.0 <= self.sim.pwt[np, 0] <= self.sim.pmaxt):
                    self.sim.ifatal += 1
                    iocode.write(f"\n     -----PWT ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
            else:
                if not (self.sim.pwt[np, i-1] <= self.sim.pwt[np, i] <= self.sim.pmaxt):
                    self.sim.ifatal += 1
                    iocode.write(f"\n     -----PWT ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if self.sim.muwt[np, i] < 0.0:
                self.sim.ifatal += 1
                iocode.write(f"\n     -----MUW ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if self.sim.bwt[np, i] < 0.0:
                self.sim.ifatal += 1
                iocode.write(f"\n     -----BW ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if self.sim.rswt[np, i] < 0.0:
                self.sim.ifatal += 1
                iocode.write(f"\n     -----RSW ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
    
    def _validate_gas_data(self, iocode: TextIO, np: int) -> None:
        """Validate gas PVT data"""
        for i in range(self.sim.mpgt[np]):
            if i == 0:
                if not (0.0 <= self.sim.pgt[np, 0] <= self.sim.pmaxt):
                    self.sim.ifatal += 1
                    iocode.write(f"\n     -----PGT ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
            else:
                if not (self.sim.pgt[np, i-1] <= self.sim.pgt[np, i] <= self.sim.pmaxt + 0.1):
                    self.sim.ifatal += 1
                    iocode.write(f"\n     -----PGT ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if self.sim.mugt[np, i] < 0.0:
                self.sim.ifatal += 1
                iocode.write(f"\n     -----MUG ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if self.sim.bgt[np, i] < 0.0:
                self.sim.ifatal += 1
                iocode.write(f"\n     -----BG ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if self.sim.crt[np, i] < 0.0:
                self.sim.ifatal += 1
                iocode.write(f"\n     -----CR ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if not (0.0 <= self.sim.prt[np, i] <= self.sim.pmaxt + 0.1):
                self.sim.ifatal += 1
                iocode.write(f"\n     -----ROCK PRES ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")
            
            if self.sim.psit[np, i] < 0.0:
                self.sim.ifatal += 1
                iocode.write(f"\n     -----PSEUDO-PRES ERROR FOR PVT REGION{np+1:5d} AND ENTRY # {i+1:5d}\n")


class Transmissibility:
    """Calculate grid block transmissibilities"""
    
    def __init__(self, simulator):
        """Initialize with reference to main simulator"""
        self.sim = simulator
    
    def trans(self, ii: int, jj: int, kk: int) -> None:
        """
        Calculate transmissibilities for all grid blocks
        
        Parameters:
        -----------
        ii, jj, kk : int
            Grid dimensions
        """
        # Calculate geometric factors A1, A2, A3
        iocode = self.sim.iocode
        
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    facx = 1.0
                    facy = 1.0
                    facz = 1.0
                    
                    # X-direction factor
                    if i > 0 and i < ii - 1:
                        xdenom = (2.0 * self.sim.dx[i, j, k] + self.sim.dx[i+1, j, k] +
                                 self.sim.dx[i-1, j, k])
                        if xdenom > 0.0:
                            facx = 4.0 * self.sim.dx[i, j, k] / xdenom
                    
                    # Y-direction factor
                    if j > 0 and j < jj - 1:
                        ydenom = (2.0 * self.sim.dy[i, j, k] + self.sim.dy[i, j+1, k] +
                                 self.sim.dy[i, j-1, k])
                        if ydenom > 0.0:
                            facy = 4.0 * self.sim.dy[i, j, k] / ydenom
                    
                    # Z-direction factor
                    if k > 0 and k < kk - 1:
                        zdenom = (2.0 * self.sim.dznet[i, j, k] + self.sim.dznet[i, j, k+1] +
                                 self.sim.dznet[i, j, k-1])
                        if zdenom > 0.0:
                            facz = 4.0 * self.sim.dznet[i, j, k] / zdenom
                    
                    # Calculate A1, A2, A3
                    self.sim.a1[i, j, k] = facx * self.sim.kx[i, j, k] * self.sim.dy[i, j, k] * self.sim.dznet[i, j, k]
                    self.sim.a2[i, j, k] = facy * self.sim.ky[i, j, k] * self.sim.dx[i, j, k] * self.sim.dznet[i, j, k]
                    self.sim.a3[i, j, k] = facz * self.sim.kz[i, j, k] * self.sim.dx[i, j, k] * self.sim.dy[i, j, k]
        
        # Calculate transmissibilities
        # X-direction
        if ii > 1:
            for k in range(kk):
                for j in range(jj):
                    for i in range(1, ii):
                        xdenom = (self.sim.dx[i-1, j, k] * self.sim.a1[i, j, k] +
                                 self.sim.dx[i, j, k] * self.sim.a1[i-1, j, k])
                        if xdenom > 0.0:
                            self.sim.tx[i, j, k] = (0.012656 * self.sim.a1[i-1, j, k] *
                                                   self.sim.a1[i, j, k] / xdenom)
        
        # Y-direction
        if jj > 1:
            for k in range(kk):
                for j in range(1, jj):
                    for i in range(ii):
                        ydenom = (self.sim.dy[i, j-1, k] * self.sim.a2[i, j, k] +
                                 self.sim.dy[i, j, k] * self.sim.a2[i, j-1, k])
                        if ydenom > 0.0:
                            self.sim.ty[i, j, k] = (0.012656 * self.sim.a2[i, j-1, k] *
                                                   self.sim.a2[i, j, k] / ydenom)
        
        # Z-direction
        if kk > 1:
            for k in range(1, kk):
                for j in range(jj):
                    for i in range(ii):
                        zdenom = (self.sim.dznet[i, j, k-1] * self.sim.a3[i, j, k] +
                                 self.sim.dznet[i, j, k] * self.sim.a3[i, j, k-1])
                        if zdenom > 0.0:
                            self.sim.tz[i, j, k] = (0.012656 * self.sim.a3[i, j, k-1] *
                                                   self.sim.a3[i, j, k] / zdenom)
        
        # Calculate pore volumes
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    self.sim.vp[i, j, k] *= (self.sim.dx[i, j, k] * self.sim.dy[i, j, k] *
                                            self.sim.dznet[i, j, k])
        
        # Pore volume calculation complete
        
        # Set zero transmissibility for zero PV blocks
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    if self.sim.vp[i, j, k] <= 0.0:
                        self.sim.tx[i, j, k] = 0.0
                        if i < ii - 1:
                            self.sim.tx[i+1, j, k] = 0.0
                        self.sim.ty[i, j, k] = 0.0
                        if j < jj - 1:
                            self.sim.ty[i, j+1, k] = 0.0
                        self.sim.tz[i, j, k] = 0.0
                        if k < kk - 1:
                            self.sim.tz[i, j, k+1] = 0.0
        
        # Read transmissibility modifications
        self._read_trans_modifications(ii, jj, kk)
    
    def _read_trans_modifications(self, ii: int, jj: int, kk: int) -> None:
        """Read transmissibility modifications from input"""
        iread = self.sim.iread
        iocode = self.sim.iocode
        
        iread.readline()  # Skip comment
        line = iread.readline()
        values = parse_fortran_line(line)
        while len(values) < 4:
            values.append(0.0)
        numtx = int(values[0])
        numty = int(values[1])
        numtz = int(values[2])
        itcode = int(values[3])
        
        # TX modifications
        if numtx > 0:
            iocode.write(f"\n\n{'':>14}**********TRANSMISSIBILITY (TX) NODE MODIFICATIONS**********\n\n")
            iocode.write(f"{'':>14}   I1  I2  J1  J2  K1  K2  NEW TX VALUE\n")
            
            for l in range(numtx):
                line = iread.readline()
                values = parse_fortran_line(line)
                while len(values) < 7:
                    values.append(0.0)
                i1, i2 = int(values[0]), int(values[1])
                j1, j2 = int(values[2]), int(values[3])
                k1, k2 = int(values[4]), int(values[5])
                regval = float(values[6])
                
                iocode.write(f"{'':>14}{i1:4d}{i2:4d}{j1:4d}{j2:4d}{k1:4d}{k2:4d}  {regval:10.4e}\n")
                
                for k in range(k1-1, k2):
                    for j in range(j1-1, j2):
                        for i in range(i1-1, i2):
                            if self.sim.vp[i, j, k] > 0.0:
                                self.sim.tx[i, j, k] = regval
        
        # TY modifications
        if numty > 0:
            iocode.write(f"\n\n{'':>14}**********TRANSMISSIBILITY (TY) NODE MODIFICATIONS**********\n\n")
            iocode.write(f"{'':>14}   I1  I2  J1  J2  K1  K2  NEW TY VALUE\n")
            
            for l in range(numty):
                line = iread.readline()
                values = parse_fortran_line(line)
                while len(values) < 7:
                    values.append(0.0)
                i1, i2 = int(values[0]), int(values[1])
                j1, j2 = int(values[2]), int(values[3])
                k1, k2 = int(values[4]), int(values[5])
                regval = float(values[6])
                
                iocode.write(f"{'':>14}{i1:4d}{i2:4d}{j1:4d}{j2:4d}{k1:4d}{k2:4d}  {regval:10.4e}\n")
                
                for k in range(k1-1, k2):
                    for j in range(j1-1, j2):
                        for i in range(i1-1, i2):
                            if self.sim.vp[i, j, k] > 0.0:
                                self.sim.ty[i, j, k] = regval
        
        # TZ modifications
        if numtz > 0:
            iocode.write(f"\n\n{'':>14}**********TRANSMISSIBILITY (TZ) NODE MODIFICATIONS**********\n\n")
            iocode.write(f"{'':>14}   I1  I2  J1  J2  K1  K2  NEW TZ VALUE\n")
            
            for l in range(numtz):
                line = iread.readline()
                values = parse_fortran_line(line)
                while len(values) < 7:
                    values.append(0.0)
                i1, i2 = int(values[0]), int(values[1])
                j1, j2 = int(values[2]), int(values[3])
                k1, k2 = int(values[4]), int(values[5])
                regval = float(values[6])
                
                iocode.write(f"{'':>14}{i1:4d}{i2:4d}{j1:4d}{j2:4d}{k1:4d}{k2:4d}  {regval:10.4e}\n")
                
                for k in range(k1-1, k2):
                    for j in range(j1-1, j2):
                        for i in range(i1-1, i2):
                            if self.sim.vp[i, j, k] > 0.0:
                                self.sim.tz[i, j, k] = regval
        
        # Validate transmissibilities
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    if self.sim.tx[i, j, k] < 0.0:
                        self.sim.ifatal += 1
                        iocode.write(f"\n     -----GRID BLOCK TX ERROR AT IJK ={i+1:5d}{j+1:5d}{k+1:5d}\n")
                    
                    if self.sim.ty[i, j, k] < 0.0:
                        self.sim.ifatal += 1
                        iocode.write(f"\n     -----GRID BLOCK TY ERROR AT IJK ={i+1:5d}{j+1:5d}{k+1:5d}\n")
                    
                    if self.sim.tz[i, j, k] < 0.0:
                        self.sim.ifatal += 1
                        iocode.write(f"\n     -----GRID BLOCK TZ ERROR AT IJK ={i+1:5d}{j+1:5d}{k+1:5d}\n")


# Utility functions are defined at module level since they may be called independently

def viscy(temprd: float, prsprd: float, spg: float, tem: float,
          cnch2s: float, cncco2: float, cncn2: float) -> Tuple[float, float, int]:
    """
    Calculate gas viscosity using correlation
    
    Parameters:
    -----------
    temprd : float
        Reduced temperature
    prsprd : float
        Reduced pressure
    spg : float
        Gas gravity
    tem : float
        Temperature (°F)
    cnch2s, cncco2, cncn2 : float
        Impurity mole percentages
        
    Returns:
    --------
    visgr : float
        Reduced viscosity ratio
    visg : float
        Gas viscosity (cp)
    ierr : int
        Error flag
    """
    # Correlation tables (simplified version)
    # In production, use full tables from Fortran code
    
    ierr = 0
    
    # Validate inputs
    if not (1.05 <= temprd <= 3.00):
        return 0.0, 0.0, 1
    if not (0.01 <= prsprd <= 20.00):
        return 0.0, 0.0, 1
    if not (0.55 <= spg <= 1.50):
        return 0.0, 0.0, 1
    if not (40.0 <= tem <= 400.0):
        return 0.0, 0.0, 1
    
    # Simplified correlation
    visgr = 1.0 + 0.05 * prsprd  # Simplified
    
    # Base viscosity
    visgu = (0.126585e-01 - 0.611823e-02 * spg + 0.164574e-02 * spg**2 +
             0.164574e-04 * tem - 0.719221e-06 * spg * tem -
             0.609046e-06 * spg**2 * tem)
    
    # Corrections for impurities
    corh2s = ((0.000113 * cnch2s * spg - 0.000038 * cnch2s + 0.000001) *
              (1.0 / (1.0 + spg)) + 0.000001)
    corco2 = ((0.000134 * cncco2 * spg - 0.000004 * cncco2 + 0.000004 * spg) *
              (1.0 / (1.0 + spg)) - 0.000003)
    corn2 = ((0.000170 * cncn2 * spg + 0.000021 * cncn2 + 0.000010 * spg) *
             (1.0 / (1.0 + spg)) - 0.000006)
    
    visga = visgu + corh2s + corco2 + corn2
    visg = visgr * visga
    
    return visgr, visg, ierr


def zandc(temp: float, tempc: float, prs: float, prspc: float,
          cnch2s: float, cncco2: float) -> Tuple[float, float, int]:
    """
    Calculate Z-factor and gas compressibility using correlations
    
    Parameters:
    -----------
    temp : float
        Temperature (°F)
    tempc : float
        Pseudo-critical temperature (°R)
    prs : float
        Pressure (psia)
    prspc : float
        Pseudo-critical pressure (psia)
    cnch2s : float
        H2S mole %
    cncco2 : float
        CO2 mole %
        
    Returns:
    --------
    zed : float
        Z-factor
    cmpg : float
        Gas compressibility (1/psi)
    ierr : int
        Error flag
    """
    # Coefficients for Z-factor calculation
    a = [0.31506237, -1.04670990, -0.57832729, 0.53530771,
         -0.61232032, -0.10488813, 0.68157001, 0.68446549]
    
    frca = (cnch2s + cncco2) / 100.0
    frcb = cnch2s / 100.0
    
    # Wichert-Aziz correction
    eps = 120.0 * (frca**0.9 - frca**1.6) + 15.0 * (frcb**0.5 - frcb**4.0)
    tempca = tempc - eps
    prspca = prspc * tempca / (tempc + frcb * (1.0 - frcb) * eps)
    
    temprd = (temp + 460.0) / tempca
    prsprd = prs / prspca
    
    ierr = 0
    
    # Validate inputs
    if not (1.05 <= temprd <= 3.0):
        return 0.0, 0.0, 1
    if not (0.0 <= prsprd <= 15.0):
        return 0.0, 0.0, 1
    if not (0.0 <= frca <= 0.85):
        return 0.0, 0.0, 1
    
    # Calculate coefficients
    t1 = a[0] * temprd + a[1] + a[2] / (temprd * temprd)
    t2 = a[3] * temprd + a[4]
    t3 = a[4] * a[5]
    t4 = a[6] / (temprd * temprd)
    t5 = a[7]
    
    # Newton-Raphson iteration for reduced density
    denrd = 1.0
    
    for iter in range(10):
        denrd2 = denrd * denrd
        denrd3 = denrd2 * denrd
        denrd4 = denrd2 * denrd2
        denrd5 = denrd2 * denrd3
        
        p = ((temprd + t1 * denrd + t2 * denrd2 + t3 * denrd5) * denrd +
             t4 * denrd3 * (1.0 + t5 * denrd2) * math.exp(-t5 * denrd2))
        
        dp = (temprd + 2.0 * t1 * denrd + 3.0 * t2 * denrd2 +
              6.0 * t3 * denrd5 + t4 * denrd2 * math.exp(-t5 * denrd2) *
              (3.0 + 3.0 * t5 * denrd2 - 2.0 * t5 * t5 * denrd4))
        
        denrd1 = denrd - (p - 0.270 * prsprd) / dp
        
        if denrd1 <= 0.0:
            denrd1 = 0.5 * denrd
        if denrd1 >= 2.2:
            denrd1 = denrd + 0.9 * (2.2 - denrd)
        
        if abs(denrd - denrd1) < 0.00001:
            denrd = denrd1
            break
        
        denrd = denrd1
    
    # Calculate Z-factor
    zed = 0.270 * prsprd / (denrd * temprd)
    
    # Calculate compressibility
    denrd2 = denrd * denrd
    denrd4 = denrd2 * denrd2
    
    dzed = (t1 / temprd + 2.0 * t2 / temprd * denrd +
            5.0 * t3 * denrd4 / temprd +
            (1.0 + t5 * denrd2 - t5 * t5 * denrd4) * 2.0 * t4 / temprd * denrd *
            math.exp(-t5 * denrd2))
    
    cmpprd = 1.0 / prsprd - 0.270 * dzed / (zed * zed * temprd *
                                            (1.0 + denrd / zed * dzed))
    cmpg = cmpprd / prspca
    
    return zed, cmpg, ierr


def fptd(a0: float, a1: float, a2: float, a3: float, x: float) -> float:
    """Calculate dimensionless pressure for Carter-Tracy aquifer"""
    return a0 + a1 * x + a2 * math.log(x) + a3 * math.log(x) * math.log(x)


def fdptd(a1: float, a2: float, a3: float, x: float) -> float:
    """Calculate derivative of dimensionless pressure for Carter-Tracy aquifer"""
    return a1 + (a2 / x) + (2.0 * a3 * math.log(x) / x)


def trikro(simulator, ireg: int, so: float, sw: float) -> float:
    """
    Calculate three-phase oil relative permeability using Stone's method
    
    Parameters:
    -----------
    simulator : BOASTSimulator
        Main simulator object
    ireg : int
        Rock region number
    so : float
        Oil saturation
    sw : float
        Water saturation
        
    Returns:
    --------
    float
        Oil relative permeability for 3-phase flow
    """
    from block2 import Interpolation
    
    sowr = 1.0 - simulator.swr[ireg]
    sl = so + sw
    sg = 1.0 - sl
    
    # If no free gas, use 2-phase oil rel perm
    if sg <= 0.001:
        rkro = Interpolation.interp(
            simulator.sat[ireg, :simulator.msat[ireg]],
            simulator.krot[ireg, :simulator.msat[ireg]], so)
        return rkro
    
    # Three-phase calculation using Stone's method
    krow = Interpolation.interp(
        simulator.sat[ireg, :simulator.msat[ireg]],
        simulator.krot[ireg, :simulator.msat[ireg]], 1.0 - sw)
    
    krog = Interpolation.interp(
        simulator.sat[ireg, :simulator.msat[ireg]],
        simulator.krogt[ireg, :simulator.msat[ireg]], sl)
    
    krowr = Interpolation.interp(
        simulator.sat[ireg, :simulator.msat[ireg]],
        simulator.krot[ireg, :simulator.msat[ireg]], sowr)
    
    krw = Interpolation.interp(
        simulator.sat[ireg, :simulator.msat[ireg]],
        simulator.krwt[ireg, :simulator.msat[ireg]], sw)
    
    krg = Interpolation.interp(
        simulator.sat[ireg, :simulator.msat[ireg]],
        simulator.krgt[ireg, :simulator.msat[ireg]], sg)
    
    # Stone's three-phase relative permeability formula
    rkro = ((krow + krw) * (krog + krg)) / krowr - (krw + krg)
    
    # Limit to [0, 1]
    if rkro > 1.0:
        rkro = 1.0
    if rkro < 0.0:
        rkro = 0.0
        
    return rkro


class Initialization:
    """Initialize pressure and saturation distributions"""
    
    def __init__(self, simulator):
        """Initialize with reference to main simulator"""
        self.sim = simulator
    
    def uinitl(self, ii: int, jj: int, kk: int):
        """
        Initialize pressure and saturation distributions
        
        Reads initialization data and sets initial conditions:
        - Equilibrium pressure initialization or user-specified pressures
        - Initial saturation distributions
        - Initial fluid in place calculations
        """
        iread = self.sim.iread
        iocode = self.sim.outfile
        
        # Read initialization options
        iread.readline()  # Skip comment "EQUILIBRIUM PRESSURE INITIALIZATION"
        line = iread.readline()
        values = parse_fortran_line(line)
        kpi = int(values[0])  # Pressure init option (0=equilibrium, 1=user-specified)
        ksi = int(values[1])  # Saturation init option (0=constant, 1=user-specified)
        pdatum = values[2]    # Datum depth (ft)
        grad = values[3]      # Pressure gradient (psi/ft)
        
        self.sim.pdatum = pdatum
        self.sim.grad = grad
        
        # Read equilibrium data or pressure array
        if kpi == 0:
            # Equilibrium initialization - read reference pressure and contacts
            line = iread.readline()
            values = parse_fortran_line(line)
            pi = values[0]     # Initial pressure at datum
            woc = values[1]    # Water-oil contact (ft)
            pgoc = values[2]   # Pressure at gas-oil contact (psi)
            goc = values[3]    # Gas-oil contact (ft)
            
            # Calculate equilibrium pressures based on depth and fluid contacts
            for k in range(kk):
                for j in range(jj):
                    for i in range(ii):
                        depth = self.sim.el[i, j, k]
                        # Simple gradient-based initialization
                        self.sim.pn[i, j, k] = pi + grad * (depth - woc) / 144.0
                        self.sim.p[i, j, k] = self.sim.pn[i, j, k]
        else:
            # User-specified pressure initialization - read full 3D array
            for k in range(kk):
                for j in range(jj):
                    line = iread.readline()
                    values = parse_fortran_line(line)
                    for i in range(ii):
                        if i < len(values):
                            self.sim.pn[i, j, k] = values[i]
                            self.sim.p[i, j, k] = values[i]
        
        # Read saturation initialization
        line = iread.readline()
        values = parse_fortran_line(line)
        
        if ksi == 0:
            # Constant saturation initialization
            soi = values[0]
            swi = values[1]
            sgi = values[2]
            
            for k in range(kk):
                for j in range(jj):
                    for i in range(ii):
                        self.sim.son[i, j, k] = soi
                        self.sim.so[i, j, k] = soi
                        self.sim.swn[i, j, k] = swi
                        self.sim.sw[i, j, k] = swi
                        self.sim.sgn[i, j, k] = sgi
                        self.sim.sg[i, j, k] = sgi
        else:
            # User-specified saturation initialization - read arrays
            # Read oil saturation
            for k in range(kk):
                for j in range(jj):
                    line = iread.readline()
                    values = parse_fortran_line(line)
                    for i in range(ii):
                        if i < len(values):
                            self.sim.so[i, j, k] = values[i]
                            self.sim.son[i, j, k] = values[i]
            
            # Read water saturation
            for k in range(kk):
                for j in range(jj):
                    line = iread.readline()
                    values = parse_fortran_line(line)
                    for i in range(ii):
                        if i < len(values):
                            self.sim.sw[i, j, k] = values[i]
                            self.sim.swn[i, j, k] = values[i]
            
            # Calculate gas saturation
            for k in range(kk):
                for j in range(jj):
                    for i in range(ii):
                        self.sim.sg[i, j, k] = 1.0 - self.sim.so[i, j, k] - self.sim.sw[i, j, k]
                        if self.sim.sg[i, j, k] < 0.0:
                            self.sim.sg[i, j, k] = 0.0
                        self.sim.sgn[i, j, k] = self.sim.sg[i, j, k]
        
        iocode.write("\n   INITIALIZATION SUCCESSFULLY COMPLETED\n")

