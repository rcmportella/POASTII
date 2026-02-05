"""
BLOCK4.FOR - Gas Properties, Repressurization, and 1D Solver

This module contains:
1. GasProperties class: Real gas property calculations using correlations
2. Repressurization class: PVT recalculation during repressurization
3. TridiagonalSolver class: 1D Gaussian elimination for tridiagonal systems

Converted from BOAST II (Release 1.2) Fortran code
"""

import numpy as np
from typing import TextIO, Tuple, List


class GasProperties:
    """
    Real gas property calculations using Standing-Katz correlations
    
    Calculates gas properties from composition:
    - Z-factor (compressibility)
    - Formation volume factor (Bg)
    - Gas compressibility (Cg)
    - Gas viscosity (μg)
    - Pseudo-pressure
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
        
        # Component critical properties
        self.ncomp = 12
        
        # Critical temperatures (°R)
        self.temci = np.array([
            672.37,  # H2S
            547.57,  # CO2
            227.27,  # N2
            343.04,  # C1
            549.76,  # C2
            665.68,  # C3
            734.65,  # iC4
            765.32,  # nC4
            828.77,  # iC5
            845.37,  # nC5
            913.37,  # C6
            1023.89  # C7+
        ])
        
        # Critical pressures (psia)
        self.prsci = np.array([
            1306.0,  # H2S
            1071.0,  # CO2
            493.00,  # N2
            667.80,  # C1
            707.80,  # C2
            616.30,  # C3
            529.10,  # iC4
            550.70,  # nC4
            490.40,  # iC5
            488.60,  # nC5
            436.90,  # C6
            360.60   # C7+
        ])
    
    def pseudo(self, np_idx: int, ifatal: int, iocode: TextIO) -> int:
        """
        Calculate real gas properties from composition
        
        Parameters:
        -----------
        np_idx : int
            PVT region index (0-based)
        ifatal : int
            Fatal error counter
        iocode : file object
            Output file handle
            
        Returns:
        --------
        ifatal : int
            Updated fatal error counter
        """
        iread = self.sim.iread
        
        # Read gas data
        line = iread.readline()
        parts = line.split()
        kodea = int(parts[0])
        mpgt = int(parts[1])
        tem = float(parts[2])  # Temperature (°F)
        spg = float(parts[3])  # Gas gravity
        
        self.sim.mpgt[np_idx] = mpgt
        
        # Read gas composition (mole fractions)
        line = iread.readline()
        frci = [float(x) for x in line.split()]
        
        # Pad to 12 components if needed
        while len(frci) < 12:
            frci.append(0.0)
        
        # Read custom C7+ properties if kodea >= 4
        if kodea >= 4:
            line = iread.readline()
            parts = line.split()
            self.prsci[11] = float(parts[0])
            self.temci[11] = float(parts[1])
            # rmwti[11] = float(parts[2])  # Molecular weight (not used here)
        
        # Setup
        ndata = mpgt
        pref = 14.7  # Reference pressure (psia)
        pmaxt = self.sim.pmaxt
        delp = (pmaxt - pref) / (mpgt - 1.0)
        
        # Extract key components
        cnch2s = frci[0] * 100.0  # H2S mole %
        cncco2 = frci[1] * 100.0  # CO2 mole %
        cncn2 = frci[2] * 100.0   # N2 mole %
        
        # Calculate pseudo-critical properties
        if kodea <= 2:
            # Standing correlation for gas gravity
            tempc = 170.491 + 307.344 * spg
            prspc = 709.604 - 58.718 * spg
        else:
            # Kay's mixing rule
            tempc = sum(self.temci[j] * frci[j] for j in range(self.ncomp))
            prspc = sum(self.prsci[j] * frci[j] for j in range(self.ncomp))
        
        # Write header
        iocode.write("\n\n\n        ****  REAL GAS PROPERTIES  ****\n\n")
        
        # Write composition
        iocode.write("\n GAS COMPOSITION (MOLE FRACTION):\n\n")
        iocode.write("   H2S    CO2     N2\n")
        iocode.write(f"  {frci[0]:5.4f}  {frci[1]:5.4f}  {frci[2]:5.4f}\n\n")
        iocode.write("    C1     C2     C3    IC4    NC4    IC5    NC5     C6    C7+\n")
        iocode.write(" " + "  ".join([f"{frci[i]:5.4f}" for i in range(3, 12)]) + "\n\n")
        
        # Write properties
        iocode.write("\n")
        iocode.write(f" RESERVOIR TEMPERATURE, DEGREES F   = {tem:10.2f}\n")
        iocode.write(f" GAS GRAVITY                        = {spg:10.4f}\n")
        iocode.write(f" PSEUDO-CRITICAL TEMP., DEGREES R   = {tempc:10.2f}\n")
        iocode.write(f" PSEUDO-CRITICAL PRESSURE, PSIA     = {prspc:10.2f}\n")
        iocode.write(f" MOLE PERCENT - HYDROGEN SULFIDE    = {cnch2s:10.2f}\n")
        iocode.write(f" MOLE PERCENT - CARBON DIOXIDE      = {cncco2:10.2f}\n")
        iocode.write(f" MOLE PERCENT - NITROGEN            = {cncn2:10.2f}\n\n")
        
        # Table header
        iocode.write(" PRESSURE    PSEUDO-PRESS   Z-FACTOR    GAS FVF    GAS COMP  GAS VISG\n")
        iocode.write("  (PSIA)     (PSIA**2/CP)              (RB/SCF)    (1/PSIA)     (CP)\n\n")
        
        # Calculate properties at each pressure
        prs = pref
        prssi = 0.0
        prssi1 = 0.0
        
        zedd = np.zeros(ndata)
        cmpgd = np.zeros(ndata)
        
        for i in range(ndata):
            if prs > 0.0:
                # Reduced properties
                temprd = (tem + 460.0) / tempc
                prsprd = prs / prspc
                
                # Calculate Z-factor and compressibility
                zed, cmpg, ierr = self.zandc(tem, tempc, prs, prspc, cnch2s, cncco2)
                ifatal += ierr
                if ierr > 0:
                    iocode.write(f"\n     -----ZANDC ERRORS (IERR = {ierr}) HAVE" +
                               " OCCURRED IN THE DEFAULT GAS PROPERTIES CALCULATION.\n" +
                               "          CHECK INPUT GAS DATA AND USER-SPECIFIED OPTIONS.\n\n")
                
                # Calculate viscosity
                visg, ierr = self.viscy(temprd, prsprd, spg, tem, cnch2s, cncco2, cncn2)
                ifatal += ierr
                if ierr > 0:
                    iocode.write(f"\n     -----VISCY ERRORS (IERR = {ierr}) HAVE" +
                               " OCCURRED IN THE DEFAULT GAS PROPERTIES CALCULATION.\n" +
                               "          CHECK INPUT GAS DATA AND USER-SPECIFIED OPTIONS.\n\n")
                
                # Pseudo-pressure integration
                prssi2 = 2.0 * prs / (visg * zed)
                prssi = (prssi1 + prssi2) / 2.0 * delp + prssi
            else:
                zed = 1.0
                cmpg = 0.0
                visg = 0.0
                prssi2 = 0.0
            
            if prs <= pref:
                prssi = 0.0
            
            # Store values
            self.sim.psit[np_idx, i] = prssi
            zedd[i] = zed
            cmpgd[i] = cmpg
            self.sim.mugt[np_idx, i] = visg
            self.sim.pgt[np_idx, i] = prs
            
            # Calculate formation volume factor
            # BGT factor 0.004675 = 14.7/((460+60)*5.6146)
            if prs > 0.0:
                bgt_rb = (tem + 460.0) * zed * 0.004675 / prs
            else:
                bgt_rb = 0.0
            
            self.sim.bgt[np_idx, i] = bgt_rb * 5.6146  # Convert to RB/SCF
            
            # Write output
            iocode.write(f" {prs:7.1f}     {prssi:13.7e}   {zed:6.3f}  " +
                        f"{bgt_rb:10.4e}  {cmpg:10.4e}  {visg:7.6f}\n")
            
            # Next pressure
            prs += delp
            prssi1 = prssi2
        
        return ifatal
    
    def zandc(self, tem: float, tempc: float, prs: float, prspc: float,
              cnch2s: float, cncco2: float) -> Tuple[float, float, int]:
        """
        Calculate Z-factor and gas compressibility using Standing-Katz correlation
        
        Parameters:
        -----------
        tem : float
            Temperature (°F)
        tempc : float
            Pseudo-critical temperature (°R)
        prs : float
            Pressure (psia)
        prspc : float
            Pseudo-critical pressure (psia)
        cnch2s : float
            H2S mole percent
        cncco2 : float
            CO2 mole percent
            
        Returns:
        --------
        zed : float
            Z-factor (compressibility factor)
        cmpg : float
            Gas compressibility (1/psi)
        ierr : int
            Error flag
        """
        ierr = 0
        
        # Reduced properties
        tr = (tem + 460.0) / tempc
        pr = prs / prspc
        
        # Standing-Katz correlation with Wichert-Aziz correction
        # Correction for H2S and CO2
        eps = 120.0 * (cnch2s/100.0 + cncco2/100.0)**0.9 - \
              (cnch2s/100.0 + cncco2/100.0)**1.6 + 15.0 * (cnch2s/100.0)**0.5 - \
              (cnch2s/100.0)**4.0
        
        # Corrected pseudo-critical temperature
        tpc_corr = tempc - eps
        tr_corr = (tem + 460.0) / tpc_corr
        
        # Corrected pseudo-critical pressure
        ppc_corr = prspc * tpc_corr / (tempc + cnch2s/100.0 * (1.0 - cnch2s/100.0) * eps)
        pr_corr = prs / ppc_corr
        
        # Standing-Katz correlation (simplified polynomial fit)
        # Z = 1 + (A + B + C) * rho_r
        # where rho_r is reduced density
        
        # Initial guess for Z
        zed = 1.0
        
        # Newton-Raphson iteration
        for _ in range(10):
            rho_r = 0.27 * pr_corr / (zed * tr_corr)
            
            # Coefficients
            a1 = 0.3265
            a2 = -1.0700
            a3 = -0.5339
            a4 = 0.01569
            a5 = -0.05165
            a6 = 0.5475
            a7 = -0.7361
            a8 = 0.1844
            
            f1 = a1 + a2/tr_corr + a3/tr_corr**3 + a4/tr_corr**4 + a5/tr_corr**5
            f2 = a6 + a7/tr_corr + a8/tr_corr**2
            f3 = a7/tr_corr + a8/tr_corr**2
            
            z_new = 1.0 + f1*rho_r + f2*rho_r**2 - f3*rho_r**5
            
            if abs(z_new - zed) < 0.0001:
                zed = z_new
                break
            
            zed = z_new
        
        # Gas compressibility
        # cg = (1/P) - (1/Z)(dZ/dP)_T
        if prs > 0.0:
            # Numerical derivative
            dprs = 0.1
            pr_plus = (prs + dprs) / ppc_corr
            rho_r_plus = 0.27 * pr_plus / (zed * tr_corr)
            z_plus = 1.0 + f1*rho_r_plus + f2*rho_r_plus**2 - f3*rho_r_plus**5
            
            dz_dp = (z_plus - zed) / dprs
            cmpg = (1.0/prs) - (1.0/zed) * dz_dp
        else:
            cmpg = 0.0
        
        return zed, cmpg, ierr
    
    def viscy(self, temprd: float, prsprd: float, spg: float, tem: float,
              cnch2s: float, cncco2: float, cncn2: float) -> Tuple[float, int]:
        """
        Calculate gas viscosity using Lee-Gonzalez-Eakin correlation
        
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
        visg : float
            Gas viscosity (cp)
        ierr : int
            Error flag
        """
        ierr = 0
        
        # Lee-Gonzalez-Eakin correlation
        temp_r = tem + 460.0  # Convert to °R
        
        # Molecular weight from gas gravity
        mw = 28.97 * spg
        
        # Gas density (lb/ft³) - approximate
        # Using ideal gas law: ρ = P*MW/(Z*R*T)
        # With corrections for real gas
        z_approx = 1.0 + 0.3 * prsprd  # Simplified
        rho = (prsprd * 14.7 * mw) / (10.73 * temp_r * z_approx) * 0.0764  # Convert to lb/ft³
        
        # Lee-Gonzalez-Eakin correlation
        k = (9.4 + 0.02*mw) * temp_r**1.5 / (209.0 + 19.0*mw + temp_r)
        x = 3.5 + 986.0/temp_r + 0.01*mw
        y = 2.4 - 0.2*x
        
        visg = 1e-4 * k * np.exp(x * (rho/62.4)**y)
        
        # Correction for impurities (simplified)
        # H2S, CO2, N2 corrections
        if cnch2s > 0.0 or cncco2 > 0.0:
            correction = 1.0 + 0.001 * (cnch2s + cncco2)
            visg *= correction
        
        return visg, ierr


class Repressurization:
    """Handles PVT recalculation during reservoir repressurization"""
    
    def __init__(self, simulator):
        """
        Initialize with reference to main simulator
        
        Parameters:
        -----------
        simulator : BOASTSimulator
            Reference to main simulator object
        """
        self.sim = simulator
    
    def reprs1(self, i: int, j: int, k: int) -> None:
        """
        Calculate repressurization PVT properties
        
        Updates bubble point pressure when pressure rises above original
        bubble point, assuming equilibration in the grid block
        
        Parameters:
        -----------
        i, j, k : int
            Grid block indices
        """
        # Get current pressure (limited to original bubble point)
        pp = self.sim.p[i, j, k]
        if pp > self.sim.pbot[i, j, k]:
            pp = self.sim.pbot[i, j, k]
        
        # Get PVT region
        ipvtr = self.sim.ipvt[i, j, k]
        
        # Interpolate properties at current pressure
        from block2 import Interpolation
        
        bbo = Interpolation.interp(
            self.sim.pot[ipvtr, :self.sim.mpot[ipvtr]],
            self.sim.bot[ipvtr, :self.sim.mpot[ipvtr]],
            pp
        )
        
        rso = Interpolation.interp(
            self.sim.pot[ipvtr, :self.sim.mpot[ipvtr]],
            self.sim.rsot[ipvtr, :self.sim.mpot[ipvtr]],
            pp
        )
        
        bbg = Interpolation.interp(
            self.sim.pgt[ipvtr, :self.sim.mpgt[ipvtr]],
            self.sim.bgt[ipvtr, :self.sim.mpgt[ipvtr]],
            pp
        )
        
        # Calculate new solution GOR
        # Gas that can dissolve = free gas converted to dissolved gas
        sgosog = 0.0
        if self.sim.so[i, j, k] > 0.0:
            sgosog = (self.sim.sg[i, j, k] * bbo) / (self.sim.so[i, j, k] * bbg)
        
        rsonew = rso + sgosog
        
        # Find new bubble point pressure from new solution GOR
        pbonew = Interpolation.interp(
            self.sim.rsot[ipvtr, :self.sim.mpot[ipvtr]],
            self.sim.pot[ipvtr, :self.sim.mpot[ipvtr]],
            rsonew
        )
        
        # Update bubble point pressure
        self.sim.pbot[i, j, k] = pbonew


class TridiagonalSolver:
    """
    Solves 1D pressure equations using Gaussian elimination
    
    For tridiagonal systems of the form:
    AZL(i)*U(i-1) + BZL(i)*U(i) + CZL(i)*U(i+1) = DZL(i)
    
    Uses Thomas algorithm for efficient solution
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
    
    def gaus1d(self, iocode: TextIO, ifatal: int, nx: int, ny: int, nz: int) -> int:
        """
        Solve 1D pressure equations
        
        Parameters:
        -----------
        iocode : file object
            Output file handle
        ifatal : int
            Fatal error counter
        nx, ny, nz : int
            Grid dimensions
            
        Returns:
        --------
        ifatal : int
            Updated fatal error counter
        """
        # Check for single block case
        if nx == 1 and ny == 1 and nz == 1:
            i, j, k = 0, 0, 0
            if self.sim.e[i, j, k] != 0.0:
                self.sim.p[i, j, k] = self.sim.b[i, j, k] / self.sim.e[i, j, k]
                return ifatal
            else:
                iocode.write("\n     -----COEF E(1,1,1) = 0 FOR" +
                           " NX=NY=NZ=1 CASE; PLEASE CHECK INPUT DATA.\n\n")
                return ifatal + 1
        
        # Check that only one dimension > 1
        ndim = sum([nx > 1, ny > 1, nz > 1])
        if ndim != 1:
            ifatal += 1
            iocode.write(f"\n     -----GAUS1D INCORRECTLY CALLED" +
                        f" WHEN NX,NY,NZ = {nx:5d}{ny:5d}{nz:5d}\n\n")
            return ifatal
        
        # Solve based on which dimension varies
        if nx > 1:
            self._solve_x_direction(nx)
        elif nz > 1:
            self._solve_z_direction(nz)
        else:  # ny > 1
            self._solve_y_direction(ny)
        
        return ifatal
    
    def _solve_x_direction(self, nx: int) -> None:
        """Solve tridiagonal system in X-direction"""
        j, k = 0, 0
        
        # Extract coefficients
        azl = np.array([self.sim.aw[i, j, k] for i in range(nx)])
        bzl = np.array([self.sim.e[i, j, k] for i in range(nx)])
        czl = np.array([self.sim.ae[i, j, k] for i in range(nx)])
        dzl = np.array([self.sim.b[i, j, k] for i in range(nx)])
        
        # Thomas algorithm
        uzl = self._thomas_algorithm(azl, bzl, czl, dzl, nx)
        
        # Store solution
        for i in range(nx):
            self.sim.p[i, j, k] = uzl[i]
    
    def _solve_y_direction(self, ny: int) -> None:
        """Solve tridiagonal system in Y-direction"""
        i, k = 0, 0
        
        # Extract coefficients
        azl = np.array([self.sim.as_[i, j, k] for j in range(ny)])
        bzl = np.array([self.sim.e[i, j, k] for j in range(ny)])
        czl = np.array([self.sim.an[i, j, k] for j in range(ny)])
        dzl = np.array([self.sim.b[i, j, k] for j in range(ny)])
        
        # Thomas algorithm
        uzl = self._thomas_algorithm(azl, bzl, czl, dzl, ny)
        
        # Store solution
        for j in range(ny):
            self.sim.p[i, j, k] = uzl[j]
    
    def _solve_z_direction(self, nz: int) -> None:
        """Solve tridiagonal system in Z-direction"""
        i, j = 0, 0
        
        # Extract coefficients
        azl = np.array([self.sim.at[i, j, k] for k in range(nz)])
        bzl = np.array([self.sim.e[i, j, k] for k in range(nz)])
        czl = np.array([self.sim.ab[i, j, k] for k in range(nz)])
        dzl = np.array([self.sim.b[i, j, k] for k in range(nz)])
        
        # Thomas algorithm
        uzl = self._thomas_algorithm(azl, bzl, czl, dzl, nz)
        
        # Store solution
        for k in range(nz):
            self.sim.p[i, j, k] = uzl[k]
    
    def _thomas_algorithm(self, a: np.ndarray, b: np.ndarray, 
                          c: np.ndarray, d: np.ndarray, n: int) -> np.ndarray:
        """
        Thomas algorithm for tridiagonal systems
        
        Solves: a[i]*u[i-1] + b[i]*u[i] + c[i]*u[i+1] = d[i]
        
        Parameters:
        -----------
        a : array
            Lower diagonal coefficients
        b : array
            Main diagonal coefficients
        c : array
            Upper diagonal coefficients
        d : array
            Right-hand side
        n : int
            System size
            
        Returns:
        --------
        u : array
            Solution vector
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
