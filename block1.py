"""
BLOCK1.FOR - Aquifer Models and Well Management
Contains:
- AQIN: Aquifer influence on pressure equation
- AQOUT: Aquifer influx rate calculations
- AQUI: Aquifer initialization
- MATBAL: Material balance calculations
- NODES: Well node management
"""

import numpy as np
from typing import Tuple, List, Optional
from dataclasses import dataclass


# Carter-Tracy dimensionless time correlation maximum values
TDCMAX = np.array([0., 0., 0.6, 5., 5., 10., 15., 30., 45., 70., 1000.])


@dataclass
class AquiferParameters:
    """Parameters for aquifer model"""
    aqcr: float = 0.0      # Aquifer rock compressibility (1/psi)
    aqcw: float = 0.0      # Aquifer water compressibility (1/psi)
    aqmuw: float = 0.0     # Aquifer water viscosity (cp)
    aqk: float = 0.0       # Aquifer permeability (md)
    aqphi: float = 0.0     # Aquifer porosity (fraction)
    aqh: float = 0.0       # Aquifer net thickness (ft)
    aqs: float = 0.0       # Aquifer/reservoir boundary interface (fraction)
    aqre: float = 0.0      # External reservoir radius (ft)
    aqkt: float = 0.0      # Carter-Tracy time parameter
    aqbeta: float = 0.0    # Carter-Tracy beta parameter


@dataclass
class WellData:
    """Well data structure"""
    welnam: str = ""       # Well name
    idwell: int = 0        # Well ID
    i: int = 0             # I-location
    j: int = 0             # J-location
    k: int = 0             # K-location (top layer)
    layer: int = 1         # Number of layers
    kip: int = 0           # Well control type
    qvo: float = 0.0       # Oil rate (STB/D)
    qvw: float = 0.0       # Water rate (RB/D)
    qvg: float = 0.0       # Gas rate (MCF/D)
    qvt: float = 0.0       # Total rate (RB/D)
    pwf: List[float] = None  # Bottom hole flowing pressure
    pid: List[float] = None  # Productivity index
    alit: float = 0.0      # A coefficient for rate-pressure relation
    blit: float = 0.0      # B coefficient for rate-pressure relation


class AquiferModel:
    """Aquifer model implementation"""
    
    def __init__(self, simulator):
        """Initialize with reference to main simulator"""
        self.sim = simulator
        
    def fptd(self, a0: float, a1: float, a2: float, a3: float, tdp: float) -> float:
        """
        Carter-Tracy dimensionless pressure function P(tD)
        P(tD) = a0 + a1*sqrt(tD) + a2*tD + a3*tD^1.5
        """
        sqrt_tdp = np.sqrt(tdp)
        return a0 + a1 * sqrt_tdp + a2 * tdp + a3 * tdp * sqrt_tdp
        
    def fdptd(self, a1: float, a2: float, a3: float, tdp: float) -> float:
        """
        Carter-Tracy dimensionless pressure derivative dP/dtD
        dP/dtD = 0.5*a1/sqrt(tD) + a2 + 1.5*a3*sqrt(tD)
        """
        sqrt_tdp = np.sqrt(tdp)
        return 0.5 * a1 / sqrt_tdp + a2 + 1.5 * a3 * sqrt_tdp
        
    def aqin(self, ii: int, jj: int, kk: int, delt: float, eti: float):
        """
        Aquifer model pressure equation terms
        Modifies matrix coefficients for aquifer influence
        """
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    iaqreg = self.sim.iaqreg[i, j, k]
                    if iaqreg < 1:
                        continue
                        
                    ppn = self.sim.pn[i, j, k]
                    ipvtr = self.sim.ipvt[i, j, k]
                    
                    # Interpolate water and gas properties
                    bbw = self.interp(ipvtr, self.sim.pwt, self.sim.bwt, ppn)
                    bbg = self.interp(ipvtr, self.sim.pgt, self.sim.bgt, ppn)
                    rsw = self.interp(ipvtr, self.sim.pwt, self.sim.rswt, ppn)
                    
                    factor = bbw - bbg * rsw
                    
                    # POT aquifer
                    if iaqreg == 1:
                        self.sim.paq[i, j, k] = ppn
                        if delt != 0:
                            self.sim.cpi1[i, j, k] = self.sim.cpiaq1[i, j, k] / delt
                        self.sim.cpi2[i, j, k] = 0.0
                        
                    # Steady-state aquifer
                    elif iaqreg == 2:
                        self.sim.cpi1[i, j, k] = self.sim.cpiaq1[i, j, k]
                        self.sim.cpi2[i, j, k] = 0.0
                        
                    # Carter-Tracy aquifer (iaqreg 3-11)
                    elif iaqreg >= 3:
                        etp = eti + delt
                        td = eti * self.sim.cpiaq1[i, j, k]
                        tdp = etp * self.sim.cpiaq1[i, j, k]
                        
                        # Select coefficients based on aquifer type
                        ptd, dptd = self._get_carter_tracy_values(iaqreg, tdp)
                        
                        # Calculate Carter-Tracy coefficients
                        denom = ptd - td * dptd
                        self.sim.cpi1[i, j, k] = (self.sim.cpiaq1[i, j, k] * 
                                                   self.sim.cpiaq2[i, j, k] / denom)
                        self.sim.cpi2[i, j, k] = (self.sim.cpiaq1[i, j, k] * 
                                                   (self.sim.cpiaq2[i, j, k] * 
                                                    (self.sim.piaq[i, j, k] - ppn) +
                                                    self.sim.cumew[i, j, k] * dptd) / denom)
                        self.sim.paq[i, j, k] = ppn
                    
                    # Modify pressure equation coefficients
                    self.sim.b[i, j, k] -= ((self.sim.cpi1[i, j, k] * 
                                             self.sim.paq[i, j, k] + 
                                             self.sim.cpi2[i, j, k]) * factor)
                    self.sim.e[i, j, k] -= self.sim.cpi1[i, j, k] * factor
                    
    def _get_carter_tracy_values(self, iaqreg: int, tdp: float) -> Tuple[float, float]:
        """Get Carter-Tracy P(tD) and dP/dtD values based on aquifer type"""
        # Ratio = 1.5
        if iaqreg == 3:
            ptd = self.fptd(0.10371, 1.66657, -0.04579, -0.01023, tdp)
            dptd = self.fdptd(1.66657, -0.04579, -0.01023, tdp)
        # Ratio = 2.0
        elif iaqreg == 4:
            ptd = self.fptd(0.30210, 0.68178, -0.01599, -0.01356, tdp)
            dptd = self.fdptd(0.68178, -0.01599, -0.01356, tdp)
        # Ratio = 3.0
        elif iaqreg == 5:
            ptd = self.fptd(0.51243, 0.29317, 0.01534, -0.06732, tdp)
            dptd = self.fdptd(0.29317, 0.01534, -0.06732, tdp)
        # Ratio = 4.0
        elif iaqreg == 6:
            ptd = self.fptd(0.63656, 0.16101, 0.15812, -0.09104, tdp)
            dptd = self.fdptd(0.16101, 0.15812, -0.09104, tdp)
        # Ratio = 5.0
        elif iaqreg == 7:
            ptd = self.fptd(0.65106, 0.10414, 0.30953, -0.11258, tdp)
            dptd = self.fdptd(0.10414, 0.30953, -0.11258, tdp)
        # Ratio = 6.0
        elif iaqreg == 8:
            ptd = self.fptd(0.63367, 0.06940, 0.41750, -0.11137, tdp)
            dptd = self.fdptd(0.06940, 0.41750, -0.11137, tdp)
        # Ratio = 8.0
        elif iaqreg == 9:
            ptd = self.fptd(0.40132, 0.04104, 0.69592, -0.14350, tdp)
            dptd = self.fdptd(0.04104, 0.69592, -0.14350, tdp)
        # Ratio = 10.0
        elif iaqreg == 10:
            ptd = self.fptd(0.14386, 0.02649, 0.89646, -0.15502, tdp)
            dptd = self.fdptd(0.02649, 0.89646, -0.15502, tdp)
        # Ratio = Infinity
        elif iaqreg == 11:
            ptd = self.fptd(0.82092, -0.000368, 0.28908, 0.28817, tdp)
            dptd = self.fdptd(-0.000368, 0.28908, 0.28817, tdp)
        else:
            ptd = dptd = 0.0
            
        return ptd, dptd
        
    def aqout(self, ii: int, jj: int, kk: int, delt: float):
        """
        Calculate aquifer influx rates
        """
        # Initialize rock region aquifer rates
        nrock = self.sim.nrock
        qwaqr = np.zeros(nrock)
        
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    iaqreg = self.sim.iaqreg[i, j, k]
                    if iaqreg < 1:
                        continue
                        
                    ppn = self.sim.pn[i, j, k]
                    pp = self.sim.p[i, j, k]
                    
                    # POT aquifer
                    if iaqreg == 1:
                        self.sim.ewaq[i, j, k] = -self.sim.cpi1[i, j, k] * (ppn - pp)
                        
                    # Steady-state aquifer
                    elif iaqreg == 2:
                        self.sim.ewaq[i, j, k] = -self.sim.cpi1[i, j, k] * (
                            self.sim.paq[i, j, k] - pp)
                        
                    # Carter-Tracy aquifer
                    else:
                        self.sim.ewaq[i, j, k] = (-self.sim.cpi2[i, j, k] + 
                                                   self.sim.cpi1[i, j, k] * (pp - ppn))
                        
    def aqcum(self, nloop: int, delt: float):
        """Calculate cumulative aquifer influx"""
        ii, jj, kk = self.sim.ii, self.sim.jj, self.sim.kk
        nrock = self.sim.nrock
        
        # Initialize
        qwaqr = np.zeros(nrock)
        cumaqr = np.zeros(nrock)
        
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    iaqreg = self.sim.iaqreg[i, j, k]
                    if iaqreg < 1:
                        continue
                        
                    # Update cumulative
                    self.sim.cumew[i, j, k] += self.sim.ewaq[i, j, k] * delt
                    
                    # Convert SCF to STB
                    self.sim.qwaq[i, j, k] = self.sim.ewaq[i, j, k] / 5.615
                    self.sim.cumaqw[i, j, k] = self.sim.cumew[i, j, k] / 5.615
                    
                    # Rock region aquifer influx rates and cumulative
                    iroc = self.sim.irock[i, j, k]
                    qwaqr[iroc] += self.sim.qwaq[i, j, k]
                    cumaqr[iroc] += delt * self.sim.qwaq[i, j, k]
                    
        # Total run summary aquifer entries
        cumrat = np.sum(qwaqr)
        cumprd = np.sum(cumaqr)
        
        self.sim.saquir[nloop] = cumrat / 1000.0
        self.sim.saquic[nloop] = cumprd / 1.0e6
        
    def aqprnt(self, outfile):
        """Print aquifer influx information"""
        for n in range(self.sim.nrock):
            outfile.write(f'\nAQUIFER MODEL FOR ROCK REGION {n+1}:\n')
            outfile.write(f'AQUIFER INFLUX RATE (STB/D) = {self.sim.qwaqr[n]:9.1f}\n')
            outfile.write(f'CUM. AQUIFER INFLUX (STB) = {self.sim.cumaqr[n]:10.4e}\n\n')
            
    def aqui(self, infile, outfile, tmax: float):
        """
        Initialize aquifer model
        Reads aquifer parameters and assigns to grid blocks
        """
        # Skip comment line
        infile.readline()
        
        # Read aquifer option
        iaqopt = int(infile.readline().split()[0])
        self.sim.iaqopt = iaqopt
        
        if iaqopt == 0:
            return
            
        # POT aquifer
        if iaqopt == 1:
            outfile.write('\n\nPOT AQUIFER PARAMETERS:\n')
            outfile.write('\n  I1  I2  J1  J2  K1  K2      POT\n')
            
            naqen = int(infile.readline().split()[0])
            for n in range(naqen):
                line = infile.readline().split()
                i1, i2, j1, j2, k1, k2 = [int(x) for x in line[:6]]
                pot = float(line[6])
                
                outfile.write(f'{i1:4d}{i2:4d}{j1:4d}{j2:4d}{k1:4d}{k2:4d}{pot:10.2f}\n')
                
                for k in range(k1-1, k2):
                    for j in range(j1-1, j2):
                        for i in range(i1-1, i2):
                            self.sim.cpiaq1[i, j, k] = pot
                            self.sim.iaqreg[i, j, k] = iaqopt
            return
            
        # Steady-state aquifer
        if iaqopt == 2:
            outfile.write('\n\nSTEADY-STATE AQUIFER PARAMETERS:\n')
            outfile.write('\n  I1  I2  J1  J2  K1  K2     SSAQ\n')
            
            naqen = int(infile.readline().split()[0])
            for n in range(naqen):
                line = infile.readline().split()
                i1, i2, j1, j2, k1, k2 = [int(x) for x in line[:6]]
                ssaq = float(line[6])
                
                outfile.write(f'{i1:4d}{i2:4d}{j1:4d}{j2:4d}{k1:4d}{k2:4d}{ssaq:10.2f}\n')
                
                for k in range(k1-1, k2):
                    for j in range(j1-1, j2):
                        for i in range(i1-1, i2):
                            self.sim.cpiaq1[i, j, k] = ssaq
                            self.sim.iaqreg[i, j, k] = iaqopt
                            self.sim.paq[i, j, k] = self.sim.pn[i, j, k]
            return
            
        # Carter-Tracy aquifer (iaqopt 3-11)
        naqreg = int(infile.readline().split()[0])
        
        for nr in range(naqreg):
            outfile.write(f'\n\nCARTER-TRACY AQUIFER PARAMETERS FOR REGION {nr+1}:\n')
            
            # Write aquifer type
            ratio_labels = {
                3: '  1.5', 4: '  2.0', 5: '  3.0', 6: '  4.0',
                7: '  5.0', 8: '  6.0', 9: '  8.0', 10: ' 10.0',
                11: 'INFINITE-ACTING AQUIFER CASE'
            }
            if iaqopt in ratio_labels:
                if iaqopt == 11:
                    outfile.write(f'                    {ratio_labels[iaqopt]}\n')
                else:
                    outfile.write(f'RE/RW                                                     {ratio_labels[iaqopt]}\n')
                    
            # Read aquifer parameters
            line = infile.readline().split()
            params = AquiferParameters(
                aqcr=float(line[0]),
                aqcw=float(line[1]),
                aqmuw=float(line[2]),
                aqk=float(line[3]),
                aqphi=float(line[4]),
                aqh=float(line[5]),
                aqs=float(line[6]),
                aqre=float(line[7])
            )
            
            # Write aquifer parameters
            outfile.write(f'AQ ROCK COMP (1/PSI)                                {params.aqcr:10.4e}\n')
            outfile.write(f'AQ WATER COMP (1/PSI)                               {params.aqcw:10.4e}\n')
            outfile.write(f'AQ WATER VISCOSITY (CP)                             {params.aqmuw:10.4f}\n')
            outfile.write(f'AQ PERMEABILITY (MD)                                {params.aqk:10.4f}\n')
            outfile.write(f'AQ POROSITY (FRACTION)                              {params.aqphi:10.4f}\n')
            outfile.write(f'AQ NET THICKNESS (FT)                               {params.aqh:10.4f}\n')
            outfile.write(f'AQ/RES BOUNDARY INTERFACE (FRACTION)                {params.aqs:10.4f}\n')
            outfile.write(f'EXTERNAL RES RADIUS (FT)                            {params.aqre:10.4f}\n\n')
            
            # Calculate Carter-Tracy parameters
            aqcomp = params.aqcr + params.aqcw
            aqkt = 0.00633 * params.aqk / (params.aqmuw * params.aqphi * aqcomp * 
                                            params.aqre * params.aqre)
            aqbeta = 6.1832 * params.aqphi * params.aqh * aqcomp * params.aqre * params.aqre * params.aqs
            
            outfile.write('\nREGION LIMITS AND C-T PARAMETERS:\n')
            outfile.write('  I1  I2  J1  J2  K1  K2  AQ PAR 1  AQ PAR 2\n')
            
            # Read region limits
            naqen = int(infile.readline().split()[0])
            for n in range(naqen):
                line = infile.readline().split()
                i1, i2, j1, j2, k1, k2 = [int(x) for x in line[:6]]
                
                outfile.write(f'{i1:4d}{i2:4d}{j1:4d}{j2:4d}{k1:4d}{k2:4d}{aqkt:10.3e}{aqbeta:10.3e}\n')
                
                for k in range(k1-1, k2):
                    for j in range(j1-1, j2):
                        for i in range(i1-1, i2):
                            self.sim.cpiaq1[i, j, k] = aqkt
                            self.sim.cpiaq2[i, j, k] = aqbeta
                            self.sim.iaqreg[i, j, k] = iaqopt
                            self.sim.piaq[i, j, k] = self.sim.pn[i, j, k]
                            
            # Check dimensionless time limits
            tdmax = tmax * aqkt
            if TDCMAX[iaqopt] < tdmax:
                outfile.write(f'\nMAX DIMENSIONLESS TIME {tdmax:10.4e} ')
                outfile.write(f'EXCEEDS CORRELATION MAX DIM TIME {TDCMAX[iaqopt]:6.1f} ')
                outfile.write(f'FOR AQUI OPTION {iaqopt}\n')
                
    def interp(self, region: int, xtable: np.ndarray, ytable: np.ndarray, 
               x: float) -> float:
        """Linear interpolation helper"""
        # This would need actual implementation based on table structure
        # Placeholder for now
        return 0.0


class MaterialBalance:
    """Material balance calculations"""
    
    def __init__(self, simulator):
        """Initialize with reference to main simulator"""
        self.sim = simulator
    
    @staticmethod
    def matbal(simulator, ii: int, jj: int, kk: int, stbo: float, stboi: float,
               stbw: float, stbwi: float, towip: float, tooip: float,
               togip: float, mcfgi: float, delt0: float, d5615: float) -> dict:
        """
        Calculate material balance errors and average reservoir pressure
        Returns dictionary with production/injection rates and material balance errors
        """
        fact = d5615 * delt0
        pavg0 = 0.0
        pavg = 0.0
        op = 0.0  # Oil production this time step
        wp = 0.0  # Water production this time step
        gp = 0.0  # Gas production this time step
        wi = 0.0  # Water injection this time step
        gi = 0.0  # Gas injection this time step
        
        resvol = 0.0
        
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    vp = simulator.vp[i, j, k]
                    resvol += vp
                    
                    pavg0 += simulator.pn[i, j, k] * vp
                    pavg += simulator.p[i, j, k] * vp
                    
                    op += simulator.qo[i, j, k] * fact
                    
                    qw = simulator.qw[i, j, k]
                    if qw > 0.0:
                        wp += qw * fact
                    else:
                        wi += qw * fact + simulator.ewaq[i, j, k] * fact
                        
                    if simulator.iaqopt != 0 and qw > 0.0:
                        wi += simulator.ewaq[i, j, k] * fact
                        
                    qg = simulator.qg[i, j, k]
                    if qg > 0.0:
                        gp += qg * delt0
                    elif qg < 0.0:
                        gi += qg * delt0
                        
        # Update cumulatives
        cop = simulator.cop + op
        cwp = simulator.cwp + wp
        cgp = simulator.cgp + gp * 0.001
        cwi = simulator.cwi + wi
        cgi = simulator.cgi + gi * 0.001
        
        # Convert SCF to MCF
        gp = gp * 0.001
        gi = gi * 0.001
        
        # Calculate rates
        div = 1.0 / delt0
        opr = op * div
        wpr = wp * div
        gpr = gp * div
        wir = wi * div
        gir = gi * div
        
        # Average pressures
        pavg = pavg / resvol if resvol > 0 else 0.0
        pavg0 = pavg0 / resvol if resvol > 0 else 0.0
        
        # Material balance errors for this time step
        mbeo = mbew = mbeg = 0.0
        
        denom1 = stboi - op
        if abs(denom1 - stbo) >= 1.0e-4:
            mbeo = (stbo / (stboi - op) - 1.0) * 100.0
            
        denom2 = stbwi - wp - wi
        if abs(denom2 - stbw) >= 1.0e-4:
            mbew = (stbw / (stbwi - wp - wi) - 1.0) * 100.0
            
        mcfgt = simulator.mcfgt
        denom3 = mcfgi - gp - gi
        if abs(denom3 - mcfgt) >= 1.0e-4:
            mbeg = (mcfgt / (mcfgi - gp - gi) - 1.0) * 100.0
            
        # Cumulative material balance errors
        cmbeo = cmbew = cmbeg = 0.0
        
        denomo = tooip * 1.0e6 - cop
        if abs(denomo - stbo) >= 1.0e-4:
            cmbeo = (stbo / (tooip * 1.0e6 - cop) - 1.0) * 100.0
            
        denomw = towip * 1.0e6 - cwp - cwi
        if abs(denomw - stbw) >= 1.0e-4:
            cmbew = (stbw / (towip * 1.0e6 - cwp - cwi) - 1.0) * 100.0
            
        denomg = togip * 1.0e6 - cgp - cgi
        if abs(denomg - mcfgt) >= 1.0e-4:
            cmbeg = (mcfgt / (togip * 1.0e6 - cgp - cgi) - 1.0) * 100.0
            
        # Gas-oil and water-oil ratios
        gor = gp / op if op != 0.0 else 0.0
        wor = wp / op if op != 0.0 else 0.0
        
        cgor = cgp / cop * 1000.0 if cop != 0.0 else 0.0
        cwor = cwp / cop if cop != 0.0 else 0.0
        
        return {
            'op': op, 'wp': wp, 'gp': gp, 'wi': wi, 'gi': gi,
            'cop': cop, 'cwp': cwp, 'cgp': cgp, 'cwi': cwi, 'cgi': cgi,
            'opr': opr, 'wpr': wpr, 'gpr': gpr, 'wir': wir, 'gir': gir,
            'pavg': pavg, 'pavg0': pavg0,
            'mbeo': mbeo, 'mbew': mbew, 'mbeg': mbeg,
            'cmbeo': cmbeo, 'cmbew': cmbew, 'cmbeg': cmbeg,
            'gor': gor, 'wor': wor, 'cgor': cgor, 'cwor': cwor,
            'resvol': resvol
        }


class WellManager:
    """Well node management"""
    
    def __init__(self, simulator):
        """Initialize with reference to main simulator"""
        self.sim = simulator
        self.wells = []
        
    def nodes(self, infile, outfile) -> int:
        """
        Read well data and establish rate-specified and pressure-specified wells
        Returns number of wells (nvqn)
        """
        # Skip comment line
        infile.readline()
        
        # Read number of new and old wells
        line = infile.readline().split()
        nwelln = int(line[0])  # New wells
        nwello = int(line[1])  # Wells with only rate/pressure updates
        
        if nwelln == 0 and nwello == 0:
            return 0
            
        nchang = nwelln + nwello
        
        # Write header
        outfile.write('\n' + '*' * 21 + '\n')
        outfile.write('*   WELL   UPDATE   *\n')
        outfile.write('*' * 21 + '\n\n')
        outfile.write('                    RESERVOIR CONTAINS THE FOLLOWING RATE NODES:\n\n')
        outfile.write(' WELL  NO.     NODE      OIL(STBD)   WATER(RBD)   GAS(MCFD)   ')
        outfile.write('TOTAL(RBD)   BHFP(PSIA)      PID      ALIT      BLIT\n')
        
        ncount = []
        idmax = 0
        
        # Read new wells
        if nwelln > 0:
            infile.readline()  # Skip comment
            
            for j in range(nwelln):
                # Read well identification
                line = infile.readline().split()
                welnam = line[0]
                i1 = int(line[1])  # Well ID
                i2 = int(line[2])  # I-location
                i3 = int(line[3])  # J-location
                i4 = int(line[4])  # K-location (top layer)
                i5 = int(line[5])  # Number of layers
                
                if idmax < i1:
                    idmax = i1
                    
                # Read pressure data for each layer
                iq3 = i4 - 1  # Convert to 0-based
                lay = iq3 + i5
                
                pid_line = infile.readline().split()
                pwf_line = infile.readline().split()
                
                pid_values = [float(x) for x in pid_line]
                pwf_values = [float(x) for x in pwf_line]
                
                # Read rate data
                rate_line = infile.readline().split()
                kip = int(rate_line[2])
                qvo = float(rate_line[3])
                qvw = float(rate_line[4])
                qvg = float(rate_line[5])
                qvt = float(rate_line[6])
                
                # Read additional coefficients if needed
                alit = blit = 0.0
                if kip == -4:
                    coef_line = infile.readline().split()
                    alit = float(coef_line[2])
                    blit = float(coef_line[3])
                    
                # Create well object
                well = WellData(
                    welnam=welnam,
                    idwell=i1,
                    i=i2-1, j=i3-1, k=iq3,  # Convert to 0-based
                    layer=i5,
                    kip=kip,
                    qvo=qvo, qvw=qvw, qvg=qvg, qvt=qvt,
                    pwf=pwf_values,
                    pid=pid_values,
                    alit=alit, blit=blit
                )
                
                self.wells.append(well)
                ncount.append(i1)
                
                # Write well information
                for k in range(iq3, lay):
                    qwv = qvw if kip != -1 else 0.0
                    outfile.write(f' {welnam:5s} {i1:5d} {i2:3d}{i3:3d}{k+1:3d}')
                    outfile.write(f'   {qvo:11.2f}{qwv:13.2f}{qvg:13.2f}{qvt:13.2f}')
                    outfile.write(f'{pwf_values[k-iq3]:13.2f}{pid_values[k-iq3]:10.4f}')
                    outfile.write(f'{alit:10.3e}{blit:10.3e}\n')
                    
        # Read updates for existing wells
        if nwello > 0:
            infile.readline()  # Skip comment
            
            for nc in range(nwello):
                line = infile.readline().split()
                welnam = line[0]
                i1 = int(line[1])
                kip = int(line[2])
                qvo = float(line[3])
                qvw = float(line[4])
                qvg = float(line[5])
                qvt = float(line[6])
                
                # Find existing well and update
                for well in self.wells:
                    if well.idwell == i1:
                        well.kip = kip
                        well.qvo = qvo
                        well.qvw = qvw
                        well.qvg = qvg
                        well.qvt = qvt
                        break
                        
                ncount.append(i1)
                
                # Write updated well information
                for well in self.wells:
                    if well.idwell == i1:
                        for k in range(well.k, well.k + well.layer):
                            outfile.write(f' {well.welnam:5s} {i1:5d} {well.i+1:3d}{well.j+1:3d}{k+1:3d}')
                            outfile.write(f'   {qvo:11.2f}{qvw:13.2f}{qvg:13.2f}{qvt:13.2f}')
                            outfile.write(f'{well.pwf[k-well.k]:13.2f}{well.pid[k-well.k]:10.4f}')
                            outfile.write(f'{well.alit:10.3e}{well.blit:10.3e}\n')
                        break
                        
        outfile.write('\n')
        
        # Write well control type descriptions
        self._write_well_descriptions(outfile, ncount)
        
        # Copy well data to simulator arrays
        for idx, well in enumerate(self.wells):
            self.sim.iqn1[idx] = well.i
            self.sim.iqn2[idx] = well.j
            self.sim.iqn3[idx] = well.k
            self.sim.layer[idx] = well.layer
            self.sim.kip[idx] = well.kip
            self.sim.idwell[idx] = well.idwell
            self.sim.qvo[idx] = well.qvo
            self.sim.qvw[idx] = well.qvw
            self.sim.qvg[idx] = well.qvg
            self.sim.qvt[idx] = well.qvt
            self.sim.alit[idx] = well.alit
            self.sim.blit[idx] = well.blit
            
            # Copy PID and PWF arrays for each layer
            for k in range(well.layer):
                if k < len(well.pid):
                    self.sim.pid[idx, well.k + k] = well.pid[k]
                if k < len(well.pwf):
                    self.sim.pwf[idx, well.k + k] = well.pwf[k]
        
        return idmax
        
    def _write_well_descriptions(self, outfile, ncount: List[int]):
        """Write descriptions of well control types"""
        for well_id in ncount:
            for well in self.wells:
                if well.idwell == well_id:
                    i, j, k = well.i + 1, well.j + 1, well.k + 1
                    
                    if well.kip == 1 and well.qvo > 0.0:
                        outfile.write(f'               BLOCK {i:3d}{j:3d}{k:3d} CONTAINS THE OIL RATE ')
                        outfile.write(f'SPECIFIED PRODUCING WELL NUMBER {well.idwell:5d}\n')
                    elif well.kip == 1 and well.qvt > 0.0:
                        outfile.write(f'               BLOCK {i:3d}{j:3d}{k:3d} CONTAINS THE TOTAL RATE ')
                        outfile.write(f'SPECIFIED PRODUCING WELL NUMBER {well.idwell:5d}\n')
                    elif well.kip == 1 and well.qvw > 0.0:
                        outfile.write(f'               BLOCK {i:3d}{j:3d}{k:3d} CONTAINS THE WATER RATE ')
                        outfile.write(f'SPECIFIED PRODUCING WELL NUMBER {well.idwell:5d}\n')
                    elif well.kip == 1 and well.qvg > 0.0:
                        outfile.write(f'               BLOCK {i:3d}{j:3d}{k:3d} CONTAINS THE GAS RATE ')
                        outfile.write(f'SPECIFIED PRODUCING WELL NUMBER {well.idwell:5d}\n')
                    elif well.kip == 2:
                        outfile.write(f'               BLOCK {i:3d}{j:3d}{k:3d} CONTAINS THE RATE ')
                        outfile.write(f'SPECIFIED WATER INJECTION WELL NUMBER {well.idwell:5d}\n')
                    elif well.kip == 3:
                        outfile.write(f'               BLOCK {i:3d}{j:3d}{k:3d} CONTAINS THE RATE ')
                        outfile.write(f'SPECIFIED GAS INJECTION WELL NUMBER {well.idwell:5d}\n')
                    elif well.kip == -1:
                        outfile.write(f'               BLOCK {i:3d}{j:3d}{k:3d} CONTAINS THE ')
                        outfile.write(f'EXPLICIT PRESSURE SPECIFIED PRODUCING WELL NUMBER {well.idwell:5d}\n')
                    elif well.kip == -2:
                        outfile.write(f'               BLOCK {i:3d}{j:3d}{k:3d} CONTAINS THE ')
                        outfile.write(f'EXPLICIT PRESSURE SPECIFIED WATER INJECTION WELL NUMBER {well.idwell:5d}\n')
                    elif well.kip == -3:
                        outfile.write(f'               BLOCK {i:3d}{j:3d}{k:3d} CONTAINS THE ')
                        outfile.write(f'EXPLICIT PRESSURE SPECIFIED GAS INJECTION WELL NUMBER {well.idwell:5d}\n')
                    elif well.kip == -4:
                        outfile.write(f'               BLOCK {i:3d}{j:3d}{k:3d} CONTAINS THE ')
                        outfile.write(f'EXPLICIT PRESSURE SPECIFIED GAS PRODUCTION WELL NUMBER {well.idwell:5d}\n')
                    elif well.kip == -11:
                        outfile.write(f'               BLOCK {i:3d}{j:3d}{k:3d} CONTAINS THE ')
                        outfile.write(f'IMPLICIT PRESSURE SPECIFIED PRODUCING WELL NUMBER {well.idwell:5d}\n')
                    elif well.kip == -12:
                        outfile.write(f'               BLOCK {i:3d}{j:3d}{k:3d} CONTAINS THE ')
                        outfile.write(f'IMPLICIT PRESSURE SPECIFIED WATER INJECTION WELL NUMBER {well.idwell:5d}\n')
                    elif well.kip == -13:
                        outfile.write(f'               BLOCK {i:3d}{j:3d}{k:3d} CONTAINS THE ')
                        outfile.write(f'IMPLICIT PRESSURE SPECIFIED GAS INJECTION WELL NUMBER {well.idwell:5d}\n')
                    break
