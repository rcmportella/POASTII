"""
BLOCK7.FOR - Flow Equation Coefficients with Single Point Upstream

This module contains:
1. FlowEquation class: Calculate flow equation coefficients using single-point
   upstream weighting for relative permeabilities

Key features:
- Single-point upstream weighting for relative permeabilities
- Gravity and capillary pressure effects
- Three-phase oil relative permeability (Stone's method) option
- Calculation of transmissibility multipliers for oil, water, and gas
- Assembly of pressure equation coefficients (7-point stencil)
- Right-hand side vector with compressibility and well terms

Converted from BOAST II (Release 1.2) Fortran code
"""

import numpy as np
from typing import TextIO


class FlowEquation:
    """
    Calculate flow equation coefficients using single-point upstream weighting
    
    This class computes the coefficients for the pressure equation by:
    1. Calculating relative permeabilities and capillary pressures at each block
    2. Computing phase mobilities for each face using upstream weighting
    3. Calculating transmissibility multipliers including gravity and capillary effects
    4. Assembling the 7-point stencil coefficients (AW, AE, AS, AN, AT, AB)
    5. Computing the right-hand side vector with accumulation and well terms
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
    
    def solone(self, ii: int, jj: int, kk: int, div1: float, d288: float,
               ksm: int, ksm1: int, n: int, nn: int, kcoff: int,
               iocode: TextIO) -> None:
        """
        Calculate flow equation coefficients with single-point upstream weighting
        
        Parameters:
        -----------
        ii, jj, kk : int
            Grid dimensions
        div1 : float
            1.0 / time step size
        d288 : float
            Gravity constant (0.001127 for field units)
        ksm : int
            Solution method code
        ksm1 : int
            Debug output flag
        n : int
            Current time step number
        nn : int
            Total number of time steps
        kcoff : int
            Coefficient output flag
        iocode : file object
            Output file handle
        """
        from block2 import Interpolation
        from block6 import trikro
        
        interp_obj = Interpolation(self.sim)
        
        # Initialize temporary arrays for relative perms and capillary pressures
        rpw = np.zeros((ii, jj, kk))      # Water rel perm
        rpg = np.zeros((ii, jj, kk))      # Gas rel perm
        rpow = np.zeros((ii, jj, kk))     # Oil rel perm (2-phase)
        rpo3 = np.zeros((ii, jj, kk))     # Oil rel perm (3-phase)
        capow = np.zeros((ii, jj, kk))    # Oil-water capillary pressure
        capgo = np.zeros((ii, jj, kk))    # Gas-oil capillary pressure
        
        # Initialize solution gas ratios
        rso1 = rso2 = rso3 = rso4 = rso5 = rso6 = 0.0
        rsw1 = rsw2 = rsw3 = rsw4 = rsw5 = rsw6 = 0.0
        
        # Calculate relative permeabilities and capillary pressures for all blocks
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    if self.sim.vp[i, j, k] <= 0.0:
                        continue
                    
                    sso = self.sim.so[i, j, k]
                    ssw = self.sim.sw[i, j, k]
                    ssg = self.sim.sg[i, j, k]
                    irockb = self.sim.irock[i, j, k]
                    
                    # CRITICAL: Convert rock region to 0-based index for Python arrays
                    # irock values are 1,2,3,... but Python arrays need 0,1,2,...
                    irockb_0 = irockb - 1
                    
                    # Interpolate relative permeabilities
                    # interp(x_table, y_table, ireg, n, xo)
                    kro1 = interp_obj.interp(self.sim.sat, self.sim.krot, irockb_0,
                                            self.sim.msat[irockb_0], sso)
                    rpow[i, j, k] = kro1
                    
                    krw1 = interp_obj.interp(self.sim.sat, self.sim.krwt, irockb_0,
                                            self.sim.msat[irockb_0], ssw)
                    rpw[i, j, k] = krw1
                    
                    krg1 = interp_obj.interp(self.sim.sat, self.sim.krgt, irockb_0,
                                            self.sim.msat[irockb_0], ssg)
                    rpg[i, j, k] = krg1
                    
                    # Three-phase oil relative permeability (still uses 1-based irockb)
                    kro3 = trikro(self.sim, irockb, sso, ssw)
                    rpo3[i, j, k] = kro3
                    
                    # Capillary pressures
                    pcow = interp_obj.interp(self.sim.sat, self.sim.pcowt, irockb_0,
                                            self.sim.msat[irockb_0], ssw)
                    capow[i, j, k] = pcow
                    
                    pcgo = interp_obj.interp(self.sim.sat, self.sim.pcgot, irockb_0,
                                            self.sim.msat[irockb_0], ssg)
                    capgo[i, j, k] = pcgo
        
        # Calculate flow coefficients for all blocks
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    if self.sim.vp[i, j, k] <= 0.0:
                        continue
                    
                    # Initialize neighbor mobilities
                    mo1 = mo2 = mo3 = mo4 = mo5 = mo6 = 0.0
                    mw1 = mw2 = mw3 = mw4 = mw5 = mw6 = 0.0
                    mg1 = mg2 = mg3 = mg4 = mg5 = mg6 = 0.0
                    
                    # Initialize gravity and transmissibility terms
                    gow1 = gow2 = gow3 = gow4 = gow5 = gow6 = 0.0
                    gww1 = gww2 = gww3 = gww4 = gww5 = gww6 = 0.0
                    ggw1 = ggw2 = ggw3 = ggw4 = ggw5 = ggw6 = 0.0
                    aow = aoe = aos = aon = aot = aob = 0.0
                    aww = awe = aws = awn = awt = awb = 0.0
                    agw = age = ags = agn = agt = agb = 0.0
                    
                    # Current block properties
                    pp = self.sim.p[i, j, k]
                    bpt = self.sim.pbot[i, j, k]
                    ipvtr = self.sim.ipvt[i, j, k]
                    ipvtr_0 = ipvtr - 1  # Convert to 0-based
                    
                    # Interpolate PVT properties
                    rso = interp_obj.intpvt(self.sim.pot, self.sim.rsot, ipvtr_0,
                                           self.sim.mpot[ipvtr_0], pp, bpt,
                                           self.sim.rslope[ipvtr_0])
                    muo = interp_obj.intpvt(self.sim.pot, self.sim.muot, ipvtr_0,
                                           self.sim.mpot[ipvtr_0], pp, bpt,
                                           self.sim.vslope[ipvtr_0])
                    rsw = interp_obj.interp(self.sim.pwt, self.sim.rswt, ipvtr_0,
                                           self.sim.mpwt[ipvtr_0], pp)
                    muw = interp_obj.interp(self.sim.pwt, self.sim.muwt, ipvtr_0,
                                           self.sim.mpwt[ipvtr_0], pp)
                    mug = interp_obj.interp(self.sim.pgt, self.sim.mugt, ipvtr_0,
                                           self.sim.mpgt[ipvtr_0], pp)
                    
                    sso = self.sim.so[i, j, k]
                    ssw = self.sim.sw[i, j, k]
                    ssg = self.sim.sg[i, j, k]
                    irockb = self.sim.irock[i, j, k]
                    
                    pcow = capow[i, j, k]
                    pcgo = capgo[i, j, k]
                    
                    # Phase densities
                    ro = (self.sim.rhosco[ipvtr] + rso * self.sim.rhoscg[ipvtr]) / self.sim.bo[i, j, k]
                    rw = (self.sim.rhoscw[ipvtr] + rsw * self.sim.rhoscg[ipvtr]) / self.sim.bw[i, j, k]
                    rg = self.sim.rhoscg[ipvtr] / self.sim.bg[i, j, k]
                    
                    # Process I-1 neighbor (West)
                    if i > 0 and self.sim.vp[i-1, j, k] > 0.0:
                        p1 = self.sim.p[i-1, j, k]
                        bpt1 = self.sim.pbot[i-1, j, k]
                        ipvtr1 = self.sim.ipvt[i-1, j, k]
                        ipvtr1_0 = ipvtr1 - 1
                        
                        rso1 = interp_obj.intpvt(self.sim.pot, self.sim.rsot, ipvtr1_0,
                                                self.sim.mpot[ipvtr1_0], p1, bpt1,
                                                self.sim.rslope[ipvtr1_0])
                        muo1 = interp_obj.intpvt(self.sim.pot, self.sim.muot, ipvtr1_0,
                                                self.sim.mpot[ipvtr1_0], p1, bpt1,
                                                self.sim.vslope[ipvtr1_0])
                        rsw1 = interp_obj.interp(self.sim.pwt, self.sim.rswt, ipvtr1_0,
                                                self.sim.mpwt[ipvtr1_0], p1)
                        muw1 = interp_obj.interp(self.sim.pwt, self.sim.muwt, ipvtr1_0,
                                                self.sim.mpwt[ipvtr1_0], p1)
                        mug1 = interp_obj.interp(self.sim.pgt, self.sim.mugt, ipvtr1_0,
                                                self.sim.mpgt[ipvtr1_0], p1)
                        
                        pcow1 = capow[i-1, j, k]
                        pcgo1 = capgo[i-1, j, k]
                        
                        ro1 = (self.sim.rhosco[ipvtr1] + rso1 * self.sim.rhoscg[ipvtr1]) / self.sim.bo[i-1, j, k]
                        rw1 = (self.sim.rhoscw[ipvtr1] + rsw1 * self.sim.rhoscg[ipvtr1]) / self.sim.bw[i-1, j, k]
                        rg1 = self.sim.rhoscg[ipvtr1] / self.sim.bg[i-1, j, k]
                        
                        # Gravity and capillary terms
                        fact = -d288 * (self.sim.el[i-1, j, k] - self.sim.el[i, j, k])
                        gow1 = (ro1 + ro) * fact
                        gww1 = (rw1 + rw) * fact + pcow - pcow1
                        ggw1 = (rg1 + rg) * fact + pcgo1 - pcgo
                        
                        # Potential differences
                        p11 = p1 - pp
                        ho1 = p11 + gow1
                        hw1 = p11 + gww1
                        hg1 = p11 + ggw1
                        
                        # Upstream weighting for relative permeabilities
                        if self.sim.ithree[irockb] == 1:
                            kro1 = rpo3[i, j, k] if ho1 < 0.0 else rpo3[i-1, j, k]
                        else:
                            kro1 = rpow[i, j, k] if ho1 < 0.0 else rpow[i-1, j, k]
                        
                        krw1 = rpw[i, j, k] if hw1 < 0.0 else rpw[i-1, j, k]
                        krg1 = rpg[i, j, k] if hg1 < 0.0 else rpg[i-1, j, k]
                        
                        # Phase mobilities
                        mo1 = 4.0 * kro1 / ((self.sim.bo[i-1, j, k] + self.sim.bo[i, j, k]) * (muo1 + muo))
                        mw1 = 4.0 * krw1 / ((self.sim.bw[i-1, j, k] + self.sim.bw[i, j, k]) * (muw1 + muw))
                        mg1 = 4.0 * krg1 / ((self.sim.bg[i-1, j, k] + self.sim.bg[i, j, k]) * (mug1 + mug))
                    
                    aow = self.sim.tx[i, j, k] * mo1
                    aww = self.sim.tx[i, j, k] * mw1
                    agw = self.sim.tx[i, j, k] * mg1
                    
                    # Process I+1 neighbor (East)
                    if i < ii - 1 and self.sim.vp[i+1, j, k] > 0.0:
                        p2 = self.sim.p[i+1, j, k]
                        bpt2 = self.sim.pbot[i+1, j, k]
                        ipvtr2 = self.sim.ipvt[i+1, j, k]
                        ipvtr2_0 = ipvtr2 - 1
                        
                        rso2 = interp_obj.intpvt(self.sim.pot, self.sim.rsot, ipvtr2_0,
                                                self.sim.mpot[ipvtr2_0], p2, bpt2,
                                                self.sim.rslope[ipvtr2_0])
                        muo2 = interp_obj.intpvt(self.sim.pot, self.sim.muot, ipvtr2_0,
                                                self.sim.mpot[ipvtr2_0], p2, bpt2,
                                                self.sim.vslope[ipvtr2_0])
                        rsw2 = interp_obj.interp(self.sim.pwt, self.sim.rswt, ipvtr2_0,
                                                self.sim.mpwt[ipvtr2_0], p2)
                        muw2 = interp_obj.interp(self.sim.pwt, self.sim.muwt, ipvtr2_0,
                                                self.sim.mpwt[ipvtr2_0], p2)
                        mug2 = interp_obj.interp(self.sim.pgt, self.sim.mugt, ipvtr2_0,
                                                self.sim.mpgt[ipvtr2_0], p2)
                        
                        pcow2 = capow[i+1, j, k]
                        pcgo2 = capgo[i+1, j, k]
                        
                        ro2 = (self.sim.rhosco[ipvtr2] + rso2 * self.sim.rhoscg[ipvtr2]) / self.sim.bo[i+1, j, k]
                        rw2 = (self.sim.rhoscw[ipvtr2] + rsw2 * self.sim.rhoscg[ipvtr2]) / self.sim.bw[i+1, j, k]
                        rg2 = self.sim.rhoscg[ipvtr2] / self.sim.bg[i+1, j, k]
                        
                        fact = -d288 * (self.sim.el[i+1, j, k] - self.sim.el[i, j, k])
                        gow2 = (ro2 + ro) * fact
                        gww2 = (rw2 + rw) * fact + pcow - pcow2
                        ggw2 = (rg2 + rg) * fact + pcgo2 - pcgo
                        
                        p22 = p2 - pp
                        ho2 = p22 + gow2
                        hw2 = p22 + gww2
                        hg2 = p22 + ggw2
                        
                        if self.sim.ithree[irockb] == 1:
                            kro2 = rpo3[i, j, k] if ho2 < 0.0 else rpo3[i+1, j, k]
                        else:
                            kro2 = rpow[i, j, k] if ho2 < 0.0 else rpow[i+1, j, k]
                        
                        krw2 = rpw[i, j, k] if hw2 < 0.0 else rpw[i+1, j, k]
                        krg2 = rpg[i, j, k] if hg2 < 0.0 else rpg[i+1, j, k]
                        
                        mo2 = 4.0 * kro2 / ((self.sim.bo[i+1, j, k] + self.sim.bo[i, j, k]) * (muo2 + muo))
                        mw2 = 4.0 * krw2 / ((self.sim.bw[i+1, j, k] + self.sim.bw[i, j, k]) * (muw2 + muw))
                        mg2 = 4.0 * krg2 / ((self.sim.bg[i+1, j, k] + self.sim.bg[i, j, k]) * (mug2 + mug))
                    
                    aoe = self.sim.tx[i+1, j, k] * mo2
                    awe = self.sim.tx[i+1, j, k] * mw2
                    age = self.sim.tx[i+1, j, k] * mg2
                    
                    # Process J-1 neighbor (South)
                    if j > 0 and self.sim.vp[i, j-1, k] > 0.0:
                        p3 = self.sim.p[i, j-1, k]
                        bpt3 = self.sim.pbot[i, j-1, k]
                        ipvtr3 = self.sim.ipvt[i, j-1, k]
                        ipvtr3_0 = ipvtr3 - 1
                        
                        rso3 = interp_obj.intpvt(self.sim.pot, self.sim.rsot, ipvtr3_0,
                                                self.sim.mpot[ipvtr3_0], p3, bpt3,
                                                self.sim.rslope[ipvtr3_0])
                        muo3 = interp_obj.intpvt(self.sim.pot, self.sim.muot, ipvtr3_0,
                                                self.sim.mpot[ipvtr3_0], p3, bpt3,
                                                self.sim.vslope[ipvtr3_0])
                        rsw3 = interp_obj.interp(self.sim.pwt, self.sim.rswt, ipvtr3_0,
                                                self.sim.mpwt[ipvtr3_0], p3)
                        muw3 = interp_obj.interp(self.sim.pwt, self.sim.muwt, ipvtr3_0,
                                                self.sim.mpwt[ipvtr3_0], p3)
                        mug3 = interp_obj.interp(self.sim.pgt, self.sim.mugt, ipvtr3_0,
                                                self.sim.mpgt[ipvtr3_0], p3)
                        
                        pcow3 = capow[i, j-1, k]
                        pcgo3 = capgo[i, j-1, k]
                        
                        ro3 = (self.sim.rhosco[ipvtr3] + rso3 * self.sim.rhoscg[ipvtr3]) / self.sim.bo[i, j-1, k]
                        rw3 = (self.sim.rhoscw[ipvtr3] + rsw3 * self.sim.rhoscg[ipvtr3]) / self.sim.bw[i, j-1, k]
                        rg3 = self.sim.rhoscg[ipvtr3] / self.sim.bg[i, j-1, k]
                        
                        fact = -d288 * (self.sim.el[i, j-1, k] - self.sim.el[i, j, k])
                        gow3 = (ro3 + ro) * fact
                        gww3 = (rw3 + rw) * fact + pcow - pcow3
                        ggw3 = (rg3 + rg) * fact + pcgo3 - pcgo
                        
                        p33 = p3 - pp
                        ho3 = p33 + gow3
                        hw3 = p33 + gww3
                        hg3 = p33 + ggw3
                        
                        if self.sim.ithree[irockb] == 1:
                            kro3 = rpo3[i, j, k] if ho3 < 0.0 else rpo3[i, j-1, k]
                        else:
                            kro3 = rpow[i, j, k] if ho3 < 0.0 else rpow[i, j-1, k]
                        
                        krw3 = rpw[i, j, k] if hw3 < 0.0 else rpw[i, j-1, k]
                        krg3 = rpg[i, j, k] if hg3 < 0.0 else rpg[i, j-1, k]
                        
                        mo3 = 4.0 * kro3 / ((self.sim.bo[i, j-1, k] + self.sim.bo[i, j, k]) * (muo3 + muo))
                        mw3 = 4.0 * krw3 / ((self.sim.bw[i, j-1, k] + self.sim.bw[i, j, k]) * (muw3 + muw))
                        mg3 = 4.0 * krg3 / ((self.sim.bg[i, j-1, k] + self.sim.bg[i, j, k]) * (mug3 + mug))
                    
                    aos = self.sim.ty[i, j, k] * mo3
                    aws = self.sim.ty[i, j, k] * mw3
                    ags = self.sim.ty[i, j, k] * mg3
                    
                    # Process J+1 neighbor (North)
                    if j < jj - 1 and self.sim.vp[i, j+1, k] > 0.0:
                        p4 = self.sim.p[i, j+1, k]
                        bpt4 = self.sim.pbot[i, j+1, k]
                        ipvtr4 = self.sim.ipvt[i, j+1, k]
                        ipvtr4_0 = ipvtr4 - 1
                        
                        rso4 = interp_obj.intpvt(self.sim.pot, self.sim.rsot, ipvtr4_0,
                                                self.sim.mpot[ipvtr4_0], p4, bpt4,
                                                self.sim.rslope[ipvtr4_0])
                        muo4 = interp_obj.intpvt(self.sim.pot, self.sim.muot, ipvtr4_0,
                                                self.sim.mpot[ipvtr4_0], p4, bpt4,
                                                self.sim.vslope[ipvtr4_0])
                        rsw4 = interp_obj.interp(self.sim.pwt, self.sim.rswt, ipvtr4_0,
                                                self.sim.mpwt[ipvtr4_0], p4)
                        muw4 = interp_obj.interp(self.sim.pwt, self.sim.muwt, ipvtr4_0,
                                                self.sim.mpwt[ipvtr4_0], p4)
                        mug4 = interp_obj.interp(self.sim.pgt, self.sim.mugt, ipvtr4_0,
                                                self.sim.mpgt[ipvtr4_0], p4)
                        
                        pcow4 = capow[i, j+1, k]
                        pcgo4 = capgo[i, j+1, k]
                        
                        ro4 = (self.sim.rhosco[ipvtr4] + rso4 * self.sim.rhoscg[ipvtr4]) / self.sim.bo[i, j+1, k]
                        rw4 = (self.sim.rhoscw[ipvtr4] + rsw4 * self.sim.rhoscg[ipvtr4]) / self.sim.bw[i, j+1, k]
                        rg4 = self.sim.rhoscg[ipvtr4] / self.sim.bg[i, j+1, k]
                        
                        fact = -d288 * (self.sim.el[i, j+1, k] - self.sim.el[i, j, k])
                        gow4 = (ro4 + ro) * fact
                        gww4 = (rw4 + rw) * fact + pcow - pcow4
                        ggw4 = (rg4 + rg) * fact + pcgo4 - pcgo
                        
                        p44 = p4 - pp
                        ho4 = p44 + gow4
                        hw4 = p44 + gww4
                        hg4 = p44 + ggw4
                        
                        if self.sim.ithree[irockb] == 1:
                            kro4 = rpo3[i, j, k] if ho4 < 0.0 else rpo3[i, j+1, k]
                        else:
                            kro4 = rpow[i, j, k] if ho4 < 0.0 else rpow[i, j+1, k]
                        
                        krw4 = rpw[i, j, k] if hw4 < 0.0 else rpw[i, j+1, k]
                        krg4 = rpg[i, j, k] if hg4 < 0.0 else rpg[i, j+1, k]
                        
                        mo4 = 4.0 * kro4 / ((self.sim.bo[i, j+1, k] + self.sim.bo[i, j, k]) * (muo4 + muo))
                        mw4 = 4.0 * krw4 / ((self.sim.bw[i, j+1, k] + self.sim.bw[i, j, k]) * (muw4 + muw))
                        mg4 = 4.0 * krg4 / ((self.sim.bg[i, j+1, k] + self.sim.bg[i, j, k]) * (mug4 + mug))
                    
                    aon = self.sim.ty[i, j+1, k] * mo4
                    awn = self.sim.ty[i, j+1, k] * mw4
                    agn = self.sim.ty[i, j+1, k] * mg4
                    
                    # Process K-1 neighbor (Top)
                    if k > 0 and self.sim.vp[i, j, k-1] > 0.0:
                        p5 = self.sim.p[i, j, k-1]
                        bpt5 = self.sim.pbot[i, j, k-1]
                        ipvtr5 = self.sim.ipvt[i, j, k-1]
                        ipvtr5_0 = ipvtr5 - 1
                        
                        rso5 = interp_obj.intpvt(self.sim.pot, self.sim.rsot, ipvtr5_0,
                                                self.sim.mpot[ipvtr5_0], p5, bpt5,
                                                self.sim.rslope[ipvtr5_0])
                        muo5 = interp_obj.intpvt(self.sim.pot, self.sim.muot, ipvtr5_0,
                                                self.sim.mpot[ipvtr5_0], p5, bpt5,
                                                self.sim.vslope[ipvtr5_0])
                        rsw5 = interp_obj.interp(self.sim.pwt, self.sim.rswt, ipvtr5_0,
                                                self.sim.mpwt[ipvtr5_0], p5)
                        muw5 = interp_obj.interp(self.sim.pwt, self.sim.muwt, ipvtr5_0,
                                                self.sim.mpwt[ipvtr5_0], p5)
                        mug5 = interp_obj.interp(self.sim.pgt, self.sim.mugt, ipvtr5_0,
                                                self.sim.mpgt[ipvtr5_0], p5)
                        
                        pcow5 = capow[i, j, k-1]
                        pcgo5 = capgo[i, j, k-1]
                        
                        ro5 = (self.sim.rhosco[ipvtr5] + rso5 * self.sim.rhoscg[ipvtr5]) / self.sim.bo[i, j, k-1]
                        rw5 = (self.sim.rhoscw[ipvtr5] + rsw5 * self.sim.rhoscg[ipvtr5]) / self.sim.bw[i, j, k-1]
                        rg5 = self.sim.rhoscg[ipvtr5] / self.sim.bg[i, j, k-1]
                        
                        fact = -d288 * (self.sim.el[i, j, k-1] - self.sim.el[i, j, k])
                        gow5 = (ro5 + ro) * fact
                        gww5 = (rw5 + rw) * fact + pcow - pcow5
                        ggw5 = (rg5 + rg) * fact + pcgo5 - pcgo
                        
                        p55 = p5 - pp
                        ho5 = p55 + gow5
                        hw5 = p55 + gww5
                        hg5 = p55 + ggw5
                        
                        if self.sim.ithree[irockb] == 1:
                            kro5 = rpo3[i, j, k] if ho5 < 0.0 else rpo3[i, j, k-1]
                        else:
                            kro5 = rpow[i, j, k] if ho5 < 0.0 else rpow[i, j, k-1]
                        
                        krw5 = rpw[i, j, k] if hw5 < 0.0 else rpw[i, j, k-1]
                        krg5 = rpg[i, j, k] if hg5 < 0.0 else rpg[i, j, k-1]
                        
                        mo5 = 4.0 * kro5 / ((self.sim.bo[i, j, k-1] + self.sim.bo[i, j, k]) * (muo5 + muo))
                        mw5 = 4.0 * krw5 / ((self.sim.bw[i, j, k-1] + self.sim.bw[i, j, k]) * (muw5 + muw))
                        mg5 = 4.0 * krg5 / ((self.sim.bg[i, j, k-1] + self.sim.bg[i, j, k]) * (mug5 + mug))
                    
                    aot = self.sim.tz[i, j, k] * mo5
                    awt = self.sim.tz[i, j, k] * mw5
                    agt = self.sim.tz[i, j, k] * mg5
                    
                    # Process K+1 neighbor (Bottom)
                    if k < kk - 1 and self.sim.vp[i, j, k+1] > 0.0:
                        p6 = self.sim.p[i, j, k+1]
                        bpt6 = self.sim.pbot[i, j, k+1]
                        ipvtr6 = self.sim.ipvt[i, j, k+1]
                        ipvtr6_0 = ipvtr6 - 1
                        
                        rso6 = interp_obj.intpvt(self.sim.pot, self.sim.rsot, ipvtr6_0,
                                                self.sim.mpot[ipvtr6_0], p6, bpt6,
                                                self.sim.rslope[ipvtr6_0])
                        muo6 = interp_obj.intpvt(self.sim.pot, self.sim.muot, ipvtr6_0,
                                                self.sim.mpot[ipvtr6_0], p6, bpt6,
                                                self.sim.vslope[ipvtr6_0])
                        rsw6 = interp_obj.interp(self.sim.pwt, self.sim.rswt, ipvtr6_0,
                                                self.sim.mpwt[ipvtr6_0], p6)
                        muw6 = interp_obj.interp(self.sim.pwt, self.sim.muwt, ipvtr6_0,
                                                self.sim.mpwt[ipvtr6_0], p6)
                        mug6 = interp_obj.interp(self.sim.pgt, self.sim.mugt, ipvtr6_0,
                                                self.sim.mpgt[ipvtr6_0], p6)
                        
                        pcow6 = capow[i, j, k+1]
                        pcgo6 = capgo[i, j, k+1]
                        
                        ro6 = (self.sim.rhosco[ipvtr6] + rso6 * self.sim.rhoscg[ipvtr6]) / self.sim.bo[i, j, k+1]
                        rw6 = (self.sim.rhoscw[ipvtr6] + rsw6 * self.sim.rhoscg[ipvtr6]) / self.sim.bw[i, j, k+1]
                        rg6 = self.sim.rhoscg[ipvtr6] / self.sim.bg[i, j, k+1]
                        
                        fact = -d288 * (self.sim.el[i, j, k+1] - self.sim.el[i, j, k])
                        gow6 = (ro6 + ro) * fact
                        gww6 = (rw6 + rw) * fact + pcow - pcow6
                        ggw6 = (rg6 + rg) * fact + pcgo6 - pcgo
                        
                        p66 = p6 - pp
                        ho6 = p66 + gow6
                        hw6 = p66 + gww6
                        hg6 = p66 + ggw6
                        
                        if self.sim.ithree[irockb] == 1:
                            kro6 = rpo3[i, j, k] if ho6 < 0.0 else rpo3[i, j, k+1]
                        else:
                            kro6 = rpow[i, j, k] if ho6 < 0.0 else rpow[i, j, k+1]
                        
                        krw6 = rpw[i, j, k] if hw6 < 0.0 else rpw[i, j, k+1]
                        krg6 = rpg[i, j, k] if hg6 < 0.0 else rpg[i, j, k+1]
                        
                        mo6 = 4.0 * kro6 / ((self.sim.bo[i, j, k+1] + self.sim.bo[i, j, k]) * (muo6 + muo))
                        mw6 = 4.0 * krw6 / ((self.sim.bw[i, j, k+1] + self.sim.bw[i, j, k]) * (muw6 + muw))
                        mg6 = 4.0 * krg6 / ((self.sim.bg[i, j, k+1] + self.sim.bg[i, j, k]) * (mug6 + mug))
                    
                    aob = self.sim.tz[i, j, k+1] * mo6
                    awb = self.sim.tz[i, j, k+1] * mw6
                    agb = self.sim.tz[i, j, k+1] * mg6
                    
                    # Average solution gas ratios for each direction
                    rso1a = 0.5 * (rso1 + rso)
                    rso2a = 0.5 * (rso2 + rso)
                    rso3a = 0.5 * (rso3 + rso)
                    rso4a = 0.5 * (rso4 + rso)
                    rso5a = 0.5 * (rso5 + rso)
                    rso6a = 0.5 * (rso6 + rso)
                    
                    rsw1a = 0.5 * (rsw1 + rsw)
                    rsw2a = 0.5 * (rsw2 + rsw)
                    rsw3a = 0.5 * (rsw3 + rsw)
                    rsw4a = 0.5 * (rsw4 + rsw)
                    rsw5a = 0.5 * (rsw5 + rsw)
                    rsw6a = 0.5 * (rsw6 + rsw)
                    
                    # Gravity terms for each phase and direction
                    ao1 = aow * gow1 if i > 0 and self.sim.vp[i-1, j, k] > 0.0 else 0.0
                    ao2 = aoe * gow2 if i < ii - 1 and self.sim.vp[i+1, j, k] > 0.0 else 0.0
                    ao3 = aos * gow3 if j > 0 and self.sim.vp[i, j-1, k] > 0.0 else 0.0
                    ao4 = aon * gow4 if j < jj - 1 and self.sim.vp[i, j+1, k] > 0.0 else 0.0
                    ao5 = aot * gow5 if k > 0 and self.sim.vp[i, j, k-1] > 0.0 else 0.0
                    ao6 = aob * gow6 if k < kk - 1 and self.sim.vp[i, j, k+1] > 0.0 else 0.0
                    
                    aw1 = aww * gww1 if i > 0 and self.sim.vp[i-1, j, k] > 0.0 else 0.0
                    aw2 = awe * gww2 if i < ii - 1 and self.sim.vp[i+1, j, k] > 0.0 else 0.0
                    aw3 = aws * gww3 if j > 0 and self.sim.vp[i, j-1, k] > 0.0 else 0.0
                    aw4 = awn * gww4 if j < jj - 1 and self.sim.vp[i, j+1, k] > 0.0 else 0.0
                    aw5 = awt * gww5 if k > 0 and self.sim.vp[i, j, k-1] > 0.0 else 0.0
                    aw6 = awb * gww6 if k < kk - 1 and self.sim.vp[i, j, k+1] > 0.0 else 0.0
                    
                    # Total gravity weight terms
                    self.sim.gowt[i, j, k] = ao1 + ao2 + ao3 + ao4 + ao5 + ao6
                    self.sim.gwwt[i, j, k] = aw1 + aw2 + aw3 + aw4 + aw5 + aw6
                    self.sim.ggwt[i, j, k] = (agw * ggw1 + age * ggw2 + ags * ggw3 +
                                              agn * ggw4 + agt * ggw5 + agb * ggw6 +
                                              rso1a * ao1 + rso2a * ao2 + rso3a * ao3 +
                                              rso4a * ao4 + rso5a * ao5 + rso6a * ao6 +
                                              rsw1a * aw1 + rsw2a * aw2 + rsw3a * aw3 +
                                              rsw4a * aw4 + rsw5a * aw5 + rsw6a * aw6)
                    
                    # Right-hand side combining all phases
                    self.sim.qowg[i, j, k] = ((self.sim.bo[i, j, k] - self.sim.bg[i, j, k] * rso) *
                                             (-self.sim.gowt[i, j, k] + self.sim.qo[i, j, k]) +
                                             (self.sim.bw[i, j, k] - self.sim.bg[i, j, k] * rsw) *
                                             (-self.sim.gwwt[i, j, k] + self.sim.qw[i, j, k]) +
                                             self.sim.bg[i, j, k] *
                                             (-self.sim.ggwt[i, j, k] + self.sim.qg[i, j, k]))
                    
                    # Pressure equation coefficients (7-point stencil)
                    self.sim.aw[i, j, k] = ((self.sim.bo[i, j, k] + 0.5 * self.sim.bg[i, j, k] * (rso1 - rso)) * aow +
                                           (self.sim.bw[i, j, k] + 0.5 * self.sim.bg[i, j, k] * (rsw1 - rsw)) * aww +
                                           self.sim.bg[i, j, k] * agw)
                    
                    self.sim.ae[i, j, k] = ((self.sim.bo[i, j, k] + 0.5 * self.sim.bg[i, j, k] * (rso2 - rso)) * aoe +
                                           (self.sim.bw[i, j, k] + 0.5 * self.sim.bg[i, j, k] * (rsw2 - rsw)) * awe +
                                           self.sim.bg[i, j, k] * age)
                    
                    self.sim.as_[i, j, k] = ((self.sim.bo[i, j, k] + 0.5 * self.sim.bg[i, j, k] * (rso3 - rso)) * aos +
                                            (self.sim.bw[i, j, k] + 0.5 * self.sim.bg[i, j, k] * (rsw3 - rsw)) * aws +
                                            self.sim.bg[i, j, k] * ags)
                    
                    self.sim.an[i, j, k] = ((self.sim.bo[i, j, k] + 0.5 * self.sim.bg[i, j, k] * (rso4 - rso)) * aon +
                                           (self.sim.bw[i, j, k] + 0.5 * self.sim.bg[i, j, k] * (rsw4 - rsw)) * awn +
                                           self.sim.bg[i, j, k] * agn)
                    
                    self.sim.at[i, j, k] = ((self.sim.bo[i, j, k] + 0.5 * self.sim.bg[i, j, k] * (rso5 - rso)) * aot +
                                           (self.sim.bw[i, j, k] + 0.5 * self.sim.bg[i, j, k] * (rsw5 - rsw)) * awt +
                                           self.sim.bg[i, j, k] * agt)
                    
                    self.sim.ab[i, j, k] = ((self.sim.bo[i, j, k] + 0.5 * self.sim.bg[i, j, k] * (rso6 - rso)) * aob +
                                           (self.sim.bw[i, j, k] + 0.5 * self.sim.bg[i, j, k] * (rsw6 - rsw)) * awb +
                                           self.sim.bg[i, j, k] * agb)
                    
                    # Store phase mobilities for output
                    self.sim.ow[i, j, k] = aow
                    self.sim.oe[i, j, k] = aoe
                    self.sim.os[i, j, k] = aos
                    self.sim.on[i, j, k] = aon
                    self.sim.ot[i, j, k] = aot
                    self.sim.ob[i, j, k] = aob
                    
                    self.sim.ww[i, j, k] = aww
                    self.sim.we[i, j, k] = awe
                    self.sim.ws[i, j, k] = aws
                    self.sim.wn[i, j, k] = awn
                    self.sim.wt[i, j, k] = awt
                    self.sim.wb[i, j, k] = awb
                    
                    # Optional debug output
                    if kcoff != 0:
                        iocode.write("\n\n")
                        iocode.write(f"{i+1:3d}{j+1:3d}{k+1:3d}{mo1:15.6e}{mo2:15.6e}{mo3:15.6e}{mo4:15.6e}{mo5:15.6e}{mo6:15.6e}\n")
                        iocode.write(f"{i+1:3d}{j+1:3d}{k+1:3d}{mw1:15.6e}{mw2:15.6e}{mw3:15.6e}{mw4:15.6e}{mw5:15.6e}{mw6:15.6e}\n")
                        iocode.write(f"{i+1:3d}{j+1:3d}{k+1:3d}{mg1:15.6e}{mg2:15.6e}{mg3:15.6e}{mg4:15.6e}{mg5:15.6e}{mg6:15.6e}\n")
                        iocode.write(f"{i+1:3d}{j+1:3d}{k+1:3d}{aow:15.6e}{aoe:15.6e}{aos:15.6e}{aon:15.6e}{aot:15.6e}{aob:15.6e}" +
                                   f"{self.sim.bo[i, j, k]:15.6e}{rso:15.6e}\n")
                        iocode.write(f"{i+1:3d}{j+1:3d}{k+1:3d}{aww:15.6e}{awe:15.6e}{aws:15.6e}{awn:15.6e}{awt:15.6e}{awb:15.6e}" +
                                   f"{self.sim.bw[i, j, k]:15.6e}{rsw:15.6e}\n")
                        iocode.write(f"{i+1:3d}{j+1:3d}{k+1:3d}{agw:15.6e}{age:15.6e}{ags:15.6e}{agn:15.6e}{agt:15.6e}{agb:15.6e}" +
                                   f"{self.sim.bg[i, j, k]:15.6e}\n")
                        iocode.write(f"{i+1:3d}{j+1:3d}{k+1:3d}{self.sim.gowt[i, j, k]:15.6e}{self.sim.qo[i, j, k]:15.6e}" +
                                   f"{self.sim.gwwt[i, j, k]:15.6e}{self.sim.qw[i, j, k]:15.6e}" +
                                   f"{self.sim.ggwt[i, j, k]:15.6e}{self.sim.qg[i, j, k]:15.6e}{self.sim.qowg[i, j, k]:15.6e}\n")
        
        # Calculate main diagonal and RHS vector
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    self.sim.sum[i, j, k] = (self.sim.aw[i, j, k] + self.sim.ae[i, j, k] +
                                            self.sim.as_[i, j, k] + self.sim.an[i, j, k] +
                                            self.sim.at[i, j, k] + self.sim.ab[i, j, k])
                    
                    self.sim.gam[i, j, k] = self.sim.vp[i, j, k] * self.sim.ct[i, j, k] * div1
                    self.sim.e[i, j, k] = -self.sim.sum[i, j, k] - self.sim.gam[i, j, k]
                    self.sim.b[i, j, k] = self.sim.qowg[i, j, k] - self.sim.gam[i, j, k] * self.sim.pn[i, j, k]
        
        # Set coefficients for zero pore volume blocks (no pressure change)
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    if self.sim.vp[i, j, k] <= 0.0:
                        self.sim.e[i, j, k] = -1.0
                        self.sim.b[i, j, k] = -self.sim.pn[i, j, k]
                        self.sim.aw[i, j, k] = 0.0
                        self.sim.ae[i, j, k] = 0.0
                        self.sim.as_[i, j, k] = 0.0
                        self.sim.an[i, j, k] = 0.0
                        self.sim.at[i, j, k] = 0.0
                        self.sim.ab[i, j, k] = 0.0
        
        # Optional debug output for coefficient matrix
        if ksm1 != 0 and (n == 1 or n == nn or n == ksm):
            iocode.write("\n\n")
            iocode.write(f"{'NODE':>9}{'AT(I,J,K)':>15}{'AS(I,J,K)':>15}{'AW(I,J,K)':>15}{'E(I,J,K)':>15}" +
                        f"{'AE(I,J,K)':>15}{'AN(I,J,K)':>15}{'AB(I,J,K)':>15}{'B(I,J,K)':>15}\n\n")
            
            for k in range(kk):
                for j in range(jj):
                    for i in range(ii):
                        iocode.write(f"{i+1:3d}{j+1:3d}{k+1:3d}{self.sim.at[i, j, k]:15.6e}{self.sim.as_[i, j, k]:15.6e}" +
                                   f"{self.sim.aw[i, j, k]:15.6e}{self.sim.e[i, j, k]:15.6e}{self.sim.ae[i, j, k]:15.6e}" +
                                   f"{self.sim.an[i, j, k]:15.6e}{self.sim.ab[i, j, k]:15.6e}{self.sim.b[i, j, k]:15.6e}\n")


def trikro(simulator, irockb: int, sso: float, ssw: float) -> float:
    """
    Calculate three-phase oil relative permeability using Stone's Model I
    
    Parameters:
    -----------
    simulator : BOASTSimulator
        Reference to main simulator
    irockb : int
        Rock region index (0-based)
    sso : float
        Oil saturation
    ssw : float
        Water saturation
        
    Returns:
    --------
    rkro : float
        Three-phase oil relative permeability
    """
    from block2 import Interpolation
    
    interp_obj = Interpolation(simulator)
    
    ssg = 1.0 - sso - ssw
    
    # Handle case with no free gas
    if ssg <= 0.0:
        kro1 = interp_obj.interp(irockb, simulator.sat, simulator.krot,
                                simulator.msat[irockb], sso)
        return kro1
    
    # Get endpoint relative permeability (KROW at SWR)
    krowr = interp_obj.interp(irockb, simulator.sat, simulator.krot,
                             simulator.msat[irockb], simulator.swr[irockb])
    
    # Interpolate two-phase relative permeabilities
    krow = interp_obj.interp(irockb, simulator.sat, simulator.krot,
                            simulator.msat[irockb], sso)
    krw = interp_obj.interp(irockb, simulator.sat, simulator.krwt,
                           simulator.msat[irockb], ssw)
    
    krog = interp_obj.interp(irockb, simulator.sat, simulator.krogt,
                            simulator.msat[irockb], sso)
    krg = interp_obj.interp(irockb, simulator.sat, simulator.krgt,
                           simulator.msat[irockb], ssg)
    
    # Stone's Model I formula
    if krowr > 0.0:
        rkro = ((krow + krw) * (krog + krg)) / krowr - (krw + krg)
    else:
        rkro = 0.0
    
    # Clip to valid range
    rkro = max(0.0, min(1.0, rkro))
    
    return rkro
