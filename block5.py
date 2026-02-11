"""
BLOCK5.FOR - Well Models and Rate Calculations

This module contains the WellRates class for calculating production and injection rates
for various well types with multiple constraints:
- Rate-controlled wells (oil, water, gas)
- Pressure-controlled wells
- GOR and WOR limits
- Maximum liquid withdrawal
- Minimum oil production
- Implicit pressure formulation

Converted from BOAST II (Release 1.2) Fortran code
"""

import numpy as np
from typing import TextIO
import math


class WellRates:
    """
    Well model calculations for production and injection rates
    
    Handles:
    - Oil producers (specified rate or BHP)
    - Water injectors (specified rate or BHP)
    - Gas injectors (specified rate or BHP)
    - Gas wells (with back-pressure equation)
    - GOR/WOR constraints with layer shutoff
    - Minimum oil rate and maximum liquid rate
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
    
    def qrate(self, ii: int, jj: int, kk: int, nvqn: int,
              gormax: float, wormax: float, eti: float) -> None:
        """
        Calculate well production and injection rates
        
        Parameters:
        -----------
        ii, jj, kk : int
            Grid dimensions
        nvqn : int
            Number of wells
        gormax : float
            Maximum producing GOR (SCF/STB)
        wormax : float
            Maximum producing WOR (STB/STB)
        eti : float
            Elapsed time (days)
        """
        from block2 import Interpolation
        
        # Initialize GOR/WOR limits
        for j in range(nvqn):
            self.sim.gort[j] = gormax
            self.sim.wort[j] = wormax
            self.sim.gorl[j] = 0.0
            self.sim.worl[j] = 0.0
            self.sim.ilimop[j] = 1
            
            # Get well location
            iq1 = self.sim.iqn1[j]
            iq2 = self.sim.iqn2[j]
            iq3 = self.sim.iqn3[j]
            
            # GOR/WOR vary with rock region if not specified
            if self.sim.gort[j] == 0.0:
                lay = iq3 + self.sim.layer[j] - 1
                irockr = self.sim.irock[iq1, iq2, iq3]
                irockr_0 = irockr - 1
                self.sim.gort[j] = self.sim.gorock[irockr_0]
                
                for k in range(iq3, lay + 1):
                    irockr = self.sim.irock[iq1, iq2, k]
                    irockr_0 = irockr - 1
                    if self.sim.gorock[irockr_0] > self.sim.gort[j]:
                        self.sim.gort[j] = self.sim.gorock[irockr_0]
            
            if self.sim.wort[j] == 0.0:
                lay = iq3 + self.sim.layer[j] - 1
                irockr = self.sim.irock[iq1, iq2, iq3]
                irockr_0 = irockr - 1
                self.sim.wort[j] = self.sim.worock[irockr_0]
                
                for k in range(iq3, lay + 1):
                    irockr = self.sim.irock[iq1, iq2, k]
                    irockr_0 = irockr - 1
                    if self.sim.worock[irockr_0] > self.sim.wort[j]:
                        self.sim.wort[j] = self.sim.worock[irockr_0]
        
        # Initialize rates
        self.sim.qo.fill(0.0)
        self.sim.qw.fill(0.0)
        self.sim.qg.fill(0.0)
        
        for m in range(nvqn):
            self.sim.qoc[m, :].fill(0.0)
            self.sim.qwc[m, :].fill(0.0)
            self.sim.qgc[m, :].fill(0.0)
        
        # Calculate mobilities for each well layer
        for j in range(nvqn):
            iq1 = self.sim.iqn1[j]
            iq2 = self.sim.iqn2[j]
            iq3 = self.sim.iqn3[j]
            
            lay = iq3 + self.sim.layer[j] - 1
            
            for k in range(iq3, lay + 1):
                self.sim.pwfc[j, k] = -1.0
                
                pp = self.sim.p[iq1, iq2, k]
                bpt = self.sim.pbot[iq1, iq2, k]
                ipvtr = self.sim.ipvt[iq1, iq2, k]
                ipvtr_0 = ipvtr - 1  # Convert to 0-based
                irockr = self.sim.irock[iq1, iq2, k]
                irockr_0 = irockr - 1  # Convert to 0-based
                
                # Interpolate fluid properties
                muo = self._intpvt(ipvtr_0, bpt, self.sim.vslope[ipvtr_0],
                                  self.sim.pot, self.sim.muot,
                                  self.sim.mpot[ipvtr_0], pp)
                
                muw = Interpolation.interp(
                    self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                    self.sim.muwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], pp)
                
                mug = Interpolation.interp(
                    self.sim.pgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]],
                    self.sim.mugt[ipvtr_0, :self.sim.mpgt[ipvtr_0]], pp)
                
                # Use OLD-TIME saturations for IMPES (explicit saturations)
                sso = self.sim.son[iq1, iq2, k]
                ssw = self.sim.swn[iq1, iq2, k]
                ssg = self.sim.sgn[iq1, iq2, k]
                
                # Interpolate relative permeabilities
                krw = Interpolation.interp(
                    self.sim.sat[irockr_0, :self.sim.msat[irockr_0]],
                    self.sim.krwt[irockr_0, :self.sim.msat[irockr_0]], ssw)
                
                if self.sim.ithree[irockr_0] == 0:
                    kro = Interpolation.interp(
                        self.sim.sat[irockr_0, :self.sim.msat[irockr_0]],
                        self.sim.krot[irockr_0, :self.sim.msat[irockr_0]], sso)
                else:
                    kro = self._trikro(irockr_0, sso, ssw)
                
                krg = Interpolation.interp(
                    self.sim.sat[irockr_0, :self.sim.msat[irockr_0]],
                    self.sim.krgt[irockr_0, :self.sim.msat[irockr_0]], ssg)
                
                # Calculate mobilities
                self.sim.gmw[j, k] = krw / muw
                self.sim.gmo[j, k] = kro / muo
                self.sim.gmg[j, k] = krg / mug
        
        # Calculate rates for rate-controlled wells
        self._calculate_rate_controlled_wells(nvqn)
        
        # Calculate rates for pressure-controlled wells
        self._calculate_pressure_controlled_wells(nvqn)
        
        # Apply rate constraints
        self._apply_rate_constraints(nvqn, eti)
        
        # Apply GOR and WOR constraints
        self._apply_gor_wor_constraints(nvqn, eti)
        
        # Calculate bottom-hole flowing pressures
        self._calculate_bhp(nvqn)
        
        # Sum rates by grid block (excluding implicit wells)
        self._sum_rates_by_block(nvqn)
    
    def _intpvt(self, ipvtr: int, bpt: float, slope: float,
                pot_table: np.ndarray, prop_table: np.ndarray,
                mpt: int, pp: float) -> float:
        """
        Interpolate PVT properties with bubble point handling
        
        Parameters:
        -----------
        ipvtr : int
            PVT region index
        bpt : float
            Bubble point pressure
        slope : float
            Slope above bubble point
        pot_table, prop_table : arrays
            Pressure and property tables
        mpt : int
            Number of table entries
        pp : float
            Current pressure
            
        Returns:
        --------
        prop : float
            Interpolated property value
        """
        from block2 import Interpolation
        
        if pp <= bpt:
            # Below bubble point - use table
            prop = Interpolation.interp(
                pot_table[ipvtr, :mpt],
                prop_table[ipvtr, :mpt], pp)
        else:
            # Above bubble point - use slope
            prop_at_pb = Interpolation.interp(
                pot_table[ipvtr, :mpt],
                prop_table[ipvtr, :mpt], bpt)
            prop = prop_at_pb + slope * (pp - bpt)
        
        return prop
    
    def _trikro(self, irockr: int, sso: float, ssw: float) -> float:
        """
        Three-phase oil relative permeability (Stone's method)
        
        Parameters:
        -----------
        irockr : int
            Rock region index
        sso, ssw : float
            Oil and water saturations
            
        Returns:
        --------
        kro : float
            Three-phase oil relative permeability
        """
        from block2 import Interpolation
        
        # Get endpoint saturations
        swr = self.sim.swr[irockr]
        
        # Normalized saturations
        sw_norm = (ssw - swr) / (1.0 - swr)
        
        # Get two-phase oil relative permeabilities
        kro_ow = Interpolation.interp(
            self.sim.sat[irockr, :self.sim.msat[irockr]],
            self.sim.krot[irockr, :self.sim.msat[irockr]], sso)
        
        kro_og = Interpolation.interp(
            self.sim.sat[irockr, :self.sim.msat[irockr]],
            self.sim.krogt[irockr, :self.sim.msat[irockr]], sso)
        
        # Stone's model (simplified)
        kro = kro_ow * kro_og * (1.0 - sw_norm)
        
        return max(kro, 0.0)
    
    def _calculate_rate_controlled_wells(self, nvqn: int) -> None:
        """Calculate rates for rate-controlled wells (KIP >= 0)"""
        from block2 import Interpolation
        
        for j in range(nvqn):
            if self.sim.kip[j] < 0:
                continue
            
            iq1 = self.sim.iqn1[j]
            iq2 = self.sim.iqn2[j]
            iq3 = self.sim.iqn3[j]
            
            lay = iq3 + self.sim.layer[j] - 1
            
            # Special case: oil injection (soluble oil)
            if self.sim.kip[j] == 1 and self.sim.qvo[j] < -0.001:
                self._calculate_oil_injection(j, iq1, iq2, iq3, j, lay)
                continue
            
            # Special case: oil production with specified oil rate (KIP=1, QVO > 0)
            if self.sim.kip[j] == 1 and self.sim.qvo[j] > 0.001:
                self._calculate_oil_production(j, iq1, iq2, iq3, j, lay)
                continue
            
            # Standard rate-controlled wells (water/gas injection)
            iterq = 0
            qdenom = 0.0
            
            # Two iterations: first to calculate qdenom, second to allocate rates
            for iterq in range(2):
                for k in range(iq3, lay + 1):
                    pp = self.sim.p[iq1, iq2, k]
                    bpt = self.sim.pbot[iq1, iq2, k]
                    ipvtr = self.sim.ipvt[iq1, iq2, k]
                    ipvtr_0 = ipvtr - 1
                    
                    bbo = self._intpvt(ipvtr_0, bpt, self.sim.bslope[ipvtr_0],
                                      self.sim.pot, self.sim.bot,
                                      self.sim.mpot[ipvtr_0], pp)
                    
                    bbw = Interpolation.interp(
                        self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                        self.sim.bwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], pp)
                    
                    bbg = Interpolation.interp(
                        self.sim.pgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]],
                        self.sim.bgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]], pp)
                    
                    if iterq == 0:
                        qdenom += self.sim.pid[j, k] * (
                            self.sim.gmo[j, k] + self.sim.gmw[j, k] + self.sim.gmg[j, k])
                    else:
                        if qdenom == 0.0:
                            continue
                        
                        # Oil injection (negative rate)
                        if self.sim.qvo[j] < -0.001:
                            self.sim.qoc[j, k] = (self.sim.qvo[j] * 5.615 *
                                self.sim.pid[j, k] * (self.sim.gmo[j, k] +
                                self.sim.gmw[j, k] + self.sim.gmg[j, k]) / qdenom)
                        
                        # Water injection
                        elif self.sim.kip[j] == 2:
                            self.sim.qwc[j, k] = (self.sim.qvw[j] * 5.615 *
                                self.sim.pid[j, k] * (self.sim.gmo[j, k] +
                                self.sim.gmw[j, k] + self.sim.gmg[j, k]) / qdenom)
                        
                        # Gas injection
                        else:
                            self.sim.qgc[j, k] = (self.sim.qvg[j] * 1000.0 *
                                self.sim.pid[j, k] * (self.sim.gmo[j, k] +
                                self.sim.gmw[j, k] + self.sim.gmg[j, k]) / qdenom)
    
    def _calculate_oil_production(self, j: int, iq1: int, iq2: int,
                                  iq3: int, ij: int, lay: int) -> None:
        """Calculate rates for oil production wells (KIP=1, QVO > 0)"""
        from block2 import Interpolation
        
        iterq = 0
        qdenom = 0.0
        
        # Two iterations: first to calculate qdenom, second to allocate rates
        for iterq in range(2):
            for k in range(iq3, lay + 1):
                pp = self.sim.p[iq1, iq2, k]
                bpt = self.sim.pbot[iq1, iq2, k]
                ipvtr = self.sim.ipvt[iq1, iq2, k]
                ipvtr_0 = ipvtr - 1
                
                bbo = self._intpvt(ipvtr_0, bpt, self.sim.bslope[ipvtr_0],
                                  self.sim.pot, self.sim.bot,
                                  self.sim.mpot[ipvtr_0], pp)
                bbw = Interpolation.interp(
                    self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                    self.sim.bwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], pp)
                bbg = Interpolation.interp(
                    self.sim.pgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]],
                    self.sim.bgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]], pp)
                rso = self._intpvt(ipvtr_0, bpt, self.sim.rslope[ipvtr_0],
                                  self.sim.pot, self.sim.rsot,
                                  self.sim.mpot[ipvtr_0], pp)
                rsw = Interpolation.interp(
                    self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                    self.sim.rswt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], pp)
                
                if iterq == 0:
                    # First iteration: calculate denominator
                    qdenom += self.sim.pid[j, k] * self.sim.gmo[j, k] / bbo
                    if self.sim.qvw[j] != 0.0:
                        qdenom += self.sim.pid[j, k] * self.sim.gmw[j, k] / bbw
                    if self.sim.qvg[j] != 0.0:
                        qdenom += self.sim.pid[j, k] * self.sim.gmg[j, k] / bbg
                else:
                    # Second iteration: calculate rates
                    if qdenom == 0.0 or self.sim.gmo[j, k] == 0.0:
                        continue
                    
                    # Oil rate (specified)
                    totor = self.sim.qvo[j]
                    self.sim.qoc[ij, k] = (totor * 5.615 * self.sim.pid[j, k] *
                                          self.sim.gmo[j, k] / (bbo * qdenom))
                    
                    # Water rate (proportional to mobility)
                    self.sim.qwc[ij, k] = (self.sim.qoc[ij, k] * self.sim.gmw[j, k] *
                                          bbo / (bbw * self.sim.gmo[j, k]))
                    
                    # Gas rate (solution gas + free gas)
                    self.sim.qgc[ij, k] = (self.sim.qoc[ij, k] *
                                          (self.sim.gmg[j, k] * bbo / (bbg * self.sim.gmo[j, k]) + rso) +
                                          rsw * self.sim.qwc[ij, k])
    
    def _calculate_oil_injection(self, j: int, iq1: int, iq2: int,
                                 iq3: int, ij: int, lay: int) -> None:
        """Calculate rates for oil injection wells with phase splitting"""
        from block2 import Interpolation
        
        iterq = 0
        qdenom = 0.0
        alphao = 0.0
        alphaw = 0.0
        alphag = 0.0
        bbosum = 0.0
        
        for iterq in range(2):
            for k in range(iq3, lay + 1):
                pp = self.sim.p[iq1, iq2, k]
                bpt = self.sim.pbot[iq1, iq2, k]
                ipvtr = self.sim.ipvt[iq1, iq2, k]
                ipvtr_0 = ipvtr - 1
                
                bbo = self._intpvt(ipvtr_0, bpt, self.sim.bslope[ipvtr_0],
                                  self.sim.pot, self.sim.bot,
                                  self.sim.mpot[ipvtr_0], pp)
                bbw = Interpolation.interp(
                    self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                    self.sim.bwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], pp)
                bbg = Interpolation.interp(
                    self.sim.pgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]],
                    self.sim.bgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]], pp)
                rso = self._intpvt(ipvtr_0, bpt, self.sim.rslope[ipvtr_0],
                                  self.sim.pot, self.sim.rsot,
                                  self.sim.mpot[ipvtr_0], pp)
                rsw = Interpolation.interp(
                    self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                    self.sim.rswt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], pp)
                
                if iterq == 0:
                    qdenom += self.sim.pid[j, k] * self.sim.gmo[j, k] / bbo
                    if self.sim.qvw[j] != 0.0:
                        qdenom += self.sim.pid[j, k] * self.sim.gmw[j, k] / bbw
                    if self.sim.qvg[j] != 0.0:
                        qdenom += self.sim.pid[j, k] * self.sim.gmg[j, k] / bbg
                    
                    gmt = self.sim.gmo[j, k] + self.sim.gmw[j, k] + self.sim.gmg[j, k]
                    alphao += self.sim.gmo[j, k] / gmt
                    alphaw += self.sim.gmw[j, k] / gmt
                    alphag += self.sim.gmg[j, k] / gmt
                    bbosum += bbo
                else:
                    # Convert total rate to oil rate
                    if self.sim.qvt[j] != 0.0:
                        bboavg = bbosum / self.sim.layer[j]
                        totor = (self.sim.qvt[j] / bboavg) * alphao / (alphao + alphaw + alphag)
                    else:
                        totor = self.sim.qvo[j]
                    
                    if qdenom == 0.0 or self.sim.gmo[j, k] == 0.0:
                        continue
                    
                    # Allocate rates by mobility
                    self.sim.qoc[ij, k] = (totor * 5.615 * self.sim.pid[j, k] *
                                          self.sim.gmo[j, k] / (bbo * qdenom))
                    self.sim.qwc[ij, k] = (self.sim.qoc[ij, k] * self.sim.gmw[j, k] *
                                          bbo / (bbw * self.sim.gmo[j, k]))
                    self.sim.qgc[ij, k] = (self.sim.qoc[ij, k] * (self.sim.gmg[j, k] *
                                          bbo / (bbg * self.sim.gmo[j, k]) + rso) +
                                          rsw * self.sim.qwc[ij, k])
    
    def _calculate_pressure_controlled_wells(self, nvqn: int) -> None:
        """Calculate rates for pressure-controlled wells (KIP < 0)"""
        from block2 import Interpolation
        
        for j in range(nvqn):
            kip = self.sim.kip[j]
            
            # Skip rate-controlled wells
            if kip >= 0:
                continue
            
            iq1 = self.sim.iqn1[j]
            iq2 = self.sim.iqn2[j]
            iq3 = self.sim.iqn3[j]
            
            lay = iq3 + self.sim.layer[j] - 1
            
            # Gas well with back-pressure equation
            if kip == -4:
                self._calculate_gas_well(j, iq1, iq2, iq3, j, lay)
                continue
            
            # Standard pressure-controlled wells
            for k in range(iq3, lay + 1):
                ppn = self.sim.pn[iq1, iq2, k]
                bpt = self.sim.pbot[iq1, iq2, k]
                ipvtr = self.sim.ipvt[iq1, iq2, k]
                ipvtr_0 = ipvtr - 1
                
                bbo = self._intpvt(ipvtr_0, bpt, self.sim.bslope[ipvtr_0],
                                  self.sim.pot, self.sim.bot,
                                  self.sim.mpot[ipvtr_0], ppn)
                bbw = Interpolation.interp(
                    self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                    self.sim.bwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], ppn)
                bbg = Interpolation.interp(
                    self.sim.pgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]],
                    self.sim.bgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]], ppn)
                rso = self._intpvt(ipvtr_0, bpt, self.sim.rslope[ipvtr_0],
                                  self.sim.pot, self.sim.rsot,
                                  self.sim.mpot[ipvtr_0], ppn)
                rsw = Interpolation.interp(
                    self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                    self.sim.rswt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], ppn)
                
                pwf_k = self.sim.pwf[j, k]
                
                # Oil producer (KIP = -1)
                if kip == -1:
                    self.sim.qoc[ij, k] = (self.sim.pid[j, k] * 5.615 *
                        self.sim.gmo[j, k] * (ppn - pwf_k) / bbo)
                    if ppn <= pwf_k:
                        self.sim.qoc[ij, k] = 0.0
                    
                    self.sim.qwc[ij, k] = (self.sim.pid[j, k] * 5.615 *
                        self.sim.gmw[j, k] * (ppn - pwf_k) / bbw)
                    if ppn <= pwf_k:
                        self.sim.qwc[ij, k] = 0.0
                    
                    if self.sim.qoc[ij, k] > 0.0:
                        qg1 = (self.sim.qoc[ij, k] * (self.sim.gmg[j, k] *
                              bbo / (bbg * self.sim.gmo[j, k]) + rso))
                    else:
                        qg1 = 0.0
                    
                    self.sim.qgc[ij, k] = qg1 + rsw * self.sim.qwc[ij, k]
                
                # Water injector (KIP = -2)
                elif kip == -2:
                    self.sim.qwc[ij, k] = (self.sim.pid[j, k] * 5.615 *
                        (self.sim.gmo[j, k] + self.sim.gmw[j, k] + self.sim.gmg[j, k]) *
                        (ppn - pwf_k) / bbw)
                    if ppn >= pwf_k:
                        self.sim.qwc[ij, k] = 0.0
                
                # Gas injector (KIP = -3)
                elif kip == -3:
                    self.sim.qgc[ij, k] = (self.sim.pid[j, k] * 5.615 *
                        (self.sim.gmo[j, k] + self.sim.gmw[j, k] + self.sim.gmg[j, k]) *
                        (ppn - pwf_k) / bbg)
                    if ppn >= pwf_k:
                        self.sim.qgc[ij, k] = 0.0
    
    def _calculate_gas_well(self, j: int, iq1: int, iq2: int,
                           iq3: int, ij: int, lay: int) -> None:
        """Calculate rates for gas well using back-pressure equation"""
        from block2 import Interpolation
        
        iterq = 0
        qdenom = 0.0
        
        for iterq in range(2):
            for k in range(iq3, lay + 1):
                pp = self.sim.p[iq1, iq2, k]
                bpt = self.sim.pbot[iq1, iq2, k]
                ipvtr = self.sim.ipvt[iq1, iq2, k]
                ipvtr_0 = ipvtr - 1
                
                bbo = self._intpvt(ipvtr_0, bpt, self.sim.bslope[ipvtr_0],
                                  self.sim.pot, self.sim.bot,
                                  self.sim.mpot[ipvtr_0], pp)
                bbw = Interpolation.interp(
                    self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                    self.sim.bwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], pp)
                bbg = Interpolation.interp(
                    self.sim.pgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]],
                    self.sim.bgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]], pp)
                rso = self._intpvt(ipvtr_0, bpt, self.sim.rslope[ipvtr_0],
                                  self.sim.pot, self.sim.rsot,
                                  self.sim.mpot[ipvtr_0], pp)
                rsw = Interpolation.interp(
                    self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                    self.sim.rswt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], pp)
                
                if iterq == 0:
                    qdenom += self.sim.pid[j, k] * self.sim.gmg[j, k] / bbg
                else:
                    pwf_k = self.sim.pwf[j, k]
                    
                    # Oil and water production
                    self.sim.qoc[ij, k] = (self.sim.pid[j, k] * 5.615 *
                        self.sim.gmo[j, k] * (pp - pwf_k) / bbo)
                    if pp <= pwf_k:
                        self.sim.qoc[ij, k] = 0.0
                    
                    self.sim.qwc[ij, k] = (self.sim.pid[j, k] * 5.615 *
                        self.sim.gmw[j, k] * (pp - pwf_k) / bbw)
                    if pp <= pwf_k:
                        self.sim.qwc[ij, k] = 0.0
                    
                    # Gas production using back-pressure equation
                    psir = Interpolation.interp(
                        self.sim.pgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]],
                        self.sim.psit[ipvtr_0, :self.sim.mpgt[ipvtr_0]], pp)
                    psiwf = Interpolation.interp(
                        self.sim.pgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]],
                        self.sim.psit[ipvtr_0, :self.sim.mpgt[ipvtr_0]], pwf_k)
                    
                    qlit = 0.0
                    qlitk = 0.0
                    
                    if psir >= psiwf and qdenom > 0.0:
                        # Back-pressure equation: Q = C*(Pr^2 - Pwf^2)^n
                        # Linearized form: Q = (-A + sqrt(A^2 + 4B(Psi_r - Psi_wf)))/(2B)
                        alit = self.sim.alit[ij]
                        blit = self.sim.blit[ij]
                        
                        qlit = 1.0e6 * (-alit + math.sqrt(alit**2 +
                               4.0 * blit * (psir - psiwf))) / (2.0 * blit)
                        
                        qlitk = (qlit * self.sim.pid[j, k] * self.sim.gmg[j, k] /
                                (qdenom * bbg))
                    
                    self.sim.qgc[ij, k] = (qlitk + rso * self.sim.qoc[ij, k] +
                                          rsw * self.sim.qwc[ij, k])
    
    def _apply_rate_constraints(self, nvqn: int, eti: float) -> None:
        """Apply minimum oil rate and maximum liquid rate constraints"""
        
        for j in range(nvqn):
            if self.sim.kip[j] != -1:
                continue
            
            iq1 = self.sim.iqn1[j]
            iq2 = self.sim.iqn2[j]
            iq3 = self.sim.iqn3[j]
            
            lay = iq3 + self.sim.layer[j] - 1
            
            # Sum rates
            qot = 0.0
            qwt = 0.0
            pidsum = 0.0
            
            for k in range(iq3, lay + 1):
                qot += self.sim.qoc[j, k]
                qwt += self.sim.qwc[j, k]
                pidsum += self.sim.pid[j, k]
            
            # Skip if well already shut-in
            if pidsum <= 0.0:
                continue
            
            # Check minimum oil rate (convert STB to SCF: 5.615)
            if qot < self.sim.qvo[j] * 5.615:
                # Shut-in well
                for k in range(iq3, lay + 1):
                    self.sim.qoc[j, k] = 0.0
                    self.sim.qwc[j, k] = 0.0
                    self.sim.qgc[j, k] = 0.0
                    self.sim.pid[j, k] = 0.0
                
                self.sim.iocode.write(
                    f"\n{'':>9}{'-'*110}\n"
                    f"{'':>9}MINIMUM OIL RATE NOT ACHIEVED BY WELL #{j+1:3d}, " +
                    f"AREAL LOCATION{iq1+1:3d},{iq2+1:3d} AFTER{eti:10.2f} DAYS OF ELAPSED TIME.\n"
                    f"{'':>9}{'-'*110}\n")
                continue
            
            # Check maximum oil rate
            fac1 = 1.0
            if self.sim.qvw[j] > 0.0 and qot > 5.615 * self.sim.qvw[j]:
                fac1 = 5.615 * self.sim.qvw[j] / qot
            
            # Check maximum liquid rate
            fac2 = 1.0
            if self.sim.qvt[j] > 0.0:
                qliqt = (qot + qwt) * fac1
                if qliqt > 5.615 * self.sim.qvt[j]:
                    fac2 = 5.615 * self.sim.qvt[j] / qliqt
            
            fac = fac1 * fac2
            
            if fac < 1.0:
                # Scale back rates
                from block2 import Interpolation
                
                for k in range(iq3, lay + 1):
                    self.sim.qoc[j, k] *= fac
                    self.sim.qwc[j, k] *= fac
                    
                    # Recalculate gas rate
                    ppn = self.sim.pn[iq1, iq2, k]
                    bpt = self.sim.pbot[iq1, iq2, k]
                    ipvtr = self.sim.ipvt[iq1, iq2, k]
                    ipvtr_0 = ipvtr - 1
                    
                    bbo = self._intpvt(ipvtr_0, bpt, self.sim.bslope[ipvtr_0],
                                      self.sim.pot, self.sim.bot,
                                      self.sim.mpot[ipvtr_0], ppn)
                    bbg = Interpolation.interp(
                        self.sim.pgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]],
                        self.sim.bgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]], ppn)
                    rso = self._intpvt(ipvtr_0, bpt, self.sim.rslope[ipvtr_0],
                                      self.sim.pot, self.sim.rsot,
                                      self.sim.mpot[ipvtr_0], ppn)
                    rsw = Interpolation.interp(
                        self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                        self.sim.rswt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], ppn)
                    
                    if self.sim.qoc[j, k] > 0.0:
                        qg1 = (self.sim.qoc[j, k] * (self.sim.gmg[j, k] *
                              bbo / (bbg * self.sim.gmo[j, k]) + rso))
                    else:
                        qg1 = 0.0
                    
                    self.sim.qgc[j, k] = qg1 + rsw * self.sim.qwc[j, k]
        
        # Apply rate constraints on injection wells
        for j in range(nvqn):
            iq1 = self.sim.iqn1[j]
            iq2 = self.sim.iqn2[j]
            iq3 = self.sim.iqn3[j]
            lay = iq3 + self.sim.layer[j] - 1
            
            facw = 1.0
            facg = 1.0
            
            # Water injector constraint
            if self.sim.kip[j] == -2:
                qwi = sum(self.sim.qwc[j, k] for k in range(iq3, lay + 1))
                if self.sim.qvw[j] < 0.0 and abs(qwi) > abs(self.sim.qvw[j]) * 5.615:
                    facw = self.sim.qvw[j] * 5.615 / qwi
            
            # Gas injector constraint
            if self.sim.kip[j] == -3:
                qgi = sum(self.sim.qgc[j, k] for k in range(iq3, lay + 1))
                if self.sim.qvg[j] < 0.0 and abs(qgi) > abs(self.sim.qvg[j]) * 1000.0:
                    facg = self.sim.qvg[j] * 1000.0 / qgi
            
            if facw < 1.0 or facg < 1.0:
                for k in range(iq3, lay + 1):
                    self.sim.qwc[j, k] *= facw
                    self.sim.qgc[j, k] *= facg
    
    def _apply_gor_wor_constraints(self, nvqn: int, eti: float) -> None:
        """Apply GOR and WOR constraints with layer shutoff"""
        
        for j in range(nvqn):
            iq1 = self.sim.iqn1[j]
            iq2 = self.sim.iqn2[j]
            iq3 = self.sim.iqn3[j]
            
            if self.sim.ilimop[j] == 0 or self.sim.kip[j] < -10:
                continue
            
            lay = iq3 + self.sim.layer[j] - 1
            
            # Iteratively shut off high GOR/WOR layers
            while True:
                # Sum rates
                qot = 0.0
                qwt = 0.0
                qgt = 0.0
                
                for k in range(iq3, lay + 1):
                    qot += self.sim.qoc[j, k]
                    qwt += self.sim.qwc[j, k]
                    qgt += self.sim.qgc[j, k]
                
                # Calculate ratios
                gor = 0.0
                wor = 0.0
                if qot > 0.0:
                    gor = qgt * 5.615 / qot
                    wor = qwt / qot
                
                # Check GOR constraint
                if gor > self.sim.gort[j]:
                    # Calculate GOR by layer
                    for k in range(iq3, lay + 1):
                        if self.sim.qoc[j, k] == 0.0:
                            self.sim.pid[j, k] = 0.0
                            self.sim.qwc[j, k] = 0.0
                            self.sim.qgc[j, k] = 0.0
                            self.sim.gorl[k] = 0.0
                        else:
                            self.sim.gorl[k] = self.sim.qgc[j, k] * 5.615 / self.sim.qoc[j, k]
                    
                    # Find layer with max GOR
                    gorsi = self.sim.gorl[iq3]
                    kmax = iq3
                    for k in range(iq3, lay + 1):
                        if self.sim.gorl[k] > gorsi:
                            gorsi = self.sim.gorl[k]
                            kmax = k
                    
                    # Shut-in layer
                    self.sim.pid[j, kmax] = 0.0
                    self.sim.qoc[j, kmax] = 0.0
                    self.sim.qwc[j, kmax] = 0.0
                    self.sim.qgc[j, kmax] = 0.0
                    
                    self.sim.iocode.write(
                        f"\n{'':>9}{'-'*110}\n"
                        f"{'':>9}GOR LIMIT EXCEEDED BY LAYER K ={kmax+1:3d}, WELL #{j+1:3d}, " +
                        f"AREAL LOCATION{iq1+1:3d},{iq2+1:3d} AFTER{eti:10.2f} DAYS OF ELAPSED TIME.\n"
                        f"{'':>9}{'-'*110}\n")
                    continue  # Repeat check
                
                # Check WOR constraint
                if wor > self.sim.wort[j]:
                    # Calculate WOR by layer
                    for k in range(iq3, lay + 1):
                        if self.sim.qoc[j, k] == 0.0:
                            self.sim.pid[j, k] = 0.0
                            self.sim.qwc[j, k] = 0.0
                            self.sim.qgc[j, k] = 0.0
                            self.sim.worl[k] = 0.0
                        else:
                            self.sim.worl[k] = self.sim.qwc[j, k] / self.sim.qoc[j, k]
                    
                    # Find layer with max WOR
                    worsi = self.sim.worl[lay]
                    kmax = lay
                    for k in range(iq3, lay + 1):
                        if self.sim.worl[k] >= worsi:
                            worsi = self.sim.worl[k]
                            kmax = k
                    
                    # Shut-in layer
                    self.sim.pid[j, kmax] = 0.0
                    self.sim.qoc[j, kmax] = 0.0
                    self.sim.qwc[j, kmax] = 0.0
                    self.sim.qgc[j, kmax] = 0.0
                    
                    self.sim.iocode.write(
                        f"\n{'':>9}{'-'*110}\n"
                        f"{'':>9}WOR LIMIT EXCEEDED BY LAYER K ={kmax+1:3d}, WELL #{j+1:3d}, " +
                        f"AREAL LOCATION{iq1+1:3d},{iq2+1:3d} AFTER{eti:10.2f} DAYS OF ELAPSED TIME.\n"
                        f"{'':>9}{'-'*110}\n")
                    continue  # Repeat check
                
                break  # No violations, exit loop
    
    def _calculate_bhp(self, nvqn: int) -> None:
        """Calculate bottom-hole flowing pressure for each layer"""
        from block2 import Interpolation
        
        for j in range(nvqn):
            iq1 = self.sim.iqn1[j]
            iq2 = self.sim.iqn2[j]
            iq3 = self.sim.iqn3[j]
            
            if self.sim.kip[j] < -10:
                continue
            
            lay = iq3 + self.sim.layer[j] - 1
            
            for k in range(iq3, lay + 1):
                self.sim.pwfc[j, k] = 0.0
                
                if self.sim.pid[j, k] <= 0.0001:
                    continue
                
                pp = self.sim.p[iq1, iq2, k]
                if pp <= 0.0:
                    continue
                
                bpt = self.sim.pbot[iq1, iq2, k]
                ipvtr = self.sim.ipvt[iq1, iq2, k]
                ipvtr_0 = ipvtr - 1
                
                bbo = self._intpvt(ipvtr_0, bpt, self.sim.bslope[ipvtr_0],
                                  self.sim.pot, self.sim.bot,
                                  self.sim.mpot[ipvtr_0], pp)
                bbw = Interpolation.interp(
                    self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                    self.sim.bwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], pp)
                bbg = Interpolation.interp(
                    self.sim.pgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]],
                    self.sim.bgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]], pp)
                
                rso = self._intpvt(ipvtr_0, bpt, self.sim.rslope[ipvtr_0],
                                  self.sim.pot, self.sim.rsot,
                                  self.sim.mpot[ipvtr_0], pp)
                rsw = Interpolation.interp(
                    self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                    self.sim.rswt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], pp)
                
                fac = self.sim.pid[j, k] * 5.615
                gmtb = (self.sim.gmo[j, k] / bbo + self.sim.gmw[j, k] / bbw +
                       self.sim.gmg[j, k] / bbg)
                soln = rso * self.sim.qoc[j, k] + rsw * self.sim.qwc[j, k]
                qt = self.sim.qoc[j, k] + self.sim.qwc[j, k] + self.sim.qgc[j, k]
                
                self.sim.pwfc[j, k] = pp - (qt - soln) / (fac * gmtb)
    
    def _sum_rates_by_block(self, nvqn: int) -> None:
        """Sum well rates by grid block (excluding implicit wells)"""
        
        for j in range(nvqn):
            if self.sim.kip[j] < -10:
                continue
            
            iq1 = self.sim.iqn1[j]
            iq2 = self.sim.iqn2[j]
            iq3 = self.sim.iqn3[j]
            
            lay = iq3 + self.sim.layer[j] - 1
            
            for k in range(iq3, lay + 1):
                self.sim.qo[iq1, iq2, k] += self.sim.qoc[j, k]
                self.sim.qw[iq1, iq2, k] += self.sim.qwc[j, k]
                self.sim.qg[iq1, iq2, k] += self.sim.qgc[j, k]
    
    def pratei(self, nvqn: int) -> None:
        """
        Add implicit pressure terms to coefficient matrix for implicit wells
        
        Parameters:
        -----------
        nvqn : int
            Number of wells
        """
        from block2 import Interpolation
        
        for j in range(nvqn):
            if self.sim.kip[j] >= -10:
                continue
            
            iq1 = self.sim.iqn1[j]
            iq2 = self.sim.iqn2[j]
            iq3 = self.sim.iqn3[j]
            lay = iq3 + self.sim.layer[j] - 1
            
            for k in range(iq3, lay + 1):
                p56 = self.sim.pid[j, k] * 5.615
                ppn = self.sim.pn[iq1, iq2, k]
                bpt = self.sim.pbot[iq1, iq2, k]
                ipvtr = self.sim.ipvt[iq1, iq2, k]
                ipvtr_0 = ipvtr - 1
                
                bbo = self._intpvt(ipvtr_0, bpt, self.sim.bslope[ipvtr_0],
                                  self.sim.pot, self.sim.bot,
                                  self.sim.mpot[ipvtr_0], ppn)
                rso = self._intpvt(ipvtr_0, bpt, self.sim.rslope[ipvtr_0],
                                  self.sim.pot, self.sim.rsot,
                                  self.sim.mpot[ipvtr_0], ppn)
                bbw = Interpolation.interp(
                    self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                    self.sim.bwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], ppn)
                bbg = Interpolation.interp(
                    self.sim.pgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]],
                    self.sim.bgt[ipvtr_0, :self.sim.mpgt[ipvtr_0]], ppn)
                rsw = Interpolation.interp(
                    self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                    self.sim.rswt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], ppn)
                
                cpio = self.sim.gmo[j, k] * p56 * (bbo - bbg * rso) / bbo
                cpiw = self.sim.gmw[j, k] * p56 * (bbw - bbg * rsw) / bbw
                cpig = self.sim.gmg[j, k] * p56
                cpi = cpio + cpiw + cpig
                
                self.sim.b[iq1, iq2, k] -= cpi * self.sim.pwf[j, k]
                self.sim.e[iq1, iq2, k] -= cpi
    
    def prateo(self, nvqn: int) -> None:
        """
        Calculate explicit rates for implicit pressure wells
        
        Parameters:
        -----------
        nvqn : int
            Number of wells
        """
        from block2 import Interpolation
        
        for j in range(nvqn):
            if self.sim.kip[j] >= -10:
                continue
            
            iq1 = self.sim.iqn1[j]
            iq2 = self.sim.iqn2[j]
            iq3 = self.sim.iqn3[j]
            
            lay = iq3 + self.sim.layer[j] - 1
            
            for k in range(iq3, lay + 1):
                pp = self.sim.p[iq1, iq2, k]
                ppn = self.sim.pn[iq1, iq2, k]
                bpt = self.sim.pbot[iq1, iq2, k]
                ipvtr = self.sim.ipvt[iq1, iq2, k]
                ipvtr_0 = ipvtr - 1
                
                rson = self._intpvt(ipvtr_0, bpt, self.sim.rslope[ipvtr_0],
                                   self.sim.pot, self.sim.rsot,
                                   self.sim.mpot[ipvtr_0], ppn)
                rso = self._intpvt(ipvtr_0, bpt, self.sim.rslope[ipvtr_0],
                                  self.sim.pot, self.sim.rsot,
                                  self.sim.mpot[ipvtr_0], pp)
                rswn = Interpolation.interp(
                    self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                    self.sim.rswt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], ppn)
                rsw = Interpolation.interp(
                    self.sim.pwt[ipvtr_0, :self.sim.mpwt[ipvtr_0]],
                    self.sim.rswt[ipvtr_0, :self.sim.mpwt[ipvtr_0]], pp)
                
                rsoav = 0.5 * (rso + rson)
                rswav = 0.5 * (rsw + rswn)
                
                factor = self.sim.pid[j, k] * 5.615 * (pp - self.sim.pwf[j, k])
                
                # Water producer/injector (KIP = -12 or -13)
                if self.sim.kip[j] == -13:
                    self.sim.qgc[j, k] = ((self.sim.gmo[j, k] + self.sim.gmw[j, k] +
                                           self.sim.gmg[j, k]) /
                                          self.sim.bg[iq1, iq2, k] * factor)
                else:
                    self.sim.qwc[j, k] = self.sim.gmw[j, k] / self.sim.bw[iq1, iq2, k] * factor
                    
                    if self.sim.kip[j] == -12:
                        self.sim.qwc[j, k] = ((self.sim.gmo[j, k] + self.sim.gmw[j, k] +
                                               self.sim.gmg[j, k]) /
                                              self.sim.bw[iq1, iq2, k] * factor)
                    else:
                        self.sim.qoc[j, k] = self.sim.gmo[j, k] / self.sim.bo[iq1, iq2, k] * factor
                        self.sim.qgc[j, k] = (self.sim.gmg[j, k] / self.sim.bg[iq1, iq2, k] * factor +
                                              rsoav * self.sim.qoc[j, k] +
                                              rswav * self.sim.qwc[j, k])
        
        # Sum rates by grid block including implicit wells
        for j in range(nvqn):
            if self.sim.kip[j] >= -10:
                continue
            
            iq1 = self.sim.iqn1[j]
            iq2 = self.sim.iqn2[j]
            iq3 = self.sim.iqn3[j]
            
            lay = iq3 + self.sim.layer[j] - 1
            
            for k in range(iq3, lay + 1):
                self.sim.qo[iq1, iq2, k] += self.sim.qoc[j, k]
                self.sim.qw[iq1, iq2, k] += self.sim.qwc[j, k]
                self.sim.qg[iq1, iq2, k] += self.sim.qgc[j, k]
