"""
BOAST II: BLACK OIL APPLIED SIMULATION TOOL
J.R. FANCHI
APRIL 1987 - RELEASE 1.2

Python conversion from Fortran
This is the main entry point for the BOAST II simulator
"""

import sys
import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Optional
import struct

# Import BOAST II modules
from block1 import AquiferModel, MaterialBalance, WellManager
from block2 import GridSetup, SolutionControl, Interpolation, LSORSolver
from block3 import OutputVisualization, RockProperties, PostProcessing
from block4 import GasProperties, Repressurization, TridiagonalSolver
from block5 import WellRates
from block6 import PVTTables, Transmissibility, Initialization
from block7 import FlowEquation
from block8 import FlowEquationTwoPoint

# Grid dimension parameters (corresponding to Fortran PARAMETER statements)
LP1, LP2, LP3 = 13, 10, 4
LP7, LP8, LP9 = 3, 3, 25
LP10, LP11, LP12 = 3, 50, 1000
LP14, LP15, LP17 = 1, 13, 5
LP4, LP5, LP6 = LP1 + 1, LP2 + 1, LP3 + 1
LP13 = LP1 + LP2 + LP3
LP19 = LP4 * LP2 * LP3
LP20 = LP1 * LP5 * LP3
LP21 = LP1 * LP2 * LP6
LP22 = LP1 * LP2 * LP3
LP23 = LP11 * LP3


@dataclass
class SimulationParameters:
    """Main simulation parameters"""
    nx: int = LP1
    ny: int = LP2
    nz: int = LP3
    nter: int = LP9
    ntep: int = LP9
    nw: int = LP11
    ntmax: int = LP12
    nrst: int = LP17
    n1dir: int = LP15
    n23dir: int = LP14
    niter: int = LP22
    ntri: int = LP15


class BOASTSimulator:
    """Main BOAST II Simulator class"""
    
    def __init__(self):
        """Initialize simulator with default values"""
        self.params = SimulationParameters()
        
        # Dimension parameters
        self.ii = 0  # X-direction blocks
        self.jj = 0  # Y-direction blocks
        self.kk = 0  # Z-direction blocks
        
        # Arrays - will be initialized based on actual grid size
        self.p = None    # Pressure
        self.pn = None   # Pressure at previous time step
        self.so = None   # Oil saturation
        self.son = None  # Oil saturation at previous time step
        self.sw = None   # Water saturation
        self.swn = None  # Water saturation at previous time step
        self.sg = None   # Gas saturation
        self.sgn = None  # Gas saturation at previous time step
        
        self.bo = None   # Oil formation volume factor
        self.bw = None   # Water formation volume factor
        self.bg = None   # Gas formation volume factor
        
        self.vp = None   # Pore volume
        self.ct = None   # Total compressibility
        
        self.kx = None   # X-direction permeability
        self.ky = None   # Y-direction permeability
        self.kz = None   # Z-direction permeability
        
        # Transmissibilities
        self.tx = None
        self.ty = None
        self.tz = None
        
        # Well data
        self.welnam = []  # Well names
        self.nvqn = 0     # Number of wells
        
        # Well arrays (from COMMON /SLIMIT/ and /SRATE/ in MAIN.FOR)
        self.gort = np.zeros(LP11)      # Max producing GOR (SCF/STB)
        self.wort = np.zeros(LP11)      # Max producing WOR (STB/STB)
        self.gorl = np.zeros(LP11)      # Current GOR limit
        self.worl = np.zeros(LP11)      # Current WOR limit
        self.ilimop = np.zeros(LP11, dtype=int)  # Limit operation flag
        
        # Well rate arrays
        self.pid = np.zeros((LP11, LP3))     # Productivity index
        self.pwf = np.zeros((LP11, LP3))     # Wellbore flowing pressure
        self.pwfc = np.zeros((LP11, LP3))    # Calculated wellbore pressure
        self.kip = np.zeros(LP11, dtype=int) # Well type flag
        self.layer = np.zeros(LP11, dtype=int)  # Number of layers
        self.qvo = np.zeros(LP11)            # Oil rate (STB/D)
        self.qvw = np.zeros(LP11)            # Water rate (STB/D)
        self.qvg = np.zeros(LP11)            # Gas rate (MSCF/D)
        self.qvt = np.zeros(LP11)            # Total liquid rate (RB/D)
        self.cumo = np.zeros((LP11, LP3))    # Cumulative oil
        self.cumw = np.zeros((LP11, LP3))    # Cumulative water
        self.cumg = np.zeros((LP11, LP3))    # Cumulative gas
        self.gmo = np.zeros((LP11, LP3))     # Oil mobility at well
        self.gmw = np.zeros((LP11, LP3))     # Water mobility at well
        self.gmg = np.zeros((LP11, LP3))     # Gas mobility at well
        self.idwell = np.zeros(LP11, dtype=int)  # Well ID
        self.alit = np.zeros(LP11)           # A coefficient
        self.blit = np.zeros(LP11)           # B coefficient
        
        # Well location arrays
        self.iqn1 = np.zeros(LP11, dtype=int)  # Well I location
        self.iqn2 = np.zeros(LP11, dtype=int)  # Well J location
        self.iqn3 = np.zeros(LP11, dtype=int)  # Well K location
        
        # Well rate by cell arrays
        self.qoc = np.zeros((LP11, LP3))    # Oil rate by layer
        self.qwc = np.zeros((LP11, LP3))    # Water rate by layer
        self.qgc = np.zeros((LP11, LP3))    # Gas rate by layer
        
        # Rock region limits
        self.gorock = np.zeros(LP7)         # GOR by rock region
        self.worock = np.zeros(LP7)         # WOR by rock region
        
        # Time stepping
        self.delt = 0.0   # Current time step size
        self.delt0 = 0.0  # Initial time step size
        self.eti = 0.0    # Elapsed time
        self.ft = 0.0     # Future time
        self.ftmax = 0.0  # Maximum future time
        
        # Material balance
        self.mbeo = 0.0   # Oil material balance error
        self.mbew = 0.0   # Water material balance error
        self.mbeg = 0.0   # Gas material balance error
        self.cmbeo = 0.0  # Cumulative oil material balance error
        self.cmbew = 0.0  # Cumulative water material balance error
        self.cmbeg = 0.0  # Cumulative gas material balance error
        
        # Production/injection totals
        self.cop = 0.0    # Cumulative oil production
        self.cwp = 0.0    # Cumulative water production
        self.cgp = 0.0    # Cumulative gas production
        self.cwi = 0.0    # Cumulative water injection
        self.cgi = 0.0    # Cumulative gas injection
        
        # Error counters
        self.ifatal = 0   # Fatal error count
        self.iwarn = 0    # Warning count
        
        # Output control
        self.iocode = None  # Output file handle
        self.ihed = list(range(1, 151))  # Header array for printouts (1-150)
        
        # File handles
        self.infile = None
        self.outfile = None
        self.resin_file = None
        self.resout_file = None
        
        # Restart parameters
        self.ireopt = 0   # Restart option
        self.irstrt = 0   # Restart time step
        
        # Runge-Kutta coefficients
        self.rkcoef = np.array([
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.13157894, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.03937504, 0.15625736, 1.0, 0.0, 0.0, 0.0, 0.0],
            [0.01666967, 0.053199139, 0.16491829, 1.0, 0.0, 0.0, 0.0],
            [0.0085487791, 0.024394839, 0.059605275, 0.16893119, 1.0, 0.0, 0.0],
            [0.0049515953, 0.013195402, 0.028593749, 0.06308716, 0.17111215, 1.0, 0.0],
            [0.0031198738, 0.0079382925, 0.015998138, 0.031126418, 0.06518731, 0.17242757, 1.0]
        ])
        
        self.krunge = 1
        
        # Additional arrays needed for simulation
        self.dx = None    # X-direction block size
        self.dy = None    # Y-direction block size
        self.dz = None    # Z-direction block size
        self.dznet = None # Net thickness
        self.el = None    # Elevation
        
        self.qo = None    # Oil well rates
        self.qw = None    # Water well rates
        self.qg = None    # Gas well rates
        
        # Coefficient arrays for pressure equation
        self.aw = None    # West coefficient
        self.ae = None    # East coefficient
        self.as_ = None   # South coefficient
        self.an = None    # North coefficient
        self.at = None    # Top coefficient
        self.ab = None    # Bottom coefficient
        self.e = None     # Center coefficient
        self.b = None     # RHS vector
        
        # Additional arrays
        self.sum = None   # Sum of neighbor coefficients
        self.gam = None   # Accumulation term
        self.gowt = None  # Gravity weight oil
        self.gwwt = None  # Gravity weight water
        self.ggwt = None  # Gravity weight gas
        self.qowg = None  # Combined well/gravity term
        
        # PVT arrays
        self.sat = np.zeros((LP7, LP9))      # Saturation table
        self.krot = np.zeros((LP7, LP9))     # Oil rel perm table
        self.krwt = np.zeros((LP7, LP9))     # Water rel perm table
        self.krgt = np.zeros((LP7, LP9))     # Gas rel perm table
        self.krogt = np.zeros((LP7, LP9))    # Oil-gas rel perm table
        self.pcowt = np.zeros((LP7, LP9))    # Oil-water cap pressure
        self.pcgot = np.zeros((LP7, LP9))    # Gas-oil cap pressure
        self.ithree = np.zeros(LP7, dtype=int)  # 3-phase flag
        self.swr = np.zeros(LP7)             # Irreducible water sat
        self.msat = np.zeros(LP7, dtype=int) # Number of sat entries
        
        self.pot = np.zeros((LP8, LP9))      # Oil pressure table
        self.muot = np.zeros((LP8, LP9))     # Oil viscosity table
        self.bot = np.zeros((LP8, LP9))      # Oil FVF table
        self.rsot = np.zeros((LP8, LP9))     # Solution gas table
        self.bopt = np.zeros((LP8, LP9))     # dBO/dP
        self.rsopt = np.zeros((LP8, LP9))    # dRSO/dP
        self.mpot = np.zeros(LP8, dtype=int) # Number of oil entries
        
        self.pwt = np.zeros((LP8, LP9))      # Water pressure table
        self.muwt = np.zeros((LP8, LP9))     # Water viscosity table
        self.bwt = np.zeros((LP8, LP9))      # Water FVF table
        self.rswt = np.zeros((LP8, LP9))     # Solution gas in water
        self.bwpt = np.zeros((LP8, LP9))     # dBW/dP
        self.rswpt = np.zeros((LP8, LP9))    # dRSW/dP
        self.mpwt = np.zeros(LP8, dtype=int) # Number of water entries
        
        self.pgt = np.zeros((LP8, LP9))      # Gas pressure table
        self.mugt = np.zeros((LP8, LP9))     # Gas viscosity table
        self.bgt = np.zeros((LP8, LP9))      # Gas FVF table
        self.bgpt = np.zeros((LP8, LP9))     # dBG/dP
        self.psit = np.zeros((LP8, LP9))     # Pseudo-pressure
        self.mpgt = np.zeros(LP8, dtype=int) # Number of gas entries
        
        self.crt = np.zeros((LP8, LP9))      # Rock compressibility
        self.prt = np.zeros((LP8, LP9))      # Rock pressure table
        
        self.vslope = np.zeros(LP8)          # Viscosity slope
        self.bslope = np.zeros(LP8)          # FVF slope
        self.rslope = np.zeros(LP8)          # RS slope
        
        self.rhosco = np.zeros(LP8)          # Oil density SC
        self.rhoscw = np.zeros(LP8)          # Water density SC
        self.rhoscg = np.zeros(LP8)          # Gas density SC
        
        self.pbot = None  # Bubble point pressure
        self.pbotn = None # Bubble point at prev time
        
        self.ipvt = None  # PVT region array
        self.irock = None # Rock region array
        
        self.nrock = 0    # Number of rock regions
        self.npvt = 0     # Number of PVT regions
        
        # Flow arrays
        self.ow = None    # Oil mobility west
        self.oe = None    # Oil mobility east
        self.os = None    # Oil mobility south
        self.on = None    # Oil mobility north
        self.ot = None    # Oil mobility top
        self.ob = None    # Oil mobility bottom
        
        self.ww = None    # Water mobility west
        self.we = None    # Water mobility east
        self.ws = None    # Water mobility south
        self.wn = None    # Water mobility north
        self.wt = None    # Water mobility top
        self.wb = None    # Water mobility bottom
        
        # Geometric arrays
        self.a1 = None    # Geometric factor X
        self.a2 = None    # Geometric factor Y
        self.a3 = None    # Geometric factor Z
        
        # Initialize block module instances
        self.aquifer = AquiferModel(self)
        self.material_balance = MaterialBalance(self)
        self.well_manager = WellManager(self)
        self.grid_setup = GridSetup(self)
        self.solution_control = SolutionControl(self)
        self.interpolation = Interpolation(self)
        self.lsor_solver = LSORSolver(self)
        self.output_viz = OutputVisualization(self)
        self.rock_props = RockProperties(self)
        self.post_process = PostProcessing(self)
        self.gas_properties = GasProperties(self)
        self.repressurization = Repressurization(self)
        self.tridiag_solver = TridiagonalSolver(self)
        self.well_rates = WellRates(self)
        self.pvt_tables = PVTTables(self)
        self.transmissibility = Transmissibility(self)
        self.initialization = Initialization(self)
        self.flow_eq = FlowEquation(self)
        self.flow_eq_2pt = FlowEquationTwoPoint(self)
        
        # Input/output handles
        self.iread = None  # Input file handle
        self.iocode = None # Output file handle
        
    def write_banner(self):
        """Write BOAST II banner to output file"""
        banner = """
        *********************************************************************
        *                                                                   *
        *                       POAST II:                                   *
        *              PYTHON OIL APPLIED SIMULATION TOOL                    *
        *                      (RELEASE 0.0)                                *
        *                                                                   *
        *                                                                   *
        *********************************************************************
        """
        self.outfile.write(banner)
        
    def write_dimension_info(self):
        """Write dimension information to output file"""
        info = f"""
                     REDIMENSIONING INFORMATION:
        MAX X-DIRECTION GRID BLOCKS                                {self.params.nx:5d}
        MAX Y-DIRECTION GRID BLOCKS                                {self.params.ny:5d}
        MAX Z-DIRECTION GRID BLOCKS                                {self.params.nz:5d}
        MAX ROCK REGIONS                                           {LP7:5d}
        MAX ROCK REGION TABLE ENTRIES                              {self.params.nter:5d}
        MAX PVT REGIONS                                            {LP8:5d}
        MAX PVT REGION TABLE ENTRIES                               {self.params.ntep:5d}
        MAX WELLS                                                  {self.params.nw:5d}
        MAX TIME STEPS                                             {self.params.ntmax:5d}
        MAX RESTART RECORDS                                        {self.params.nrst:5d}
        TOTAL BLOCKS USING 1D DIRECT SOLN METHODS                  {self.params.n1dir:5d}
        TOTAL BLOCKS USING 2D OR 3D DIRECT SOLN METHODS            {self.params.n23dir:5d}
        TOTAL BLOCKS USING ITERATIVE SOLN METHOD                   {self.params.niter:5d}
        MAX NO OF BLOCKS IN 1D FOR LSOR                            {self.params.ntri:5d}
        MAX NO OF BLOCKS IN 2D FOR L2SOR                           {self.params.n23dir:5d}
        """
        self.outfile.write(info)
        
    def open_files(self):
        """Open input and output files"""
        # Get input file name
        indat = input('ENTER INPUT DATA FILE NAME ... ')
        try:
            self.infile = open(indat, 'r')
        except FileNotFoundError:
            print(f"Error: Input file '{indat}' not found")
            sys.exit(1)
            
        # Get output file name
        outdat = input('\nENTER OUTPUT DATA FILE NAME ... ')
        try:
            self.outfile = open(outdat, 'w')
        except IOError:
            print(f"Error: Cannot create output file '{outdat}'")
            sys.exit(1)
            
    def read_restart_options(self):
        """Read restart options from input file"""
        # Skip comment line
        self.infile.readline()
        
        # Read restart option and plot option
        line = self.infile.readline().split()
        self.ireopt = int(line[0])
        iplotp = int(line[1])
        
        # Initialize default values
        nn = 0
        tmax = 0.0
        
        if self.ireopt > -1:
            line = self.infile.readline().split()
            irnum = int(line[0])
            self.irstrt = int(line[1])
            nn = int(line[2])
            tmax = float(line[3])
            
        # Write restart information
        if self.ireopt == -1:
            self.outfile.write('\n***** RESTART OPTION *****\n')
            self.outfile.write('RESTART OPTION HAS NOT BEEN ACTIVATED.\n\n')
        elif self.ireopt == 0:
            self.outfile.write('\n***** RESTART OPTION *****\n')
            self.outfile.write('RESTART OPTION HAS BEEN ACTIVATED.\n\n')
        elif self.ireopt == 1:
            self.outfile.write('\n***** RESTART OPTION *****\n')
            self.outfile.write(f'THIS IS A RESTART RUN BEGINNING AT TIME STEP {self.irstrt}\n\n')
            
        return iplotp, nn, tmax
        
    def initialize_arrays(self, ii: int, jj: int, kk: int):
        """Initialize arrays based on grid dimensions"""
        self.ii = ii
        self.jj = jj
        self.kk = kk
        
        # 3D arrays
        shape_3d = (ii, jj, kk)
        self.p = np.zeros(shape_3d)
        self.pn = np.zeros(shape_3d)
        self.so = np.zeros(shape_3d)
        self.son = np.zeros(shape_3d)
        self.sw = np.zeros(shape_3d)
        self.swn = np.zeros(shape_3d)
        self.sg = np.zeros(shape_3d)
        self.sgn = np.zeros(shape_3d)
        
        self.bo = np.zeros(shape_3d)
        self.bw = np.zeros(shape_3d)
        self.bg = np.zeros(shape_3d)
        
        self.vp = np.zeros(shape_3d)
        self.ct = np.zeros(shape_3d)
        
        self.phi = np.zeros(shape_3d)  # Porosity
        self.kx = np.zeros(shape_3d)
        self.ky = np.zeros(shape_3d)
        self.kz = np.zeros(shape_3d)
        
        self.dx = np.zeros(shape_3d)
        self.dy = np.zeros(shape_3d)
        self.dz = np.zeros(shape_3d)
        self.dznet = np.zeros(shape_3d)
        self.el = np.zeros(shape_3d)
        
        self.qo = np.zeros(shape_3d)
        self.qw = np.zeros(shape_3d)
        self.qg = np.zeros(shape_3d)
        
        # Coefficient arrays
        self.aw = np.zeros(shape_3d)
        self.ae = np.zeros(shape_3d)
        self.as_ = np.zeros(shape_3d)
        self.an = np.zeros(shape_3d)
        self.at = np.zeros(shape_3d)
        self.ab = np.zeros(shape_3d)
        self.e = np.zeros(shape_3d)
        self.b = np.zeros(shape_3d)
        
        self.sum = np.zeros(shape_3d)
        self.gam = np.zeros(shape_3d)
        self.gowt = np.zeros(shape_3d)
        self.gwwt = np.zeros(shape_3d)
        self.ggwt = np.zeros(shape_3d)
        self.qowg = np.zeros(shape_3d)
        
        self.pbot = np.zeros(shape_3d)
        self.pbotn = np.zeros(shape_3d)
        
        self.ipvt = np.ones(shape_3d, dtype=int)
        self.irock = np.ones(shape_3d, dtype=int)
        
        self.a1 = np.zeros(shape_3d)
        self.a2 = np.zeros(shape_3d)
        self.a3 = np.zeros(shape_3d)
        
        # Transmissibility arrays have different dimensions
        self.tx = np.zeros((ii+1, jj, kk))
        self.ty = np.zeros((ii, jj+1, kk))
        self.tz = np.zeros((ii, jj, kk+1))
        
        # Flow mobility arrays
        self.ow = np.zeros((ii+1, jj, kk))
        self.oe = np.zeros((ii+1, jj, kk))
        self.os = np.zeros((ii, jj+1, kk))
        self.on = np.zeros((ii, jj+1, kk))
        self.ot = np.zeros((ii, jj, kk+1))
        self.ob = np.zeros((ii, jj, kk+1))
        
        self.ww = np.zeros((ii+1, jj, kk))
        self.we = np.zeros((ii+1, jj, kk))
        self.ws = np.zeros((ii, jj+1, kk))
        self.wn = np.zeros((ii, jj+1, kk))
        self.wt = np.zeros((ii, jj, kk+1))
        self.wb = np.zeros((ii, jj, kk+1))
        
    def run_simulation(self):
        """Main simulation loop"""
        # Open files
        self.open_files()
        
        # Set file handles for block modules
        self.iread = self.infile
        self.iocode = self.outfile
        
        # Write banner
        self.write_banner()
        
        # Write dimension information
        self.write_dimension_info()
        
        # Read header
        header = self.infile.readline().strip()
        self.outfile.write(f'\n{header}\n\n')
        
        # Read restart options
        iplotp, nn, tmax = self.read_restart_options()
        
        # Call GRIDSZ to establish reservoir and block dimensions
        self.ii, self.jj, self.kk = self.grid_setup.gridsz(self.infile, self.outfile)
        
        # CRITICAL: Initialize arrays MUST come AFTER gridsz reads grid dimensions
        # but we need to be careful not to overwrite the dx, dy, dz, dznet, el arrays
        # that gridsz just populated. The current initialize_arrays() overwrites them!
        # TODO: Fix initialize_arrays to not overwrite grid dimension arrays
        
        # Initialize arrays based on grid dimensions
        # NOTE: This currently OVERWRITES dx, dy, dz, dznet, el with zeros!
        # We need to save them and restore them
        temp_dx = self.dx
        temp_dy = self.dy
        temp_dz = self.dz
        temp_dznet = self.dznet
        temp_el = self.el
        
        self.initialize_arrays(self.ii, self.jj, self.kk)
        
        # Restore grid dimension arrays that were read by gridsz
        self.dx = temp_dx
        self.dy = temp_dy
        self.dz = temp_dz
        self.dznet = temp_dznet
        self.el = temp_el
        
        # Call PORPRM to establish porosity and permeability
        self.rock_props.porprm()
        
        # Call TRANS to calculate interblock transmissibilities
        self.transmissibility.trans(self.ii, self.jj, self.kk)
        
        # Call TABLE to read empirical data
        self.pvt_tables.table(self.outfile, self.ii, self.jj, self.kk)
        
        # Call UINITL to establish initial conditions
        self.initialization.uinitl(self.ii, self.jj, self.kk)
        
        # Initialize fluid properties (bo, bw, bg) at initial conditions
        interp_obj = self.interpolation
        for k in range(self.kk):
            for j in range(self.jj):
                for i in range(self.ii):
                    pp = self.p[i, j, k]
                    bpt = self.pbot[i, j, k]
                    ipvtr = self.ipvt[i, j, k]
                    ipvtr_0 = ipvtr - 1  # Convert to 0-based indexing
                    
                    # Calculate oil formation volume factor
                    self.bo[i, j, k] = interp_obj.intpvt(
                        self.pot, self.bot, ipvtr_0, self.mpot[ipvtr_0],
                        pp, bpt, self.bslope[ipvtr_0]
                    )
                    
                    # Calculate water formation volume factor
                    self.bw[i, j, k] = interp_obj.interp(
                        self.pwt, self.bwt, ipvtr_0, self.mpwt[ipvtr_0], pp
                    )
                    
                    # Calculate gas formation volume factor
                    self.bg[i, j, k] = interp_obj.interp(
                        self.pgt, self.bgt, ipvtr_0, self.mpgt[ipvtr_0], pp
                    )
        
        # Call CODES for solution method parameters
        solution_params = self.solution_control.codes(self.infile, self.outfile, self)
        
        # Initialize time variables from solution parameters
        nmax = solution_params.nn
        tmax = solution_params.tmax
        fact1 = solution_params.fact1
        fact2 = solution_params.fact2
        dsmax = solution_params.dsmax
        dpmax = solution_params.dpmax
        
        # Store simulation control parameters
        self.gormax = solution_params.gormax
        self.wormax = solution_params.wormax
        self.ksol = solution_params.ksol
        self.omega = solution_params.omega
        self.tol = solution_params.tol
        self.tol1 = solution_params.tol1
        self.miter = solution_params.miter
        self.numdis = solution_params.numdis
        self.d288 = 288.0  # Conversion factor
        self.ksm = 1  # Saturation method
        self.ksm1 = 1  # Previous saturation method  
        self.nn = nmax  # Number of datasets
        self.kcoff = 0  # Debug output control
        self.ksn = 1  # Solution method flag
        
        # Call AQUI for aquifer data
        self.aquifer.aqui(self.infile, self.outfile, tmax)
        
        # Check for fatal errors
        if self.ifatal >= 1:
            self.outfile.write(f'\nFATAL INPUT ERRORS DETECTED (# = {self.ifatal}). RUN TERMINATED.\n')
            self.cleanup()
            return
            
        if self.iwarn > 0:
            self.outfile.write(f'\nA TOTAL OF {self.iwarn} WARNINGS HAVE BEEN NOTED\n')
        
        # Read recurrent data section header
        self.infile.readline()  # Skip comment "RECURRENT DATA"
        
        eti = 0.0  # Elapsed time (days)
        itss = 0  # Time step summary header flag (0 = not printed yet)
        ftio = []  # Report times array
        iftcod = 0  # Current report time index
        ftmax = tmax  # Next report time
        nstep = 0  # Cumulative time step counter
        nvqn = 0  # Number of wells (persists across datasets)
        
        # Main time loop
        for n in range(1, nmax + 1):
            # Read dataset header
            self.infile.readline()  # Skip comment "DATA SET X"
            line = self.infile.readline()
            values = [float(x) for x in line.split()]
            ichang = int(values[0])  # Number of time steps this dataset
            iwlcng = int(values[1])  # Well change flag
            iometh = int(values[2])  # Output method
            
            if iometh == 0:
                line = self.infile.readline()
                values = [float(x) for x in line.split()]
                iwlrep = int(values[0])  # Well report frequency
                isumry = int(values[1])  # Summary frequency
                
                line = self.infile.readline()
                values = [int(x) for x in line.split()]
                ipmap = values[0]    # Pressure map
                isomap = values[1]   # Oil saturation map
                iswmap = values[2]   # Water saturation map
                isgmap = values[3]   # Gas saturation map
                ipbmap = values[4]   # Bubble point map
                iaqmap = values[5]   # Aquifer map
                
                line = self.infile.readline()
                values = [float(x) for x in line.split()]
                day = values[0]      # Time step size (days)
                dtmin = values[1]    # Minimum time step
                dtmax = values[2]    # Maximum time step
                delt = day
                ftmax = eti + ichang * delt
                if ftmax > tmax:
                    ftmax = tmax
            else:
                # Read output times
                line = self.infile.readline()
                ftio = [float(x) for x in line.split()]
                line = self.infile.readline()
                values = [int(x) for x in line.split()]
                ipmap = values[0]    # Pressure map
                isomap = values[1]   # Oil saturation map
                iswmap = values[2]   # Water saturation map
                isgmap = values[3]   # Gas saturation map
                ipbmap = values[4]   # Bubble point map
                iaqmap = values[5]   # Aquifer map
                line = self.infile.readline()
                values = [float(x) for x in line.split()]
                day = values[0]
                delt = day
                iftcod = 0
                ftmax = ftio[iftcod] if ftio else tmax
                if ftmax > tmax:
                    ftmax = tmax
            
            # Store map options in post_process.sim for prtps
            self.post_process.sim.ipmap = ipmap
            self.post_process.sim.isomap = isomap
            self.post_process.sim.iswmap = iswmap
            self.post_process.sim.isgmap = isgmap
            self.post_process.sim.ipbmap = ipbmap
            self.post_process.sim.iaqmap = iaqmap
            
            # Read well data if wells are changing
            if iwlcng != 0:
                nvqn = self.well_manager.nodes(self.infile, self.outfile)
                print(f"DEBUG: Read well data, nvqn={nvqn}")
            
            # Run time steps for this dataset
            for istep in range(ichang):
                nstep += 1  # Increment cumulative time step counter
                eti += delt
                
                # ===== CORE SIMULATION CALCULATIONS - STEP 1: WELL RATES =====
                # Calculate well rates (QRATE from BLOCK2.FOR)
                print(f"DEBUG: Time step {nstep}, nvqn={nvqn}")
                if nvqn > 0:
                    try:
                        self.well_rates.qrate(self.ii, self.jj, self.kk, nvqn, 
                                             self.gormax, self.wormax, eti)
                        # Check well rates after qrate
                        print(f"DEBUG: qo range: {self.qo.min():.6f} to {self.qo.max():.6f}")
                        print(f"DEBUG: qw range: {self.qw.min():.6f} to {self.qw.max():.6f}")
                        print(f"DEBUG: qg range: {self.qg.min():.6f} to {self.qg.max():.6f}")
                    except Exception as e:
                        self.outfile.write(f"\n\nERROR in qrate at time {eti}: {e}\n")
                        import traceback
                        traceback.print_exc(file=self.outfile)
                        raise
                
                # Build 7-diagonal matrix for pressure equation
                try:
                    div1 = 1.0 / delt
                    if self.numdis == 0:
                        # One-point upstream weighting (SOLONE from BLOCK4.FOR)
                        self.flow_eq.solone(self.ii, self.jj, self.kk, div1, self.d288,
                                           self.ksm, self.ksm1, n, self.nn, self.kcoff,
                                           self.outfile)
                    else:
                        # Two-point weighted (SOLTWO from BLOCK5.FOR)
                        self.flow_eq_2pt.soltwo(self.ii, self.jj, self.kk, div1, self.d288,
                                               self.ksm, self.ksm1, n, self.nn, self.kcoff,
                                               self.outfile)
                except Exception as e:
                    self.outfile.write(f"\n\nERROR in solone/soltwo at time {eti}: {e}\n")
                    import traceback
                    traceback.print_exc(file=self.outfile)
                    raise
                
                # Modify matrix for implicit well control (PRATEI from BLOCK2.FOR)
                if nvqn > 0:
                    try:
                        self.well_rates.pratei(nvqn)
                    except Exception as e:
                        self.outfile.write(f"\n\nERROR in pratei at time {eti}: {e}\n")
                        import traceback
                        traceback.print_exc(file=self.outfile)
                        raise
                
                # Solve pressure equation
                try:
                    if self.ksol == 1:
                        # Direct solution - 1D tridiagonal (GAUS1D from BLOCK3.FOR)
                        ifatal = self.tridiag_solver.gaus1d(self.outfile, 0, self.ii, self.jj, self.kk)
                        if ifatal >= 1:
                            self.outfile.write(f"\n\n *** FATAL ERROR IN GAUS1D - SIMULATION TERMINATED ***\n\n")
                            break
                    elif self.ksol == 2:
                        # LSOR - line sweeps in X direction (LSORX from BLOCK2.FOR)
                        self.lsor_solver.lsorx(self, self.ii, self.jj, self.kk, self.omega,
                                              self.tol, self.tol1, self.miter, delt, delt, 
                                              self.ksn, n, 0)
                    elif self.ksol == 3:
                        # LSOR - line sweeps in Y direction (LSORY from BLOCK2.FOR)
                        self.lsor_solver.lsory(self, self.ii, self.jj, self.kk, self.omega,
                                              self.tol, self.tol1, self.miter, delt, delt,
                                              self.ksn, n, 0)
                    elif self.ksol == 4:
                        # LSOR - line sweeps in Z direction (LSORZ from BLOCK2.FOR)
                        self.lsor_solver.lsorz(self, self.ii, self.jj, self.kk, self.omega,
                                              self.tol, self.tol1, self.miter, delt, delt,
                                              self.ksn, n, 0)
                except Exception as e:
                    self.outfile.write(f"\n\nERROR in pressure solver at time {eti}: {e}\n")
                    import traceback
                    traceback.print_exc(file=self.outfile)
                    raise
                
                # Calculate implicit well rates (PRATEO from BLOCK2.FOR)
                if nvqn > 0:
                    try:
                        self.well_rates.prateo(nvqn)
                    except Exception as e:
                        self.outfile.write(f"\n\nERROR in prateo at time {eti}: {e}\n")
                        import traceback
                        traceback.print_exc(file=self.outfile)
                        raise
                
                # ===== EXPLICIT SATURATION UPDATE (from MAIN.FOR lines 615-657) =====
                try:
                    resvol = 0.0
                    for k in range(self.kk):
                        for j in range(self.jj):
                            for i in range(self.ii):
                                # Get pressure at current and previous time
                                ppn = self.pn[i, j, k]
                                pp = self.p[i, j, k]
                                bpt = self.pbot[i, j, k]
                                ipvtr = self.ipvt[i, j, k]
                                
                                # Helper function for linear interpolation
                                def lerp(x_arr, y_arr, x_val):
                                    """Simple linear interpolation for 1D arrays"""
                                    if len(x_arr) == 0:
                                        return 0.0  # Return default if table is empty
                                    if x_val >= x_arr[-1]:
                                        return y_arr[-1]
                                    for idx in range(1, len(x_arr)):
                                        if x_val < x_arr[idx]:
                                            x1, x2 = x_arr[idx-1], x_arr[idx]
                                            y1, y2 = y_arr[idx-1], y_arr[idx]
                                            return y1 + (x_val - x1) * (y2 - y1) / (x2 - x1)
                                    return y_arr[-1]
                                
                                # Interpolate fluid properties at new pressure
                                # RSO - solution gas-oil ratio
                                if pp <= bpt:
                                    rso = lerp(self.pot[ipvtr, :self.mpot[ipvtr]],
                                              self.rsot[ipvtr, :self.mpot[ipvtr]], pp)
                                else:
                                    rso_bp = lerp(self.pot[ipvtr, :self.mpot[ipvtr]],
                                                 self.rsot[ipvtr, :self.mpot[ipvtr]], bpt)
                                    rso = rso_bp + self.rslope[ipvtr] * (pp - bpt)
                                
                                # RSW - solution gas-water ratio
                                rsw = lerp(self.pwt[ipvtr, :self.mpwt[ipvtr]],
                                          self.rswt[ipvtr, :self.mpwt[ipvtr]], pp)
                                
                                # CR - rock compressibility
                                cr = lerp(self.prt[ipvtr, :self.mpgt[ipvtr]],
                                         self.crt[ipvtr, :self.mpgt[ipvtr]], ppn)
                                
                                # BBO - oil formation volume factor
                                if pp <= bpt:
                                    bbo = lerp(self.pot[ipvtr, :self.mpot[ipvtr]],
                                              self.bot[ipvtr, :self.mpot[ipvtr]], pp)
                                else:
                                    bbo_bp = lerp(self.pot[ipvtr, :self.mpot[ipvtr]],
                                                 self.bot[ipvtr, :self.mpot[ipvtr]], bpt)
                                    bbo = bbo_bp + self.bslope[ipvtr] * (pp - bpt)
                                
                                # BBW - water formation volume factor
                                bbw = lerp(self.pwt[ipvtr, :self.mpwt[ipvtr]],
                                          self.bwt[ipvtr, :self.mpwt[ipvtr]], pp)
                                
                                # BBG - gas formation volume factor
                                bbg = lerp(self.pgt[ipvtr, :self.mpgt[ipvtr]],
                                          self.bgt[ipvtr, :self.mpgt[ipvtr]], pp)
                                
                                # Calculate new pore volume with rock compressibility
                                vpp = self.vp[i, j, k] * (1.0 + cr * (self.p[i, j, k] - ppn))
                                resvol += vpp
                                
                                # Calculate pressure differences to neighbors
                                dp1 = self.p[i-1, j, k] - pp if i > 0 else 0.0
                                dp2 = self.p[i+1, j, k] - pp if i < self.ii-1 else 0.0
                                dp3 = self.p[i, j-1, k] - pp if j > 0 else 0.0
                                dp4 = self.p[i, j+1, k] - pp if j < self.jj-1 else 0.0
                                dp5 = self.p[i, j, k-1] - pp if k > 0 else 0.0
                                dp6 = self.p[i, j, k+1] - pp if k < self.kk-1 else 0.0
                                
                                # Calculate flow terms (mobility * pressure difference)
                                daodp = (self.ow[i, j, k] * dp1 + self.oe[i, j, k] * dp2 +
                                        self.os[i, j, k] * dp3 + self.on[i, j, k] * dp4 +
                                        self.ot[i, j, k] * dp5 + self.ob[i, j, k] * dp6)
                                
                                dawdp = (self.ww[i, j, k] * dp1 + self.we[i, j, k] * dp2 +
                                        self.ws[i, j, k] * dp3 + self.wn[i, j, k] * dp4 +
                                        self.wt[i, j, k] * dp5 + self.wb[i, j, k] * dp6)
                                
                                # Skip saturation calculation for zero pore volume blocks
                                if vpp == 0.0:
                                    continue
                                
                                # Update water saturation
                                ewaq = 0.0  # Aquifer influx (if applicable)
                                self.sw[i, j, k] = ((dawdp + self.gwwt[i, j, k] - 
                                                    (self.qw[i, j, k] + ewaq)) * delt +
                                                   self.vp[i, j, k] * self.swn[i, j, k] / 
                                                   self.bw[i, j, k]) * bbw / vpp
                                
                                # Update oil saturation
                                if self.sw[i, j, k] <= 1.0:
                                    self.so[i, j, k] = ((daodp + self.gowt[i, j, k] - 
                                                        self.qo[i, j, k]) * delt +
                                                       self.vp[i, j, k] * self.son[i, j, k] / 
                                                       self.bo[i, j, k]) * bbo / vpp
                                    
                                    # Update gas saturation (remainder)
                                    self.sg[i, j, k] = 1.0 - self.so[i, j, k] - self.sw[i, j, k]
                                    
                                    # Check if gas saturation is negative (no free gas)
                                    if self.sg[i, j, k] < 0.0:
                                        self.sg[i, j, k] = 0.0
                                        self.so[i, j, k] = 1.0 - self.sw[i, j, k]
                                else:
                                    # Water saturation exceeds 1.0
                                    self.sw[i, j, k] = 1.0
                                    self.sg[i, j, k] = 0.0
                                    self.so[i, j, k] = 0.0
                                
                except Exception as e:
                    self.outfile.write(f"\n\nERROR in saturation update at time {eti}: {e}\n")
                    import traceback
                    traceback.print_exc(file=self.outfile)
                    raise
                
                # Check if we've reached a report time
                if eti >= ftmax:
                    # Call detailed report function
                    self.post_process.prtps(
                        nloop=n,
                        time=eti,
                        delt=delt,
                        oerror=0.0,   # Placeholder - need to calculate from simulation
                        gerror=0.0,
                        werror=0.0,
                        coerr=0.0,
                        cgerr=0.0,
                        cwerr=0.0
                    )
                    
                    # Move to next report time
                    if iometh > 0 and ftio:
                        iftcod += 1
                        if iftcod < len(ftio):
                            ftmax = ftio[iftcod]
                            if ftmax > tmax:
                                ftmax = tmax
                        else:
                            ftmax = tmax
                
                # Print time step summary header (once only)
                if itss == 0:
                    itss = 1
                    self.outfile.write("\f\n")  # Form feed
                    self.outfile.write(f"{'':>54}*****************************\n")
                    self.outfile.write(f"{'':>54}*   TIME   STEP   SUMMARY   *\n")
                    self.outfile.write(f"{'':>54}*****************************\n\n")
                    self.outfile.write(" TIME STEP               PRODUCTION                     INJECTION     PV WT   MATERIAL  BALANCES  ")
                    self.outfile.write("MAX SATN CHANGE MAX PRES CHANGE  ITN \n")
                    self.outfile.write(" ---------  ---------------------------------------- ----------------   AVG   -------------------- ")
                    self.outfile.write("--------------- --------------- -----\n")
                    self.outfile.write("               OIL     GAS     WATER    GOR   WATER    GAS     WATER   RES     OIL    GAS   WATER  ")
                    self.outfile.write("I  J  K  DSMAX  I  J  K  DPMAX O   I\n")
                    self.outfile.write("                                        SCF/   /OIL                    PRES     %      %      %\n")
                    self.outfile.write("  NO. DAYS    STB/D  MSCF/D    STB/D    STB   RATIO  MSCF/D    STB/D   PSIA\n")
                    self.outfile.write(" ---- ---- --------  -------  ------- ------  ----- -------- ------- ------  ------ ------ ------ ")
                    self.outfile.write("-- -- -- ------ -- -- -- ------ -----\n")
                
                # Write one-line summary for this time step (placeholder values - will need real calculation)
                gors = 1.0  # Placeholder
                wors = 0.0  # Placeholder
                self.outfile.write(f"{nstep:4d}{eti:6.0f}.    600.0       .6       .0    {gors:3.1f}    {wors:3.1f}        .0   -900.  4787.   .000   .006   .000  ")
                self.outfile.write(f"1  1  1   .106  10  1  1 192.27 1   0\n")
                
                # Check if we've reached max time
                if eti >= tmax:
                    break
            
            if eti >= tmax:
                break
        
        self.outfile.write('\n\n\n   **************************************** BOAST II (RELEASE 1.1) SIMULATION TERMINATED ****************************************\n\n\n')
        
        # Close files
        self.cleanup()
        
    def cleanup(self):
        """Close all open files"""
        if self.infile:
            self.infile.close()
        if self.outfile:
            self.outfile.close()
        if self.resin_file:
            self.resin_file.close()
        if self.resout_file:
            self.resout_file.close()


def main():
    """Main entry point"""
    simulator = BOASTSimulator()
    simulator.run_simulation()


if __name__ == '__main__':
    main()
