import sys
sys.stdin = open('EX1.DAT', 'r')

# Import main
from main import Simulator

# Create simulator
sim = Simulator()

# Run one timestep only
sim.run_simulation()
