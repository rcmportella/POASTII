"""Test runner with automatic input"""
import sys
from main import BOASTSimulator

# Override input to provide automatic responses
original_input = __builtins__.input
def mock_input(prompt):
    if 'INPUT DATA' in prompt:
        return 'EX1.DAT'
    elif 'OUTPUT' in prompt:
        return 'OUT1_14.DAT'
    return original_input(prompt)

__builtins__.input = mock_input

# Run simulation
try:
    simulator = BOASTSimulator()
    simulator.run_simulation()
    print("\n=== SIMULATION COMPLETED SUCCESSFULLY ===")
except Exception as e:
    print(f"\n=== ERROR: {type(e).__name__}: {e} ===")
    import traceback
    traceback.print_exc()
