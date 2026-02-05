"""Quick test to see what error occurs"""
import sys
sys.stdout = open('test_output.txt', 'w', buffering=1)
sys.stderr = sys.stdout

try:
    print("Starting import...")
    from main import main
    
    print("Running simulation...")
    main()
    
    print("Completed successfully!")
except Exception as e:
    print(f"ERROR: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
finally:
    sys.stdout.close()
