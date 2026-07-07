from pathlib import Path
import subprocess
import sys

exe = sys.argv[1] if len(sys.argv) > 1 else 'FEMaster'
root = Path(__file__).resolve().parent

for inp in sorted(root.glob('**/*.inp')):
    rel = inp.relative_to(root)
    print(f'Running {rel}')
    subprocess.run([exe, str(inp)], check=False)
