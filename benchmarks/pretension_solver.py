"""Generate the same n*n*3 C3D8 bars as pretension_prepare.cpp for CLI solves.

Usage: python3 benchmarks/pretension_solver.py /tmp/pretension-models
Timing uses the CLI's separate 'solving constraint system' timer (integer ms).
"""
from pathlib import Path
import sys

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
for n in (4, 8, 16, 32, 64):
    def node(x, y, z):
        return 1 + z*(n+1)*(n+1) + y*(n+1) + x
    lines = ['*NODE, NSET=NALL']
    for z in range(4):
        for y in range(n+1):
            for x in range(n+1):
                lines.append(f'{node(x,y,z)}, {x/n}, {y/n}, {z}')
    lines.append('*ELEMENT, TYPE=C3D8, ELSET=SOLID')
    e = 0
    cut = []
    for z in range(3):
        for y in range(n):
            for x in range(n):
                e += 1
                ids = [node(x,y,z), node(x+1,y,z), node(x+1,y+1,z), node(x,y+1,z),
                       node(x,y,z+1), node(x+1,y,z+1), node(x+1,y+1,z+1), node(x,y+1,z+1)]
                lines.append(', '.join(map(str, [e]+ids)))
                if z == 1:
                    cut.append(f'{len(cut)+1}, {e}, 1')
    lines += ['*SURFACE, SFSET=CUT'] + cut
    lines += ['*PRETENSIONSECTION, NAME=PT, SURFACE_A=CUT', '*MATERIAL, NAME=STEEL',
              '*ELASTIC, TYPE=ISOTROPIC', '210000, 0', '*SOLIDSECTION, ELSET=SOLID, MATERIAL=STEEL',
              '*NSET, NAME=ENDS']
    ends = [node(x,y,z) for z in (0,3) for y in range(n+1) for x in range(n+1)]
    lines += [', '.join(map(str, ends[i:i+16])) for i in range(0,len(ends),16)]
    lines += ['*SUPPORT, SUPPORT_COLLECTOR=BC', 'ENDS, 0, 0, 0, NAN, NAN, NAN',
              '*LOADCASE, TYPE=LINEARSTATIC', '*SUPPORTS', 'BC',
              '*PRETENSION, SECTION=PT, ACTION=LOAD, CONTROL=FORCE, VALUE=1000',
              '*SOLVER, DEVICE=CPU, METHOD=DIRECT', '*CONSTRAINTMETHOD, TYPE=ELIMINATION', '*END']
    (out/f'bar-{n}.inp').write_text('\n'.join(lines)+'\n')
