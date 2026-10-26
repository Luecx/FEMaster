"""Generate the same n*n*3 C3D8 bars as pretension_prepare.cpp for CLI solves.

Usage: python3 benchmarks/pretension_solver.py /tmp/pretension-models
Timing uses the CLI's separate 'solving constraint system' timer (integer ms).
"""
from pathlib import Path
import sys

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
for n in (4, 8, 16, 32, 64):
    layer_size = (n + 1) * (n + 1)
    def node(x, y, z, upper=False):
        return 1 + z*layer_size + y*(n+1) + x + (4*layer_size if upper and z == 1 else 0)
    lines = ['*NODE, NSET=NALL']
    for z in range(4):
        for y in range(n+1):
            for x in range(n+1):
                lines.append(f'{node(x,y,z)}, {x/n}, {y/n}, {z}')
    for y in range(n+1):
        for x in range(n+1):
            lines.append(f'{node(x,y,1,True)}, {x/n}, {y/n}, 1')
    lines.append('*ELEMENT, TYPE=C3D8, ELSET=SOLID')
    e = 0
    cut_a = []
    cut_b = []
    for z in range(3):
        for y in range(n):
            for x in range(n):
                e += 1
                ids = [node(x,y,z,z==1), node(x+1,y,z,z==1), node(x+1,y+1,z,z==1), node(x,y+1,z,z==1),
                       node(x,y,z+1), node(x+1,y,z+1), node(x+1,y+1,z+1), node(x,y+1,z+1)]
                lines.append(', '.join(map(str, [e]+ids)))
                if z == 0:
                    cut_a.append(f'{len(cut_a)+1}, {e}, 2')
                if z == 1:
                    cut_b.append(f'{n*n+len(cut_b)+1}, {e}, 1')
    lines += ['*SURFACE, SFSET=CUT_A'] + cut_a + ['*SURFACE, SFSET=CUT_B'] + cut_b
    lines += ['*PRETENSIONSECTION, NAME=PT, SURFACE_A=CUT_A, SURFACE_B=CUT_B', '*MATERIAL, NAME=STEEL',
              '*ELASTIC, TYPE=ISOTROPIC', '210000, 0', '*SOLIDSECTION, ELSET=SOLID, MATERIAL=STEEL',
              '*NSET, NAME=ENDS']
    ends = [node(x,y,z) for z in (0,3) for y in range(n+1) for x in range(n+1)]
    lines += [', '.join(map(str, ends[i:i+16])) for i in range(0,len(ends),16)]
    lines += ['*SUPPORT, SUPPORT_COLLECTOR=BC', 'ENDS, 0, 0, 0, NAN, NAN, NAN',
              '*LOADCASE, TYPE=LINEARSTATIC', '*SUPPORTS', 'BC',
              '*PRETENSION, SECTION=PT, MODE=FORCE, VALUE=1000',
              '*SOLVER, DEVICE=CPU, METHOD=DIRECT', '*CONSTRAINTMETHOD, TYPE=ELIMINATION', '*END']
    (out/f'bar-{n}.inp').write_text('\n'.join(lines)+'\n')
