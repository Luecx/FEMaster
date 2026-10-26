"""Compare separate equilibrium-solver times; one warmup plus five samples.
Usage: python3 benchmarks/pretension_compare.py BASELINE_EXE OPTIMIZED_EXE MODEL_DIR CSV
Generate decks with pretension_solver.py first. Keeps logs and RES for inspection.
"""
from pathlib import Path
import subprocess
import os
import re
import csv
import sys
baseline, optimized, model_directory, csv_path = map(Path, sys.argv[1:])
# Equilibrium timer excludes preparation/compile/output; the CLI reports integer ms.
with csv_path.open('w') as f:
    writer=csv.writer(f, lineterminator="\n");writer.writerow(['revision','n','repeat','solver_ms','nodes','elements','q','force'])
    for rev,exe in [('baseline',baseline.resolve()),('optimized',optimized.resolve())]:
        for n in [4,8,16,32,64]:
            for repeat in range(6):
                output=model_directory / f'{rev}-{n}-{repeat}'
                result=subprocess.run([str(exe),str(model_directory / f'bar-{n}.inp'),'--output',str(output),'--output-format','res'],env={**os.environ,'OMP_NUM_THREADS':'1'},capture_output=True,text=True,check=True)
                (output.with_suffix('.log')).write_text(result.stdout+result.stderr)
                timer=re.search(r'solving constraint system.*?\[\s*([\d.]+)\s*ms\]',result.stdout)
                if not timer: raise ValueError(result.stdout[-12000:])
                text=output.with_suffix('.res').read_text()
                fields={}
                for name in ['PTGAP','PTFORC']:
                    rows=text.split('FIELD, NAME='+name+',')[1].split('END FIELD')[0].splitlines()[1:]
                    fields[name]={int(parts[0]): [float(value) for value in parts[1:]]
                                  for row in rows if (parts := row.split())}
                first_b=1+5*(n+1)**2
                side_b=range(first_b, first_b+(n+1)**2)
                q=fields['PTGAP'][first_b][0]
                force=sum(fields['PTFORC'][node][2] for node in side_b)
                assert abs(q-1000*3/210000)<1e-6,(n,q)
                assert abs(force-1000)<1e-2,(n,force)
                writer.writerow([rev,n,repeat,timer[1],5*(n+1)**2,3*n*n,q,force]);f.flush()
                print(rev,n,repeat,timer[1],flush=True)

# Compare every exported field row, including semantic IDs, on identical meshes.
# Numerical zero values can differ in their last bits between process runs.
import math
for n in [4, 8, 16, 32, 64]:
    paths = [model_directory / f'{rev}-{n}-1.res' for rev in ['baseline', 'optimized']]
    texts = [path.read_text().splitlines() for path in paths]
    assert len(texts[0]) == len(texts[1])
    largest = 0.0
    for before, after in zip(*texts):
        if before == after:
            continue
        a, b = before.split(), after.split()
        assert len(a) == len(b) and a[0] == b[0], (n, before, after)
        for x, y in zip(a[1:], b[1:]):
            x, y = float(x), float(y)
            if math.isnan(x) and math.isnan(y):
                continue  # Inactive result components retain their NaN markers.
            assert math.isclose(x, y, rel_tol=1e-6, abs_tol=1e-5), (n, x, y)
            largest = max(largest, abs(x-y))
    print(f'RES comparison n={n}: maximum absolute field difference {largest:.6g}', flush=True)
