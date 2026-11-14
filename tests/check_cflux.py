"""Check CFLUX against analytical unit-cube conduction solutions.

Run with: python3 tests/check_cflux.py /path/to/FEMaster
Only the Python standard library and an existing solver binary are required.
"""

from pathlib import Path
import argparse
import math
import subprocess
import tempfile

cli = argparse.ArgumentParser(description=__doc__)
cli.add_argument("solver", type=Path)
args = cli.parse_args()
workspace = Path(__file__).resolve().parents[1]

solver = args.solver.resolve()
base = (workspace / 'examples/13_thermal/131_steady_state_thermal_cflux_c3d8.inp').read_text()
model=base.split('*STEP')[0]
step='*STEP\n*HEAT TRANSFER,STEADY STATE\n*BOUNDARY\nCOLD,11,,100\n'

def temperatures(text):
    fields=[]
    for block in text.split('FIELD, ')[1:]:
        header,_,tail=block.partition('\n')
        if 'NAME=TEMPERATURE' in header:
            rows=[list(map(float,row.split())) for row in tail.split('END FIELD')[0].splitlines() if row.strip()]
            fields.append(rows)
    return fields

def run(name,deck,expected=None,reject=None):
    path = scratch / ('cflux_' + name + '.inp')
    path.write_text(deck)
    p=subprocess.run([str(solver),str(path),'--output-format','res'],capture_output=True,text=True)
    path.with_suffix('.log').write_text(p.stdout + p.stderr)
    if reject:
        assert p.returncode != 0, (name, 'invalid input was accepted')
        assert reject in p.stdout+p.stderr,(name,p.returncode,p.stdout[-1800:],p.stderr)
    else:
        assert p.returncode==0,(name,p.stdout[-1800:],p.stderr)
        result_text=path.with_suffix('.res').read_text()
        data=temperatures(result_text)
        assert len(data)==len(expected),(name,len(data),expected)
        for rows,(low,high) in zip(data,expected):
            # Temperature field rows may include explicit node IDs.
            vals=[row[-1] for row in rows]
            assert abs(min(vals)-low)<2e-4 and abs(max(vals)-high)<2e-4,(name,min(vals),max(vals))
        flux_fields=[]
        for block in result_text.split('FIELD, ')[1:]:
            header,_,tail=block.partition('\n')
            if 'NAME=HEAT_FLUX' in header:
                flux_fields.append([list(map(float,row.split()))[-3:] for row in tail.split('END FIELD')[0].splitlines() if row.strip()])
        assert len(flux_fields)==len(expected),(name, 'missing heat flux')
        for rows,(low,high) in zip(flux_fields,expected):
            expected_flux=10*(high-low)
            assert all(abs(math.sqrt(sum(v*v for v in row))-expected_flux)<2e-3 for row in rows),(name,'heat flux mismatch',rows)
    print('PASS',name,flush=True)

with tempfile.TemporaryDirectory(prefix='femaster_cflux_') as directory:
    scratch = Path(directory)
    run('history',base,[(100,110),(90,100),(100,100)])
    run('amplitude',model+'*AMPLITUDE,NAME=A\n0,0,1,2\n'+step+'*CFLUX,AMPLITUDE=A\nHOT,11,25\n*END STEP\n',[(100,120)])
    run('root_and_persistence',model+'*CFLUX\nHOT,11,25\n'+step+'*END STEP\n*STEP\n*HEAT TRANSFER,STEADY STATE\n*END STEP\n',[(100,110),(100,110)])
    run('superposition',model+step+'*CFLUX\nHOT,11,25\n2,11,5\n3,11,5\n6,11,5\n7,11,5\n*END STEP\n',[(100,112)])
    run('bad_dof',model+step+'*CFLUX\nHOT,1,25\n*END STEP\n',reject='only temperature DOF 11')
    run('unknown_amplitude',model+step+'*CFLUX,AMPLITUDE=MISSING\nHOT,11,25\n*END STEP\n',reject='does not exist')
