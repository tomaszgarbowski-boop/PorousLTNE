"""Regenerate independent fixtures; run from any directory with Python/SciPy.

The reference uses SciPy BDF, not MATLAB/Octave ode15s. Agreement checks a
translation of the same reduced equations, not independent physical validity.
"""
from pathlib import Path
import platform
import json
import numpy as np
import scipy
from scipy.io import savemat
from reference_python import Model, Params

root = Path(__file__).resolve().parents[1]
destination = root / 'tests' / 'fixtures'
destination.mkdir(exist_ok=True)
m = Model(Params(N=6))
y = m.initial()
z = (np.arange(6) + .5) / 6
for j in range(3):
    v = y.copy()
    v[:6] = 10 + j*5 + 3*z
    v[6:12] = 8 + j*4 + 2*z
    v[18:24] = 2 + .02*z
    v[12:18] *= 1 + .02*j
    v[24:30] *= 1 - .01*z
    c = np.array([30-.5*j, .95-.1*j, .03+.005*j])
    savemat(destination / f'rhs_{j}.mat',
            dict(y=v[:,None], control=c[None,:], rhs=m.rhs(0,v,c)[:,None]))
for n in [10,20,40,60]:
    r = Model(Params(N=n)).run(horizon=24)
    savemat(destination / f'python_grid{n}.mat',
            dict(t=r['t'][:,None]*3600, Ts=r['y'][:,n:2*n],
                 activationHours=r['metrics']['t95'], conditioningKWh=r['metrics']['E']))
    print(n, r['metrics']['t95'], flush=True)
(destination / 'environment.json').write_text(json.dumps(
    dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
    indent=2))
