# View one snapshot of the hydrodynamic jet (files jet_001.txt ... written by jet2d.f90)
# Usage: python jet2d_view.py        (last snapshot)
#        python jet2d_view.py 12     (snapshot number 12)
import sys, glob, io, re
import numpy as np
import matplotlib.pyplot as plt

def load(f):
    # repairs old files where 3-digit exponents were written without the letter E (1.0-100)
    with open(f) as fh:
        txt = re.sub(r'(\d)([+-]\d{3})\b', r'\1E\2', fh.read())
    return np.loadtxt(io.StringIO(txt), unpack=True)

files = sorted(glob.glob('jet_*.txt'))
if not files:
    sys.exit('No jet_*.txt file found: run jet2d first')
n = int(sys.argv[1]) if len(sys.argv) > 1 else len(files)
f = files[n - 1]
t = n * 2.4 / len(files)

x, y, rho, vx, vy, p = load(f)
ux, uy = np.unique(x), np.unique(y)
ext = [x.min(), x.max(), y.min(), y.max()]

fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
for a, q, name in zip(ax, [np.log10(rho), np.log10(p)], ['log$_{10}$ density', 'log$_{10}$ pressure']):
    im = a.imshow(q.reshape(len(ux), len(uy)).T, origin='lower', extent=ext, cmap='inferno')
    a.set_title(name); a.set_xlabel('x'); a.set_ylabel('y')
    fig.colorbar(im, ax=a, fraction=0.046, pad=0.02)
fig.suptitle(f'Hydrodynamic jet (γ = 1.4, η = 0.1, Mach 5), t = {t:.1f}')
fig.tight_layout()
fig.savefig('jet_hydro.png', dpi=200)
plt.show()
