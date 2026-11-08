# Position of the end of the supersonic beam on the axis versus time (hydrodynamic jet)
import glob, io, re, sys
import numpy as np
import matplotlib.pyplot as plt

tmax, yc, gamma = 2.4, 4.0, 1.4
vj = 5 * np.sqrt(gamma / 0.1)                 # jet speed: Mach 5, eta = 0.1
v1d = vj * np.sqrt(0.1) / (1 + np.sqrt(0.1))  # 1D momentum-balance estimate

def load(f):
    with open(f) as fh:
        txt = re.sub(r'(\d)([+-]\d{3})\b', r'\1E\2', fh.read())
    return np.loadtxt(io.StringIO(txt), unpack=True)

files = sorted(glob.glob('jet_*.txt'))
if not files:
    sys.exit('No jet_*.txt file found: run jet2d first')
t, xh = [], []
for n, f in enumerate(files, start=1):
    x, y, rho, vx, vy, p = load(f)
    ux, uy = np.unique(x), np.unique(y)
    VX = vx.reshape(len(ux), len(uy))
    j = np.argmin(abs(uy - yc))
    idx = np.where(VX[:, j] > 0.5 * vj)[0]     # end of the beam on the axis
    xh.append(ux[idx.max()] if len(idx) else 0.0)
    t.append(n * tmax / len(files))

t, xh = np.array(t), np.array(xh)
m = t > 1.0
slope = np.polyfit(t[m], xh[m], 1)[0]
print(f'advance speed of the beam end (t > 1): {slope:.2f}   (1D estimate: {v1d:.2f})')

fig, a = plt.subplots(figsize=(6.4, 4.6))
a.plot(t, xh, 'o-', label='simulation: beam end on the axis')
a.plot(t, v1d * t, '--', label=f'1D estimate  ({v1d:.2f} t)')
a.set_xlabel('t'); a.set_ylabel('position on the axis'); a.grid(alpha=0.25); a.legend(frameon=False)
fig.tight_layout()
fig.savefig('jet_hydro_head.png', dpi=200)
plt.show()
