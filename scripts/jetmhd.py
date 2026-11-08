import sys, glob
import numpy as np
import matplotlib.pyplot as plt

tag = sys.argv[1] if len(sys.argv) > 1 else '1.0'
num = int(sys.argv[2]) if len(sys.argv) > 2 else None
files = sorted(glob.glob(f'jetmhd_{tag}_*.txt'))
if not files:
    sys.exit(f"Aucun fichier jetmhd_{tag}_*.txt : lance d'abord  .\\jetmhd {tag}")
f = files[num - 1] if num else files[-1]

x, y, rho, vx, vy, p, bx, by = np.loadtxt(f, unpack=True)
ux, uy = np.unique(x), np.unique(y)
nx, ny = len(ux), len(uy)
ext = [x.min(), x.max(), y.min(), y.max()]
bmag = np.sqrt(bx**2 + by**2)

fig, ax = plt.subplots(3, 1, figsize=(9, 11))
for a, q, name in zip(ax, [np.log10(rho), np.log10(p), bmag],
                      ['log10 densité', 'log10 pression', '|B| (plan)']):
    im = a.imshow(q.reshape(nx, ny).T, origin='lower', extent=ext, cmap='inferno')
    a.set_title(f'{name}  ({f})')
    a.set_xlabel('x'); a.set_ylabel('y')
    plt.colorbar(im, ax=a)
plt.tight_layout()
plt.savefig(f'jetmhd_{tag}.png', dpi=130)
plt.show()
