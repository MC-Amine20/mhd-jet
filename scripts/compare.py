import sys, glob
import numpy as np
import matplotlib.pyplot as plt

# Usage : python compare.py 0 1.0 10      (une valeur de beta par argument ; 0 = sans champ)
tags = sys.argv[1:] if len(sys.argv) > 1 else ['0', '1.0']
tmax, yc, gam = 2.4, 4.0, 5.0/3.0
vj = 5.0 * np.sqrt(gam / 0.1)

def load_run(tag):
    files = sorted(glob.glob(f'jetmhd_{tag}_*.txt'))
    if not files:
        sys.exit(f"Aucun fichier jetmhd_{tag}_*.txt : lance d'abord  .\\jetmhd {tag}")
    return files

fig, ax = plt.subplots(len(tags), 1, figsize=(8, 3.3 * len(tags)), squeeze=False)
plt.figure(2, figsize=(7, 5))
for k, tag in enumerate(tags):
    files = load_run(tag)
    T, XB = [], []
    for n, f in enumerate(files, 1):
        x, y, rho, vx, vy, p, bx, by = np.loadtxt(f, unpack=True)
        ux, uy = np.unique(x), np.unique(y)
        V = vx.reshape(len(ux), len(uy))
        j = np.argmin(abs(uy - yc))
        idx = np.where(V[:, j] > 0.5 * vj)[0]
        XB.append(ux[idx.max()] if len(idx) else 0.0)
        T.append(n * tmax / len(files))
    # carte de densite du dernier snapshot
    R = rho.reshape(len(ux), len(uy)).T
    im = ax[k, 0].imshow(np.log10(R), origin='lower', cmap='inferno',
                         extent=[ux.min(), ux.max(), uy.min(), uy.max()])
    ax[k, 0].set_title(f'log10 densité, t = {T[-1]:.1f}, beta = {tag}' + ('  (sans champ)' if tag in ('0', '0.0') else ''))
    ax[k, 0].set_ylabel('y')
    fig.colorbar(im, ax=ax[k, 0])
    plt.figure(2)
    plt.plot(T, XB, 'o-', ms=4, label=f'beta = {tag}')
ax[-1, 0].set_xlabel('x')
fig.tight_layout()
fig.savefig('compare_densite.png', dpi=130)
plt.figure(2)
plt.xlabel('t'); plt.ylabel('fin du faisceau sur l\'axe (vx > 0.5 v_jet)')
plt.legend(); plt.grid(alpha=0.3)
plt.savefig('compare_tete.png', dpi=130)
plt.show()
