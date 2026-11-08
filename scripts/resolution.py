import sys, glob, os
import numpy as np
import matplotlib.pyplot as plt

# Usage : python resolution.py            (compare le dossier courant et le dossier res360)
#         python resolution.py . res360   (ou toute liste de dossiers)
folders = sys.argv[1:] if len(sys.argv) > 1 else ['.', 'res360']
tags = ['0', '1.0']
tmax, yc, gam = 2.4, 4.0, 5.0/3.0
vj = 5.0 * np.sqrt(gam / 0.1)

def analyse(folder, tag):
    files = sorted(glob.glob(os.path.join(folder, f'jetmhd_{tag}_*.txt')))
    if not files:
        return None
    T, XB, XF = [], [], []
    for n, f in enumerate(files, 1):
        x, y, rho, vx, vy, p, bx, by = np.loadtxt(f, unpack=True)
        ux, uy = np.unique(x), np.unique(y)
        V = vx.reshape(len(ux), len(uy)); R = rho.reshape(len(ux), len(uy))
        j = np.argmin(abs(uy - yc))
        idx = np.where(V[:, j] > 0.5 * vj)[0]
        XB.append(ux[idx.max()] if len(idx) else 0.0)         # fin du faisceau
        idx = np.where(R[:, j] > 1.3)[0]
        XF.append(ux[idx.min()] if len(idx) else np.nan)      # front du gaz choque
        T.append(n * tmax / len(files))
    i6 = np.argmin(abs(ux - 6.0)); col = R[i6, :]; w = np.where(col < 0.5)[0]
    width = (uy[w.max()] - uy[w.min()]) if len(w) else 0.0    # cocon a x=6, t final
    return np.array(T), np.array(XB), np.array(XF), width, len(ux), len(uy)

fig, ax = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
print(f'{"beta":5s} {"grille":9s} {"fin faisceau (t=2.4)":22s} {"front choque (t=2.0)":22s} {"cocon x=6"}')
for k, tag in enumerate(tags):
    vus = []
    for folder in folders:
        res = analyse(folder, tag)
        if res is None:
            print(f'{tag:5s} (aucun fichier dans le dossier {folder!r})')
            continue
        T, XB, XF, width, nx, ny = res
        if (nx, ny) in vus:
            print(f'{tag:5s} ATTENTION : le dossier {folder!r} a deja la grille {nx}x{ny}, courbe ignoree')
            continue
        vus.append((nx, ny))
        i20 = np.argmin(abs(T - 2.0))
        print(f'{tag:5s} {nx}x{ny:<5d} {XB[-1]:<22.2f} {XF[i20]:<22.2f} {width:.2f}   ({len(T)} snapshots, dossier {folder!r})')
        ax[k].plot(T, XB, 'o-', ms=4, label=f'{nx}x{ny}')
    ax[k].set_title(f'beta = {tag}' + (' (sans champ)' if tag == '0' else ''))
    ax[k].set_xlabel('t'); ax[k].grid(alpha=0.3)
    if vus: ax[k].legend()
ax[0].set_ylabel("fin du faisceau sur l'axe (vx > 0.5 v_jet)")
plt.tight_layout(); plt.savefig('resolution.png', dpi=130); plt.show()