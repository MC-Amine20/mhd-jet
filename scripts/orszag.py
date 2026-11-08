import numpy as np
import matplotlib.pyplot as plt

x, y, rho, vx, vy, bx, by, p, psi = np.loadtxt('orszag.txt', unpack=True)
n = int(round(np.sqrt(len(x))))
R = rho.reshape(n, n)       # indices [j, i]
P = p.reshape(n, n)
print('rho min/max', rho.min(), rho.max(), ' p min/max', p.min(), p.max())

fig, ax = plt.subplots(1, 2, figsize=(11, 5))
for a, q, name in zip(ax, [R, P], ['densité', 'pression']):
    im = a.imshow(q, origin='lower', extent=[0, 1, 0, 1], cmap='viridis')
    a.set_title(f'{name} (t = 0.5)'); a.set_xlabel('x'); a.set_ylabel('y')
    plt.colorbar(im, ax=a, fraction=0.046)
plt.tight_layout()
plt.savefig('orszag.png', dpi=130)
