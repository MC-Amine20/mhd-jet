# Plot the Brio-Wu shock tube produced by mhd.f90 (file mhd.txt: x, rho, vx, vy, vz, By, Bz, p)
import numpy as np
import matplotlib.pyplot as plt

x, rho, vx, vy, vz, by, bz, p = np.loadtxt('mhd.txt', unpack=True)

fig, ax = plt.subplots(2, 2, figsize=(9.4, 5.8))
for a, y, name in zip(ax.flat, [rho, vx, by, p],
                      ['Density  ρ', 'Velocity  v$_x$', 'Magnetic field  B$_y$', 'Pressure  p']):
    a.plot(x, y, '.', ms=3)
    a.set_title(name); a.set_xlabel('x'); a.grid(alpha=0.25)
fig.suptitle('Brio–Wu shock tube, t = 0.1, N = 800 cells')
fig.tight_layout()
fig.savefig('briowu.png', dpi=200)
plt.show()
