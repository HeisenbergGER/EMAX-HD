from fdtd import *
from grid import *
from constants import *
from pml import *
import numpy as np
import cython
import matplotlib.pyplot as plt
dt = 1/(200*np.sqrt(3)*c)
grid = Grid(L=1,du=0.01,t_end=1e-7)
grid.precision = np.float64
pml = Pml(grid,noc=10)
pml.configure()
print(np.max(pml.PML_z))
print(np.max(pml.PML_y))
print(np.max(pml.PML_x))
grid.load_phantom("phantom.png")
#plt.imshow(grid.inv_permittivity[:,:,0,0])
#plt.imshow(pml.PML_y[:,:,0])
#plt.colorbar()
plt.show()
print("Grid size: {:.0f}".format(grid.grid_size))
print("Courant number: {:.3f}".format(grid.sc))
sim = Fdtd(grid,pml)
sim.run()
