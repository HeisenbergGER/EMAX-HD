import numpy as np
import matplotlib.pyplot as plt
from numba import jit
@jit(nopython=True, parallel=True, cache=True, fastmath=True)
def curl_E(E,precision="float64"):
        curl_E = np.zeros(E.shape,dtype=precision)
        curl_E[:,:-1,:,0] += E[:,1:,:,2] - E[:,:-1,:,2]
        curl_E[:,:,:-1,0] -= E[:,:,1:,1] - E[:,:,:-1,1]

        curl_E[:,:,:-1,1] += E[:,:,1:,0] - E[:,:,:-1,0]
        curl_E[:-1,:,:,1] -= E[1:,:,:,2] - E[:-1,:,:,2]

        curl_E[:-1,:,:,2] += E[1:,:,:,1] - E[:-1,:,:,1]
        curl_E[:,:-1,:,2] -= E[:,1:,:,0] - E[:,:-1,:,0]
        
        return curl_E
@jit(nopython=True, parallel=True, cache=True, fastmath=True)
def curl_H(H,precision="float64"):
    curl_H = np.zeros(H.shape,dtype=precision)
    curl_H[:,1:,:,0] += H[:,1:,:,2] - H[:,:-1,:,2]
    curl_H[:,:,1:,0] -= H[:,:,1:,1] - H[:,:,:-1,1]

    curl_H[:,:,1:,1] += H[:,:,1:,0] - H[:,:,:-1,0]
    curl_H[1:,:,:,1] -= H[1:,:,:,2] - H[:-1,:,:,2]

    curl_H[1:,:,:,2] += H[1:,:,:,1] - H[:-1,:,:,1]
    curl_H[:,1:,:,2] -= H[:,1:,:,0] - H[:,:-1,:,0]
        
    return curl_H

@jit(nopython=True, parallel=True, cache=True, fastmath=True)
def compute_magnitude(x):
    return np.sqrt(x[:,:,0,0]**2 + x[:,:,0,1]**2 + x[:,:,0,2]**2)

def plot_E(grid,path):
    xx,yy = np.meshgrid(np.linspace(-grid.L/2,grid.L/2,int(grid.L/grid.du)),np.linspace(-grid.L/2,grid.L/2,int(grid.L/grid.du)))
    #E = np.sqrt(grid.E[:,:,0,1]**2)
    E = compute_magnitude(grid.E)
    plt.figure(figsize=(7,7))
    plt.contourf(xx,yy,E,np.linspace(0,100,100),cmap="plasma")
    plt.colorbar()
    ax = plt.gca()
    ax.set_aspect("equal", adjustable="box")
    plt.savefig(path,format="PNG")
    plt.close()

