from constants import *

class Pml:
    def __init__(self,grid,noc=10):
        print(mu0,e0)
        self.grid = grid
        self.noc = noc # number of PML cells
        self.kappa_max = 10 # max conductivity
        self.sigma_max = (self.kappa_max*e0*mu0)**0.5/(2*grid.du) # max conductivity in S/m
        self.alpha_max = self.sigma_max/(e0)
        self.thickness = self.noc*grid.du # thickness of PML layer in m
        # init arrays
        self.sigma_x = np.zeros((grid.grid_size,grid.grid_size,grid.grid_size), dtype=grid.precision)
        self.sigma_y = np.zeros((grid.grid_size,grid.grid_size,grid.grid_size), dtype=grid.precision)
        self.sigma_z = np.zeros((grid.grid_size,grid.grid_size,grid.grid_size), dtype=grid.precision)
        self.alpha_x = np.zeros((grid.grid_size,grid.grid_size,grid.grid_size), dtype=grid.precision)
        self.alpha_y = np.zeros((grid.grid_size,grid.grid_size,grid.grid_size), dtype=grid.precision)
        self.alpha_z = np.zeros((grid.grid_size,grid.grid_size,grid.grid_size), dtype=grid.precision)
        
        self.PML_x = np.zeros((grid.grid_size,grid.grid_size,grid.grid_size),dtype=grid.precision)
        self.PML_y = np.zeros((grid.grid_size,grid.grid_size,grid.grid_size),dtype=grid.precision)
        self.PML_z = np.zeros((grid.grid_size,grid.grid_size,grid.grid_size),dtype=grid.precision)

    def configure(self):
        for i in range(self.noc):
            x = (i+0.5)*self.grid.du
            self.sigma_x[i,:,:] = self.sigma_max*((self.thickness-x)/self.thickness)**2
            self.alpha_x[i,:,:] = self.alpha_max*((self.thickness-x)/self.thickness)**3
            self.PML_x[i,:,:] = np.exp(-(self.sigma_x[i,:,:]/e0 + self.alpha_x[i,:,:])*self.grid.dt/e0)
            self.sigma_x[self.grid.grid_size-i-1,:,:] = self.sigma_max*((self.thickness-x)/self.thickness)**2
            self.alpha_x[self.grid.grid_size-i-1,:,:] = self.alpha_max*((self.thickness-x)/self.thickness)**3
            self.PML_x[self.grid.grid_size-i-1,:,:] = np.exp(-(self.sigma_x[self.grid.grid_size-i-1,:,:]/e0 + self.alpha_x[self.grid.grid_size-i-1,:,:])*self.grid.dt/e0)
            y = (i+0.5)*self.grid.du
            self.sigma_y[:,i,:] = self.sigma_max*((self.thickness-y)/self.thickness)**2
            self.alpha_y[:,i,:] = self.alpha_max*((self.thickness-y)/self.thickness)**3
            self.PML_y[:,i,:] = np.exp(-(self.sigma_y[:,i,:]/e0 + self.alpha_y[:,i,:])*self.grid.dt/e0)
            self.sigma_y[:,self.grid.grid_size-i-1,:] = self.sigma_max*((self.thickness-y)/self.thickness)**2
            self.alpha_y[:,self.grid.grid_size-i-1,:] = self.alpha_max*((self.thickness-y)/self.thickness)**3
            self.PML_y[:,self.grid.grid_size-i-1,:] = np.exp(-(self.sigma_y[:,self.grid.grid_size-i-1,:]/e0 + self.alpha_y[:,self.grid.grid_size-i-1,:])*self.grid.dt/e0)
            z = (i+0.5)*self.grid.du
            self.sigma_z[:,:,i] = self.sigma_max*((self.thickness-z)/self.thickness)**2
            self.alpha_y[:,:,i] = self.alpha_max*((self.thickness-z)/self.thickness)**3
            self.PML_z[:,:,i] = np.exp(-(self.sigma_z[:,:,i]/e0 + self.alpha_z[:,:,i])*self.grid.dt/e0)
            self.sigma_z[:,:,self.grid.grid_size-i-1] = self.sigma_max*((self.thickness-z)/self.thickness)**2
            self.alpha_z[:,:,self.grid.grid_size-i-1] = self.alpha_max*((self.thickness-z)/self.thickness)**3
            self.PML_z[:,:,self.grid.grid_size-i-1] = np.exp(-(self.sigma_z[:,:,self.grid.grid_size-i-1]/e0 + self.alpha_z[:,:,self.grid.grid_size-i-1])*self.grid.dt/e0)





