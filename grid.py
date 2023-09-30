from functions import *
from constants import *
from PIL import Image
class Grid:
    def __init__(self, L=1, du=0.005, t_end=1, precision=np.float64):
        # define grid parameters
        self.precision = precision
        self.L = L
        self.du = du
        self.t_end = t_end
        self.dt = self.du/(np.sqrt(3)*c)
        self.grid_size = int(self.L/self.du)
        self.time_steps = int(self.t_end/self.dt)
        self.d = 3 # dimension of the grid
        self.sc = c*self.dt/self.du # sc should be smaller than 1/sqrt(dim)
        # initialize fields
        self.E = np.zeros((self.grid_size,self.grid_size,self.grid_size,self.d),dtype=self.precision)
        self.H = np.zeros((self.grid_size,self.grid_size,self.grid_size,self.d),dtype=self.precision)
        self.inv_permittivity = 1
        self.inv_permeability = 1

    def load_phantom(self,filepath):
        image = Image.open(filepath)
        if ((image.size[0] != self.grid_size) or (image.size[1] != self.grid_size)):
            print("Phantom does not match grid resolution. Resizing...")
            image = image.resize((self.grid_size,self.grid_size), Image.Resampling.LANCZOS)
        frame = np.array(image)
        self.inv_permittivity = np.ones((self.grid_size,self.grid_size,self.grid_size,3))
        fat_mask = np.where(frame[:,:,0]==255)
        bone_mask = np.where(frame[:,:,1] == 255)
        self.inv_permittivity[fat_mask[0],fat_mask[1],0,:] = 1.27e1
        self.inv_permittivity[bone_mask[0],bone_mask[1],0,:] = 2.76e1

    def update_E(self):
        self.E += self.sc * self.inv_permittivity**(-1) * curl_H(self.H,self.precision)

    def update_H(self):
        self.H -= self.sc * self.inv_permeability**(-1) * curl_E(self.E,self.precision)

    def do_step(self):
        self.update_E()
        self.update_H()

    def do_pml_step(self,pml):
        self.E[:pml.noc,:,:,0] *= pml.PML_x[:pml.noc,:,:]
        self.E[self.grid_size - pml.noc:,:,:,0] *= pml.PML_x[self.grid_size - pml.noc:,:,:]
        self.E[:,:pml.noc,:,1] *= pml.PML_y[:,:pml.noc,:] 
        self.E[:,self.grid_size - pml.noc:,:,1] *= pml.PML_y[:,self.grid_size - pml.noc:,:]
        self.E[:,:,:pml.noc,2] *= pml.PML_z[:,:,:pml.noc] 
        self.E[:,:,self.grid_size - pml.noc:,2] *= pml.PML_z[:,:,self.grid_size - pml.noc:]
        self.update_E()
        self.H[:pml.noc,:,:,0] *= pml.PML_x[:pml.noc,:,:] 
        self.H[self.grid_size - pml.noc:,:,:,0] *= pml.PML_x[self.grid_size - pml.noc:,:,:]
        self.H[:,:pml.noc,:,1] *= pml.PML_y[:,:pml.noc,:] 
        self.H[:,self.grid_size - pml.noc:,:,1] *= pml.PML_y[:,self.grid_size - pml.noc:,:]
        self.H[:,:,:pml.noc,2] *= pml.PML_z[:,:,:pml.noc] 
        self.H[:,:,self.grid_size - pml.noc:,2] *= pml.PML_z[:,:,self.grid_size - pml.noc:]
        self.update_H()



