from functions import *
from constants import *
class Fdtd:
    def __init__(self, grid, pml=None):
        self.grid = grid
        self.pml = pml

    def run(self):
        for time_step in range(self.grid.time_steps):
            print("\rtime step {:.0f}/{:.0f}".format(time_step+1,self.grid.time_steps),end='',flush=True)
            # apply plane wave source
            self.grid.E[int(self.grid.grid_size/2 - 40),int(self.grid.grid_size/2),0,1] = np.sin(2*np.pi*time_step/10)
            #self.grid.H[int(self.grid.grid_size/2),int(self.grid.grid_size/2),0,0] = np.cos(2*np.pi*time_step/50)
            # update fields
            if self.pml == None:
                self.grid.do_step()
            if self.pml != None:
                self.grid.do_pml_step(self.pml)
            # plot snapshot at current time step
            plot_E(self.grid,"output/output-"+str(time_step+1)+".png")
