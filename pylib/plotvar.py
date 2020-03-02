#=================================================================#
#================= Plotting UEDGE 2D data on mesh  ===============#
#
#-Usage:
# execfile("../pylib/plotvar.py") 
# plotvar(var)
# plotvar(var, iso=True)
#
# Input arguments:
#   var(0:com.nx+2, 0:com.ny+2)
#
#-Optional arguments:
#   iso (True/False) - True for equal aspect ratio
#
#-Expects imported modules:
# import matplotlib.pyplot as plt; import numpy as np
#
# Usage example:
# fig1 = plt.figure(); plotvar(com.rm[:,:,0])
#
# First coding: MVU, 22-jul-17
# log,axis,cmap keywords: BD, 28-Feb-20
#=================================================================#

from uedge import com

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

def plotvar(var, iso=False, title="UEDGE data", show=True, axis=None, cmap=None, log=False):
    """Plot UEDGE 2D data on mesh
    
    Input arguments:
      var(0:com.nx+2, 0:com.ny+2)

    Optional arguments:
      iso (True/False)  - True for equal aspect ratio
      title             - The figure title
      show (True/False) - Call plt.show() at the end?
      axis              - Plot on given axis.
      cmap              - Specify a colormap
      log (True/False)  - If true, use a logarithmic scale

    Returns
      axis     The axis plotted on

    Usage example:
       plotvar(com.rm[:,:,0])

       plotvar(bbb.te/bbb.ev, log=True, show=False, iso=True)
       plt.savefig("logte.pdf")
       plt.show()
    """
    patches = []

    for iy in np.arange(0,com.ny+2):
        for ix in np.arange(0,com.nx+2):
            rcol=com.rm[ix,iy,[1,2,4,3]]
            zcol=com.zm[ix,iy,[1,2,4,3]]
            rcol.shape=(4,1)
            zcol.shape=(4,1)
            polygon = Polygon(np.column_stack((rcol,zcol)), True)
            patches.append(polygon)

    #-is there a better way to cast input data into 2D array?
    vals=np.zeros((com.nx+2)*(com.ny+2))

    for iy in np.arange(0,com.ny+2):
        for ix in np.arange(0,com.nx+2):
            k=ix+(com.nx+2)*iy
            vals[k] = var[ix,iy]

    if log:
        p = PatchCollection(patches, norm=matplotlib.colors.LogNorm())
    else:
        p = PatchCollection(patches)
    p.set_array(np.array(vals))

    if axis is None:
        fig,axis = plt.subplots(1)
        
    axis.add_collection(p)
    axis.autoscale_view()
    plt.colorbar(p)
    
    axis.set_title(title)
    axis.set_xlabel('R [m]')
    axis.set_ylabel('Z [m]')
    axis.grid(True)

    if iso:
        axis.set_aspect('equal')

    if show:
        plt.show()
    return axis

#=================================================================#
