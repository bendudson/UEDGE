# Convert a half-plane solution to a double null solution
#
from uedge import *

import numpy as np

def dnbot_to_dnull():
    """
    Convert a dnbot geometry to dnull

    This will modify, by reflecting:
    
    Plasma state:
      ni, ti, te, up, ng, phi

    Diffusion coefficients:
      kye_use, kyi_use, dif_use
   
    If bbb.isimpon == 2 then impurity fraction:
      afracs

    If bbb.iflcore == 1 then double:
      pcore, pcorei
    
    """

    # Check that the geometry is bottom half of a double null
    # Note that com.geometry can be an array with a padded string
    assert "dnbot" in str(com.geometry)
    
    # Grid sizes for half domain
    # Note: working array sizes are [half_nx + 2, half_ny + 2, com.nisp]
    half_nx = com.nx
    half_ny = com.ny

    half_inner = com.nxleg[0,0] + com.nxcore[0,0] + 1  # Includes boundary cell
    half_outer = com.nxleg[0,1] + com.nxcore[0,1] + 1
    
    full_inner = 2 * half_inner - 2
    full_outer = 2 * half_outer - 2
    
    # Note: Numbering in lower single null goes:
    # Inner target -> inboard midplane, then outboard midplane -> outer target
    #
    # In double null, ordered from lower inboard target to upper inboard, then
    # upper outboard to lower outboard
    
    def reflect_centred(lower_var):
        """
        Reflect a cell-centred variable (i.e. not up)
        Input is the half-domain array
        Returns the full-domain array
        """
        # Shape is mostly the same as the lower-half variable, except first index
        full_shape = (2 * half_nx,) + lower_var.shape[1:]
        
        full_var = np.zeros(full_shape)
        
        # Lower inner
        full_var[:half_inner, :] = lower_var[:half_inner, :]
        
        # Upper inner
        full_var[half_inner-1:2*half_inner-2, :] = lower_var[:half_inner-1, :][::-1,:]

        # Lower outer
        full_var[full_inner + half_outer-1:, :] = lower_var[half_inner+1:, :]

        # Upper outer
        full_var[full_inner:full_inner + half_outer-1, :] = lower_var[half_inner+1:, :][::-1,:]

        return full_var


    # Plasma state except parallel flow (ups)
    niss = reflect_centred(bbb.nis)  # Ion density
    tess = reflect_centred(bbb.tes)  # Electron temperature
    tiss = reflect_centred(bbb.tis)  # Ion temperature
    ngss = reflect_centred(bbb.ngs)  # Ion density
    phiss = reflect_centred(bbb.phis)  # Electrostatic potential

    # Diffusion coefficients
    kye_use = reflect_centred(bbb.kye_use)
    kyi_use = reflect_centred(bbb.kyi_use)
    dif_use = reflect_centred(bbb.dif_use)

    if bbb.isimpon == 2:
        # Impurity fraction
        afracs = reflect_centred(bbb.afracs)
    
    ##############################################
    # Staggered parallel flow velocity
    # The ups variable is staggered, so that ups[i,j,:] is ups[i+1/2,j,:]
    # i.e. at upper edge of cell
    # -> First index ups[0,:,:] is value on lower boundary
    # -> Last index ups[-1,:,:] is not used
    
    upss = np.zeros((2 * half_nx, half_ny + 2, com.nusp))  # Ion parallel flow
    
    # Lower inner
    upss[:half_inner, :, :] = bbb.ups[:half_inner, :, :com.nusp]

    # Upper inner
    upss[half_inner-2:2*half_inner-3, :, :] = -bbb.ups[:half_inner-1, :, :com.nusp][::-1,:,:]  # Note: Reversed sign
    
    # Lower outer
    upss[full_inner + half_outer-1:, :, :] = bbb.ups[half_inner+1:, :, :com.nusp]
    
    # Upper outer
    upss[full_inner-1:full_inner + half_outer-2, :, :] = -bbb.ups[half_inner+1:, :, :com.nusp][::-1,:,:] # Reversed sign
    
    # Correct flow velocity in upper inner boundary cell
    upss[2*half_inner-3, :, :] = upss[2*half_inner-4, :, :]  # Upper inner boundary cell
    
    # Set switches
    com.nxomit = 0
    bbb.isfixlb = 0
    com.isudsym = 0  # Solution not necessarily up-down symmetric
    
    # Reflect the lower grid to make an exactly up-down symmetric grid
    com.geometry="dnull"
    
    # Generate new grid, allocate arrays
    bbb.restart=1
    bbb.newgeo=1
    bbb.icntnunk=0
    bbb.ftol=1e10; bbb.dtreal = 1e-6; bbb.exmain()
    bbb.newgeo = 0
    bbb.ftol=1e-8;
    
    # Replace state variables
    bbb.nis = niss
    bbb.ups = upss
    bbb.tes = tess
    bbb.tis = tiss
    bbb.ngs = ngss
    bbb.phis = phiss

    # Diffusion coefficients
    bbb.kye_use = kye_use
    bbb.kyi_use = kyi_use
    bbb.dif_use = dif_use

    if bbb.isimpon == 2:
        # Fixed impurity fraction
        bbb.afracs = afracs
    
    # If core power is fixed, need to double
    if bbb.iflcore == 1:
        bbb.pcoree *= 2
        bbb.pcorei *= 2
