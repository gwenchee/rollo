import os, sys 
import shutil
import openmc
import numpy as np
from numpy import sin, cos, tan, pi
sys.path.insert(1, '../scripts/')
from constants import *

def run_openmc(name,total_pf,poly_coeff):

    path = 'openmc_' + name
    os.mkdir(path)
    os.chdir(path)

    mats.export_to_xml()

    # shape of reactor
    outer = openmc.Cell(fill=graphite)
    z_height = 25
    outer.region = +openmc.XPlane(x0=0,boundary_type='reflective') & -openmc.XPlane(x0=0.4,boundary_type='reflective') \
                & +openmc.YPlane(y0=0,boundary_type='reflective') & -openmc.YPlane(y0=0.4,boundary_type='reflective') \
                & +openmc.ZPlane(z0=0,boundary_type='reflective') & -openmc.ZPlane(z0=0.2+z_height,boundary_type='reflective') 

    # triso PF distribution
    total_core_vol = 1
    vol_triso = 4/3*pi*T_r5**3
    total_trisos = round(total_pf*total_core_vol/vol_triso)
    dz = 10
    z_vals = np.arange(1,dz+1)
    z = poly_coeff[0]*z_vals**3 + poly_coeff[1]*z_vals**2 + poly_coeff[2]*z_vals + poly_coeff[3]
    z_trisos = z/(sum(z))*total_trisos
    pf_z = z_trisos*vol_triso/(total_core_vol/dz)
    z_thick = z_height/dz

    # core
    all_prism_univ = openmc.Universe()
    small_prism = +openmc.XPlane(x0=0.1) & -openmc.XPlane(x0=0.3) \
            & +openmc.YPlane(y0=0.1) & -openmc.YPlane(y0=0.3) \
            & +openmc.ZPlane(z0=0.1) & -openmc.ZPlane(z0=0.1+z_thick) 
    all_prism_regions = small_prism
    for i in range(dz): 
        prism_region = small_prism 
        #prism_cell = openmc.Cell(fill=lm_graphite,)
        #print('PACKFRAC',pf_z[i])
        try: 
            centers = openmc.model.pack_spheres(radius=T_r5, region=prism_region, pf=pf_z[i])
        except ZeroDivisionError: 
            centers = []
            #print('here')
        trisos = [openmc.model.TRISO(T_r5, triso_univ, c) for c in centers]
        #print(pf_z[i],len(trisos)*4/3*pi*T_r5**3/(z_thick*0.2*0.2))
        prism_cell = openmc.Cell(region=prism_region)
        lower_left, upper_right = prism_cell.region.bounding_box
        shape = (1,1,1)
        pitch = (upper_right - lower_left)/shape
        lattice = openmc.model.create_triso_lattice(trisos, lower_left, pitch, shape, lm_graphite)
        prism_cell.fill = lattice
        prism_univ = openmc.Universe(cells=(prism_cell,))
        z_trans = i*z_thick
        prism_region_new = prism_region.translate((0,0,z_trans))
        prism_cell_new = openmc.Cell(fill=prism_univ,region=prism_region_new)
        prism_cell_new.translation = (0,0,z_trans)
        all_prism_univ.add_cell(prism_cell_new)
        all_prism_regions |= prism_region_new 
    prism_areas = openmc.Cell(fill=all_prism_univ, region=all_prism_regions)
    outer.region &= ~all_prism_regions
    #print('out')
    # geometry
    univ = openmc.Universe(cells=[outer,prism_areas])
    geom = openmc.Geometry(univ)
    geom.export_to_xml()

    # settings 
    point = openmc.stats.Point((0.2, 0.2, 12.5))
    src = openmc.Source(space=point)
    settings = openmc.Settings()
    settings.source = src
    settings.batches = 20
    settings.inactive = 5
    settings.particles = 500
    settings.temperature = {'multipole':True,'method':'interpolation'}
    settings.export_to_xml()

    openmc.run()

    sp = openmc.StatePoint('statepoint.20.h5')
    keff = sp.k_combined.nominal_value
    os.chdir('../')
    #shutil.rmtree(path)
    return keff, total_pf  

#run_openmc(1,0.01,[1,1,1,1])

# openmc input + output
"""
def evalOpenMC(ind):
    # return keff, total PF  
    #keff, total_pf = run_openmc(name=str(ind.gen)+'_'+str(ind.num), 
    #                            total_pf=ind[0],
    #                            poly_coeff=[ind[1],ind[2],ind[3],ind[4]])
    return (random.randint(-2,3),ind[0])"""