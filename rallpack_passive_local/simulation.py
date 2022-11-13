import steps.interface
from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *
from constants import *
from matplotlib import pyplot as plt
import numpy as np
import math
import sys


SEED = 1234
MESH_FILES = "local_mesh/local_16cores_v2_dist"


# # # # # # # # # # # # # # # # SIMULATION CONTROLS # # # # # # # # # # # # # #

# Sim end time (seconds)
SIM_END = 0.25

# The current injection in amps
Iinj = 0.1e-9

EF_DT = 1e-6
SAVE_DT = 5e-5  #

# # # # # # # # # # # # # # # # PARAMETERS # # # # # # # # # # # # # #

# Leak conductance
L_G = 0.25 #Siemens/m^2


# Leak reveral potential, V
leak_rev = -65.0e-3

# Ohm.m
Ra = 1.0


mesh=DistMesh(MESH_FILES,1e-6)


########################### BIOCHEMICAL MODEL ###############################
mdl= Model()
r=ReactionManager()
with mdl:
    vsys= VolumeSystem.Create()
    ssys=SurfaceSystem.Create()
    
    #Leak
    leaksus = SubUnitState.Create()
    Leak = Channel.Create([leaksus])
    


########### MESH & COMPARTMENTALIZATION & Data collection #################

with mesh:
    __MESH__ = Compartment.Create(vsys)
    __MESH_BOUNDARY__ = Patch.Create(__MESH__, None, ssys)
    memb = Membrane.Create([__MESH_BOUNDARY__],capacitance=0.01)
    __MESH__.Conductivity = 1 / Ra
    
    
    #Record Point .(for current injection)
    record_point1 = [
        [-0.875033654923e-6, 2.37318054762e-6, -1.48584303865e-6], #inject point
    ]
    record_tets1 = TetList(mesh.tets[point] for point in record_point1)
    record_tris1 = TriList([])
    for tet in record_tets1:
        tris = tet.faces & mesh.surface
        if len(tris) > 0:
            record_tris1.append(tris[0])
            
    injverts=record_tris1.verts

    #record point(for record the potential)
    record_point2 = [
        [-1.11999521668e-6, 3.0658738637e-6, -0.107680016989e-6] #potential record
    ]
    record_tets2 = TetList(mesh.tets[point] for point in record_point2)
   
#set single channel conductance
with ssys:
    L_G_tot= L_G * memb.Area
    g_leak_sc=L_G_tot/len(memb.tris)
    OC_L=OhmicCurr.Create(Leak[leaksus],g_leak_sc, leak_rev)


# # # # # # # # # # # SIMULATION  # # # # # # # # # # #
rng = RNG('mt19937', 512, SEED)

sim = Simulation('DistTetOpSplit', mdl, mesh, rng, searchMethod=NextEventSearchMethod.GIBSON_BRUCK)
print("Rank %i: Job Done!" % (MPI.rank))

#Data saving
rs=ResultSelector(sim)
Vrs=rs.TETS(record_tets2).V
sim.toSave(Vrs,dt=SAVE_DT)
######################################
sim.newRun()
# Inject channels
sim.TRIS(memb.tris).Leak[leaksus].Count = 1

sim.memb.Potential=-65e-3
sim.VERTS(injverts).IClamp=Iinj/len(injverts)
sim.EfieldDT = EF_DT
sim.run(SIM_END)

#Analatical solution



if MPI.rank ==0:
    plt.figure(figsize=(10,7))
    plt.plot(Vrs.time[0]*1e3,Vrs.data[0]*1e3)
    plt.xlabel('time[ms]')
    plt.ylabel('potential(mV)')
    plt.show()
   


