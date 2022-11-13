########################################################################                                                                          
import steps.interface
                                                                                                                                                  
from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *
from constants import *

import sys

SEED = int(sys.argv[1])
MESH_FILES = "local_mesh/local_16cores"


########################### BIOCHEMICAL MODEL ###############################

mdl = Model()
r = ReactionManager()
with mdl:
    # Species
    Ca, Pump, CaPump, iCBsf, iCBsCa, iCBCaf, iCBCaCa, CBsf, CBsCa, CBCaf, CBCaCa, PV, PVMg, PVCa, Mg = Species.Create()

    # Vol/surface systems
    vsys = VolumeSystem.Create()
    ssys = SurfaceSystem.Create()

    with vsys:
        diff_Ca =     Diffusion.Create(Ca, DCST)
        diff_CBsf =   Diffusion.Create(CBsf, DCB)
        diff_CBsCa =  Diffusion.Create(CBsCa, DCB)
        diff_CBCaf =  Diffusion.Create(CBCaf, DCB)
        diff_CBCaCa = Diffusion.Create(CBCaCa, DCB)
        diff_PV =     Diffusion.Create(PV, DPV)
        diff_PVCa =   Diffusion.Create(PVCa, DPV)
        diff_PVMg =   Diffusion.Create(PVMg, DPV)

        # iCBsf fast and slow
        (iCBsf + Ca <r[1]> iCBsCa) + Ca <r[2]> iCBCaCa
        (iCBsf + Ca <r[3]> iCBCaf) + Ca <r[4]> iCBCaCa
        r[1].K = iCBsf1_f_kcst, iCBsf1_b_kcst
        r[2].K = iCBsCa_f_kcst, iCBsCa_b_kcst
        r[3].K = iCBsf2_f_kcst, iCBsf2_b_kcst
        r[4].K = iCBCaf_f_kcst, iCBCaf_b_kcst

        # CBsf fast and slow
        (CBsf + Ca <r[1]> CBsCa) + Ca <r[2]> CBCaCa
        (CBsf + Ca <r[3]> CBCaf) + Ca <r[4]> CBCaCa
        r[1].K = CBsf1_f_kcst, CBsf1_b_kcst
        r[2].K = CBsCa_f_kcst, CBsCa_b_kcst
        r[3].K = CBsf2_f_kcst, CBsf2_b_kcst
        r[4].K = CBCaf_f_kcst, CBCaf_b_kcst

        # PVCa
        PV + Ca <r[1]> PVCa
        r[1].K = PVca_f_kcst, PVca_b_kcst

        # PVMg
        PV + Mg <r[1]> PVMg
        r[1].K = PVmg_f_kcst, PVmg_b_kcst

        # Ca Influx converted from P Type current
        None >r['CaInflux']> Ca
        r['CaInflux'].K = 0.0

    with ssys:
        # Ca Pump
        Pump.s + Ca.i <r[1]> CaPump.s >r[2]> Pump.s
        r[1].K = P_f_kcst, P_b_kcst
        r[2].K = P_k_kcst
    
########### MESH & COMPARTMENTALIZATION #################

mesh = DistMesh(MESH_FILES, 1e-6)

with mesh:
    __MESH__ = Compartment.Create(vsys)
    __MESH_BOUNDARY__ = Patch.Create(__MESH__, None, ssys)
    memb = Membrane.Create([__MESH_BOUNDARY__])
    __MESH__.Conductivity = 1 / Ra #may need change for Purkinje cell

# # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

rng = RNG('mt19937', 512, SEED)

sim = Simulation('DistTetOpSplit', mdl, mesh, rng, searchMethod=NextEventSearchMethod.GIBSON_BRUCK)

print("Rank %i: Job Done!" % (MPI.rank))
