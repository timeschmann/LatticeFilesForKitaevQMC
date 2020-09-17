# -------------------------------------------------------------------
# JULIA TOOL FOR EXPORTING INTERACTION MATRICES AS C++-HEADER FILES
# AS USED IN KITAEV_QMC AND KITAEV_QMC_KPM
# Author: Tim Eschmann
# Date: 2017-01-26
# Modified for new version of LatticePhysics: 2019-01-31
#
# This tool requires the Julia package "LatticePhysics.jl"
# (Copyright (c) 2017: Jan Attig) which is available under
# https://github.com/janattig/LatticePhysics.jl
# -------------------------------------------------------------------

include("matrix_creator.jl")

# use the modules
using LatticePhysics
using MatrixCreator

#####################################################################
# 1. Define System: 

# Lattice:
lat = "(10,3)a"

# Size:
Lx = 8
Ly = 8
Lz = 8

# Periodicity:
per_dir = "periodic"
whichone = "Any"

# Gauge field:
gauge = "All"

# Loop length:
l = 10

# Coupling Constants
Jx = 1/3.
Jy = 1/3.
Jz = 1/3.

#####################################################################
# 2. Choose unit cell:

if (lat == "(6,3)a")
	uc = getUnitcellHoneycomb(4)
	filename = "honeycomb_$(Lx)_$(Ly)_Jz_$(Jz)_disorder.hpp"
elseif (lat == "(8,3)a")
	uc = getUnitcell_8_3_a(4)
	filename = "8_3_a_$(Lx)x$(Ly)x$(Lz)_$(per_dir)_$(l)flux_Jz_$(Jz).hpp"
elseif (lat == "(8,3)b")
	uc = getUnitcell_8_3_b(4)
	filename = "8_3_b_$(Lx)x$(Ly)x$(Lz)_$(per_dir)_$(l)flux_Jz_$(Jz).hpp"
elseif (lat == "(8,3)c")
	uc = getUnitcell_8_3_c(4)
	filename = "8_3_c_$(Lx)x$(Ly)x$(Lz)_$(per_dir)_$(l)flux_Jz_$(Jz).hpp"
elseif (lat == "(8,3)n")
	uc = getUnitcell_8_3_n(4)
	filename = "8_3_n_$(Lx)x$(Ly)x$(Lz)_$(per_dir)_$(l)flux_Jz_$(Jz).hpp"
elseif (lat == "(9,3)a")
	uc = getUnitcell_9_3_a(4)
	filename = "9_3_a_$(Lx)x$(Ly)x$(Lz)_$(per_dir)_$(l)flux_Jz_$(Jz).hpp"
elseif (lat == "(10,3)a")
	uc = getUnitcell_10_3_a(4)
	filename = "10_3_a_$(Lx)x$(Ly)x$(Lz)_$(per_dir)_$(l)flux_Jz_$(Jz).hpp"
elseif (lat == "(10,3)b")
	uc = getUnitcell_10_3_b(4)
	filename = "10_3_b_$(Lx)x$(Ly)x$(Lz)_$(per_dir)_$(l)flux_Jz_$(Jz).hpp"
elseif (lat == "(10,3)c")
	uc = getUnitcell_10_3_c(4)
	filename = "10_3_c_$(Lx)x$(Ly)x$(Lz)_$(per_dir)_$(l)flux_Jz_$(Jz).hpp"
elseif (lat == "(10,3)d")
	uc = getUnitCell_10_3_d(4)
	filename = "10_3_d_$(Lx)x$(Ly)x$(Lz)_$(per_dir)_$(l)flux_Jz_$(Jz).hpp"
else
	println("ERROR: Choose unitcell")
end

#####################################################################
# 3. Generate Lattice:

# Boundary conditions:
if (lat == "(6,3)a" && per_dir == "periodic")
	lt = getLatticePeriodic2D(uc, [Lx, Ly])
elseif (lat == "(6,3)a" && per_dir == "open")
	lt = getLatticeOpen2D(uc, [Lx, Ly])
elseif (per_dir == "open")
	lt = getLatticeOpen3D(uc, [Lx, Ly, Lz])
#elseif (per_dir == "poo")
#	lt = getStripLattice3D(uc, -Lx, Ly, Lz, whichone)
#elseif (per_dir == "opo")
#	lt = getStripLattice3D(uc, Lx, -Ly, Lz, whichone)
#elseif (per_dir == "oop")
#	lt = getStripLattice3D(uc, Lx, Ly, -Lz, whichone)
#elseif (per_dir == "ppo")
#	lt = getStripLattice3D(uc, -Lx, -Ly, Lz, whichone)
#elseif (per_dir == "pop")
#	lt = getStripLattice3D(uc, -Lx, Ly, -Lz, whichone)
#elseif (per_dir == "opp")
#	lt = getStripLattice3D(uc, Lx, -Ly, -Lz, whichone)
elseif (per_dir == "periodic")
	lt = getLatticePeriodic3D(uc, [Lx,Ly,Lz])
end

#####################################################################
# 4. Create Kitaev interaction matrix A from lattice:

matrix = getInteractionMatrixCouplings(lt, Jx, Jy, Jz, gauge)

# Disordered version:
#p = 0.25
#dJ = 0.8
#matrix = getInteractionMatrixDisorder(lt, Jx, Jy, Jz, dJ, p)

#####################################################################
# 5. Find closed (even-length) loops inside the lattice:
netto_plaquettes = getPlaquettesOfLattice(lt,l)
    
#####################################################################
# 6. Write Header-File for C++ Monte Carlo Program:

writehpp(lt, netto_plaquettes, l, matrix, filename)

# Include only zz-bonds
writehpp_zz(lt, netto_plaquettes, l, matrix, filename)
#####################################################################

exit()

