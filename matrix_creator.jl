# -------------------------------------------------------------------
# JULIA MODULE FOR EXPORTING INTERACTION MATRICES AS C++-HEADER FILES 
# AS USED IN KITAEV_QMC AND KITAEV_QMC_KPM
# Author: Tim Eschmann
# Date: 2017-01-26
# Modified for new version of LatticePhysics: 2019-01-31
#
# This module requires the Julia package "LatticePhysics.jl"
# (Copyright (c) 2017: Jan Attig) which is available under
# https://github.com/janattig/LatticePhysics.jl
# -------------------------------------------------------------------

module MatrixCreator

using JLD
using LatticePhysics


#####################################################################
# Generate Hamiltonian for all-bond sampling in QMC
function writehpp(lat, plq, plaq_length, matrix, filename)
	N1 = size(matrix,1)
	N2 = size(matrix,2)	
	epsil = 0.0001

	file = open(filename, "w")
	
	# write the headerstring
	write(file, getHeaderString())
	write(file, declare_def_matrix(matrix))
	
	for j1 in (1:N1)
		for j2 in (j1:N2)
			if abs(matrix[j1,j2]) > epsil
				#if (lat.positions_indices[j1] == 2.0)
				#	write(file, giveMinusEntries(matrix, j1, j2))
				#else
				write(file, giveEntries(matrix, j1, j2))
				#end			
			end		
		end	
	end

	write(file, endofDec()) 
	write(file, declare_coord())

	for j1 in (1:N1)
	    for j2 in (j1:N2)
			if matrix[j1,j2] < 0
			   write(file, giveCoordinates(N1, j1, j2))			
		    end		
		end	
	end

	write(file, endofCooDec())

	write(file, declare_plaquettes(length(plq), plaq_length))

	for j1 in (1:length(plq))
		for j2 in (1:plaq_length)
			write(file, givePlaquetteCoords(N1, plq, j1, j2, Int(plq[j1][j2]), Int(plq[j1][j2+1])))
		end
	end

	write(file, endof_plaquettes())

	close(file)
end

#####################################################################
# Header string for .hpp-file:
function getHeaderString()
	#include <iostream>
	headerstring = """
	#include <armadillo>
	#include <complex>

	using namespace arma;
	
	"""

	return headerstring
end

# Declaration for matrix generating function:
function declare_def_matrix(matrix)
	N1 = size(matrix,1)
	N2 = size(matrix,2)
	str = """cx_mat def_matrix()
	{
	   cx_mat A($(N1),$(N2), fill::zeros);
	"""
	return str
end

# Fill matrix A with entries:
function giveEntries(matrix, i1, i2)
	str = """
             A($(i1-1), $(i2-1)) = std::complex<double>(0.0, 2.0*$(matrix[i1, i2]));
             A($(i2-1), $(i1-1)) = -std::complex<double>(0.0, 2.0*$(matrix[i1, i2]));
          """
	return str;
end

function giveMinusEntries(matrix, i1, i2)
	str = """
             A($(i1-1), $(i2-1)) = -std::complex<double>(0.0, 2.0*$(matrix[i1, i2]));
             A($(i2-1), $(i1-1)) = std::complex<double>(0.0, 2.0*$(matrix[i1, i2]));
          """
	return str;
end



# Close first function:
function endofDec()
	str = """
             return A;
          }
          
          """
end

# Second function: Create vector with matrix coordinates != 0:
function declare_coord()
	str = """std::vector<int> non_zeros()
{
    std::vector<int> nz;
    int coo;
"""
end

# Fill vector:
function giveCoordinates(N, i1, i2)
	str = """
              coo = $(N)*$(i1-1) + $(i2-1);
              nz.push_back(coo);
          """

end

# Close second function:
function endofCooDec()
	str = """
              return nz;
          }
          
          """
end

# Third function: Create matrix with coordinates of elementary plaquettes:

# Declaration for matrix generating function:
function declare_plaquettes(size1, size2)
	

	str = """Mat<int> create_plaquettes()
	{
	   Mat<int> P($(size1),$(size2));
	"""
	return str

end

# Fill matrix with encoded coordinates:
function givePlaquetteCoords(N, plq, j1, j2, z1, z2)
	str = """
             P($(j1-1), $(j2-1)) = $(N)*$(z1-1)+$(z2-1);
          """
	return str
end

# Close third function:
function endof_plaquettes()
	str = """
			 return P;
		  }
		  """
end


#####################################################################

export writehpp
export writehpp_zz
export getHeaderString

end
