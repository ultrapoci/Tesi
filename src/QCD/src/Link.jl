using LinearAlgebra: I

export Link

struct Link{D}
	position::NTuple{D, Integer}
	direction::Integer
	m::Matrix{Complex} # TODO: see if there's a better way to store Sp(2) in memory

	# TODO: cold and hot start
	function Link(position::NTuple{D, Integer}, direction::Integer) where D
		if direction > D 
			throw(ArgumentError("direction = $direction of given link is bigger than dimensions = $D of the lattice"))
		end
		m = Matrix{Complex}(I, 4, 4) # TODO: see if Hermitian(m) is useful 
		new{D}(position, direction, m)
	end            
end