using LinearAlgebra: I

export Link

struct Link{D}
	position::NTuple{D, Integer}
	direction::Integer
	m::Sp2Element

	# TODO: cold and hot start
	function Link(position::NTuple{D, Integer}, direction::Integer) where D
		if direction > D 
			throw(ArgumentError("direction = $direction of given link is bigger than dimensions = $D of the lattice"))
		end
		m = Sp2Element() 
		new{D}(position, direction, m)
	end            
end