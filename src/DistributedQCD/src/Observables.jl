export averageplaquette

function averageplaquette(L::Lattice{D}, inds::Indices{D}) where D
	current_workers = [w for w in procs(L)]
	n_current_workers = length(current_workers)
	# initialize a vector with nworkers() elements, of which each worker owns one cell
	partial = dzeros((n_current_workers,), current_workers, [n_current_workers])
	with_workers() do _	
		tot = 0.0
		for x in CartesianIndices(L[:L]), u in 1:D-1, v in u+1:D
			tot += tr(L[:L][x][u] * staple(L, v, u, inds[:L][x]))
		end
		partial[:L][begin] = tot
	end
	res = sum(partial)
	V = length(L) # lattice's volume
	np = D * (D - 1) ÷ 2 # number of plaquettes per site

	# the 4 term is due to 2N for N=2
	res / (4 * np * V) 
end
averageplaquette(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}) where D = averageplaquette(T[1], T[3])
averageplaquette(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}) where D = averageplaquette(T.lattice, T.inds)