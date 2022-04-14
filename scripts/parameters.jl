params = @dict(
	dims = (8, 8, 8),
	lattice_start = :hot,
	β = 1.0,

	nterm = 50, # iterations of termalization
	norm_every = 5, # after how many cycle to normalize lattice
	nover = 3 # how many cycles of overrelaxation to do
)
