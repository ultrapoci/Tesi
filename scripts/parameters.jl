params = @dict(
	dims = (8, 8, 8),
	lattice_start = :cold,
	β = 2.0,
	sp2type = Sp2ElementB,

	nterm = 200, # iterations of termalization
	norm_every = 10, # after how many cycle to normalize lattice
	nover = 3 # how many cycles of overrelaxation to do
)
