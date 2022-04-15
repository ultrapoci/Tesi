params = @dict(
	dims = (8, 8, 8),
	lattice_start = :cold,
	β = 1.0,
	sp2type = Sp2ElementB,

	nterm = 100, # iterations of termalization
	norm_every = 5, # after how many cycle to normalize lattice
	nover = 3 # how many cycles of overrelaxation to do
)
