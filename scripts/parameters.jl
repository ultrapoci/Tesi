allparams = @dict(
	dims = (8, 8, 8),
	latticestart = :cold,
	β = [1.0, 2.0],
	sp2type = Sp2ElementB,

	nterm = 20, # iterations of termalization
	nnorm = 10, # after how many cycle to normalize lattice
	nover = 3, # how many cycles of overrelaxation to do

	observables = (averageplaquette,), # function taking Lattice as input
	to_plot = true
)
