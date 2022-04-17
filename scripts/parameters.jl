allparams = @strdict(
	dims = (8, 8, 8),
	latticestart = :cold,
	β = [i for i in 1.0:15.0],
	sp2type = Sp2ElementB,

	nterm = 100, # iterations of termalization
	nnorm = 10, # after how many cycle to normalize lattice
	nover = 3, # how many cycles of overrelaxation to do
)

obsparams = @strdict(
	observables = averageplaquette, # function taking Lattice as input
	to_plot = true,
	meanoffset = 1, # from which iteration to being calculating the mean
)
