using DrWatson
@quickactivate "Tesi"

# QCD and DistributedQCD are local packages
import QCD as Q
import DistributedQCD as DQ
import Term.progress
import Distributed: rmprocs

function one_termalization!(L::Q.Lattice)
	Q.lattice_overrelaxation!(L, 3)
	Q.lattice_heatbath!(L, 1.0)
	Q.lattice_normalization!(L)
end

function testQ!(L)
	progress.@track for i in 1:20
		one_termalization!(L)
	end
end

function testQsleep!(L)
	progress.@track for i in 1:20
		sleep(0.01)
		one_termalization!(L)
	end
end

function testDQ!(L)
	progress.@track for i in 1:20
		DQ.one_termalization!(L, 3, 1.0, true)
	end
end

function testDQsleep!(L)
	progress.@track for i in 1:20
		sleep(0.01)
		DQ.one_termalization!(L, 3, 1.0, true)
	end
end

## Old algorithm

LQ = Q.Lattice(4, 4, 4)
@info "Without sleep"
testQ!(LQ) # doesn't show the progress bar

@info "With sleep"
testQsleep!(LQ) # shows the progress bar

## New algorithm

DQ.nworkers() == 1 && DQ.initprocs(4)
DQ.@everywhere using DrWatson
DQ.@everywhere @quickactivate "Tesi"
DQ.@everywhere import DistributedQCD

LDQ = DQ.newlattice(4, 4, 4)

@info "Without sleep"
testDQ!(LDQ) # shows the progress bar

@info "With sleep"
testDQsleep!(LDQ) # shows the progress bar

rmprocs(DQ.workers())