import Statistics

"""
	autocorrelation_function(x)
Returns the autocorrelation function: it accepts an integer `t` which represents the time at which evaluate the \
autocorrelation. `t` goes from 0 to N-1, where N is the length of `x`.
"""
function autocorrelation_function(x)
	x̄ = Statistics.mean(x)
	N = length(x)
	t -> sum((x[i] - x̄) * (x[i+t] - x̄) for i in 1:N-t) / (N - t)
end

"""
	autocorrelations(x)
Returns a vector containing the values of the autocorrelation function of `x` calculate at all times `t`. \
The first element is t=0, the second t=1, ...
"""
function autocorrelations(x)
	x̄ = Statistics.mean(x)
	N = length(x)
	[sum(((x[i] - x̄) * (x[i+t] - x̄)) for i in 1:N-t) / (N - t) for t in 0:N-1]
end

"""
	autocorrelation_time(x)
Returns the *integrated* autocorrelation time multiplied by two (see Dalla Brida's notes). Binning the values with a bin \
size bigger than the autocorrelation time guarantees the errors between measurements are independent.
"""
function autocorrelation_time(x, cutoff)
	N = cutoff
	Γ = autocorrelations(x)
	1 + 2*sum(Γ[t] / Γ[1] for t in 2:N)
end