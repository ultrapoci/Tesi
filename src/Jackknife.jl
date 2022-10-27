import Statistics, Measurements

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

"""
	binsamples(v, binsize)
Divide a vector `v` into bins of size `binsize`. If `v` cannot be divided perfectly, increase the binsize of the first bins \
by one, until the vector is completely divided. This implies that each bin is identical to each other or at most has one more element.
"""
function binsamples(v, binsize)
	r = length(v) % binsize
	if r == 0
		vcat(convert.(Array, Iterators.partition(v, binsize)))
	else
		n = length(v) ÷ binsize # number of bins 
		new_binsize = binsize + (r ÷ n)
		new_r = length(v) % new_binsize
		# place the remaining elements in each bin, starting from the first one
		split_at = new_r * (new_binsize + 1)
		head = convert.(Array, Iterators.partition(v[begin:split_at], new_binsize + 1))
		tail = convert.(Array, Iterators.partition(v[split_at+1:end], new_binsize))
		vcat(head, tail)
	end
end

"""
	jackknife_sample(v, provided_mean = nothing)
Given a list of samples `v`, returns the *jackknife samples* of such list. A jackknife sample is the mean of all samples \
removing the considered sample from the list. 
"""
function jackknife_sample(v, provided_mean = nothing)
	v̄ = isnothing(provided_mean) ? Statistics.mean(v) : provided_mean
	N = length(v)
	[v̄ - (x - v̄) / (N - 1) for x in v]
end

"""
	jackknife(f, args...)
Return the expected value and the error of the function `f` calculated over `args`, using the jackknife method. 
	
`args` is a list of vectors of the the same length, representing the samples over which to calculate mean and error.

The returned value is a `measurement` object from Measurements.jl.

# Example
```jldoctest
julia> v = rand(100);

julia> w = rand(100);

julia> f(x, y) = x + y;

julia> jackknife(f, v, w)
0.995 ± 0.039
```
"""
function jackknife(f, args...)
	N = length(first(args))
	all(x -> length(x) == N, args) || throw(ArgumentError("args must have the same length"))

	args_mean = Statistics.mean.(args)
	jackknife_samples = jackknife_sample.(args, args_mean)

	f_mean = f(args_mean...)
	f_jackknife_samples = f.(jackknife_samples...)
	f_error = (N - 1) * sum((f_j - f_mean)^2 for f_j in f_jackknife_samples) / N

	Measurements.measurement(f_mean, sqrt(f_error))
end