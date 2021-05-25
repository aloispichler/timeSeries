"""
	visualize exemplary time series: ARMA, Yule–Walker, constant ACF
	created: 2021, April
	author©: Alois Pichler
"""

using Gnuplot; Gnuplot.options.gpviewer= true	# external viewer
using Distributions, StatsBase

struct timeSeries
	Xt::Vector{Float64}	# hold the time series' realization
	meta::Any			# meta information
end

#	│	ARMA time series
#	╰────────────────────────────────────────────────────
function ARMA(T::Int64; φ= [0.8, 0.1, -0.1], θ= [1.7, 1.4, .2, 1.], Z::UnivariateDistribution= Normal(0,1))::timeSeries
	p= length(φ); X= rand(Z, T)		# AR: initialize the time series
	q= length(θ); Z= rand(Z, T)		# MA: random noise
	for t= max(p,q)+1: T			# apply the linear time series
		X[t]= reverse(φ)' * X[t-p:t-1] + reverse([1; θ])' * Z[t-q:t]
	end
	return timeSeries(X, "ARMA(Φ=$φ, θ=$θ)")
end

#	│	this time series has constant autocovariance function
#	╰────────────────────────────────────────────────────────────
function constantACF(T::Int64; ρ= 0.9, Z::UnivariateDistribution= Normal(0,1))::timeSeries
	X̅= 0.0; X= Vector{Float64}(undef, T)
	Z= rand(Z, T)			# random noise
	for t= 1:T
		ρt= (t-1)*ρ/(1+ (t-2)*ρ)
		X[t]= ρt* X̅+ sqrt(1- ρt* ρ)* Z[t]
		X̅= (X̅*(t-1) + X[t])/ t	# update the running average
	end
	return timeSeries(X, "constant acf: ρ= $ρ")
end

#	│	Yule–Walker: the autocorrlation ρ is given
#	╰────────────────────────────────────────────────────
function YuleWalker(T::Int64; ρ= [0.1, -0.1], Z::UnivariateDistribution= Normal(0,1))::timeSeries
	σ²= var(Z); γ= σ² * [1; ρ]			# auto covariance vector
	Γ= Matrix{Float64}(undef, length(γ)-1, length(γ)-1)
	for i= 1:length(γ)-1, j= 1:i
		Γ[i,j]= Γ[j,i]= γ[i-j+1]	# symmetric Töplitz matrix
	end
	φ= Γ\ γ[2:end]			# solve Yule–Walker equations
	ψ²= σ² - φ'* γ[2:end]; ψ² ≥ 0 || @error "ψ² is negative."
	X= rand(Z, T); Z= X		# make some noise
	for t= length(γ): T		# apply the linear time series
		X[t]= reverse(φ)' * X[t-length(γ)+1:t-1] + sqrt(ψ²)* Z[t]
	end
	return timeSeries(X, "Yule–Walker: Φ=$(round.(φ, sigdigits=3)), ψ=$(round(sqrt(ψ²), sigdigits=3)), γ=$(round.(γ, sigdigits=3))")
end


#	│	main: visualize the time series
#	╰────────────────────────────────────────────────────
maxTime= 1000; copies= 3; ts= timeSeries
distZ= Uniform(-1, 1)		# Normal(0.0, 1.), Uniform(-1, 1), Cauchy(0, 1), Exponential(0.3)
@gp "reset; set encoding utf8; set multiplot layout 2,1" :-
@gp :- 1 "set border 0; set xlabel 'time t'; set zeroaxis linetype 1 linecolor 'black'; set xtics axis add ('' 0)" :-
for i= 1:copies
	ts= YuleWalker(maxTime; ρ= [0.9, 0.8], Z= distZ)	# ARMA(maxTime; Z= distZ), YuleWalker(maxTime; ρ= [0.9, 0.8]), constantACF(maxTime; ρ= 0.99, Z= distZ)
	@gp :- 1 ts.Xt "with lines title 'realization $i'" :-
	@gp :- 2 autocov(ts.Xt) "with boxes title ''" :-
end
@gp :- 1 "set title '$(ts.meta)'" :-
@gp :- 2 "set title 'autocovariance'; set border 0; set xlabel 'lag ℓ'"
time()