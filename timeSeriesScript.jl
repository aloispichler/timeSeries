"""
	time series: ARMA, Yule–Walker, constant ACF
	created: 2021, April
    author©: Alois Pichler
"""

using Gnuplot; Gnuplot.options.gpviewer= true;	# external viewer
using Distributions, StatsBase

struct timeSeries
	Xt::Vector{Float64}	# holds the time series' elements
	meta::Any			# meta information
end

#	│	this time series has constant autocovariance function
#	╰────────────────────────────────────────────────────────
function constantACF(T; ρ= 0.9, distZ::UnivariateDistribution= Normal(0,1))::timeSeries
	Z= rand(distZ, T)			# random noise
	X̅= 0.0; X= Vector{Float64}(undef, T)
	for t= 1:T
		ρt= (t-1)*ρ/(1+ (t-2)*ρ)
		X[t]= ρt* X̅+ sqrt(1- ρt* ρ)* Z[t]
		X̅= X̅*(t-1)/t + X[t]/ t	# update the mean
	end
	return timeSeries(X, "constant acf: ρ= $ρ")
end

#	│	ARMA time series
#	╰────────────────────────────────────────────────────
function ARMA(T; φ= [0.8, 0.1, -0.1], θ= [1.7, 1.4, .2, 1.], distZ::UnivariateDistribution= Normal(0,1))
	p= length(φ); X= rand(distZ, T)		# AR: to start the time series
	q= length(θ); Z= rand(distZ, T)		# MA: random noise
	for t= max(p,q)+1: T
		X[t]= reverse(φ)' * X[t-p:t-1] + reverse([1; θ])' * Z[t-q:t]
	end
	return timeSeries(X, "ARMA(Φ=$φ, θ=$θ)")
end

#	│	Yule–Walker: the autocovariance ρ is given
#	╰────────────────────────────────────────────────────
function YuleWalker(T; ρ= [0.1, -0.1], distZ::UnivariateDistribution= Normal(0,1))
	σ²= var(distZ); γ= σ² * [1; ρ]			# auto covariance
	Γ= σ²* Matrix{Float64}(undef, length(γ)-1, length(γ)-1)
	for i= 1:length(γ)-1		# Töplitz matrix
		for j= 1:length(γ)-1
			Γ[i,j]= γ[abs(i-j)+1]
	end	end
	φ= Γ\ γ[2:end]				# solve Yule–Walker equations
	ψ²= σ² - φ'* γ[2:end]
	ψ² ≥ 0 || @error "ψ² is negative."
	X= rand(distZ, T); Z= X		# make some noise
	for t= length(γ): T
		X[t]= reverse(φ)' * X[t-length(γ)+1:t-1] + sqrt(ψ²)* Z[t]
	end
	return timeSeries(X, "Yule–Walker: Φ=$(round.(φ, sigdigits=3)), γ=$γ")
end


#	│	main
#	╰────────────────────────────────────────────────────
T= 1000; copies= 5; ts= timeSeries
distZ= Normal(0.0, 1.)		# Normal(0.0, 1.), Uniform(-1, 1), Cauchy(0, 1)
@gp "reset; set encoding utf8"
@gp :- "set multiplot layout 2,1" :-
@gp :- 1 "set border 0; set xlabel 'time t'; set zeroaxis linetype 1 linecolor 'black'; set xtics axis add ('' 0)" :-
for i= 1:copies
	ts= YuleWalker(T; ρ= [0.9, 0.8, 0.7])	# ARMA(T; distZ= distZ), YuleWalker(T; ρ= [0.9, 0.8]), constantACF(T; ρ= 0.9)
	@gp :- 1 ts.Xt "with lines title 'realization $i'" :-
	@gp :- 2 autocov(ts.Xt) "with boxes title ''" :-
end
@gp :- 1 "set title '$(ts.meta)'" :-
@gp :- 2 "set title 'autocovariance'; set border 0; set xlabel 'lag ℓ'"
time()