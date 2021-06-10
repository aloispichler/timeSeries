"""
	Temperature in Germany
	created: 2021, April
	author©: Alois Pichler
"""

using Gnuplot; Gnuplot.options.gpviewer= true;	# external viewer
using CSV, DataFrames, Distributions, Dates

function basisFunction(t)	#	basis for linear regression
	return [1., t-2020, (t-2020)^2, sin(2π*t), cos(2π*t)]
end

println("────────────────── Temperature")
println("──────── https://de.wikipedia.org/wiki/Zeitreihe_der_Lufttemperatur_in_Deutschland ──────────")

df= CSV.read("HistoricTemperatureGermany.csv", DataFrame; delim=';', decimal=',', dateformat="m/d/yyyy")
sort!(df, [:date])							# dates ascending
insertcols!(df, 2, :time=> (1970).+Dates.datetime2unix.(DateTime.(df.date))./ 60/60/24/365.2422)

A= Array{Float64,2}(undef, size(df,1), 5)	# regression matrix
[A[i,:]= basisFunction(df.time[i]) for i= 1:size(df,1)]
weight= A\ df.temperature					# regression parameters, weights
insertcols!(df, 4, :regression=> A*weight)	# push the regression to the data frame

maℓag= 12		# moving average, 12 months
ma= copy(df.temperature); ma[1:maℓag-1].= NaN	# nothing to see here
for i= length(ma):-1:maℓag
	ma[i]= mean(ma[i-11:i])					# compute the moving average
end

@gp "reset; set multiplot layout 2,1; set border 0" :-
@gp :- 1 "set title 'historic temperatures, Germany'" :-
@gp :- "set xdata time; set timefmt '\"%Y-%m-%d\"'" :-
@gp :- "set format x '%m-%Y'; set xtics rotate by -30" :-
@gp :- """set xrange ['"1749-01-01"':'"2025-01-01"']""" :-
@gp :- "set style line 1 lc rgb 'blue' lt 1 lw 8 pt 7 ps 2.0"
@gp :- string.(df.date) df.regression  "using 1:2 ls -1 lt rgb'red' title 'regression' with linespoints" :-
@gp :- string.(df.date) df.temperature "using 1:2 ls -1 title 'temperature' with linespoints" :-
@gp :- string.(df.date) ma "using 1:2 ls -1 lw 3 lt rgb'blue' title 'annual moving average' with linespoints" :-
@gp :- 2 "unset title" :-
@gp :- 2 "set xzeroaxis linetype 1 linecolor 'black'" :-
@gp :- 2 string.(df.date) df.temperature-df.regression "using 1:2 ls -1 lt rgb'violet' title 'residual' with linespoints"

df	# display the dataframe