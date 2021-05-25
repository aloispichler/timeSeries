"""
	Temperature in Germany
	created: 2021, April
    author©: Alois Pichler
"""

using CSV, DataFrames, Distributions 
using Gnuplot, Dates

function basisFunction(t)	#	basis for linear regression
	return [1., t, sin(2π*t), cos(2π*t)]
end

println("________________________ Temperature")
df= CSV.read("C:\\Users\\alopi\\Dropbox\\Julia\\StochasticProcess\\HistoricTemperatureGermany.csv", DataFrame; delim=';', decimal=',', dateformat="m/d/yyyy")
insertcols!(df, 2, :time=> (1970).+Dates.datetime2unix.(DateTime.(df.date))./ 60/60/24/365.2422)

A= Array{Float64,2}(undef, size(df,1), 4)
[A[i,:]= basisFunction(df.time[i]) for i= 1:size(df,1)]
weight= A\ df.temperature
insertcols!(df, 4, :regression=> A*weight)

ma= copy(df.temperature)
for i= length(ma):-1:12
	for j=1:11
		ma[i]+= ma[i-j]
	end
	ma[i]/= 12
end
ma[1:11].= NaN

@gp "reset"
@gp :- "set title 'historic temperatures, Germany'"
@gp :- "set xdata time" "set timefmt '\"%Y-%m-%d\"'"
@gp :- "set format x '%m-%Y'" "set xtics rotate by -30"
@gp :- """set xrange ['"1749-01-01"':'"2025-01-01"']"""
@gp :- "set xzeroaxis linetype 1 linecolor 'black'; set border 0"
@gp :- "set style line 1 lc rgb 'blue' lt 1 lw 8 pt 7 ps 2.0"
@gp :- string.(df.date) df.regression  "using 1:2 ls -1 lt rgb'red' title 'regression' with linespoints"
@gp :- string.(df.date) df.temperature "using 1:2 ls -1 title 'temperature' with linespoints"
@gp :- string.(df.date) ma "using 1:2 ls -1 lw 3 lt rgb'blue' title 'annual moving average' with linespoints"
@gp :- string.(df.date) df.temperature-df.regression.-15 "using 1:2 ls -1 lt rgb'violet' title 'residual' with linespoints"

df