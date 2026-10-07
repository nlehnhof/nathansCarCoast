using CSV, DataFrames, LsqFit, Plots

# Clean up data
df = CSV.read("LexusData.csv", DataFrame)
d1 = dropmissing(df[:, [:t1, :v1]])
d2 = dropmissing(df[:, [:t2, :v2]])

# set up parameters
m  = 1760 + 150
g  = 9.81
θ  = 1.03                      # degrees
F  = m*g*sind(θ)

# Define the model for velocity as a function of time, initial velocity, slope direction, and drag coefficient
# v(t) = fs*F/b + (v0 - fs*F/b) e^(-b t/m), fs = + downhill, - uphill
vmodel(t, v0, fs, b) = fs*F/b .+ (v0 - fs*F/b) .* exp.(-b/m .* t)

# Stack both runs into one fit; the first n1 points are run 1
n1 = nrow(d1)
tall = vcat(d1.t1, d2.t2)
vall = vcat(d1.v1, d2.v2)
model(x, p) = vcat(vmodel(x[1:n1],     d1.v1[1], +1, p[1]),
                   vmodel(x[n1+1:end], d2.v2[1], -1, p[1]))

fit = curve_fit(model, tall, vall, [10.0])
b = coef(fit)[1]
println(b)
# b = 20.6252

plot(d1.t1, d1.v1, label="South (data)")
plot!(d1.t1, vmodel(d1.t1, d1.v1[1], +1, b), label="South (model)")
plot!(d2.t2, d2.v2, label="North (data)")
plot!(d2.t2, vmodel(d2.t2, d2.v2[1], -1, b), label="North (model)")
xlabel!("t (s)"); ylabel!("v (m/s)")

savefig("linear_model_vs_data.png")