using MAT
using CSV
using DataFrames
using JuMP
using Gurobi
using DifferentialEquations
using Plots

# Cargar modelo desde .mat (ajusta el nombre del archivo y variable)
matfile = matopen("yeast9-gem.mat")
model_data = read(matfile, "model")
close(matfile)

# Cargar parámetros desde CSV
params = CSV.read("R6_params.csv", DataFrame)
parameters = Matrix(params)[:]
y0 = parameters[1:6]
theta = parameters[7:end]

# Definir índices de reacciones (ajusta según tu modelo)
index_rxns = [1,2,3,4] # Reemplaza con los índices correctos

# Crear modelo JuMP desde datos (debes adaptar esto a tu estructura)
model = Model(Gurobi.Optimizer)
# ... aquí debes construir el modelo JuMP a partir de model_data ...
# Por ejemplo: variables, restricciones, función objetivo, etc.

# Parámetros para la función ODE
p = (theta, model, index_rxns)

# Configuración de integración
prob = ODEProblem(model_dFBA2fit_sanchez!, y0, (0.0, 11.0), p)
sol = solve(prob, Tsit5(), abstol=1e-6, reltol=1e-6)

# Cargar datos experimentales
exp_data = CSV.read("Aerobic_data_sanchez_1.csv", DataFrame)

# Graficar resultados
plot(sol.t, sol[1,:], label="Biomass")
scatter!(exp_data.TIME, exp_data.X, label="Exp. Biomass")
plot!(sol.t, sol[2,:], label="Glucose")
scatter!(exp_data.TIME, exp_data.GLC, label="Exp. Glucose")
plot!(sol.t, sol[3,:], label="Ethanol")
scatter!(exp_data.TIME, exp_data.ETOH, label="Exp. Ethanol")
# ... agrega más subplots según tus variables ...

# Guardar resultados
results = DataFrame(Time=sol.t, Biomass=sol[1,:], Glucose=sol[2,:], Ethanol=sol[3,:])
CSV.write("DFBA_results.csv", results)
