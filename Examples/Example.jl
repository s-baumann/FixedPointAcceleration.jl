x = [1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6]
y = log(x) + sqrt(x)
gradients = "Not-Supplied"

using SchumakerSpline


########################
# Linear Extrapolation #
extrapolation = "Linear"
Spline, DerivSpline, SecondDerivSpline = schumaker(x,y, gradients, extrapolation)

#Pkg.add("PyPlot")
using PyPlot
#Pkg.update()
xArray = linspace(-5, 10, 100)
Spp = map(x -> DerivSpline(x), xArray)
plt = plot(xArray, Spp, color="green", linewidth=2.0, linestyle="--")
plt
Spp = map(x -> Spline(x), xArray)
plt = plot(x, y, color="blue", linewidth=2.0, linestyle="-")
plt
plt = plot(xArray, Spp, color="red", linewidth=2.0, linestyle="--")
plt
Spp = map(x -> DerivSpline(x), xArray)
plt = plot(xArray, Spp, color="green", linewidth=2.0, linestyle="--")
plt
Spp = map(x -> SecondDerivSpline(x), xArray)
plt = plot(xArray, Spp, color="Black", linewidth=2.0, linestyle="--")
plt


##########################
# Constant Extrapolation #
extrapolation = "Constant"
Spline, DerivSpline, SecondDerivSpline = schumaker(x,y, gradients, extrapolation)

#Pkg.add("PyPlot")
using PyPlot
#Pkg.update()
xArray = linspace(-5, 10, 100)
Spp = map(x -> Spline(x), xArray)
plt = plot(x, y, color="blue", linewidth=2.0, linestyle="-")
plt
plt = plot(xArray, Spp, color="red", linewidth=2.0, linestyle="--")
plt

##########################
# Curve    Extrapolation #
Spline, DerivSpline, SecondDerivSpline = schumaker(x,y, gradients)

#Pkg.add("PyPlot")
using PyPlot
#Pkg.update()
xArray = linspace(-5, 10, 100)
Spp = map(x -> Spline(x), xArray)
plt = plot(x, y, color="blue", linewidth=2.0, linestyle="-")
plt
plt = plot(xArray, Spp, color="red", linewidth=2.0, linestyle="--")
plt
