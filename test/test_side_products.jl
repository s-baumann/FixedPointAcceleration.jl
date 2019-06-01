# Generating data from two two-dimensional gaussian processes
using Distributions
using FixedPointAcceleration
using Random
using DataFrames
using Plots
using LinearAlgebra
true_tau = 0.6
nobs_1 = 400
nobs_2 = convert(Int, round(nobs_1 * ((1-true_tau)/true_tau)))
Random.seed!(1234)
mu_1 = [0.0,8.0]
cov_1 = [2.0,0.5,2.0]
covar_1 = Symmetric([cov_1[1] cov_1[2]; cov_1[2] cov_1[3]])
md_1 = MultivariateNormal(mu_1,covar_1)
mu_2 = [-4.0,10.0]
cov_2 = [2.0,-0.75,12.0]
covar_2 = Symmetric([cov_2[1] cov_2[2]; cov_2[2] cov_2[3]])
md_2 = MultivariateNormal(mu_2,covar_2)

rands_from_1 = transpose(rand(md_1, nobs_1))
rands_from_2 = transpose(rand(md_2, nobs_2))
data1 = DataFrame([rands_from_1[:,1], rands_from_1[:,2]], [:x1, :x2])
data2 = DataFrame([rands_from_2[:,1], rands_from_2[:,2]], [:x1, :x2])
dd  = vcat(data1,data2)
# Plotting it:
plot(data1.x1, data1.x2,seriestype=:scatter)
plot!(data2.x1, data2.x2,seriestype=:scatter)



function z_estimate_given_theta(x::Array{Float64,1}, md_1::MultivariateNormal, md_2::MultivariateNormal, tau::Float64)
    pdf_1 = pdf(md_1, x)
    pdf_2 = pdf(md_2, x)
    return tau*pdf_1 / (tau*pdf_1 + (1-tau)*pdf_2)
end

function update_tau(Z::Array{Float64,1})
    return mean(Z)
end

function update_mu(dd::DataFrame, Z::Array{Float64,1})
    X = convert(Array{Float64,2}, dd[[:x1, :x2]])
    sum_Z = sum(Z)
    updated_mu = (transpose(Z) * X) ./sum_Z
    return vec(updated_mu)
end

function update_cov(dd::DataFrame, updated_mu::Array{Float64,1}, Z::Array{Float64,1})
    X_minus_mu = convert(Array{Float64,2}, dd[[:x1, :x2]]) .- transpose(updated_mu)
    sum_Z = sum(Z)
    updated_cov = (transpose(Z .* X_minus_mu) * X_minus_mu) ./sum_Z
    return [updated_cov[1,1], updated_cov[1,2], updated_cov[2,2]]
end

function update_theta(theta::Array{Float64,1}, dd::DataFrame)
    # We will use the convention that theta's 11 entries are (mu_1, cov_1, mu_2, cov_2, tau). First unpacking theta:
    mu_1    = theta[[1,2]]
    cov_1   = theta[[3,4,5]]
    covar_1 = Symmetric([cov_1[1] cov_1[2]; cov_1[2] cov_1[3]])
    md_1 = MultivariateNormal(mu_1,covar_1)
    mu_2    = theta[[6,7]]
    cov_2   = theta[[8,9,10]]
    covar_2 = Symmetric([cov_2[1] cov_2[2]; cov_2[2] cov_2[3]])
    md_2 = MultivariateNormal(mu_2,covar_2)
    tau     = theta[11]
    # Getting Z
    Z = Array{Float64,1}(undef,size(dd)[1])
    for i in 1:size(dd)[1]
        Z[i] = z_estimate_given_theta([dd[i,:x1], dd[i,:x2]], md_1, md_2, tau)
    end

    # Updating Tau
    updated_tau = update_tau(Z)
    # Updating mu1
    updated_mu_1 = update_mu(dd,Z)
    updated_mu_2 = update_mu(dd, 1 .- Z)
    # Updating Cov
    updated_cov_1 = update_cov(dd, updated_mu_1, Z)
    updated_cov_2 = update_cov(dd, updated_mu_2, 1 .- Z)
    # Returning theta
    updated_theta = vcat(updated_mu_1, updated_cov_1, updated_mu_2, updated_cov_2, updated_tau)
    side_product = (z = Z,)
    return (updated_theta, side_product)
end



InitialGuess = [0.5, 7.5, 2.0, 0.0, 2.0, -5.0, 7.5, 2.0, 0.0, 10.0, 0.5]
fp_anderson = fixed_point(x -> update_theta(x,dd), InitialGuess; Algorithm = Anderson, PrintReports = true)
isa(fp_anderson.Other_Output_, NamedTuple)
