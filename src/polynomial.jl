# Intervals used for dispatch in interval

export Interval, Confidence, Prediction

abstract type Interval end

struct Confidence <: Interval
    α
end

struct Prediction <: Interval
    α
end


"""
    PolynomialCalibrationModel{T,F}(k,β,c::Calibration{T,F})

A CalibrationModel representing a polynomial regression.
"""
struct PolynomialCalibrationModel{T<:Quantity, F<:Quantity} <: CalibrationModel{T,F}
    k::Int
    β::Vector{Float64}
    c::Calibration{T,F}
end

StatsBase.coef(m::PolynomialCalibrationModel) = m.β

StatsBase.nobs(m::PolynomialCalibrationModel) = length(to_quantity(m.c))

"""
    design_polynomial(Q,k)

Create the design matrix for a polynomial regression
of highest order k
"""
function design_polynomial(Q,k)
    X = ones(length(Q),k+1)
    for i in 1:k+1
        X[:,i] = Q.^(i-1)
    end
    X
end

function design_polynomial(Q::Quantity,k)
    design_polynomial(quantity(Q),k)
end

function design_polynomial(c::Calibration,k)
    qi = interpolatecal(c)
    design_polynomial(qi,k)
end

design_polynomial(m::PolynomialCalibrationModel,k) = design_polynomial(m.c,k)

function StatsBase.fit(::Type{PolynomialCalibrationModel},c::Calibration,k)
    X = design_polynomial(c,k)
    y = quantity(to_quantity(c))
    β = X\y
    PolynomialCalibrationModel(k,β,c)
end

function StatsBase.predict{T,F}(m::PolynomialCalibrationModel{T,F},q::F)
    X = design_polynomial(q,m.k)
    T(times(q),X*m.β)
end

function StatsBase.predict{T,F}(m::PolynomialCalibrationModel{T,F},q)
    X = design_polynomial(q,m.k)
    X*m.β
end

"""
Prediction from a PolynomialCalibrationModel with an Interval
"""
function StatsBase.predict(m::PolynomialCalibrationModel{T,F},q,interval::Interval) where {T,F}
    X = design_polynomial(q,m.k)
    p = X*m.β
    i = interval(m,X)
    p, i
end

function StatsBase.predict(m::PolynomialCalibrationModel{T,F},q::F,interval::Interval) where {T,F}
    X = design_polynomial(q,m.k)
    p = X*m.β
    i = interval(m,X)
    T(times(q),p), i
end

StatsBase.residuals(m::PolynomialCalibrationModel) = quantity(to_quantity(m.c)) - predict(m,m.c)

function resstderr(m::PolynomialCalibrationModel)
    e = residuals(m)
    n = nobs(m)
    k = length(coef(m))
    sqrt(e'e/(n-k))
end    

function StatsBase.vcov(m::PolynomialCalibrationModel)
    s = resstderr(m)
    X = design_polynomial(m.c,m.k)
    Q = X'X
    s*inv(Q)
end

# This is a confidence interval for the coefficients
function StatsBase.confint(m::PolynomialCalibrationModel,α)
    q = quantile(Normal(),1-α/2)
    s = stderr(m)
    m.β.+q*s,m.β+q*s
end

# This is a confidence interval for predictions on the matrix X
function StatsBase.confint(m::PolynomialCalibrationModel,X,α)
    q = quantile(Normal(),1-α/2)
    V = vcov(m)
    q*sqrt.([X[i,:]'V*X[i,:] for i in 1:size(X,1)])
end

"""
Calculate the prediction interval width for a PolynomialCalibrationModel
"""
function predint(m::PolynomialCalibrationModel,X,α)
    q = quantile(Normal(),1-α/2)
    s2 = resstderr(m)
    V = vcov(m)
    q*sqrt.([X[i,:]'V*X[i,:] + s2 for i in 1:size(X,1)])
end

function (p::Confidence)(m::PolynomialCalibrationModel,X)
    confint(m,X,p.α)
end

function (p::Prediction)(m::PolynomialCalibrationModel,X)
    predint(m,X,p.α)
end
