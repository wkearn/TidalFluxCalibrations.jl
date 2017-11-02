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
