module TidalFluxCalibrations

export calibratePolynomial, calibrateData

using DischargeData, Interpolations

"""
    interpolateCalibration(c::Calibration)

Takes a Calibration and returns the from quantity
that has been interpolated to the times of the 
to quantity
"""
function interpolateCalibration{T}(cal::Calibration{T})
    qi = interpolate((times(from_quantity(cal)),),
                     quantity(from_quantity(cal)),
                     Gridded(Linear())
                     )
    ts = times(to_quantity(cal))
    qs = qi[ts]
    T(ts,qs)
end

"""
Constructs a polynomial regression of order k to calibrate 
ADCP discharge to the true discharge
"""
function calibratePolynomial(cal::Calibration,k)
    qs = interpolateCalibration(cal)
    X = ones(length(qs),k+1)
    Q = quantity(qs)
    for i in 1:k+1
        X[:,i] = Q.^(i-1)
    end
    X\quantity(to_quantity(cal))
end

"""
Performs a global calibration for multiple Calibrations
"""
function calibratePolynomial{T}(cals::Vector{Calibration{T}},k)
    qs = vcat(quantity.(interpolateCalibration.(cals))...)
    X = ones(length(qs),k+1)
    for i in 1:k+1
        X[:,i] = qs.^(i-1)
    end
    qq = vcat(quantity.(to_quantity.(cals))...)
    X\qq
end

function calibrateData{T<:Quantity}(q::T,β::Vector{Float64})
    k = length(β)-1
    X = zeros(length(q),k+1)
    Q = quantity(q)
    for i in 1:k+1
        X[:,i] = Q.^(i-1)
    end
    Discharge(times(q),X*β)
end

function calibrateData{T<:Quantity}(cal::Calibration{T},q::T,k::Int)
    β = calibratePolynomial(cal,k)
    calibrateData(q,β)
end

function calibrateData{T<:Quantity}(cals::Vector{Calibration{T}},q::T,k::Int)
    β = calibratePolynomial(cals,k)
    calibrateData(q,β)
end

end # module
