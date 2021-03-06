module TidalFluxCalibrations

using Reexport

@reexport using TidalFluxQuantities, StatsBase, Interpolations

export CalibrationModel,
    PolynomialCalibrationModel,
    Calibration,
    to_quantity,
    from_quantity

abstract type CalibrationModel{T<:Quantity, F<:Quantity} <: RegressionModel end

include("calibrations.jl")

"""
    interpolateCalibration(c::Calibration)

Takes a Calibration and returns the from quantity
that has been interpolated to the times of the 
to quantity
"""
function interpolatecal{T,F}(cal::Calibration{T,F})
    qi = interpolate((times(from_quantity(cal)),),
                     quantity(from_quantity(cal)),
                     Gridded(Linear())
                     )
    ts = times(to_quantity(cal))
    qs = qi[ts]
    T(ts,qs)
end


include("polynomial.jl")

end # module
