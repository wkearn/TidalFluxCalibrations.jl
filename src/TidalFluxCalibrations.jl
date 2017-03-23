module TidalFluxCalibrations

export calibratePolynomial, calibrateData

using DischargeData, Interpolations

function interpolateCalibration(cal::CalibrationData)
    Qi = interpolate((cal.dd.ts,),cal.dd.Q,Gridded(Linear()))
    Qs = Qi[cal.t]
end

"""
Constructs a polynomial regression of order k to calibrate 
ADCP discharge to the true discharge
"""
function calibratePolynomial(cal::CalibrationData,k)
    Qs = interpolateCalibration(cal)
    X = ones(length(Qs),k+1)
    for i in 1:k+1
        X[:,i] = Qs.^(i-1)
    end
    X\cal.Q
end

"""
Performs a global calibration for multiple CalibrationDatas
"""
function calibratePolynomial(cals::Vector{CalibrationData},k)
    Qs = vcat(interpolateCalibration.(cals)...)
    Qq = vcat([cals[i].Q for i in 1:length(cals)]...)
    X = ones(length(Qs),k+1)
    for i in 1:k+1
        X[:,i] = Qs.^(i-1)
    end
    X\Qq
end

function calibrateData(dd::Discharge,β::Vector{Float64})
    k = length(β)-1
    X = zeros(length(dd.Q),k+1)
    for i in 1:k+1
        X[:,i] = dd.Q.^(i-1)
    end
    Discharge(dd.cp,dd.ts,dd.vs,dd.A,X*β)
end

function calibrateData(cal::CalibrationData,dd::Discharge,k::Int)
    β = calibratePolynomial(cal,k)
    calibrateData(dd,β)
end

function calibrateData(cals::Vector{CalibrationData},dd::Discharge,k::Int)
    β = calibratePolynomial(cals,k)
    calibrateData(dd,β)
end

end # module
