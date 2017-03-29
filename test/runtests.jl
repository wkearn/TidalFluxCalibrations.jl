using DischargeData, TidalFluxCalibrations, Base.Dates
using Base.Test
TFC = TidalFluxCalibrations

# Generate an hourly time series
tst = DateTime(2016,01,01):Hour(1):DateTime(2016,02,01)
nt = length(tst)

# Generate a high-frequency time series
tsn = DateTime(2016,01,01):Minute(10):DateTime(2016,02,01)
nn = length(tsn)

# Create random data
qtn = randn(nn)
qt = qtn[1:6:end]
qn = 0.5*qtn+0.1*randn(nn)

# Create quantities out of it
dt = Discharge(tst,qt)
dn = Discharge(tsn,qn)

# Make a calibration
cal = Calibration(dt,dn)

# Interpolate the calibration to the to times
dc = TFC.interpolateCalibration(cal)

# Calibrate linearly
β1 = calibratePolynomial(cal,1)

# Calibrate on more than one Calibration
# This should be roughly equal to the calibration
# with just one Calibration
cals = [cal;cal]
β2 = calibratePolynomial(cals,1)

@test_approx_eq β1 β2

dn1 = calibrateData(dn,β1)
dn2 = calibrateData(cal,dn,1)
dn3 = calibrateData(cals,dn,1)

@test dn1 == dn2
@test_approx_eq quantity(dn2) quantity(dn3)

###################################################
# Now do it all over again, but with a different
# Quantity to make sure the generic interface works

# Create quantities out of it
dt = Stage(tst,qt)
dn = Stage(tsn,qn)

# Make a calibration
cal = Calibration(dt,dn)

# Interpolate the calibration to the to times
dc = TFC.interpolateCalibration(cal)

# Calibrate linearly
β1 = calibratePolynomial(cal,1)

# Calibrate on more than one Calibration
# This should be roughly equal to the calibration
# with just one Calibration
cals = [cal;cal]
β2 = calibratePolynomial(cals,1)

@test_approx_eq β1 β2

dn1 = calibrateData(dn,β1)
dn2 = calibrateData(cal,dn,1)
dn3 = calibrateData(cals,dn,1)

@test dn1 == dn2
@test_approx_eq quantity(dn2) quantity(dn3)
