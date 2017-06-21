ts = collect(now():Minute(10):(99*Minute(10)+now()))
N = length(ts)

xf = randn(N)
xt = 2*xf + 1 + 0.1*randn(N)

qf = Stage(ts,xf)
qt = Temperature(ts,xt)

c = Calibration(qt,qf)

m1 = fit(PolynomialCalibrationModel,c,1)

xf2 = randn(N)
qf2 = Stage(ts,xf2)

qt2 = predict(m1,qf2)

@testset "Polynomial calibration models" begin
    
    @test quantity(qt2) == m1.β[1] + m1.β[2]*xf2
    @test typeof(qt2) == Temperature{Float64}
    
end
