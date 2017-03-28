using DischargeData, PIEMetData, TidalFluxExampleData, ADCPDataProcessing, TidalFluxCalibrations
using Base.Test

setADCPdatadir!(Pkg.dir("TidalFluxExampleData","data","adcp"))
setmetdatadir!(Pkg.dir("TidalFluxExampleData","data","met"))

creek = Creek{:sweeney}()
deps = parse_deps(creek)
adata = load_data.(deps)

cs = parse_cs(creek)
csd = load_data(cs)

dd = Discharge(adata[4],csd)

cals = parse_cals(creek)
cc = load_data(cals[1])

calibrateData([cc],dd,1)
