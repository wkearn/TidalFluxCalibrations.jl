# TidalFluxCalibrations

[![Build Status](https://travis-ci.org/wkearn/TidalFluxCalibrations.jl.svg?branch=master)](https://travis-ci.org/wkearn/TidalFluxCalibrations.jl)

Part of [TidalFluxes.jl](https://github.com/wkearn/TidalFluxes.jl).

This package holds a set of routines to perform calibrations of discharge and other concentration data from tidal channels. The basic idea is that you have an instrument in the channel which records a time series and you take point measurements of a related quantity from the channel. You want to convert your time series into a time series of that quantity. The simplest way to do this is to fit a regression of the quantity of interest on the time series quantity. Currently, this package implements a polynomial regression.

# Roadmap

- [ ] Generic calibration routines for tidal flux quantities
- [ ] Integration with Julian statistical APIs
- [ ] More complex models (especially errors-in-variables)
