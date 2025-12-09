# Calibration of Reanalysis to GB wind power generation
Started as a Global Wind Atlas calibration to UK power, however to keep working with spatiotemporal properties now we're calibrating ERA5.

## Next Steps

- **Figures and tables**
  - [x] Turbine distribution 
    - height
    - capacity
    - number of turbines per site
  - [x] turbine height by type
  - [x] Comparison with interactive map
  - [ ] elexon time variability
  - [ ] elexon panel data characteristics
  - [ ] overall evolution generation
  - [x] Curtailment
    - [x] Seasonality
    - [x] Potential yield
    - [x] Double check curtailment days 
    - [x] Curtailment map
    - [x] Wind farms with most curtailment
  - [x] Revised capacity factors
  - [x] Wind speed seasonality
  - [ ] Power curve *
    - [x] Scatter for select BMUs (curtailment adjusted)
    - [x] Comparison with scaled generic curve
    - [ ] Curve matching
  - [x] outages statistics
    - [x] number of events
    - [x] duration
    - [x] capacity impact
    - [x] spatial distribution
  - [x] Power curve scatter by season
  - [x] Full history plots
  - [ ] Time series of a segment before, after calibration

- **Inputs and EDA**
  - [x] cleaning data code
  - [x] Turbine characteristics 
  - [x] Get generic power curves
  - [x] Wind speed vertical interpolation 
  - [x] function to extract GWA wind speed
  - [x] add GWA wind speed at multiple heights
  - [x] interpolate wind speed to hub height
  - [x] Elexon capacity factors check
  - [x] ERA 5 time series for wind farm locations
  - [ ] Question to Dan: GWA to MIDAS (Onshore) then *
  - [ ] Get land elevation 
  - [x] Get curtailment (DISPTAV)
  - [x] Data for previous years
  - [x] prepare wind speed for power conversion *
  - [x] Data capping at 0 and max capacity instead of removing *
  - [x] Combine BMUs with same coordinates *
  - [x] remit messages download
  - [x] unwind remit outage profile
  - [x] merge with clean generation data
  - [x] wind vertical interpolation


- **Power curve fitting**
  - Using generic curves
    - [x] IEC class
    - [ ] allocate n_turb_bmu based on revised capacity
    - [x] function to rescale pc
    - [x] add zeroes at the end of power curve
    - map curves
      - [x] fill in capacity_turb
      - [x] n_turb per BMU
      - [x] hub height imputation
    - [x] convert GWA to power with mean
    - [x] convert GWA to power using weibull
    - [x] compare with generation
  - [ ] Using turbine characteristics
    - [ ] find power curves for some models
    - [x] interpolate wind speed to turbine height  
  - [ ] Using historical power data  
  - [ ] Bayesian approach (PC + data)

- **Slides**
  - [x] New set with era 5
  - [x] Power curve updates 
  - [x] Curtailment and Maintenance
  - [x] New comparison
  - [x] lm, 
  - [ ] QM, QR
  - [ ] Summary Research review
    - [x] Power curve
    - [x] Cannon/Brayshaw
    - [x] Potisomporn relevant
    - [ ] Thornton relevant

- **Research topics**
  -[ ] power curve estimation
  -[ ] calibration of power curve estimates

- **Loss factors**
  - Estimate using history

- **Linear regression**
  - [x]
- **Quantile mapping**
  - [ ] Uncertainty around bias correction
  - [ ] Varying parameters through time or space

- **Quantile regression**
- **Spatial properties**
  - [x] Correlation decay
  - [x] Non-stationary spatial patterns

- **Model validation**
  - [ ] out-of-sample validation
  - [ ] UQ and error metrics
  - [ ] Extreme value properties?
  - [ ] Sequences of low wind power

- **Expansion ideas**
  - Model QM through hierarchical model
  - Model vertical interp. / PC / calibration jointly
  - Add spatial component
  - add loss due to outages model
  - Estimate generation from locations without elexon data
  - UQ of Power estimates for locations without elexon data
  - WTG efficiency loss: type / time / weather / errors / maintenance events
  - Extreme draughts
  - Extreme low wind speed scenarios
  - Extreme low wind power scenarios
  - Solar energy calibration
  - Renewable energy defict scenarios



