# GWAcalibration
Global Wind Atlas calibration to UK power

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
  - [ ] Curtailment
    - [ ] Seasonality
    - [ ] Potential yield
    - [ ] Double check curtailment days 
    - [ ] Curtailment map
    - [ ] Wind farms with most curtailment
  - [ ] Power curve
    - [ ] Scatter for select BMUs (curtailment adjusted)
    - [ ] Comparison with scaled generic curve
    - [ ] Curve matching

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
  - [ ] Data for previous years


- **Power curve fitting**
  - Using generic curves
    - [x] IEC class
    - [ ] allocate n_turb_bmu based on elexon capacity
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
    - [ ] find power curves for each model  
  - [ ] Using historical power data  
  - [ ] Bayesian approach (PC + data)

- **Loss factors**
  - Estimate using history

- **Linear regression**

- **Quantile mapping**

- **Quantile regression**

