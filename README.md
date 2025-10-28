# GWAcalibration
Global Wind Atlas calibration to UK power

## Next Steps

- **Figures and tables**
  - [x] Turbine distribution 
    - height
    - capacity
    - number of turbines per site
  - [x] turbine height by type
  - [ ] elexon time variability
  - [ ] elexon panel data characteristics
  - [ ] overall evolution generation

- **Inputs and EDA**
  - [x] cleaning data code
  - [x] Turbine characteristics 
  - [x] Get generic power curves
  - [x] Wind speed vertical interpolation 
  - [x] function to extract GWA wind speed
  - [x] add GWA wind speed at multiple heights
  - [x] interpolate wind speed to hub height
  - [ ] Get land elevation 
  - [ ] Get curtailment *

- **Power curve fitting**
  - Using generic curves
    - [x] IEC class
    - [x] function to rescale pc
    - [x] add zeroes at the end of power curve
    - map curves
      - [x] fill in capacity_turb
      - [x] n_turb per BMU
      - [x] hub height imputation
    - [x] convert GWA to power with mean
    - [x] convert GWA to power using weibull
    - [x] compare with generation
  - Using turbine characteristics  
  - Using historical power data  

- **Linear regression**

- **Quantile mapping**

- **Quantile regression**

