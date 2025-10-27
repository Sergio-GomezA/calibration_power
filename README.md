# GWAcalibration
Global Wind Atlas calibration to UK power

## Next Steps

- **Figures and tables**
  - [x] Turbine distribution 
    - height
    - capacity
    - number of turbines per site
  - [ ] turbine height by type
  - [ ] elexon time variability
  - [ ] elexon panel data characteristics
  - [ ] overall evolution generation

- **Inputs and EDA**
  - [ ] Turbine characteristics 
  - [x] Get generic power curves
  - [ ] Wind speed vertical interpolation 
  - [x] function to extract GWA wind speed
  - [x] add GWA wind speed at multiple heights
  - [x] interpolate wind speed to hub height
  - [ ] Get land elevation 
  - [ ] Get curtailment *

- **Power curve fitting**
  - Using generic curves
    - [x] IEC class
    - [ ] function to rescale pc
    - map curves
      - [x] fill in capacity_turb
      - [ ] n_turb per BMU
      - [x] hub height imputation
    - convert GWA to power
    - compare with generation
  - Using turbine characteristics  
  - Using historical power data  

- **Linear regression**

- **Quantile mapping**

- **Quantile regression**

