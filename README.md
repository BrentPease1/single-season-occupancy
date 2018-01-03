# Single-species, Single-season Occupancy Model using package `unmarked`

---

#### This is a repository to conduct single-species, single-season occupancy models using maximum likelihood methods in package `unmarked`. 

---

### Repository contents:  
  
  **1.** An R script `single-season_occupancy.R` for conducting the analyses  
  **2.** Example datasets. Camera trap data were originally downloaded from eMammal's [browse data tab]((http://emammal.si.edu/analysis/data-download)).  
     **2a** `deer_detection_history.csv` is real camera trap data from 87 sites that operated for ~3 weeks. Continous data were binned into daily intervals ("visits")  
     **2b** `detection_covariates.csv` is a file containing two simulated detection covariates.  
  **3.** An R script `correlation_pvalue.R` for visualizing correlation among covariates. The script is sourced automatically within the `single-season_occupancy` script, thus it is unlikely this file will need to be used unless you would like to make edits to the output. 