
File for analysis:

fulldata_k.csv   - here k refers to the number of QUS acquisitions that were averaged
This file is created by running the preprocessing scripts: 

1. "preprocess-qus.qmd" (set the number of acquisitions to average in this script)
2. "preprocess-clinical.qmd" - merge output file form the first script with "fulldata.csv" and write the file "fulldata_k.csv".


This data set includes the clinical information for each participant, the birth outcome (in two fields), and the QUS averaged measurements.

For birth outcome:

sPTB = Yes  -> spontaneous preterm birth
sPTB = No  -> not a spontaneous preterm birth
mPTB = Yes -> medical preterm birth
mPTB = No -> not a medical preterm birth

We will analyze only the data with mPTB=No. Within this subset the response of interest is sPTB (Yes measn a preterm birth, No means a full term birth).

Other variables:
 please see the reference paper McFarlin et al. 2024, "Enhanced identification...cohort study," for descriptions of the variables. The variable names have changed slightly here. 

Conventions:

GAV1 - gestational age at research visit 1 (around 20 weeks)
GAV2 - gestational age at research visit 2 (around 24 weeks)
CLV1 - cervical length at research visit 1
CLV2 - cervical length at research visit 2

AC_V1, AC_V2 - average attenuation coefficient at visits 1 and 2 respectively
LF_Intercept_V1, LF_Intercept_V2 - Lizzi-Feleppa Intercept at visits 1 and 2 respectively
...
SWS_V1, SWS_V2 - shearwave speed at visits 1 and 2 respectively




Other files:

- fulldata.csv - original data with clinical and outcome data along with first acquisition (out of 10) QUS data. We use this to get the clinical data.


- UIC_QUS_8.8.22_C425updates.xlsx - excel file with all QUS acquisitions using two different methods. 
We are using data from the "per-image (SLD) sheet.  See the QMD file "preprocess-qus.qmd" for how we process these data.

- QUS_means_k.csv  - intermediate file created by the preprocess-qus script and read by the preprocess-clinical script.




