# RSEModels
Functions to compare race and pooling models of the redundant signal effect
Data should always be arranged in a four-column format: participant number, auditory condition, visual condition, audiovisual condition unless otherwise noted.
Empirical Data
rseviolation - Calculates the violation for each subject and returns a table
•	Utilized sampleDown, getGrice, getViolation - written by Otto (2018)

rsegain - Calculates the empirical gain for each subject and returns a table
•	Utilized sampleDown, getGrice, getCP, GetGain – written by Otto (2018)

rtcdfs – creates CDF plots for all participants

subject – pulls the 3 conditions for an individual subject

block – Returns a table with the average mean, SD and gain of each experimental block
•	Utilized getGain, sampleDown, - Written by Otto (2018)

Modelling
Pooling
poolfittting – Fits the pooling model to the data for one set of starting values
•	Utilized functions to randomly select start values and calculate the model quantiles – written by Zeheitleitner et al. (2015)
•	Utilized exwaldcdf and fminsearchbnd function – written by Palmer et al (2012)

fitexwald – Fits the pooling model to the data for multiple start values

poolmodelmean – Input must be a table with the participants organized vertically and the best fitting parameter values organized horizontally. Returns the mean, SD and gain for each participant
•	Utilized sampleDown and getGain – written by Otto (2018)
•	Utilized exwaldcdf and fminsearchbnd function – written by Palmer et al (2012)

Race
raceq – Fits the full six parameter race model to the data
•	Utilized function to calculate the model quantiles – written by Zeheitleitner et al. (2015)
•	Utilized maxNormCDF – Written by Otto (2018)
•	Utilized fminsearchbnd written by Palmer et al (2012)

raceqeta – Fits the five-parameter race model with additional parameter eta
•	Utilized function to calculate the model quantiles – written by Zeheitleitner et al. (2015)
•	Utilized maxNormCDF – Written by Otto (2018)
•	Utilized fminsearchbnd written by Palmer et al (2012)

raceqrho - Fits the five-parameter race model with additional parameter eho
•	Utilized function to calculate the model quantiles – written by Zeheitleitner et al. (2015)
•	Utilized maxNormCDF – Written by Otto (2018)
•	Utilized fminsearchbnd written by Palmer et al (2012)

raceqraab – Fits the independent race model
•	Utilized function to calculate the model quantiles – written by Zeheitleitner et al. (2015)
•	Utilized maxNormCDF – Written by Otto (2018)
•	Utilized fminsearchbnd written by Palmer et al (2012)

racemodelmean - Input must be a table with the participants organized vertically and the best fitting parameter values organized horizontally. Returns the mean, SD and gain for each participant
•	Utilized laterCDF, raceCDF, getGain and sampleDown – written by Otto (2018)

Other
 comparesingle – Input is a single column vector of RTs from a single signal condition. Outputs the fitting results for the ex-Wald model versus the LATER model
•	Utilized exwaldcdf and fminsearchbnd written by Palmer et al (2012)
