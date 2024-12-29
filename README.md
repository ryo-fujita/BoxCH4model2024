# BoxCH4model2024
Simple atmospheric one-box model for CH₄ and its isotopes with particle filter used in Fujita et al. 2025, JGR-Atmosphere.

Coded in Matlab (Tested in versions 2023a)

This repository contains the main programs and text files necessary for running the model and generating main outputs.

## Procedure to Run the Model
1. Install  
Dowload ZIP of this repository and unzip it. 

2. Preparation  
Ensure all necessary input files are prepared and placed in the correct directory.

3. Run the model  
You can run the model in two different modes, depending on your requirements:
  - Basic Mode (Without Plots and Smoothing):  
Use this mode to quickly check the code functionality without additional processing.
  Run the following command in the MATLAB command window:  
    `Model().main`
  - Full Mode (With All Processes):  
Use this mode for a complete and proper run of the model, including plotting and smoothing processes.
Run the following command in the MATLAB command window:  
    `RunModel().main`

4. Check the output  
The output files will be generated in the current directory.


## File configuration
### Input Files
The following input files are required to run the model:
1. `./ModelInputParameter.info`  
Parameter file that controls the number of ensemble members and specifies key inputs.
2. `./parameter/targetAtmosPara.txt`   
Target file of atmospheric observations and CH₄ source and sink parameters  (see "README" included in the original Excel sheet).
3. `./emission/*.txt`  
A priori CH₄ emission scenario files.  
4. `./loss/*.txt`  
A priori CH₄ loss scenario files.

### Key Output Files  
The model generates the following key output files:
1. `objects_SmoothSummaryPostFull_atmos.mat`  
Timeserires of simulated atmospheric CH₄ and isotopic histories based on the posterior scenario.
2. `objects_SmoothSummaryPostFull_para.mat`  
Timeserires of posterior CH₄ source and sink parameters.
3. `outStatics_*.txt`  
Statistical summaries of posterior parameters, emissions, source fractions, and losses.

### Key Programs  
Below is a brief description of the key MATLAB programs:
1. `RunModel.m`  
- The primary program for controlling model, smoothing processes, and generating outputs.
- Can be executed independently (e.g., `RunModel().main`).
- Also called by `RunScript.m` when batch processing is utilized.

2. `Model.m`  
- The core program called by RunModel.m.
- Contains basic properties and functions for the atmospheric box model.
- Can also be executed independently (e.g., `Model().main`).

3. `Atmosphere.m`  
Manages atmospheric CH₄ concentration, δ¹³C, δD, and Δ¹⁴C.

4. `ParameterSourceSink.m`  
Handles calculations for CH₄ source and sink parameters.

5. `Emission.m`  
Manages emissions for total CH₄, ¹³CH₄, CH₃D, and ¹⁴CH₄.

6. `Isotopes.m`  
Calculates isotopic emissions for ¹³CH₄, CH₃D, and ¹⁴CH₄.

7. `Loss.m`  
Calculates loss rates for total CH₄, ¹³CH₄, CH₃D, and ¹⁴CH₄.

8. `JudgeAtmosToAtmosTarget.m`  
Implements particle filter code for filtering. It compares ensembles of simulated atmospheric CH₄, δ¹³C, δD, and Δ¹⁴C with observational targets and resamples the outputs based on observation likelihoods.

9. `Smoothing.m`  
Implements particle filter logic for smoothing. It extracts timeseries of posterior distributions for atmospheric components or parameters from filtered distributions.

10. `MyCalc.m`  
Aggregates various data processing functions.

11. `Plot.m`  
An abstract superclass that consolidates properties and functions for figure processing. 
