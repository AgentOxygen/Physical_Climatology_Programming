# 2MIP
### 2 Dimensional-Orographic-Model-Inter-Comparison-Project

## Overview
The goal os 2MIP is to produce an ensemble of toy models for simulation orographic lift of an arbitrary air parcel over a mountain. Variation is captured in differences between the equations by each model and the implementation of these equations into the code. 
All models are initialized with the same starting conditions run over identical topography. The code is divided into three Python files: 
1. `orographic_models.py` - Contains individual model time step functions (all update processes for each model are bundeled into a simple I/O function that takes old state variables as input and produce new, updated state variables as output).
2. `equations.py` - Contains the equations called by the individual models. Some models use different equations, but they all rely on a few common equations.
3. `run_ensemble.py` - Executable script for producing the results of running all models indepedently.
