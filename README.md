
To run the codes, users will need MATLAB (R2021a or later), MATPOWER (version 7.1 or later), and MATLAB's Deep Learning Toolbox.
For the electrical power flow simulation, open and run the script Power_15nodes.m, which uses the case15nbr.m file to simulate an AC power flow on a 15-bus network. To simulate the gas network, run TaintrustGasModel.m. This program calculates nodal pressures, gas flow through each branch, compressor performance, and power equivalents at coupling points using the gas calorific value
For ML, the scripts TestMLCHP.m and TestMLP2G.m can be used
