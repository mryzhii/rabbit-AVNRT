# rabbit-AVNRT
This folder contains MATLAB source code for AVNRT simulation from the paper:
"Atrioventricular nodal reentrant tachycardia onset, 
 sustainability, and spontaneous termination in rabbit atrioventricular 
 node model with autonomic nervous system control"
Front. Physiol. Volume 15 - 2024 | doi: 10.3389/fphys.2024.1529426<br>

by M.Ryzhii (University of Aizu) and 
E.Ryzhii (Fukushima Medical University).<br>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tested with MATLAB R2023b

rabbit_AVNRT.m - main file
              Set 'Mode' to the corresponding value (see inside)<br>

avn_data_avnrt_v1.m      - rabbit AVN model structure and parameters (first variant)<br>
avn_data_avnrt_v1_scen.m - pacing scenarios for the first model variant<br>
avn_rabbit_fun_avnrt.m   - function for MATLAB ODE solver<br>
avn_plot_avnrt.m         - analyses output and plots ladder diagrams<br>
To start simulations, run the main rabbit_ANVRT.m file. <br>
