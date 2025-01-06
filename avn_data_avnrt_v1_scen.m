%%% AVNRT SCENARIO for VARIANT 2 
% AVNRT in 1D multi-functional model of rabbit AV node with dual pathways
%     (c) 2024
%     Maxim Ryzhii, University of Aizu, Japan
%     Elena Ryzhii, Fukushima Medical University, Japan
%
% Code for the paper "Atrioventricular nodal reentrant tachycardia onset, 
% sustainability, and spontaneous termination in rabbit atrioventricular 
% node model with autonomic nervous system control"
% 
% ANS_time:  Moments of ANS state (gamma) changing
% ANS_gamma: Value of gamma corresponding to ANS_time
% Stim_time: Stimulation moments 
%
% Tested with MATLAB R2023b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of the rabbit AVN model with ANS control

function [ANS_time,ANS_gamma,Stim_time]= avn_data_avnrt_v1_scen()

stack = dbstack('-completenames');
if numel(stack) <2
   fprintf('Wrong call. Start with ''rabbit_AVNRT''\n');
   return
end


ANS_time.PAC =  [  0   2.0  3.0  4.0  5.0  6.0  7.0  8.0 ];
ANS_gamma.PAC = [1.5   0.8  0.9  1.1  1.2  1.3  1.5  1.7 ];
Stim_time.PAC = [0.900 2.620];

ANS_time.PJC =   [  0  1.5 2.5  3.0  4.0  5.0  6.0  7.0  8.0 ];
ANS_gamma.PJC =  [1.5  .65 .57  0.7  1.0  1.2  1.3  1.5  1.7 ];
Stim_time.PJC =  [11 13];

ANS_time.PVC =   [  0   1.5   2.5   3.5   5.5   7.5 ];
ANS_gamma.PVC =  [0.8   1.0   1.2   1.5   1.7   1.8 ];
Stim_time.PVC =  [1.000 1.970 2.930 4.440 6.140 7.930];

