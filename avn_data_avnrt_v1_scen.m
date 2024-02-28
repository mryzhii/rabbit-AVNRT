%%% VARIANT 1 SCENARIO
% Parameters of the rabbit AVN model with ANS control
function [ANS_time,ANS_gamma,Stim_time]= avn_data_avnrt_v1_scen()

stack = dbstack('-completenames');
if numel(stack) <2
   fprintf('Wrong call. Start with ''rabbit_AVNRT''\n');
   return
end


ANS_time.PAC =  [ 0    2.0  3.0  4.0  5.0  6.0  7.0  8.0 ]; % PAC
ANS_gamma.PAC = [1.5  0.8  0.9  1.1  1.2  1.3  1.5  1.7 ];
Stim_time.PAC = [0.900 2.620];

ANS_time.PJC =   [0    1.5 2.5  3.0  4.0  5.0  6.0  7.0  8.0 ]; % PJC
ANS_gamma.PJC =  [1.5  .65 .57  0.7  1.0  1.2  1.3  1.5  1.7 ];
Stim_time.PJC =  [11 13];

ANS_time.PVC =   [0    2.0  3.5  5.0  6.5  7.5  8.5 ];  % PVC
ANS_gamma.PVC =  [2.0  1.7  1.5  1.3  1.1  1.0  0.8 ];
Stim_time.PVC =  [1.120 2.643 4.390 5.977 7.121 8.112 8.988];

