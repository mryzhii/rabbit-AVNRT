%%% VARIANT 1 
% Parameters of the rabbit AVN model with ANS control

function [Str,ap]= avn_data_avnrt_v1()

stack = dbstack('-completenames');
if numel(stack) <2
   fprintf('Wrong call. Start with ''rabbit_AVNRT''\n');
   return
end
% Indexes of the model cells
Str.iSN    = 1;  
Str.iAM3   = 7;  % Last atrium 
Str.iSlow1 = 7;  % First element of SP
Str.ipSlow = 13; % First SP pacemaker 
Str.iSlow2 = 16; % Last element of SP 
Str.iPB    = 16; % PB Pacemaker
Str.ipHB   = 20; % 
Str.iHB2   = 22; % Last element of His
Str.pointA = Str.iAM3;  % Measurement point A
Str.pointB = Str.iPB;   % Measurement point B
Str.imp_Antero = Str.iAM3-1; % pulse application point for Antegrade
Str.imp_Retro  = Str.iPB+1;  % pulse application point for Retrograde
Str.Fcut_cellS = 11;
Str.Fcut_cellF = 11;

TC0 = 1000.0/4.0; % Time constant for physical time
TC  = zeros(Str.iHB2,2);
TC(:,1) = TC0; TC(Str.iAM3-2:Str.iPB,2) = TC0;

a_base = 0.13;
bAP = a_base.*ones(Str.iHB2,2);

bAP_normal = -0.1100; %-0.0694;
bAP(Str.iSN,1) = bAP_normal; 
bAP(Str.iSN+1,1) = 0.08;

bAP(Str.iPB-1,1)    = -0.0368; %
bAP(Str.iPB,1)      = -0.0368; %
bAP(Str.ipSlow,2)   = -0.0369; % 
bAP(Str.ipSlow+1,2) = -0.0369; % 
bAP(Str.ipSlow+2,2) = -0.0368; %
bAP(Str.ipSlow+3,2) = -0.0368; %

aa = a_base.*ones(Str.iHB2,2);
aa(Str.iSN,1) = 0.20;
aa(Str.iSN+1,1) = 0.16;

aa(Str.iPB-1,1)    = 0.28; 
aa(Str.iPB,1)      = 0.28; 
aa(Str.ipSlow,2)   = 0.28; 
aa(Str.ipSlow+1,2) = 0.28; 
aa(Str.ipSlow+2,2) = 0.28; 
aa(Str.ipSlow+3,2) = 0.28; 

mu1 = [0.16 0.18 0.20 0.20 0.22 0.22 0.22 ...                         % SN, PSN, AM
       0.15 0.15 0.15 0.15 0.15 0.14 0.14 0.065 0.060 ...             % Fast
       0.02 0.02 0.02 0.02 0.02 0.02; ...                             % HB
       0    0    0     0     0  0.22  ...                             % AM* 
       0.23 0.22 0.215 0.205 0.190  0.175 0.07 0.070 0.060 0.060 ...  % Slow
       0  0  0  0  0  0 ];
       
mu2 = [4.20 4.20 3.30 2.50 1.60 1.60  1.60 ...                 % SN, PSN, AM
       1.70 1.80 1.85 1.95 2.05 2.15  2.25 4.60  4.60 ...      % Fast
       0.35 0.10 0.10 0.10 0.10 0.10 ; ...                     % HB
       0    0    0    0    0    1.60 ...                       % AM*
       1.80 1.80 1.80 1.80 1.80 1.80  2.50 3.30 3.70 3.90 ...  % Slow
       0   0   0   0   0   0 ];

k1 = [9  9  10 14 18 18 18  ...      % SN, PSN, AM
      18 18 18 18 18 15 15  8  8 ... % Fast
      10 10 10 10 10 10;  ...        % His
      0  0  0  0  0  18  ...         % AM*
      11 11 11 11 11 11 ...          % Slow
      8  8  8  8 ...                % Slow pacemakers
      0  0  0  0  0  0 ];
k2 = k1; 

eps0 = [ 0.040 0.036 0.025 0.025 0.018 0.018 0.018 ...
         0.011 0.011 0.011 0.011 0.011 0.011 0.011  0.042 0.042 ... % Fast
         0.042 0.050 0.050 0.050 0.050 0.050; ...                   % HB
         0     0     0     0     0     0.018  ...                   % AM*
         0.060 0.060 0.060 0.060 0.060 0.060 ...                    % Slow
         0.044 0.044   0.042  0.042 ...                             % Slow
         0  0  0  0  0  0 ];

% ----------- COUPLING
dx = [ 43 47 49 56 68 95 ...   
       95 84 70 64 62 60 59 59 55 ... 
       80 50 50 70 90 100; ...
       0  0  0  0  0  0 ...
       59 57 57 56 56   ...
       53 32 32 32 ...
       0  0  0  0  0  0 ];

dy = [   0  0  0  0  0  0  69 ...       % AM1, Slow
         0  0  0  0  0  0  0  0 42 ...
         0  0  0  0  0  0];

adx = [1    1    1    1    0.75 0.53 ...
       0.62 0.75 0.80 0.85 0.95 1.00 1.25 1.40 1.50 ...  
       0.65 0.70 1    1    1    1;
       0    0    0    0    0    1 ...
       1    1    1    1    1    1    1    1    1 ... 
       1    0    0    0    0    0];

ady = [  0 0 0 0 0 1 1 ...         
         0 0 0 0 0 0 0 0 1 ... 
         0 0 0 0 0 0];

drnd = [1;  0.075; 0.125]; % Parameters for random generator

ap.pulse_width = 0.001;  % [s]
ap.pulse_amp =  -280.0; 
ap.Tc = TC;
ap.a = aa; 
ap.b = bAP;
ap.mu1 = mu1';
ap.mu2 = mu2';
ap.k1 = k1';
ap.k2 = k2';
ap.eps0 = eps0';
ap.dx = dx';
ap.dy = dy'; 
ap.adx = adx'; 
ap.ady = ady';
ap.drnd  = drnd;

