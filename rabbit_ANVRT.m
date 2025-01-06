% Rabbit AVNRT simulation with LADDER DIAGRAMS
% AVNRT in 1D multi-functional model of rabbit AV node with dual pathways
%     (c) 2024
%     Maxim Ryzhii, University of Aizu, Japan
%     Elena Ryzhii, Fukushima Medical University, Japan
%
% Code for the paper "Atrioventricular nodal reentrant tachycardia onset, 
% sustainability, and spontaneous termination in rabbit atrioventricular 
% node model with autonomic nervous system control"
% 
%
% Tested with MATLAB R2023b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
clear all
warning ('off','all');
global tN; tN = 1;
Figure = 1;  
Tictoc = 1;

%  Mode=81  => Initiation of AVNRT with PAC and termination by vagal tone
%  Mode=82  => Initiation of AVNRT with PJC and termination by vagal tone
%  Mode=83  => Initiation of AVNRT with PVC and termination by vagal tone
Mode = 83; 
Variant = 3; % 1-3

model_name = ['avn_data_avnrt_v' num2str(Variant,'%d')];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = str2func(model_name);      
[Str,ap]=f();
 
V_cut = 0.25;       % Cut waveforms above this voltage
min_interval(1:3) = 0.050;  % Minimal interval between points
max_interval(1:2) = 0.120;  % 

NFP_cells = Str.iPB-Str.iAM3-1;
NSP_cells = Str.iSlow2-Str.iSlow1+1;
DN_cells = NSP_cells-NFP_cells;

tmp_distrib = NaN(1,[]);

Last_ANS_point = Str.iPB;  %
Basic_SN_rate = ap.b(Str.iSN,1); 

ANSdata.time =  0;
ANSdata.gamma = 0;
ANSdata.var = Variant;

add_name = '_Unknown';
switch Mode
    case 81
        add_name = '_PAC';
    case 82
        add_name = '_PJC';
    case 83
        add_name = '_PVC';
end
fprintf("== DataSet: %s, Mode=%d %s ==\n",model_name(1:end),Mode,add_name(2:end));


switch Mode
%%%%%%%%%%%%%%%%%%%%%%%%%% 81 & 82 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {81,82,83} % 81=PAC;  82=PJC;  83=PVC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    min_interval(2) = 0.050; 
    min_interval(3) = 0.050;
    plot_st_time   = 0.5;
    plot_end_time  = 9.5; % 
    End_time = plot_end_time +0.3;

    scen_name = [model_name(1:end) '_scen'];
    sc = str2func(scen_name); 
    [ANS_time,ANS_gamma,Stim_time]=sc(); % Read scenario (PAC, PJC, PVC)

    if (Mode==81)      % PAC
        ANSdata.time =  ANS_time.PAC;
        ANSdata.gamma = ANS_gamma.PAC;
        tmp_distrib =   Stim_time.PAC;
        ap.pulse_width = 0.001; 

    elseif (Mode==82) % PJC
        ANSdata.time =  ANS_time.PJC;
        ANSdata.gamma = ANS_gamma.PJC;
        tmp_distrib =   Stim_time.PJC;

    elseif (Mode==83) % PVC
        ANSdata.time =  ANS_time.PVC;
        ANSdata.gamma = ANS_gamma.PVC;
        tmp_distrib =   Stim_time.PVC;
        ap.pulse_width = 0.002;  

    end
    ap.var = ANSdata.gamma(1);
    ap.mu2(1:Last_ANS_point,:) =  ap.mu2(1:Last_ANS_point,:).*ap.var;
    ap.mu1(1:Last_ANS_point,:) =  ap.mu1(1:Last_ANS_point,:)./ap.var;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% END CASE MODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
otherwise   
    fprintf('Wrong case number\n');
    return
end

t_step = 0.1*1.0e-3;  
opts = odeset('Reltol',1e-6,'AbsTol',1e-8, 'Stats','off', ...
    'NormControl','on','InitialStep',1e-6); 
initcond = 0.0001.*ones(Str.iHB2*4,1);
tspan = [0 t_step End_time];
t1 = 0:t_step:End_time; 


if Tictoc, tic; end
fprintf('= Solving ODEs ...\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sol = ode23(@(t,y) avn_rabbit_fun_avnrt(t,y,Mode,Str,ap,tmp_distrib, ...
                   ANSdata,Last_ANS_point), ...
                   tspan, initcond, opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y11 = sol.y';
y11 = reshape(y11,[length(y11),Str.iHB2,4]);
if Tictoc, toc; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 y1 = 0.001.*ones(Str.iHB2,length(t1));
 y2 = 0.001.*ones(Str.iHB2,length(t1));
 for i =1:Str.iHB2
    y1(i,:) = interp1(sol.x,y11(:,i,1),t1); % u in 1st row  (FP)
    y2(i,:) = interp1(sol.x,y11(:,i,3),t1); % u in 2nd row  (SP)
 end
 y1c = zeros(Str.iHB2,length(t1)); % "Cut" arrays
 y2c = zeros(Str.iHB2,length(t1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of excitation latencies (Fast/Slow)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
locs_F = NaN(Str.iHB2,[]);
locs2  = NaN(Str.iHB2,[]);
locs_S = NaN(Str.iHB2+DN_cells,[]);
Peaks1 = NaN(Str.iHB2,1);
Peaks2 = NaN(Str.iHB2,1);
alen = NaN(Str.iHB2,1);
blen = NaN(Str.iHB2,1);

locs_indexF = 1:Str.iHB2;
locs_indexS = 1:Str.iHB2+DN_cells;

y1c(:,:) = min(V_cut,y1(:,:));    % Cut arrays, to get max V point
y2c(:,:) = min(V_cut,y2(:,:));

for i = 1:Str.iHB2
    MPH = V_cut-0.01; 
    MPP = V_cut-0.01;
    if i == Str.iSN,  MPP = V_cut*2/3; end 

    [pks,lcs] = findpeaks(y1c(i,1:end),t1(1:end),'MinPeakHeight',MPH, ...
            'MinPeakProminence',MPP,'MinPeakDistance',0.03);  %0.03);
    if ~isnan(lcs)
        Peaks1(i) = pks(end);
        alen(i) = length(lcs);
        locs_F(i,1:length(lcs)) = lcs();
    else
        % fprintf('No peaks (Fast) found\n');
        % return 
    end
    [pks2,lcs2] = findpeaks(y2c(i,1:end),t1(1:end),'MinPeakHeight',MPH, ...
            'MinPeakProminence',MPP,'MinPeakDistance',0.03); %0.03);
    if ~(lcs2)
        Peaks2(i) = pks2(end);
    %else
    %    fprintf('No peaks (Slow) found\n');
    %    return
    end
    blen(i) = length(lcs2);
    locs2(i,1:length(lcs2)) = lcs2();
end

al = max(alen);
bl = max(blen); 
locs_S(1:Str.iAM3,1:al)                = locs_F(1:Str.iAM3,1:al);
locs_S(Str.iSlow1+1:Str.iSlow2+1,1:bl) = locs2(Str.iAM3:Str.iPB,1:bl);
locs_S(Str.iSlow2+2:Str.iHB2+DN_cells,1:al)   = locs_F(Str.iPB:Str.iHB2,1:al);

for i = Str.iSlow1+1 : Str.iSlow2+2
    locs_indexS(i) = Str.iSlow1 +(i-Str.iSlow1)/(Str.iPB-Str.iAM3+DN_cells)*(Str.iPB-Str.iAM3);
end
for i = Str.iPB+1+DN_cells : Str.iHB2+DN_cells
    locs_indexS(i) = locs_indexF(i-DN_cells);
end


%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avn_plot_avnrt(Mode,Str,ap,al,bl,t1,y1,y2, ...
               locs_F,locs_S,locs_indexS, min_interval,max_interval, ...
               ANSdata,tmp_distrib,plot_st_time,plot_end_time);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


