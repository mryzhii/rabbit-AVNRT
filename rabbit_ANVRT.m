% Rabbit AVNRT  1-3 Variant LADDER DIAGRAMS
% 
clear all
clear all
warning ('off','all');
global tN; tN = 1;
Figure = 1;  
Tictoc = 1;

model_name = 'avn_data_avnrt_v1';

 %  81 = Initiation of AVNRT with PAC and termination by vagal tone
 %  82 = Initiation of AVNRT with PJC and termination by vagal tone
 %  83 = Initiation of AVNRT with PVC and termination by vagal tone
Mode = 81; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = str2func(model_name);      
[Str,ap]=f();
 
nS1 = 9; % Number of S1 pulses  
nS2 = 1;
nS = nS1+nS2;
PathwayS1 = 0; 
PathwayS2 = 0;
V_cut = 0.25;       % Cut waveforms above this voltage
ph_delayS1 = 0.360;
ph_delayS2 = 0.360;
min_interval(1:3) = 0.030;  % Minimal interval between points
max_interval(1:2) = 0.120;  % 

NFP_cells = Str.iPB-Str.iAM3-1;
NSP_cells = Str.iSlow2-Str.iSlow1+1;
DN_cells = NSP_cells-NFP_cells;

tmp_distrib = NaN(1,[]);
Pathway11 = NaN(1,nS1+nS2+1); 
Pathway22 = NaN(1,nS1+nS2+1);

dx1_12 = ap.dx(Str.Fcut_cellF,1);  % Store FP coupling at ablation point (11th cell)
dx2_12 = ap.dx(Str.Fcut_cellS,2);  % Store SP coupling at ablation point (11th cell)

Last_ANS_point = Str.iPB;  %
Basic_SN_rate = ap.b(Str.iSN,1); 

%%%%%%%%%%%%%%%%%%%%%%%%
%%%ap.dx(Str.Fcut_cellF,1)=0;  % Ablation Fast PW at 12th cell
%%%ap.dx(Str.Fcut_cellS,2)=0;  % Ablation Slow PW at 12th cell
%%%%%%%%%%%%%%%%%%%%%%%%

add_name = '_Control';
if ap.dx(Str.Fcut_cellF,1)==0, add_name = '_FPabl'; end 
if ap.dx(Str.Fcut_cellS,2)==0, add_name = '_SPabl'; end 
if (Mode==99)
    add_name = [add_name '_Retro']; 
elseif (Mode==81)
    add_name = '_PAC';
elseif (Mode==82)
    add_name = '_PJC';
elseif (Mode==83)
    add_name = '_PVC';
end
fprintf("== Set %s, Mode=%d %s ==\n",model_name(1:end),Mode,add_name(2:end));

switch Mode
case {8}   %;  
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

    elseif (Mode==82) % PJC
        ANSdata.time =  ANS_time.PJC;
        ANSdata.gamma = ANS_gamma.PJC;
        tmp_distrib =   Stim_time.PJC;

    elseif (Mode==83) % PVC
        ANSdata.time =  ANS_time.PVC;
        ANSdata.gamma = ANS_gamma.PVC;
        tmp_distrib =   Stim_time.PVC;

    end
    ap.var = ANSdata.gamma(1);
    ap.mu2(1:Last_ANS_point,:) =  ap.mu2(1:Last_ANS_point,:).*ap.var;
    ap.mu1(1:Last_ANS_point,:) =  ap.mu1(1:Last_ANS_point,:)./ap.var;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 9  %%%%%%%%%%%%%%% Auto S1S2S3 Anterograde  <<<<<<<<<
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    fname1 = [model_name '_Mode9H' add_name '.txt'];
    fileID1=fopen(fname1,"w");
    fprintf(fileID1,"S1S2_[ms] H1H2_[ms] H1A2_[ms] A2H2_[ms] %0.0f %0.0f \n",ph_delayS1*1000,ph_delayS2*1000);
    fclose(fileID1);

    t_shift = 0.0; % % May be negative!
    dx = 0;
    if(ap.dx(Str.Fcut_cellS,2)==0), dx = 0.001; end
    if ap.dx(Str.Fcut_cellS,2) ~= 0
    %     ph_delayS3 = [0.094 0.095 0.096  0.098 0.100 0.105 0.11 0.115 ...
    %         0.120 0.125 0.128 0.130 0.135 0.14  ...
    %         0.15 0.16  0.18   0.20 0.30 0.36];
        ph_delayS3 = [0.094 0.095 0.0955 0.096 0.097 0.098 0.099 0.100 0.102 0.105 0.110 0.115 ...
                 0.120 0.125 0.127 0.128 0.130 0.135  0.140 0.150 0.160 0.170 ...
                 0.18 0.19 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34 0.36];
    else
      %   ph_delayS3 = [0.124  0.125 0.126 0.127 ...
      %     0.128-dx 0.129-dx 0.132-dx 0.135 0.140  0.15 0.16  ...
      %    0.18 0.20 0.30 0.36];
         ph_delayS3 = [0.124  0.125 0.1255 0.126 0.127 ...
           0.128 0.129 0.132 0.135 0.140 0.145 0.15 0.16 0.17   ...
           0.18 0.19 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34 0.36];
    end
    ph_delayS3(:)=ph_delayS3(:)+dx; 
    Pathway_array = NaN(1,length(ph_delayS3));
    HH_array = NaN(1,length(ph_delayS3));
    HA_array = NaN(1,length(ph_delayS3));
    AH_array = NaN(1,length(ph_delayS3));
    tmp_distrib = NaN(1,[]);
    %fprintf('Mode9: S2=%0.1f S3=%0.1f\n', ph_delayS2*1.e3,ph_delayS3*1.e3);

    tt = 0;
    for i=0:nS
         if i <= nS1-1
            tt = t_shift+i*ph_delayS1;
         elseif i<= nS-1
            tt = t_shift+nS1*ph_delayS1 +(i-nS1)*ph_delayS2;
        else 
            tt = tt+ph_delayS3(1);
        end
        tmp_distrib(i+1) = tt; 
    end

    End_time=  tmp_distrib(end) + 1.2;
    plot_st_time   = t_shift; % + tmp_distrib(end-1) -0.040; 
    plot_end_time  = End_time -0.1;  
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 case 99  %%%%%%%%%%%% Auto Retrograde
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    fname1 = [model_name '_Mode99H' add_name '.txt'];
    fileID1=fopen(fname1,"w");
    fprintf(fileID1,"S1S2_[ms] H1H2_[ms] A1H2_[ms] A2H2_[ms] %0.0f %0.0f \n",ph_delayS1*1000,ph_delayS2*1000);
    fclose(fileID1);

    dx = 0;
    if(ap.dx(Str.Fcut_cellS,2)==0), dx = 0.002; end
        ph_delayS3 = [0.155 0.16 0.165 0.17 0.18 0.19 ...
                  0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34 0.36]; 

    ph_delayS3(1:end-1)=ph_delayS3(1:end-1)+dx; 
    Pathway_array = NaN(1,length(ph_delayS3));
    HH_array = NaN(1,length(ph_delayS3));
    HA_array = NaN(1,length(ph_delayS3));
    AH_array = NaN(1,length(ph_delayS3));
    tmp_distrib = NaN(1,[]);
    tt = 0;
    for i=0:nS
         if i <= nS1-1
            tt = t_shift+i*ph_delayS1;
         elseif i<= nS-1
            tt = t_shift+nS1*ph_delayS1 +(i-nS1)*ph_delayS2;
        else 
            tt = tt+ph_delayS3(1);
        end
        tmp_distrib(i+1) = tt; 
    end

    End_time = tmp_distrib(end)+2*max(ph_delayS1,ph_delayS2(end))+0.5;  % Total simulation time [s]
    plot_st_time   = t_shift; % + tmp_distrib(end-1) -0.040;    
    plot_end_time  = End_time;   

otherwise   
    fprintf('Wrong case number\n');
    return
end

%%%%%%%%%%%%%%%%%%%%%%% END CASE MODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_step = 0.1*1.0e-3;  %   0.005*1.0e-3; <<<<<<<<<<<<<<<<<<<<<<<<<<<
opts = odeset('Reltol',1e-6,'AbsTol',1e-8, 'Stats','off', ...
    'NormControl','on','InitialStep',1e-6); 
initcond = 0.0001.*ones(Str.iHB2*4,1);
tspan = [0 t_step End_time];
t1 = 0:t_step:End_time; 

NN = 1;
if (Mode==9||Mode==99), NN = length(ph_delayS3); end

%%%%%%%% LOOP %%%%%%%%%%%%%%%%%%
for nn = 1:NN

if (Mode==9||Mode==99)
    tmp_distrib(end) = tmp_distrib(end-1) + ph_delayS3(nn);
    tN =1;
end


if Tictoc, tic; end
fprintf('= Solving ODEs ...\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (Mode==81 || Mode==82 || Mode==83) % PAC, PJC, PVC
sol = ode23(@(t,y) avn_rabbit_fun_ladder(t,y,Mode,Str,ap,tmp_distrib, ...
                   ANSdata,Last_ANS_point), ...
                   tspan, initcond, opts);
else
sol = ode23(@(t,y) avn_rabbit_fun(t,y,Mode,nn,Str,ap,tmp_distrib, ...
                   term_time0,term_time1), ...
                   tspan, initcond, opts);
end
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

% Only for images of S1S2 stimulation pulses 
if (Mode==8 ||Mode==81||Mode==82|| Mode==9)   
    ii = Str.imp_Antero;
elseif (Mode==83|| Mode==88||Mode==99)
    ii = Str.imp_Retro;   
end 
 
for i=1:length(t1)
    for j=1:length(tmp_distrib)  
        if(t1(i) >=  tmp_distrib(j) && t1(i) < tmp_distrib(j) +ap.pulse_width*2)
              y2(ii,i) = 0.8;
        end
    end
end    


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------- Check conduction times ------------
iTime.pointA  = 0;
iTime.pointB  = 0;
PathwayS2 =0;
PathwayS3 = 0;

iend = min(size(locs_F,2),size(locs2,2))-1;
if(Mode==8||Mode==9)
      for i=2:size(locs_F,2)
          if abs(locs_F(Str.imp_Antero,i) - tmp_distrib(end)) < min_interval
              iend = i;
              break;
          end
      end
elseif (Mode==88||Mode==99)
      for i=2:min(size(locs_F,2),size(locs2,2))
          if abs(locs_F(Str.imp_Retro,i) - tmp_distrib(end)) < min_interval
              iend = i;
              break;
          end
      end
end


iTime.pointA  = 1.e3*(locs_F(Str.pointA,iend)  - locs_F(Str.iSN,iend) );
iTime.PB  = 1.e3*(locs_F(Str.iPB,iend)  - locs_F(Str.iSN,iend) );
iTime.pointB   = 1.e3*(locs_F(Str.pointB,iend)   - locs_F(Str.iSN,iend) );
iTime.HB2  = 1.e3*(locs_F(Str.iHB2,iend)  - locs_F(Str.iSN,iend) );
iTime.Slow_max = 1.e3*(max(locs2(Str.iSlow1:Str.iSlow2,iend)) - locs_F(Str.iSN,iend) );

PathwayS1 = 1.e3*(locs_F(Str.pointB,iend-2) - locs_F(Str.pointA,iend-2));
PathwayS2 = 1.e3*(locs_F(Str.pointB,iend-1) - locs_F(Str.pointA,iend-1));
PathwayS3 = 1.e3*(locs_F(Str.pointB,iend)   - locs_F(Str.pointA,iend));
Delta = PathwayS3-PathwayS2;
if (Mode>82) 
      PathwayS1=-PathwayS1; PathwayS2=-PathwayS2; PathwayS3=-PathwayS3; 
      Delta=-Delta;
      if(Delta<-2.0), PathwayS3=NaN; end
end
if (PathwayS3>200 || PathwayS3<40), PathwayS3=NaN; Delta=NaN; end
  

%if Tictoc, tic; end
%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avn_plot3(nn,Mode,Str,ap,al,bl,t1,y1,y2, ...
    locs_F,locs_S,locs_indexS, min_interval,max_interval, ...
    ANSdata,tmp_distrib,iTime, ...
    nS,plot_st_time,plot_end_time);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if Tictoc, toc; end


if(Mode==9 || Mode == 99 || Mode==8 || Mode==88 )  % H1H2, H1A2, A2H2
    if ~isnan(PathwayS2) && PathwayS1~=0 && ~isnan(PathwayS3)
        fileID1 = 1; 
        if(Mode~=8 && Mode~=9 )
            fileID1=fopen(fname2,"a");
        end
        d = length(tmp_distrib);
        if(d<2) d=2; end
        ta1 = locs_F(Str.pointA,d-1);   
        ta2 = locs_F(Str.pointA,d);    
        th1 = locs_F(Str.pointB,d-1);   
        th2 = locs_F(Str.pointB,d); 
   
        Pathway_array(nn)=PathwayS3;
        HH_array(nn)=1.0e3*abs(th2-th1); 
        HA_array(nn)=1.0e3*abs(ta2-th1);
        AH_array(nn)=1.0e3*abs(th2-ta2);
        if(Mode==99)
           HA_array(nn)=1.0e3*abs(ta1-th2);
        end
        if(Mode~=8 && Mode~=9 )
            fclose(fileID1);
        end
    end
end
end %%%%%%%%%%% Loop for Mode

if(Mode==9 || Mode == 99)
    figure(9);
    plot(ph_delayS3*1.0e3,Pathway_array,'-o');
    axis tight
    figure(10);
    hold on
    plot(ph_delayS3*1.0e3,HH_array,'-o');
    plot(ph_delayS3*1.0e3,HA_array,'-*');
    plot(ph_delayS3*1.0e3,AH_array,'-s');
    axis tight
    hold off
end
