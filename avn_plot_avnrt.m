% PLOTTING LADDER DIAGRAMS
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
function avn_plot_avnrt(Mode,S,ap,al,bl,t1,y1,y2, ...
    locs_F,locs_S,locs_indexS,min_interval,max_interval, ...
    ANSdata,tmp_distrib,plot_st_time,plot_end_time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stack = dbstack('-completenames');
if numel(stack) <2
   fprintf('Wrong call. Start with ''rabbit_AVNRT''\n');
   return
end

fprintf('= Analysing laddergram and plotting...  \n');

CLtime = 0;
CL = 0;
locccF=0;
locccS=0;

plot_height = 330; 
plot_width =  1200;

plot_scale = 1.0;  
pscale2 = 0.9;
pscale3 = 0.9;


%%%%%%% FIGURES %%%%%%%%%%%%%%%
fig1 = figure(1);
clf; 
axis tight

NFP_cells = S.iPB-S.iAM3-1; 
NSP_cells = S.iSlow2-S.iSlow1+1; 
DN_cells = NSP_cells - NFP_cells; 

fs1 = 12;
fs2 = 14; 
fs_xlab = 17; 
fn1 = "Helvetica";

main_col  = "#333333"; 
ap_col    = "#333333";
apS_col   = [0     0.447 0.741]; 
apF_col   = "#E50200";
pulse_col = "#000000"; 
osc_col   = ap_col;
point_col = ap_col; 
point_col1= ap_col;
coltest   = "#009000";
set(groot, 'defaultAxesTickLabelInterpreter','tex'); 

lw1 = 0.5;  
ms = 0.1; % 1.5; % normal cell marker
ms2= 5;         % Stimulus marker
msp= 1;         % Pacemaker marker 
mcol = '#555555';
lw2 = lw1;      % pacemakers 
lw3 = 0.3;
lw_top = 1.5;
lw_ap = 1.0;
ls = '-';

nskip0 = 5;
switch Mode
    case {81, 82, 83}
        x_pos = -0.065;
        xlim1 = plot_st_time;
        xlim2 = plot_end_time;
        xtk = 1.0;
        fig1.Position =  [30 277 plot_width plot_height];
        nskip = 50;         

    otherwise
        x_pos = -0.065;
        xlim1 = plot_st_time;
        xlim2 = plot_end_time;
        xtk = 0.1;
        nskip  = 50; %25;
        fig1.Position =      [30 277 450 420]; 
end



startt = min(find(t1 > xlim1));
endt =   max(find(t1 < xlim2));

ysub1 = fix((S.iHB2+1)/1.0); 
ysub2 = fix((S.iHB2+1)*plot_scale); 
ysub3 = fix((S.iHB2+1+DN_cells)*plot_scale);
ysub0 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Top panel, ladder diagrams %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = subplot(ysub0+ysub1+ysub2+ysub3,1,ysub0+1:ysub0+ysub1); 
vs1 = 0.000; vs2 = 0.005;  %0.01;
h1.InnerPosition = [0.03 0.6622+vs1 0.940 0.2618];
hold on
set(gcf,'color','w');

n_start = 1;      
n_end = S.iHB2-1; 

%%%%%%%% Measuring points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_1 = [0 xlim2];
y_1 = [-S.pointA -S.pointA];
line(x_1,y_1,'Color',point_col,'LineWidth',0.5,'LineStyle',':'); % --'); 
y_1 = [-S.pointB -S.pointB];
line(x_1,y_1,'Color',point_col,'LineWidth',0.5,'LineStyle',':'); % --');

hh1 = gca;
xtk = 0.5; 
set(gca,'XTick',[0:xtk:25]); 
xlim([xlim1 xlim2]);
ylim([-(n_end+2) n_start]);
hh1.XAxis.TickLength=[0.005 0.01];

set(hh1,'YTickMode','manual');
hh1.YAxis.TickLength=[0 0];
set(hh1,'YTick',[-19 -14 -10 -6 -2]) 
set(hh1,'YTickLabel',{'HB','\color[rgb]{.850 .325 .098}FP','\color[rgb]{0 .447 .741}SP','AM','SN'})
set(hh1,'XTickLabel',[]);

%%%% ------------------------------------ FAST PW
pos_F = zeros([],1);
end_F = NaN([],1);
k_endF = 0;
Lw = lw_top;

index_PW = ones(al,1);
index = 0;

for  n=n_start:n_end %%%%%%%%%%%%%%%%%%%%%%%%%

    for j= 1:al     
        for  i= 1:al 
            if(n<S.iAM3)
                min_int = min_interval(3); 
            else
                min_int = min_interval(2); 
            end
            if abs(locs_F(n,j)-locs_F(n+1,i)) < min_int
                 x = [locs_F(n,j) locs_F(n+1,i)];
                 y = [-n -n-1];

                  if n < S.iAM3
                    line(x,y,'Color',main_col,'LineWidth',Lw,'Marker','o','MarkerSize',ms,...
                                          'MarkerFaceColor',main_col, ...
                                          'MarkerEdgeColor',main_col);
                    if (n==S.iSN && ap.b(S.iSN,1)<0)
                        plot(locs_F(n,j),-S.iSN,'o','MarkerFaceColor',mcol, ...
                            'MarkerEdgeColor',mcol,'MarkerSize',msp); 
                    end
                  elseif n >= S.iAM3 && n<=S.iPB
                    line(x,y,'Color',apF_col,'LineWidth',Lw,'Marker','o', ...
                        'MarkerSize',ms,'LineStyle',ls, ...
                        'MarkerFaceColor',apF_col,'MarkerEdgeColor',apF_col);
                    if (n==S.iPB && ap.b(S.iPB,1)<0) ||  (n==S.iPB-1 && ap.b(S.iPB-1,1)<0)
                      plot(locs_F(n,j),-n,'o','MarkerFaceColor',apF_col, ...
                          'MarkerEdgeColor',apF_col,'MarkerSize',msp,'LineStyle',ls); 
                    end
                  else
                    col = main_col;
                    line(x,y,'Color',col,'LineWidth',Lw,'Marker','o','MarkerSize',ms,'LineStyle',ls, ...
                        'MarkerFaceColor',col,'MarkerEdgeColor',col);
                    if (n==S.ipHB && ap.b(S.ipHB,1)<0) ||  (n==S.ipHB+1 && ap.b(S.ipHB+1,1)<0)
                      plot(locs_F(n,j),-n,'o','MarkerFaceColor', col, ...
                          'MarkerEdgeColor', col,'MarkerSize',msp,'LineStyle',ls); 
                    end
                  end

                if n == S.iPB-1 
                    if x(1)~=0 && x(2)-x(1)>0
                        index= index +1;
                        pos_F(index)=x(1);
                    end
                end
            end
        end
    end
end


%%%% ------------------------------------ SLOW PW
pos_S = zeros([],1);
end_S = NaN([],1);
k_endS = 0;
index = 0; 

nnn=0; 

for n=n_start : n_end + DN_cells   % larger number in SP
    for j=1:al
        for i=1:bl 
            if (abs(locs_S(n,j)-locs_S(n+1,i)) < min_interval(1)) 
                x = [locs_S(n,j) locs_S(n+1,i)];
                y = [-locs_indexS(n) -locs_indexS(n+1)];
                if (n>=S.iAM3 && n<S.iSlow2+2  &&nnn==0)  
                      line(x,y,'Color',apS_col,'LineWidth',Lw,'Marker','o', ...
                          'MarkerSize',ms,'LineStyle',ls, ...
                          'MarkerFaceColor',apS_col,'MarkerEdgeColor',apS_col);  
                      if ( (n==S.ipSlow+1 && ap.b(S.ipSlow,  2)<0) || ...
                           (n==S.ipSlow+2 && ap.b(S.ipSlow+1,2)<0) || ...
                           (n==S.ipSlow+3 && ap.b(S.ipSlow+2,2)<0) || ...
                           (n==S.ipSlow+4 && ap.b(S.ipSlow+3,2)<0) ) 
                          plot(locs_S(n,j),-locs_indexS(n),'o', ...
                              'MarkerFaceColor',apS_col, 'LineStyle',ls, ...
                              'MarkerEdgeColor',apS_col,'MarkerSize',msp); 
                      end
                end
                if n==S.iSlow2+1
                    if x(1) < x(2) 
                       index1 = 0;
                       for k = 1:length(pos_F)
                           a = abs(x(1)-pos_F(k));
                           if abs(x(1)-pos_F(k)) < max_interval(1) 
                               if x(1)<pos_F(k)
                                  index = index +1;
                                  pos_S(index) = x(1);
                                  index1 = k+1;
                                  break;
                               else
                                   index1 = 1;
                                   break
                               end
                           end
                       end
                       if index1 ==length(pos_F) || index1==0
                          index = index +1;
                          pos_S(index) = x(1);
                       end
                    end
                end
               

                if n>= S.iSlow2+2 && x(1)~=0 
                    col = main_col;
                    for k = 1:length(nonzeros(pos_S))
                       if x(1)-pos_S(k)>0 && x(1)-pos_S(k) < max_interval(1)
                           if(nnn==0)
                             col = main_col;
                             line(x,y,'Color',col,'LineWidth',Lw,'Marker','o', ...
                                 'MarkerSize',ms,'LineStyle',ls, ...
                                 'MarkerFaceColor',col,'MarkerEdgeColor',col); 
                           end
                           if (n == S.iHB2+DN_cells-1)
                               k_endS = k_endS+1;
                               end_S(k_endS) = locs_S(n+1,j); % Store His traces with SP leading
                               index_PW(j) = 2;
                           end
                      end
                    end
                end


            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%% Mark external signal input place %%%%%%%%%%%%%%%%%%
dt=0.030;
switch Mode
    case {81,82}
      plot(tmp_distrib-dt,-S.imp_Antero,'>',...
        'MarkerFaceColor',pulse_col,'MarkerEdgeColor',pulse_col,'MarkerSize',ms2);
    case {83}
      plot(tmp_distrib-dt,-S.imp_Retro,'>',...
        'MarkerFaceColor',pulse_col,'MarkerEdgeColor',pulse_col,'MarkerSize',ms2);
    otherwise
end

set(hh1,'fontname',fn1) 
set(hh1,'fontsize',fs1) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BPM marks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch Mode
case {81, 82, 83} %%%%%%%%%%%%%%%%  81, 82, 83 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ypos = 2.5;
    xposs= ANSdata.time;
    xposs(1) = xposs(1) + plot_st_time;
    xposs(1:end-1) = xposs(1:end-1) + (xposs(2:end)-xposs(1:end-1))/2;
    xposs(end) = xposs(end)+(plot_end_time-xposs(end))/2;

    ya_min = 0.16; 
    ya_max = 0.96;
    xa0 = h1.InnerPosition(1);

    for i=1:length(xposs)
        if      ANSdata.gamma(i) ==0.56, ph_del = 0.1074; 
        elseif  ANSdata.gamma(i) ==0.57, ph_del = 0.1117; 
        elseif  ANSdata.gamma(i) ==0.65, ph_del = 0.1366; 
        elseif  ANSdata.gamma(i) == 0.7, ph_del = 0.1520; 
        elseif  ANSdata.gamma(i) == 0.8, ph_del = 0.1812; 
        elseif  ANSdata.gamma(i) ==0.85, ph_del = 0.1951;      
        elseif  ANSdata.gamma(i) == 0.9, ph_del = 0.2085; 
        elseif  ANSdata.gamma(i) == 1.0, ph_del = 0.2338; 
        elseif  ANSdata.gamma(i) == 1.1, ph_del = 0.2572; 
        elseif  ANSdata.gamma(i) == 1.2, ph_del = 0.2786; 
        elseif  ANSdata.gamma(i) == 1.3, ph_del = 0.2982; 
        elseif  ANSdata.gamma(i) == 1.4, ph_del = 0.3154;           
        elseif  ANSdata.gamma(i) == 1.5, ph_del = 0.3323; 
        elseif  ANSdata.gamma(i) == 1.6, ph_del = 0.3471; 
        elseif  ANSdata.gamma(i) == 1.7, ph_del = 0.3605;           
        elseif  ANSdata.gamma(i) == 1.8, ph_del = 0.3727;   
        elseif  ANSdata.gamma(i) == 1.9, ph_del = 0.3838; 
        elseif  ANSdata.gamma(i) == 2.0, ph_del = 0.3939;  
        else
            ph_del = 0.000; 
        end
        txts = num2str(ANSdata.gamma(i),'%0.2f');
        d = 0;
        if txts(end) == '0', d = 1; end
        txts = [txts(1:end-d) ' (' num2str(60.0/ph_del,'%.0f') ')'];
        text('Units','data','Position',[xposs(i) ypos],'String',txts, ...
          'HorizontalAlignment','center','Color',[0 0 0],'fontsize',fs2-1,'fontname',fn1);    
    end
    for i=2:length(ANSdata.time)
        xa = xa0 +  ...
        (ANSdata.time(i)-xlim1)/xlim1*(h1.InnerPosition(3)-h1.InnerPosition(1))*1.033/(xlim2-xlim1)/2.0;
        annotation('line',[xa xa],[ya_min ya_max],'Color','k','LineStyle','--');
    end
    switch Mode
     case 82  %----------------  82 PJC -------------------------------
          stim_string = 'PJC';
          text('Units','data','Position',[plot_st_time-0.15 ypos],'String',stim_string, ...
          'Color',[0.9 0 0],'fontsize',fs2,'fontname',fn1);  
     case 81  %----------------  81 PAC -------------------------------
          stim_string = 'PAC';
          text('Units','data','Position',[plot_st_time-0.15 ypos],'String',stim_string, ...
          'Color',[0.9 0 0],'fontsize',fs2,'fontname',fn1);  
     case 83  %----------------  83 PVC -------------------------------
          stim_string = 'PVC';
          text('Units','data','Position',[plot_st_time-0.15 ypos],'String',stim_string, ...
          'Color',[0.9 0 0],'fontsize',fs2,'fontname',fn1);  
    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Middle panel, action potentials via Fast %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
del = 0.30; % Distance between lines with APs

h2 = subplot(ysub0+ysub1+ysub2+ysub3,1,ysub0+ysub1+1:ysub0+ysub1+ysub2);
h2.InnerPosition = [0.03 0.3981+vs1-vs2+0.2618*(1-pscale2) ...
                    0.940 0.2618*pscale2];
cla(h2);
lw = lw_ap; 
hold on
ylim([-(n_end+1)*del 1]);  


for i = n_start : n_end+1 
                        lw = lw1; col = apF_col; %ap_col;
    if (i<=S.iAM3),     lw = lw1; col = ap_col;
    end
    if (i==S.iSN && ap.b(S.iSN,1)<0),           lw = lw2; col = osc_col;
    elseif (i==S.iPB && ap.b(S.iPB,1)<0),       lw = lw2; col = osc_col;  
    elseif (i==S.pointA || i==S.pointB),    lw = lw2; col = point_col1;
    elseif (i>S.iPB), lw = lw2;  col = ap_col;
    end 
    plot(t1(startt:nskip:endt),y1(i,startt:nskip:endt)-del*(i-1), ...
        'LineWidth',lw,'Color',col);
end

hh1 = gca;
xlim([xlim1 xlim2]);
set(hh1,'XTick',[0:xtk:25]);

hh1.YAxis.TickLength=[0 0];
set(hh1,'YTickMode','manual');
set(hh1,'YTick',[-5.0 -3.0 -1 0.5]) 
set(hh1,'YTickLabel',{'HB','\color[rgb]{.850 .325 .098}FP','AM','SN'})

hh1.XAxis.TickLength=[0.005 0.01];
set(hh1,'XTickLabel',[]);
set(hh1,'FontName',fn1);
set(hh1,'FontSize',fs1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Bottom panel, action potentials via Slow %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h3 = subplot(ysub0+ysub1+ysub2+ysub3,1,ysub0+ysub1+ysub2+1:ysub0+ysub1+ysub2+ysub3);
h3.InnerPosition = [0.03 0.1110+vs1-2*vs2+0.2618*(1-pscale2)+0.2849*(1-pscale3) ...
                    0.940 0.2849*pscale3];
cla(h3);
hold on

for i = n_start : S.iAM3 
                          lw = lw1; col = ap_col;
    if (i<=S.iAM3),       lw = lw1; col = ap_col; end
    if (i==S.iSN  && ap.b(S.iSN,1)<0), lw = lw2; col = osc_col;
    elseif (i==S.pointA),              lw = lw2; col = point_col1;
    end 
    plot(t1(startt:nskip:endt),y1(i,startt:nskip:endt)-del*(i-1), ...
        'LineWidth',lw,'Color',col);
end
j1 = i;
for i=S.iSlow1:S.iSlow2
                       lw = lw1; col = apS_col; %ap_col;
    if ((i==S.ipSlow || i==S.ipSlow+1 || i==S.ipSlow+2 || i==S.ipSlow+3  ) ...
            && ap.b(S.ipSlow,2)<0 ), lw = lw2; col = apS_col; %osc_col;  
    end 
    plot(t1(startt:nskip:endt),y2(i,startt:nskip:endt)-del*(j1+i-S.iSlow1), ...
        'LineWidth',lw,'Color',col);
end

j2 = i;
for i = S.iPB : n_end+1
                               lw = lw1; col = ap_col;
    if     (i==S.pointB)       lw = lw2; col = point_col1;
    elseif ( i==S.iPB  && ap.b(S.iPB,1)<0) lw = lw2; col = osc_col;
    end
    plot(t1(startt:nskip:endt),y1(i,startt:nskip:endt)-del*(j2+i-S.iPB+1), ...
        'LineWidth',lw,'Color',col);
end

hh3 = gca;
ylim([-(n_end+3)*del 1]);

xlim([xlim1 xlim2]);
set(hh3,'XTick',[0:xtk:25]);
hh3.XAxis.TickLength=[0.005 0.01];

set(hh3,'FontName',fn1);
set(hh3,'FontSize',fs1);

hh3.YAxis.TickLength=[0 0];
set(hh3,'YTickMode','manual');
set(hh3,'YTick',[-5.8 -3.0  -1 0.5]) 
set(hh3,'YTickLabel',{'HB','\color[rgb]{0 .447 .741}SP','AM','SN'})
hh3.XAxis.FontSize = fs2; ax.XAxis.FontName = fn1;
hh3.YAxis.FontSize = fs1; ax.YAxis.FontName = fn1;

xlabel('Time (s)','FontSize',fs_xlab);

hold off
fprintf('= Exporting to file... \n');
fname_export = ['Laddergram_var' num2str(ANSdata.var,'%d_') stim_string];
sppi = get(groot,'ScreenPixelsPerInch');
exportgraphics(fig1,[fname_export '.png'],'Resolution',sppi);