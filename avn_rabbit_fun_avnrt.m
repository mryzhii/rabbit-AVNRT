%%%%%%%%%%%%%%% Differential equations functions  %%%%%%%%%%%%%%%%%%%%%%%%
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
function dydt = avn_rabbit_fun_avnrt(t,y,Mode,S,ap,tmp_distrib, ...
                                      ANSdata, ANS_point)
global tN

stack = dbstack('-completenames');
if numel(stack) <2
   fprintf('Wrong call. Start with ''rabbit_AVNRT''\n');
   return
end

y1 = reshape(y,[S.iHB2,4]);
Cx = zeros(S.iHB2,2); 
Cy = zeros(S.iHB2,1);
Cyp1 = zeros(S.iHB2,1);
Cyp2 = zeros(S.iHB2,1);
dydt = zeros(S.iHB2,4);

tdur=ap.pulse_width;     % pulse duration (s)
tau=0.00002;             %pulse time constant
pulse = 0;

% Set S1-S2 stimulus
r = tmp_distrib(tN);
if (t < r+tdur*10 )
    tmod=mod(t-r/2,r)-r/2;       %modulo time
    yup=1./(1+exp(-tmod/tau));   %this part goes up, from 0 to 1, as tmod crosses 0
    ydown=1./(1+exp((tmod-tdur)/tau)); 
    pulse = ap.pulse_amp*yup*ydown;
elseif(tN<length(tmp_distrib))
    tN = tN+1;
end



%%%% First row
Cx(    S.iSN,1) = ap.dx(S.iSN,1).*(y1(S.iSN+1,1) - ...
                   y1(S.iSN,1).*ap.adx(S.iSN,1)); 
Cx(S.iSN+1:S.iHB2-1,1)= ap.dx(S.iSN:S.iHB2-2,1).*( y1( S.iSN:S.iHB2-2,1) - ...   % Left
                        y1(S.iSN+1:S.iHB2-1,1).*ap.adx(S.iSN:S.iHB2-2,1)) + ...  % Left
                        ap.dx(S.iSN+1:S.iHB2-1,1).*(-y1(S.iSN+1:S.iHB2-1,1) + ...
                        y1(S.iSN+2:S.iHB2,  1).*ap.adx(S.iSN+1:S.iHB2-1,1));     % Right
Cx(    S.iHB2,1) = 1.0.*ap.dx(S.iHB2-1,1).*(y1(S.iHB2-1,1) - ...
                         y1(S.iHB2  ,1).*ap.adx(S.iHB2-1,1));  
%%%% Second row
Cx(S.iSN+1:S.iHB2-1,2)= ap.dx(S.iSN:S.iHB2-2,2).*( y1(S.iSN:S.iHB2-2,3) - ...    % Left
                        y1(S.iSN+1:S.iHB2-1,3).*ap.adx(S.iSN:S.iHB2-2,2)) + ...  % Left
                        ap.dx(S.iSN+1:S.iHB2-1,2).*(-y1(S.iSN+1:S.iHB2-1,3) + ...
                        y1(S.iSN+2:S.iHB2,  3).*ap.adx(S.iSN+1:S.iHB2-1,2));     % Right

%%%%  Vertical coupling between FP and SP (with asymmetry)
Cy(S.iAM3) = ap.dy(S.iAM3).*(y1(S.iAM3,1) - y1(S.iAM3,3).*ap.ady(S.iAM3 ));  
Cy(S.iPB)  = ap.dy(S.iPB ).*(y1(S.iPB, 1).*ap.ady(S.iPB) - y1(S.iPB,3));    

switch Mode
    case {81,82}, Cyp1(S.imp_Antero) = pulse; 
    case  83,    Cyp1(S.imp_Retro) = pulse; 
end

gamma = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%% Stimulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (length(ANSdata.time)>=9 &&  t> ANSdata.time(9))    %  9
    gamma = ANSdata.gamma(9)/ap.var;
    ap.mu2(1:ANS_point,:) =  ap.mu2(1:ANS_point,:).*gamma;
    ap.mu1(1:ANS_point,:) =  ap.mu1(1:ANS_point,:)./gamma; 
elseif (length(ANSdata.time)>=8 &&  t> ANSdata.time(8)) % 8
    gamma = ANSdata.gamma(8)/ap.var;
    ap.mu2(1:ANS_point,:) =  ap.mu2(1:ANS_point,:).*gamma; 
    ap.mu1(1:ANS_point,:) =  ap.mu1(1:ANS_point,:)./gamma; 
elseif (length(ANSdata.time)>=7 &&  t> ANSdata.time(7)) % 7
    gamma = ANSdata.gamma(7)/ap.var;
    ap.mu2(1:ANS_point,:) =  ap.mu2(1:ANS_point,:).*gamma; 
    ap.mu1(1:ANS_point,:) =  ap.mu1(1:ANS_point,:)./gamma; 
elseif (length(ANSdata.time)>=6 &&  t> ANSdata.time(6)) % 6
    gamma = ANSdata.gamma(6)/ap.var;
    ap.mu2(1:ANS_point,:) =  ap.mu2(1:ANS_point,:).*gamma;
    ap.mu1(1:ANS_point,:) =  ap.mu1(1:ANS_point,:)./gamma; 
elseif (length(ANSdata.time)>=5 &&  t> ANSdata.time(5)) %  5
    gamma = ANSdata.gamma(5)/ap.var;
    ap.mu2(1:ANS_point,:) =  ap.mu2(1:ANS_point,:).*gamma; 
    ap.mu1(1:ANS_point,:) =  ap.mu1(1:ANS_point,:)./gamma;
elseif (length(ANSdata.time)>=4 && t> ANSdata.time(4)) % 4
    gamma = ANSdata.gamma(4)/ap.var;
    ap.mu2(1:ANS_point,:) =  ap.mu2(1:ANS_point,:).*gamma;  
    ap.mu1(1:ANS_point,:) =  ap.mu1(1:ANS_point,:)./gamma; 
elseif (length(ANSdata.time)>=3 &&  t> ANSdata.time(3)) % 3
    gamma = ANSdata.gamma(3)/ap.var; 
    ap.mu2(1:ANS_point,:) =  ap.mu2(1:ANS_point,:).*gamma;  
    ap.mu1(1:ANS_point,:) =  ap.mu1(1:ANS_point,:)./gamma; 
elseif (length(ANSdata.time)>=2 &&  t> ANSdata.time(2)) % 2
    gamma = ANSdata.gamma(2)/ap.var; 
    ap.mu2(1:ANS_point,:) =  ap.mu2(1:ANS_point,:).*gamma; 
    ap.mu1(1:ANS_point,:) =  ap.mu1(1:ANS_point,:)./gamma; 
end

% Calculate derivatives
% First raw (Fast)   
dydt(:,1)=5*0.5e2*(-ap.k1(:,1).*y1(:,1).*(y1(:,1)-ap.b(:,1)).*(y1(:,1)-1)  ...
    -y1(:,1).*y1(:,2)) +Cx(:,1) -Cy(:)  -Cyp1(:);
dydt(:,2)=5*0.5e2*(ap.eps0(:,1)+ap.mu1(:,1).*y1(:,2)./(ap.mu2(:,1)+y1(:,1))).*(-y1(:,2)-ap.k2(:,1).*y1(:,1).*(y1(:,1)-ap.a(:,1)-1.0));
% Second raw (AM, Slow)
dydt(:,3)=5*0.5e2*(-ap.k1(:,2).*y1(:,3).*(y1(:,3)-ap.b(:,2)).*(y1(:,3)-1)  ...
    -y1(:,3).*y1(:,4)) +Cx(:,2) +Cy(:) -Cyp2(:);
dydt(:,4)=5*0.5e2*(ap.eps0(:,2)+ap.mu1(:,2).*y1(:,4)./(ap.mu2(:,2)+y1(:,3))).*(-y1(:,4)-ap.k2(:,2).*y1(:,3).*(y1(:,3)-ap.a(:,2)-1.0));

if (length(ANSdata.time)>=2 && t> ANSdata.time(2))
    ap.mu2(1:ANS_point,:) =  ap.mu2(1:ANS_point,:)./gamma;
    ap.mu1(1:ANS_point,:) =  ap.mu1(1:ANS_point,:).*gamma; 
end

dydt = reshape(dydt,[S.iHB2*4,1]);

end % function end





