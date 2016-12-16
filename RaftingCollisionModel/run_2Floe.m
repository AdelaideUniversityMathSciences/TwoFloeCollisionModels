% function run_2Floe(sim_case,runno,sep,trans_e,Cd)
%
% LJ YIEW
% Created on  Mar 2016
% Last edited Dec 2016
%
% Runs the two-floe model using slope-sliding theory.
% Predicts collisions AND rafting.
% Runs simulations using data from two-floe rafting experiments (DUT)
%
% INPUTS:
% sim_case = 1:mooring tests, 2:rafting experiments
% runno    = run number (see DUT_Runsheet.xlsx in
%                         ../TwoFloeExperiments/RaftingExperiments)
% sep      = initial separation
% trans_e  = transient end time
% Cd       = drag coefficient
%
% OUTPUTS:
% Plots of x-displacement vs time
%
% FILES NEEDED:
% fn_SS_ode.m
% Param_DUT.m
% wavefield.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_2Floe(sim_case,runno,sep,trans_e,Cd)

close all
tic

if ~exist('sim_case','var'), sim_case = 2; end % ###
if ~exist('runno','var'), runno = 100; end % ###
if ~exist('trans_e','var'), trans_e = 5; end % ###
if ~exist('Cd','var'), Cd = 0.025; end % ###
if ~exist('sep','var'), sep = 30; end % ###
sep = sep*1e-3;


% SELECT SIMULATION CASE
switch sim_case
 case 1
  % PLOT MOORING TESTS - TO TUNE K & C
  no_floes   = 1;     % number of floes
  do_trans   = 0;     % run transients in SS model
  trans_e    = 0;
  
  K          = 0.6;   % ### adjust spring constant
  C          = 0.35;  % ### adjust damping coeff
  Cm         = 0.1;
  Cd         = 0.025; % ### adjust drag coeff
  
  H = 0;
  omega = 0;
  k = 0;
   
  if runno == 156
   InitialX = 455*1e-3; % initial displacement [m]
   phase = 7.8; % to match phases of model and experiment data
  elseif runno == 157
   InitialX = 649.3*1e-3;
   phase = 9.4;
  elseif runno == 158
   InitialX = 789.6*1e-3;
   phase = 10;
  elseif runno == 159
   InitialX = 946.3*1e-3;
   phase = 12.97;
  elseif runno == 160
   InitialX = 956.5*1e-3;
   phase = 15.87;
  end
  
  
 case 2
  % PLOT COLLISION/RAFTING TESTS
  no_floes   = 2;
  do_trans   = 1;
  do_simul   = 0;
  
  addpath('../TwoFloeExperiments/RaftingExperiments')
  dat = xlsread('DUT_Runsheet','E122:L177');
  
  for j = 1:length(dat)
   if dat(j,1) == runno
    H = dat(j,8)*1e-3;
    T = dat(j,2);
   end
  end
  
%   % ENTER INPUTS MANUALLY (WITHOUT SPECIFYING RUN NUMBER):
%   H  = 100*1e-3; % ### wave height [m]
%   T  = 1.5;     % ### wave period [s]
%   Cd = 0.025;
%   sep = 60e-3;

  freq       = 1/T;  % wave frequency [Hz]
  
  [field] = wavefield('f',freq,0.5);
  lambda  = cell2mat(field(4,2)); % corresponding wavelength
  WLn     = lambda/0.4;
  omega   = 2*pi*freq;
  k       = 2*pi/lambda;
  
  K          = 0.6;
  C          = 0.35;
  Cm         = 0.1;
  Cr         = 1;     % coefficient of restitution (0-1:perfectly elastic)
  n_cc       = 10;    % number of collisions to extract
  
  cprintf('blue',['WAVE PERIOD: ' num2str(T) '\n'])
  cprintf('blue',['WAVELENGTH*: ' num2str(WLn) '\n'])
  cprintf('blue',['WAVE HEIGHT: ' num2str(H) '\n'])
   
end


% FLAGS:
%  no_floes = number of floes (1/2)
%  do_trans = include transient waves
%  do_simul = show animation


 
%% 

display(['TUNING COEFFICIENTS:'])
display([' Added Mass = ',num2str(Cm)])
display([' Drag Coeff = ',num2str(Cd)])

display(['MOORING:'])
display([' Spring Constant     = ',num2str(K)])
display([' Damping Coefficient = ',num2str(C)])



%% 
% OTHER CONSTANTS

% SIMULATION PARAMETERS
Param = Param_DUT;
r     = Param.L;          % radius
h     = Param.h;          % water depth
dr    = Param.d;          % draft
D     = Param.D;          % thickness
g     = Param.g;          % gravity
A     = pi*r^2+2*pi*r*dr; % wetted surface area
rho_b = 650;              % floe density
rho   = 1000;             % fluid density
m     = rho_b*D*pi*r^2;   % floe mass [kg]
a     = H/2;

% INPUTS FOR SS MODEL
%  WAVE PARAMETERS
WaveParam.H     = H;
WaveParam.omega = omega;
WaveParam.k     = k;
WaveParam.rho   = rho;
WaveParam.h     = h;
%  FLOE PARAMETERS
FloeParam.m = m;
FloeParam.A = A;
%  DAMPING & ADDED MASS
Coeff.Cd = Cd;
Coeff.Cm = Cm; 
%  MOORING PARAMETERS
Mooring.K = K;
Mooring.C = C;
%  TRANSIENT WAVE AMPLITUDE PARAMETERS
%   A(transient) = f.a*1e-3*t^f.b*sin(-omega*t + f.c)
f.c = 0;
f.b = 3;
f.a = H*1e3/2/(trans_e^f.b);
if isnan(f.a)
 f.a = 0;
end
Trans.f = f;

%  PLOT TRANSIENT WAVE PROFILE
if sim_case == 2
 tt = linspace(0,trans_e,10000);
 A_ft = 1e-3*f.a.*tt.^f.b.*sin(-omega.*tt + f.c);
 ts = linspace(trans_e,60,10000);
 A_fs = H/2.*sin(-omega.*ts + f.c);
 figure(1)
 hold on
 plot(tt,A_ft*1e3,'r')
 plot(ts,A_fs*1e3,'k')
 xlim([0 60])
 ylim([min(A_fs)*1e3*1.2 max(A_fs)*1e3*1.2])
 ylabel('Wave Amplitude [mm]')
 xlabel('t [s]')
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COLLISION MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RANGE OF TIME
step = 0.005;
tspan_t = 0:step:trans_e;  % transient period
tspan_s = trans_e:step:60; % steady state period


% INITIAL CONDITIONS [DISPLACEMENT, VELOCITY]
IC1  = [0,0];
if no_floes == 2
 dist = 2*r + sep;
 IC2  = [dist,0];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SINGLE FLOE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if no_floes == 1
   
   if do_trans == 1 % run transients
    
    % TRANSIENT SOLUTION
    Trans.t = 1; % flag for transient solution (time dependent amplitude)
    [t1_t,X1_t] = ode45(@(t,X) ...
                fn_SS_ode(t,X,0,WaveParam,FloeParam,Coeff,Mooring,Trans), ...
                tspan_t,IC1);
    % STEADY STATE SOLUTION
    IC1s    = X1_t(end,:);
    Trans.t = 0;
    WaveParam.H = H;
    [t1_s,X1_s] = ode45(@(t,X) ...
                fn_SS_ode(t,X,0,WaveParam,FloeParam,Coeff,Mooring,Trans), ...
                tspan_s,IC1s);
    % STITCH TRANSIENT + STEADY STATE 
    t1 = [tspan_t,tspan_s]';
    X1 = [X1_t;X1_s];

   elseif do_trans == 0 % no transients
    
    IC1 = [InitialX,0];
    tspan = 0:step:60;
    Trans.t = 0;
    WaveParam.H = H;
    [t1,X1] = ode45(@(t,X) ...
                fn_SS_ode(t,X,0,WaveParam,FloeParam,Coeff,Mooring,Trans), ...
                tspan,IC1);
   
   end
              
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TWO FLOES   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  elseif no_floes == 2

   % COLLISION ALGORITHM
   
   cc   = n_cc; % collision counter
   X1c = [];
   X2c = [];
   t1c = [0];
   t2c = [0];
   
%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   while cc > 0
    
    if do_trans == 1

     if t1c(end) < trans_e
     
      tspan_t = t1c(end):step:trans_e;
      tspan_s = trans_e:step:trans_e+60;
     
      % SOLUTION FOR [TIME, [DISPLACEMENT,VELOCITY]]
      %  FLOE F
      %   TRANSIENT SOLUTION
      Trans.t = 1;
      Trans.f.c = f.c;
      [t1_t,X1_t] = ode45(@(t,X) ...
                 fn_SS_ode(t,X,0,WaveParam,FloeParam,Coeff,Mooring,Trans), ...
                 tspan_t,IC1);
      %   STEADY STATE SOLUTION
      IC1s    = X1_t(end,:);
      Trans.t = 0;
      Trans.f.c = f.c;
      WaveParam.H = H;
      [t1_s,X1_s] = ode45(@(t,X) ...
                 fn_SS_ode(t,X,0,WaveParam,FloeParam,Coeff,Mooring,Trans), ...
                 tspan_s,IC1s);
      % STITCH SOLUTION
      t1 = [tspan_t,tspan_s]';
      X1 = [X1_t;X1_s];         

      %  FLOE R
      %   TRANSIENT SOLUTION
      Trans.t = 1;
      [t2_t,X2_t] = ode45(@(t,X) ... 
                 fn_SS_ode(t,X,dist,WaveParam,FloeParam,Coeff,Mooring,Trans), ...
                 tspan_t,IC2);
      %   STEADY STATE SOLUTION          
      IC2s    = X2_t(end,:);
      Trans.t = 0;
      [t2_s,X2_s] = ode45(@(t,X) ...
                 fn_SS_ode(t,X,dist,WaveParam,FloeParam,Coeff,Mooring,Trans), ...
                 tspan_s,IC2s);
      % STITCH SOLUTION
      t2 = [tspan_t,tspan_s]';
      X2 = [X2_t;X2_s]; 
      
    elseif t1c(end) >= trans_e
     
     tspan_s = t1c(end):step:t1c(end)+60;

     % SOLUTION FOR [TIME, [DISPLACEMENT,VELOCITY]]
     %  FLOE F
     %   STEADY STATE SOLUTION
     Trans.t = 0;
     Trans.f.c = f.c;
     WaveParam.H = H;
     [t1,X1] = ode45(@(t,X) ...
                fn_SS_ode(t,X,0,WaveParam,FloeParam,Coeff,Mooring,Trans), ...
                tspan_s,IC1);         
     %  FLOE R
     %   STEADY STATE SOLUTION
     Trans.t = 0; 
%      Trans.f.c = f.c + k*dist; % ### NO DRIFT
     [t2,X2] = ode45(@(t,X) ...
                fn_SS_ode(t,X,dist,WaveParam,FloeParam,Coeff,Mooring,Trans), ...
                tspan_s,IC2);
    end
    
   elseif do_trans == 0
      
      % SOLUTION FOR [TIME, [DISPLACEMENT,VELOCITY]]
      %  FLOE F
      tspan = t1c(end):step:t1c(end)+30;
%       tspan = linspace(t1c(end),t1c(end)+30,10000);
      Trans.t = 0;
      Trans.f.c = f.c;
      WaveParam.H = H;
      [t1,X1] = ode45(@(t,X) ...
                 fn_SS_ode(t,X,0,WaveParam,FloeParam,Coeff,Mooring,Trans), ...
                 tspan,IC1);    

      %  FLOE R
      Trans.f.c = f.c + k*dist;
      [t2,X2] = ode45(@(t,X) ...
                 fn_SS_ode(t,X,dist,WaveParam,FloeParam,Coeff,Mooring,Trans), ...
                 tspan,IC2);
      
   end
    
      
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COLLISION CRITERIA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cn = 0; % counter for no collisions
    
    resolution = 200:1:length(t1);
    
    for j = resolution
     
     [FRL,FRU,FLL,FLU,RLL,RLU,RRL,RRU] = FloePosition(k,a,r,D,omega,X1(j,1),X2(j,1),t1(j),t2(j),trans_e,f);
          
     lineFB = [FLL(1) FLL(2); FRL(1) FRL(2)]; % bottom surface of floe F
     lineFT = [FLU(1) FLU(2); FRU(1) FRU(2)]; % top    surface of floe F
     lineFR = [FRL(1) FRL(2); FRU(1) FRU(2)]; % right surface of floe F
     lineRB = [RLL(1) RLL(2); RRL(1) RRL(2)]; % likewise for floe R
     lineRT = [RLU(1) RLU(2); RRU(1) RRU(2)]; % 
     lineRL = [RLL(1) RLL(2); RLU(1) RLU(2)]; % left surface

     slope = @(line) (line(2,2) - line(1,2))/(line(2,1) - line(1,1));
     mFB = slope(lineFB);
     mFT = slope(lineFT);
     mFR = slope(lineFR);
     mRB = slope(lineRB);
     mRT = slope(lineRT);
     mRL = slope(lineRL);
     
     intercept = @(line,m) line(1,2) - m*line(1,1);
     bFB = intercept(lineFB,mFB);
     bFT = intercept(lineFT,mFT);
     bFR = intercept(lineFR,mFR);
     bRB = intercept(lineRB,mRB);
     bRT = intercept(lineRT,mRT);
     bRL = intercept(lineRL,mRL);
     xintersectFRRL = (bFR-bRL)/(mRL-mFR);       % side collision
     yintersectFRRL = mFR*xintersectFRRL + bFR; 
     xintersectFRRB = (bFR-bRB)/(mRB-mFR);       
     yintersectFRRB = mFR*xintersectFRRB + bFR;  
     xintersectFRRT = (bFR-bRT)/(mRT-mFR);
     yintersectFRRT = mFR*xintersectFRRT + bFR;  
     xintersectRLFB = (bRL-bFB)/(mFB-mRL);
     yintersectRLFB = mRL*xintersectRLFB + bRL;  
     xintersectRLFT = (bRL-bFT)/(mFT-mRL);
     yintersectRLFT = mRL*xintersectRLFT + bRL; 
     
     isPointInside = @(xint,myline) ...
                      (xint >= myline(1,1) && xint <= myline(2,1)) || ...
                      (xint >= myline(2,1) && xint <= myline(1,1));
     clear insideFRRL insideFRRB insideFRRT insideRLFB insideRLFT
     insideFRRL = isPointInside(xintersectFRRL,lineFR) && ...
                  isPointInside(xintersectFRRL,lineRL);
     insideFRRB = isPointInside(xintersectFRRB,lineFR) && ...
                  isPointInside(xintersectFRRB,lineRB);
     insideFRRT = isPointInside(xintersectFRRT,lineFR) && ...
                  isPointInside(xintersectFRRT,lineRT);
     insideRLFB = isPointInside(xintersectRLFB,lineRL) && ...
                  isPointInside(xintersectRLFB,lineFB);
     insideRLFT = isPointInside(xintersectRLFT,lineRL) && ...
                  isPointInside(xintersectRLFT,lineFT);

     
     % IF INTERSECT
     if insideFRRL == 1 || ...
        insideFRRB == 1 || insideFRRT == 1 || ...
        insideRLFB == 1 || insideRLFT == 1
       
      if FRU(2) < RLL(2) || FRL(2) > RLU(2)
       display('Rafting!')
       cc = 0;
      else
       display('Collision!')
      end
      VFf = ((1-Cr)*X1(j,2) + (1+Cr)*X2(j,2))/2; 
      VBf = ((1+Cr)*X1(j,2) + (1-Cr)*X2(j,2))/2; 
      IC1 = [X1(j,1),VFf];
      IC2 = [X2(j,1),VBf]; 
      % STITCH RESULTS
      X1c = [X1c;X1(1:j,:)];
      X2c = [X2c;X2(1:j,:)];
      t1c = [t1c(:);t1(1:j)];
      t2c = t1c;
      % COUNTERS
      cc = cc-1; % collision counter
      cy = 1;    % flag for collisions
      ccj(-cc+n_cc) = j; % counter for collision position (in time)
      break
     
     % NO INTERSECTION
     else
      cn = cn+1;
     end

    end

    c_num = -cc+n_cc; % total number of collisions
    
    % IF NO MORE COLLISIONS OVER TSPAN
    if cn == length(resolution)
 %     display('No (more) collisions')
     cc = cc-n_cc;
     % STITCH RESULTS
     X1c = [X1c;X1(1:end,:)];
     X2c = [X2c;X2(1:end,:)];
     t1c = [t1c(:);t1(1:end)];
     t2c = t1c;
    end

   end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%
   

%    display(['Total Collisions: ',num2str(c_num)])

   % ARRANGE COLLISION DATA
   if exist('cy')
    X1 = X1c;
    X2 = X2c;
    t1 = t1c(2:end);
    t2 = t2c(2:end);
    for j = 2:length(ccj)
     ccj(j) = ccj(j) + ccj(j-1);
    end
   end 
   
   % MEAN COLLISION VELOCITY
   if exist('ccj')  
    Vcol   = (X1(ccj,2) - X2(ccj,2))/2;
    Mod.Vc = [mean(Vcol) median(Vcol)];
   else
    Mod.Vc = [0 0];
   end
   
   % CALCULATE COLLISION FREQUENCY 
   if exist('ccj') && length(ccj)>1
    col_time = t1(ccj); % time at each collision
    for j = 1:length(ccj)-1
     df(j) = col_time(j+1) - col_time(j);
    end
    Mod.Tc = [mean(df),median(df)]; % collision period
    Mod.Fc = 1./Mod.Tc;             % collision frequency
    % CALCULATE NUMBER OF COLLISIONS (IN 60 SECONDS)
    cc_60 = 0;
    for j = 1:length(ccj)
     if col_time(j) <= 60
      cc_60 = cc_60 + 1;
     end
    end
    Mod.Nc = cc_60; % for direct comparison with Exp.Nc

   elseif exist('ccj') && length(ccj)==1
    Mod.Nc = 1; % only 1 collision
    Mod.Fc = [0 0];   
   else
    Mod.Nc = 0; % no collisions
    Mod.Fc = [0 0];
    
   end
    cprintf('blue',['NO. OF COLLISIONS / FREQUENCY [Hz] (MEAN,MEDIAN) \n'])
    display([' 2-Floe Model: ',num2str(Mod.Nc),' / ',...
                                num2str(Mod.Fc(1),'%0.4f'),' ',...
                                num2str(Mod.Fc(2),'%0.4f')])
                               
    cprintf('blue',['COLLISION VELOCITY [mm/s] (MEAN,MEDIAN) \n'])
    display([' 2-Floe Model: ',num2str(Mod.Vc(1)*1e3,'%0.1f'),' ',...
                                num2str(Mod.Vc(2)*1e3,'%0.1f')])
                               
  
  end
  
  

%% PLOT FIGURES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if no_floes == 1
    % PLOT X VS TIME
    % plot experimental data
    addpath('../TwoFloeExperiments/RaftingExperiments/Figures')
    openfig(['Run' num2str(runno) '.fig'])
    hold on
    model = plot(t1,X1(:,1)*1e3,'k-');
    if exist('lambda')
     title(['Wavelength = ',num2str(lambda),' m, Wave Amplitude = ',num2str(a),' m'])
    end
    ylabel('x [mm]')
    xlabel('t [s]')
  %   ylim([1.5e3*min(X1(:,1)) 1.5e3*max(X2(:,1))])
    ylim('auto')
    box on
    grid on
    hold off
    
    % adjust phase
    figure(1)
    hold on
    cprintf('magenta','Enter Phase Difference [s]:')
    if exist('runno')
    else
     phase = input(' ');
    end
    delete(model)
    plot(t1+phase,X1(:,1)*1e3,'k--')
    
    
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  elseif no_floes == 2
   % PLOT X VS TIME
   [FRL,FRU,FLL,FLU,RLL,RLU,RRL,RRU] = FloePosition(k,a,r,D,omega,X1(:,1),X2(:,1),t1,t2,trans_e,f);
   R1x = FRL(:,1);
   L2x = RLL(:,1);
   
   addpath('../TwoFloeExperiments/RaftingExperiments/Figures')
   openfig(['Run' num2str(runno) '.fig'])
   
   figure(2)
   hold on
   plot(t1,-(R1x-r)*1e3,'b--') % right edge of floe F 
   plot(t2,-(L2x-r)*1e3,'k--') % left edge of floe R
   %  nb. axis reversed to from experiments - positive x direction points to the wavemaker
   % PLOT COLLISIONS
   if exist('ccj')
    for j = 1:length(ccj)
     plot(t1(ccj(j)),-(R1x(ccj(j))-r)*1e3,'m.','MarkerSize',20)
    end
   end
%    title(['Wavelength = ',num2str(lambda),' m, Wave Amplitude = ',num2str(a),' m'])
   ylabel('X [mm]')
   xlabel('t [s]')
   xlim([0 60])
   box on
   grid on
   hold off
   

%%
   % SIMULATE WAVES AND FLOE IN TIME / CHECK IF COLLISION OR RAFTING
   if do_simul == 1;
    
    xspan = linspace(-0.8,1.2,1000);
    n = 1;   % counter for collisions ccj(n)
    jjj = 1; % counter for simulation step
    for j = 1:5:length(t1)
     t = t1(j);
     for jj = 1:length(xspan)
      x = xspan(jj);
      if t < trans_e
       eta(jj) = 1e-3*f.a.*t.^f.b*sin(k*x-omega*t);
      elseif t >= trans_e
       eta(jj) = a*sin(k*x-omega*t);
      end
     end
     figure(5)
     clf
     hold on
     % wave profile
     plot(-xspan,eta,'k-') 
     % floe
     [FRL,FRU,FLL,FLU,RLL,RLU,RRL,RRU] = FloePosition(k,a,r,D,omega,X1(j,1),X2(j,1),t,t,trans_e,f);
     %
     lineFB = [FLL(1) FLL(2); FRL(1) FRL(2)]; % lower surface of floe F
     lineFT = [FLU(1) FLU(2); FRU(1) FRU(2)]; % upper surface of floe F
     lineFR = [FRL(1) FRL(2); FRU(1) FRU(2)]; % right surface of floe F
     lineFL = [FLL(1) FLL(2); FLU(1) FLU(2)]; % left  surface of floe F
     lineRB = [RLL(1) RLL(2); RRL(1) RRL(2)]; % likewise for floe R
     lineRT = [RLU(1) RLU(2); RRU(1) RRU(2)]; % 
     lineRL = [RLL(1) RLL(2); RLU(1) RLU(2)]; % 
     lineRR = [RRL(1) RRL(2); RRU(1) RRU(2)]; % 
     plot(-lineFB(:,1),lineFB(:,2),'b');
     plot(-lineFT(:,1),lineFT(:,2),'b');
     plot(-lineFR(:,1),lineFR(:,2),'b');
     plot(-lineFL(:,1),lineFL(:,2),'b');
     plot(-lineRB(:,1),lineRB(:,2),'g');
     plot(-lineRT(:,1),lineRT(:,2),'g');
     plot(-lineRR(:,1),lineRR(:,2),'g');
     plot(-lineRL(:,1),lineRL(:,2),'g');
     

     grid on
     axis equal
     box on
     xlim(-[xspan(end) xspan(1)]), ylim([-0.3 0.3]) %ylim([-50*a 50*a])
     xlabel('x [m]');ylabel('z [m]')
     % show time
     annotation('textbox', [0.15,0.78,0.1,0.1],...
             'String', ['time = ',num2str(t,'%0.4f'),' [s]'],'LineStyle','none');
%      % indicate collisions
%      if exist('ccj')
% %       if j == ccj(n)
%       if ismember(j,ccj)
%        annotation('textbox', [0.42,0.3,0.1,0.1],...
%               'String', ['COLLISION! / RAFTING!'],...
%               'LineStyle','none','Color','r','FontSize',20);
% %        pause      
%        if n == length(ccj)
%         n = n;
%        else
%         n = n+1;
%        end
%        pause
%       end
%      end
     hold off
     
     mov(jjj) = getframe(gcf);

     jjj = jjj+1;
    end
    
   end
 
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



toc

return


function [FRL,FRU,FLL,FLU,RLL,RLU,RRL,RRU] = FloePosition(k,a,r,D,omega,X1,X2,t1,t2,trans_e,f)
% nb. t1 = t2

% FLU(x,y) => o-------------o <= FRU(x,y)  RLU(x,y) => o-------------o <= RRU(x,y)
%             |   FLOE F    |                          |   FLOE R    |
% FLL(x,y) => o-------------o <= FRL(x,y)  RLL(x,y) => o-------------o <= RRL(x,y)


% calculates floe elevation (see research journal pg 16)
for j = 1:length(t1)
 if t1(j) < trans_e
  FG(j,1) = 1e-3*f.a.*t1(j).^f.b.*k.*cos(k.*X1(j)-omega.*t1(j)); % gradient at centre of lower surface of floe F #######
  BG(j,1) = 1e-3*f.a.*t2(j).^f.b.*k.*cos(k.*X2(j)-omega.*t2(j)); % gradient at centre of lower surface of floe R
 elseif t1(j) >= trans_e
  FG(j,1) = k*a.*cos(k.*X1(j)-omega.*t1(j)); % gradient at centre of lower surface of floe F #######
  BG(j,1) = k*a.*cos(k.*X2(j)-omega.*t2(j)); % gradient at centre of lower surface of floe R
 end
end
  
 FO = r.*FG./sqrt(1+FG.^2);
 BO = r.*BG./sqrt(1+BG.^2);
 FA = FO./FG;
 BA = BO./BG;
 for j = 1:length(t1)
  if t1(j) < trans_e
   eta_F(j,1) = 1e-3*f.a.*t1(j).^f.b.*sin(k.*X1(j)-omega.*t1(j)); % wave profile (y coordinate) at centre of floe F
   eta_B(j,1) = 1e-3*f.a.*t2(j).^f.b.*sin(k.*X2(j)-omega.*t2(j)); % wave profile (y coordinate) at centre of floe R   
  elseif t1(j) >= trans_e
   eta_F(j,1) = a.*sin(k.*X1(j)-omega.*t1(j)); % wave profile (y coordinate) at centre of floe F
   eta_B(j,1) = a.*sin(k.*X2(j)-omega.*t2(j)); % wave profile (y coordinate) at centre of floe R
  end
 end
 FRLx = X1+FA;      
 FRLy = FO+eta_F;   
 RLLx = X2-BA;      
 RLLy = -BO+eta_B; 
 %
 thetaF = -atan(FG);
 thetaB = -atan(BG);
%  thetaF = -atan(k.*a.*cos(k.*X1-omega.*t1));
%  thetaB = -atan(k.*a.*cos(k.*X2-omega.*t2));
 FRUx = FRLx+D.*sin(thetaF);
 FRUy = FRLy+D.*cos(thetaF);
 RLUx = RLLx+D.*sin(thetaB);
 RLUy = RLLy+D.*cos(thetaB);
 %
 FLLx = X1-FA;
 FLLy = eta_F-FO;
 FLUx = FLLx+D.*sin(thetaF);
 FLUy = FLLy+D.*cos(thetaF);
 %
 RRLx = X2+BA;
 RRLy = eta_B+BO;
 RRUx = RRLx+D.*sin(thetaB);
 RRUy = RRLy+D.*cos(thetaB);
 %
 FRL = [FRLx,FRLy];
 FRU = [FRUx,FRUy];
 FLL = [FLLx,FLLy];
 FLU = [FLUx,FLUy];
 RLL = [RLLx,RLLy];
 RLU = [RLUx,RLUy];
 RRL = [RRLx,RRLy];
 RRU = [RRUx,RRUy];
 
 
return
