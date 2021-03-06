% function run_2Floe(runno)
%
% LJ YIEW
% Created on  May 2015
% Last edited Oct 2016
%
% Runs the two-floe model using slope-sliding theory.
% Runs simulations using data from two-floe non-rafting experiments (AMC)
%
% INPUTS:
% runno = run number (see Data_ExperimentVsModel.xlsx)
%
% OUTPUTS:
% Plots of x-displacement vs time
%
% FILES NEEDED:
% fn_SS_ode.m
% fn_WaveFit.m
% Param_AMC.m
% wavefield
% cprintf.m
% Data_ExperimentVsModel.xlsx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_2Floe(runno,Cd)

close all
tic

%% MODEL INPUTS

 if ~exist('runno','var');
 runno = 94;  end   % run number (good examples: 94,120,121)
 K     = 0.32;      % spring constant  (mooring)
 C     = 0.50;      % damping constant (mooring)
 Cm    = 0.1;       % added mass coefficient
 if ~exist('Cd','var');
 Cd    = 0.005; end % drag coefficient ###adjust accordingly
 Cr    = 1;         % coefficient of restitution
 n_cc  = 10;        % number of collisions to extract
 sep   = 20e-3;     % initial separation [m]

 cprintf('blue',['RUN ' num2str(runno) '\n'])
 display(['MOORING COEFFS:'])
 display([' Spring     = ',num2str(K)])
 display([' Damping    = ',num2str(C)])
 display(['HYDRODYNAMIC COEFFS:'])
 display([' Added Mass = ',num2str(Cm)])
 display([' Drag       = ',num2str(Cd)])
  
% FLAGS:
%  do_trans = include transient wave signals from experiments
%  do_colex = read data experiments (for comparison)
%  do_simul = show animation
 do_colex   = 1;
 do_trans   = 1;
 do_simul   = 0;

%%
% EXTRACT EXPERIMENT DATA (FREQUENCY)
dat = xlsread('Data_ExperimentVsModel.xlsx','B6:P67');
run = dat(:,1); % run numbers

for j = 1:length(run)
 if run(j) == runno
  freq    = dat(j,2);      % frequency
  H       = dat(j,3)*1e-3; % wave height
  [field] = wavefield('f',freq,0.831);
  lambda  = cell2mat(field(4,2)); % corresponding wavelength
  omega   = 2*pi*freq;
  k       = 2*pi/lambda;
  trans   = dat(j,12);
 end
end
 

%% 
% EXTRACT COLLISION DATA FROM EXPERIMENTS (FOR TUNING)

if do_colex == 1
 cd('../TwoFloeExperiments/NonRaftingExperiments')
 tt.s=22;  % ###
 tt.e=40;  % ###
 d_factor=0.5;
 dosg = 0; % plot surge?
 [res,out] = fn_ColsExp(['find(run==',num2str(runno),')'],tt,d_factor,dosg,1,0);
 Exp.Tc = out.Tc; % period of collisions
 Exp.Nc = out.Nc; % number of collisions
 Exp.Vc = out.Vc; % collision velocity (mean/median)
 cd('../../NonRaftingCollisionModel')
 fig_sg = gcf;
end

%%
% TUNE TRANSIENTS - BASED ON EXPERIMENT WAVE PROBE DATA

cd('../TwoFloeExperiments/NonRaftingExperiments')
[time,~,~,~,Ap] = WaveProbes(runno,0);
cd('../../NonRaftingCollisionModel')

%  SHOW WAVE PROFILE
figure
set(gcf,'position',[100 400 1000 400]);
set(gca,'FontSize',14)
hold on
plot(time,Ap,'b')
xlabel('t [s]')
ylabel('Wave Amplitude [mm]')
title('Wave Profile from Non-Rafting Experiment')
grid on
box on

%  SELECT APPROX TRANSIENT END TIME
% cprintf('magenta',['Enter Approx Transient End Time:'])
% trans_e = input(' '); % transient end time
trans_e = trans;
close(gcf)

%  FITTED WAVE PARAMETERS
[f] = fn_WaveFit(runno,0,trans_e,1,[0.1 0 0]); % ### if error: change +/- 0.1
% f.a=H/2;f.b=3;f.c=0;

%  FIND CLOSEST TIME TO STEADY STATE AMPLITUDE (H/2)
tt = linspace(0,trans_e,10000);
A_f = 1e-3*f.a.*tt.^f.b.*sin(-omega.*tt + f.c);
for j = 1:length(tt)
 if A_f(j) >= H/2
  trans_e = tt(j);
  break
 end
end
cprintf('magenta',['(Calculated transient end time = ',num2str(trans_e,'%0.1f'),'s)','\n'])

%  STEADY STATE AMPLITUDE [m] = TRANSIENT AMPLITUDE AT t = trans_e
A_ss = 1e-3*f.a*trans_e^f.b; 

%  PLOT ADJUSTED AMPLITUDE
tt = linspace(0,trans_e,10000);
A_ft = 1e-3*f.a.*tt.^f.b.*sin(-omega.*tt + f.c);
ts = linspace(trans_e,60,10000);
A_fs = A_ss.*sin(-omega.*ts + f.c);
fig_wp = gcf;
figure(fig_wp)
hold on
plot(tt,A_ft*1e3,'r--')
plot(ts,A_fs*1e3,'k--')
xlim([0 60])


%% 
% OTHER CONSTANTS

% PHYSICAL PARAMETERS
Param = Param_AMC;
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
Trans.f = f;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COLLISION ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for count = 1:length(lambda)
 
  k     = 2*pi/lambda(count);
  omega = sqrt(g*k*tanh(k*h));
  c     = omega/k;	% wave celerity

  % RANGE OF TIME
  tspan_t = linspace(0,trans_e,1000);   % transient period
  tspan_s = linspace(trans_e,60,1000); % steady state period
   

  %% SOLVE NUMERICAL SOLUTION
  
  % INITIAL CONDITIONS [DISPLACEMENT, VELOCITY]
  IC1  = [0,0];
  dist = 2*r + sep;
  IC2  = [dist,0];
   
  cc   = n_cc; % collision counter
  X1c = [];
  X2c = [];
  t1c = [0];
  t2c = [0];
   
%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   while cc > 0
    
    if do_trans == 1

     if t1c(end) < trans_e
     
      tspan_t = linspace(t1c(end),trans_e,10000);
      tspan_s = linspace(trans_e,trans_e+60,10000);
     
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
      WaveParam.H = A_ss*2;
      [t1_s,X1_s] = ode45(@(t,X) ...
                 fn_SS_ode(t,X,0,WaveParam,FloeParam,Coeff,Mooring,Trans), ...
                 tspan_s,IC1s);
      % STITCH SOLUTION
      t1 = [tspan_t,tspan_s]';
      X1 = [X1_t;X1_s];         

      %  FLOE R
      %   TRANSIENT SOLUTION
      Trans.t = 1;
%       Trans.f.c = f.c + k*dist; % ### NO DRIFT
      [t2_t,X2_t] = ode45(@(t,X) ... 
                 fn_SS_ode(t,X,dist,WaveParam,FloeParam,Coeff,Mooring,Trans), ...
                 tspan_t,IC2);
      %   STEADY STATE SOLUTION          
      IC2s    = X2_t(end,:);
      Trans.t = 0;
%       Trans.f.c = f.c + k*dist; % ### NO DRIFT
      [t2_s,X2_s] = ode45(@(t,X) ...
                 fn_SS_ode(t,X,dist,WaveParam,FloeParam,Coeff,Mooring,Trans), ...
                 tspan_s,IC2s);
      % STITCH SOLUTION
      t2 = [tspan_t,tspan_s]';
      X2 = [X2_t;X2_s]; 
      
    elseif t1c(end) >= trans_e
     
     tspan_s = linspace(t1c(end),t1c(end)+60,10000);

     % SOLUTION FOR [TIME, [DISPLACEMENT,VELOCITY]]
     %  FLOE F
     %   STEADY STATE SOLUTION
     Trans.t = 0;
     Trans.f.c = f.c;
     WaveParam.H = A_ss*2;
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
      tspan = linspace(t1c(end),t1c(end)+30,10000);
      Trans.t = 0;
      Trans.f.c = f.c;
      WaveParam.H = A_ss*2;
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
    
    for j = 2:length(t1)
     
     [R1x,~,L2x,~] = FloePosition(k,a,r,omega,X1(j,1),X2(j,1),t1(j),t2(j));
     
     % IF OVERLAP
     if L2x - R1x <= 0 && j > 1000
 %      display('Collision!')
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
     
     % NO OVERLAP
     else
      cn = cn+1;
     end

    end

    c_num = -cc+n_cc; % total number of collisions
    
    % IF NO MORE COLLISIONS OVER TSPAN
    if cn == length(t1)-1
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
    display([' 2-Floe Model : ',num2str(Mod.Nc),' / ',...
                                num2str(Mod.Fc(1),'%0.4f'),' ',...
                                num2str(Mod.Fc(2),'%0.4f')])
    Exp.Fc = 1./Exp.Tc;             % collision frequency (experiments)
    display([' NR Experiment: ',num2str(Exp.Nc),' / ',...
                                num2str(Exp.Fc(1),'%0.4f'),' ',...
                                num2str(Exp.Fc(2),'%0.4f')])
                               
    cprintf('blue',['COLLISION VELOCITY [mm/s] (MEAN,MEDIAN) \n'])
    display([' 2-Floe Model : ',num2str(Mod.Vc(1)*1e3,'%0.1f'),' ',...
                                num2str(Mod.Vc(2)*1e3,'%0.1f')])
    display([' NR Experiment: ',num2str(Exp.Vc(1),'%0.1f'),' ',...
                                num2str(Exp.Vc(2),'%0.1f')])
                               
  
  
  

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

   % PLOT X VS TIME
   [R1x,~,L2x,~] = FloePosition(k,a,r,omega,X1(:,1),X2(:,1),t1,t2);
   figure
   set(gcf,'position',[100 400 1000 400]);
   set(gca,'FontSize',14)
   hold on
   plot(t1,R1x*1e3,'b-') % right edge of floe F
   plot(t2,L2x*1e3,'g-') % left edge of floe B
   % PLOT COLLISIONS
   if exist('ccj')
    for j = 1:length(ccj)
     plot(t1(ccj(j)),R1x(ccj(j))*1e3,'r.','MarkerSize',20)
    end
   end
   title(['Wavelength = ',num2str(lambda(count)),' m, Wave Amplitude = ',num2str(a),' m'])
   ylabel('x [mm]')
   xlabel('t [s]')
   ylim([min(R1x)-0.02 max(L2x)+0.02]*1e3)
   xlim([0 60])
   box on
   grid on
   hold off
   
   % overlay above figure on experimental data for surge
   figure(fig_sg)
   set(gcf,'position',[100 400 1000 400]);
   set(gca,'FontSize',14)
   hold on
   modelF = plot(t1,-(R1x-r)*1e3,'k--');
   modelB = plot(t2,-(L2x-r)*1e3,'m--');
   cprintf('magenta','Enter Phase Difference [s]:')
   phase = input(' ');
   delete(modelF)
   delete(modelB)
   plot(t1+phase,-(R1x-r)*1e3,'k--')
   plot(t2+phase,-(L2x-r)*1e3,'m--')
   if exist('ccj')
    plot(t1(ccj)+phase,-(R1x(ccj)-r)*1e3,'ro','MarkerSize',4)
   end

   % SIMULATE WAVES AND FLOE IN TIME
   if do_simul == 1;
    xspan = linspace(-0.8,1.2,1000);
    n = 1; % counter for collisions ccj(n)
    for j = 1:length(t1)
     t = t1(j);
     for jj = 1:length(xspan)
      x = xspan(jj);
      eta(jj) = a*sin(k*x-omega*t);
     end
     figure(5)
     clf
     hold on
     set(gcf,'position',[100 50 700 300]);
     set(gca,'FontSize',14)
     % wave profile
     plot(xspan,eta,'k-') 
     % floe
     G1 = k*a*cos(k*X1(j,1)-omega*t);
     G2 = k*a*cos(k*X2(j,1)-omega*t);
     O1 = r*G1/sqrt(1+G1^2);
     O2 = r*G2/sqrt(1+G2^2);
     A1 = O1/G1;
     A2 = O2/G2;
     % heave
     eta_X1 = a*sin(k*X1(j,1)-omega*t);
     eta_X2 = a*sin(k*X2(j,1)-omega*t);
     % angle
     theta1 = -atan(k*a*cos(k*X1(j,1)-omega*t));
     theta2 = -atan(k*a*cos(k*X2(j,1)-omega*t));
     % lower surface
     plot([X1(j,1)-A1 X1(j,1)+A1],[-O1+eta_X1 O1+eta_X1],'b')
     plot([X2(j,1)-A2 X2(j,1)+A2],[-O2+eta_X2 O2+eta_X2],'g')
     % upper surface
     plot([X1(j,1)-A1+D*sin(theta1) X1(j,1)+A1+D*sin(theta1)],...
          [-O1+eta_X1+D*cos(theta1)  O1+eta_X1+D*cos(theta1)],'b')
     plot([X2(j,1)-A2+D*sin(theta2) X2(j,1)+A2+D*sin(theta2)],...
          [-O2+eta_X2+D*cos(theta2)  O2+eta_X2+D*cos(theta2)],'g')
     % vertical surfaces
     plot([X1(j,1)-A1 X1(j,1)-A1+D*sin(theta1)],[-O1+eta_X1 -O1+eta_X1+D*cos(theta1)],'b')
     plot([X1(j,1)+A1 X1(j,1)+A1+D*sin(theta1)],[ O1+eta_X1  O1+eta_X1+D*cos(theta1)],'b')
     plot([X2(j,1)-A2 X2(j,1)-A2+D*sin(theta2)],[-O2+eta_X2 -O2+eta_X2+D*cos(theta2)],'g')
     plot([X2(j,1)+A2 X2(j,1)+A2+D*sin(theta2)],[ O2+eta_X2  O2+eta_X2+D*cos(theta2)],'g')

     grid on
     axis equal
     box on
     xlim([xspan(1) xspan(end)]), ylim([-0.3 0.3]) %ylim([-50*a 50*a])
     xlabel('x [m]');ylabel('z [m]')
     % show time
     annotation('textbox', [0.15,0.78,0.1,0.1],...
             'String', ['time = ',num2str(t,'%0.4f'),' [s]'],'LineStyle','none');
     % indicate collisions
     if exist('ccj')
      if j == ccj(n)
       annotation('textbox', [0.42,0.3,0.1,0.1],...
              'String', ['COLLISION!!!'],...
              'LineStyle','none','Color','r','FontSize',20);
       if n == length(ccj)
        n = n;
       else
        n = n+1;
       end
       pause
      end
     end
     hold off
     mov(j) = getframe(gcf);
  %    pause
    end
    movie2avi(mov,['Run',num2str(runno)])
   end
 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

toc

return


function [R1x,R1y,L2x,L2y] = FloePosition(k,a,r,omega,X1,X2,t1,t2)
% calculates floe elevation
 G1 = k*a.*cos(k.*X1-omega.*t1);
 G2 = k*a.*cos(k.*X2-omega.*t2);
 O1 = r.*G1./sqrt(1+G1.^2);
 O2 = r.*G2./sqrt(1+G2.^2);
 A1 = O1./G1;
 A2 = O2./G2;
 eta_X1 = a.*sin(k.*X1-omega.*t1);
 eta_X2 = a.*sin(k.*X2-omega.*t2);
 R1x = X1+A1; 
 R1y = O1+eta_X1; % right edge of floe F
 L2x = X2-A2;
 L2y = -O2+eta_X2; % left  edge of floe B
 
 
return
