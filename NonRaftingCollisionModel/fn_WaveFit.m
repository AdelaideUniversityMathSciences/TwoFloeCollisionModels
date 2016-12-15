% function [f] = fn_WaveFit(runno,ts,te,doplot,f_start)
%
% LJ YIEW
% Created on  May 2015
% Last edited Oct 2016
%
% Plots the wave profiles from the two-floe non-rafting experiments.
% Fits an equation to account for transient wave amplitudes.
%
% INPUTS:
% runno   = select which test (run number)
% ts      = transient start time
% te      = transient end time
% doplot  = flag for wave profile plot
% f_start = starting values for fitted parameters
%
% OUTPUTS:
% f = structure for f.a, f.b, f.c for the equation: f.a.*t.^f.b.*sin(-omega.*t + f.c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = fn_WaveFit(runno,ts,te,doplot,f_start)

cd('../TwoFloeExperiments/NonRaftingExperiments')

% Wave profile:
%  Ai = incident
%  Ap = phase
%  Af = floe
if ~exist('runno','var'); runno = 110; end
[time,~,~,~,Ap] = WaveProbes(runno,0);
Ap = -Ap; % change direction of x axis (to match SS model)
% close all

dat = xlsread('AMC_DataRAO','B6:L67');
run = dat(:,1); % run numbers

for j = 1:length(run)
 if run(j) == runno
  freq = dat(j,2); % frequency
  H = dat(j,3)*1e-3; % wave height
  [field] = wavefield('f',freq,0.831);
  lambda  = cell2mat(field(4,2)); % corresponding wavelength
  k       = 2*pi/lambda;
  omega   = 2*pi*freq;
 end
end


%%
if ~exist('ts','var'); ts = 0;  end  % start
if ~exist('te','var'); te = 12.32; end  % end
ts_c = ts/60*length(time)+1;
te_c = round(te/60*length(time));

% fit equation to curve
if ~exist('f_start','var'); f_start = [0.1 0 0]; end
f = fit(time(ts_c:te_c),Ap(ts_c:te_c),...
        ['a*x^b*sin(',num2str(-omega),'*x + c)'],...
        'Start',f_start); % ### +/-

warning('off','all')

% plot fitted curve
t = linspace(0,te,1000);
Af = f.a.*t.^f.b.*sin(-omega.*t + f.c);
ts = linspace(te,60,1000);
Afs = f.a*te^f.b*sin(-omega.*ts + f.c);

if ~exist('doplot','var'); doplot = 1; end
if doplot == 1
 figure
 set(gcf,'position',[100 400 1000 400]);
 set(gca,'FontSize',14)
 hold on
 plot(time,Ap,'b')
%  plot(t,Af,'r')
%  plot(ts,Afs,'g')
 xlabel('t [s]')
 ylabel('Wave Amplitude [mm]')
 title('Wave Profile from Non-Rafting Experiments')
 grid on
 box on
end

cd('../../NonRaftingCollisionModel')

return