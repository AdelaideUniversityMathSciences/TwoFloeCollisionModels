% function plot_ColsExp
%
% LJ YIEW
% Created on  Jun 2015
% Last edited Dec 2016
%
% Compares results from two-floe model and non-rafting experiments.
% Plots:
%  Mean collision frequency vs wavelength
%  Mean collision velocity vs wavelength
%  Empirical drag coefficient vs wavelength
%
% FILES NEEDED:
%  Data_ExperimentVsModel.xlsx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_ModVsExp

close all
clear all
clc


dat = xlsread('Data_ExperimentVsModel','B6:P67');

%%
% REMOVE DUD RUNS
jj = 1;
for j = 1:length(dat)
 if ~isnan(dat(j,8))
  gdat(jj,:) = dat(j,:);
  jj = jj+1;
 end
end
dat = gdat;

% ORGANIZE DATA
freq    = dat(:,2);        % wave frequency [Hz]
A       = dat(:,3)/2*1e-3; % wave amplitude [m]
At      = dat(:,4)/2*1e-3; % target wave amplitude [m]
ncol_e  = dat(:,6);        % number of collisions (qualisys)
ncol_m  = dat(:,7);        % number of collisions (model)
cfreq_e = dat(:,8);        % collision frequency (experiments)
cfreq_m = dat(:,9);        % collision frequency (model)
cvel_e  = dat(:,10);       % collision velocity (experiments)
cvel_m  = dat(:,11);       % collision velocity (model)
drag_e  = dat(:,13);       % drag coefficient (empirical)
pitch   = dat(:,14);       % pitch rao (from PF model)
c_rest  = dat(:,15);       % coefficient of restitution
 for j = 1:length(freq)    % wavelength [m]
  [field]   = wavefield('f',freq(j),0.831);
  lambda(j) = cell2mat(field(4,2));
 end
ka = 2*pi./lambda'.*A;
 
%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOT COLLISION FREQUENCY VS WAVELENGTH
%  EXPERIMENTS
figure(1)
hold on
set(gcf,'position',[100 400 600 400]);
box on
for j = 1:length(At)
 if At(j) == 10*1e-3
  plot(freq(j),cfreq_e(j),'kv','MarkerSize',14)
 elseif At(j) == 20*1e-3
  plot(freq(j),cfreq_e(j),'bo','MarkerSize',14)
 elseif At(j) == 40*1e-3
  plot(freq(j),cfreq_e(j),'r^','MarkerSize',14)
 end
end

%  MODEL
for j = 1:length(At)
 if At(j) == 10*1e-3
  plot(freq(j),cfreq_m(j),'k.','MarkerSize',20)
 elseif At(j) == 20*1e-3
  plot(freq(j),cfreq_m(j),'b.','MarkerSize',20)
 elseif At(j) == 40*1e-3
  plot(freq(j),cfreq_m(j),'r.','MarkerSize',20)
 end
end

plot([1.1405 1.1405],[0 0.2],'k--')

xlim([0.4 1.6])
ax1 = gca;
set(ax1,'XTick',[0.4 0.6 0.8 1 1.2 1.4 1.6])

% SECOND X SCALE FOR NONDIMENSIONAL WAVELENGTH
xlabels{1} = 'Wave Frequency [Hz]';
xlabels{2} = 'Wavelength/Floe Diameter';
ylabels{1} = 'Collision Frequency [Hz]';
ylabels{2} = '';
hold on
plotxx(0,0,0,0,xlabels,ylabels);
xlim([0.4 1.6])
ylim([0 0.2])
ax = gca;
set(ax,'XTick',[0.4 0.6 0.8 1 1.2 1.4 1.6])
set(ax,'YTick',[0 0.05 0.1 0.15 0.2])
set(ax,'XTickLabel',{'','9.5','5.9','3.9','2.7','2','1.5'})
set(ax,'YTickLabel',{'','','','',''})



%%

% PLOT COLLISION VELOCITY VS WAVELENGTH

figure(2)
hold on
set(gcf,'position',[300 100 600 400]);
box on
% EXPERIMENTS
for j = 1:length(At)
 if ~isnan(cvel_e(j))
  if At(j) == 10*1e-3
   plot(freq(j),cvel_e(j),'kv','MarkerSize',14)
  elseif At(j) == 20*1e-3
   plot(freq(j),cvel_e(j),'bo','MarkerSize',14)
  elseif At(j) == 40*1e-3
   plot(freq(j),cvel_e(j),'r^','MarkerSize',14)
  end
 end
end
% MODEL
for j = 1:length(At)
 if ~isnan(cvel_m(j))
  if At(j) == 10*1e-3
   plot(freq(j),cvel_m(j),'k.','MarkerSize',20)
  elseif At(j) == 20*1e-3
   plot(freq(j),cvel_m(j),'b.','MarkerSize',20)
  elseif At(j) == 40*1e-3
   plot(freq(j),cvel_m(j),'r.','MarkerSize',20)
  end
 end
end

plot([1.1405 1.1405],[0 200],'k--')
xlim([0.4 1.6])
ylim([0 200])
ax1 = gca;
set(ax1,'XTick',[0.4 0.6 0.8 1 1.2 1.4 1.6])
set(ax1,'XTickLabel',{'0.4','0.6','0.8','1','1.2','1.4','1.6'})


% SECOND X SCALE FOR NONDIMENSIONAL WAVELENGTH
xlabels{1} = 'Wave Frequency [Hz]';
xlabels{2} = 'Wavelength/Floe Diameter';
ylabels{1} = 'Collision Velocity [mm/s]';
ylabels{2} = '';
hold on
plotxx(0,0,0,0,xlabels,ylabels);
xlim([0.4 1.6])
ylim([0 0.2])
ax = gca;
set(ax,'XTick',[0.4 0.6 0.8 1 1.2 1.4 1.6])
set(ax,'XTickLabel',{'','9.5','5.9','3.9','2.7','2','1.5'})
set(ax,'YTickLabel',{'','','','',''})



%%

% % PLOT EMPIRICAL DRAG COEFFICIENT VS FREQUENCY
% 
% figure(3)
% hold on
% set(gcf,'position',[300 200 500 300]);
% box on
% % MODEL
% for j = 1:length(At)
%  if ~isnan(drag_e(j))
%   if At(j) == 10*1e-3
%    plot(freq(j),drag_e(j),'kv','MarkerSize',14)
%   elseif At(j) == 20*1e-3
%    plot(freq(j),drag_e(j),'bo','MarkerSize',14)
%   elseif At(j) == 40*1e-3
%    plot(freq(j),drag_e(j),'r^','MarkerSize',14)
%   end
%  end
% end
% xlim([0.4 1.6])
% ylim([0 0.05])
% ax1 = gca;
% set(ax1,'XTick',[0.4 0.6 0.8 1 1.2 1.4 1.6])
% set(ax1,'XTickLabel',{'0.4','0.6','0.8','1','1.2','1.4','1.6'})
% 
% 
% % SECOND X SCALE FOR NONDIMENSIONAL WAVELENGTH
% xlabels{1} = 'Wave Frequency [Hz]';
% xlabels{2} = 'Wavelength/Floe Diameter';
% ylabels{1} = 'Drag Coefficient';
% hold on
% plotxx(0,0,0,0,xlabels,ylabels);
% xlim([0.4 1.6])
% ylim([0 0.05])
% ax = gca;
% set(ax,'XTick',[0.4 0.6 0.8 1 1.2 1.4 1.6])
% set(ax,'XTickLabel',{'16.2','9.5','5.9','3.9','2.7','2','1.5'})
% set(ax,'YTickLabel',{'','','','',''})







