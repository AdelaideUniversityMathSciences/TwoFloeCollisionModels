% function ExtractWP(condition)
%
% LJ YIEW
% Created on  Feb 2016
% Last edited Dec 2016
%
% Extracts the wave probe data (DUT rafting experiments)
% Plots the time series of wave profiles at
%  incident probe
%  in front of floes
%  behind floes
% Calculates wave amplitudes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ExtractWP(condition)
%%
clear all
close all
% clc

if ~exist('condition','var')
 condition = 131415;
end

switch condition
 case 2 % mooring tests
  runnos = [156;158;159;160]';
 case 10
  runnos = [95;96;97;98;99;100;101;102;103;104]';
 case 11
  runnos = [105;106;107;108;109;110;111;112;113;114]';
 case 12
  runnos = [115;116;117;118;119;120;121;122]';
 case 13
  runnos = [123;124;125;126;127;128;129;130;131;132]';
 case 14
  runnos = [133;134;135;136;137;138;139;140;141;142]';
 case 15
  runnos = [143;144;145;146;147;148;149;150]';
 case 101112
  runnos = 95:122;
 case 131415
  runnos = 123:150;
end

c1 = 1; % counter for run number

for run = 123%runnos 
 
 run = run

 close all

 runno  = num2str(run);

 % READ DATA
 filename = ['RawData/WaveProbes/',runno,'.txt'];
 wp       = fopen(filename,'rt');
 data     = textscan(wp,'%f','HeaderLines',2);
 data     = data{1};
 fclose(wp);

 % ARRANGE AND EXTRACT DATA
 n = 4;   % (corresponds to no. of columns)
 if run == 2
  n = 3;
 end
 time = linspace(0,40.96,2048);
 % extract data
 WPI    = data(1:n:end); % units: cm
 WPF    = data(2:n:end); % 
 WPB    = data(3:n:end); % 
 WP2    = data(4:n:end); % 
 % set units to mm and degrees
 WPI = WPI(1:length(WPI)-1)*10;
 WPF = WPF(1:length(WPF)-1)*10;
 WPB = WPB(1:length(WPB)-1)*10;
 WP2 = WP2(1:length(WP2)-1)*10;


 %%
 % PLOT RAW DATA VS TIME
 %
 f1 = figure(1);
 title(['Run ',runno])
 set(gcf,'position',[400 200 1200 800]);
 subplot(4,1,1)
 plot(time,WPI,'b')
 ylabel('Incident [mm]')
 xlim([0 time(end)])
 ylim([min([min(WPI) min(WPF) min(WP2) min(WPB)]) max([max(WPI) max(WPF) max(WP2) max(WPB)])]*1.5)
 %
 subplot(4,1,2)
 plot(time,WPF,'b')
 ylabel('Front 1 [mm]')
 xlim([0 time(end)])
 ylim([min([min(WPI) min(WPF) min(WP2) min(WPB)]) max([max(WPI) max(WPF) max(WP2) max(WPB)])]*1.5)
 %
 subplot(4,1,3)
 plot(time,WP2,'b')
 ylabel('Front 2 [mm]')
 xlim([0 time(end)])
 ylim([min([min(WPI) min(WPF) min(WP2) min(WPB)]) max([max(WPI) max(WPF) max(WP2) max(WPB)])]*1.5)
 %
 subplot(4,1,4)
 plot(time,WPB,'b')
 ylabel('Back [mm]');xlabel('t [s]')
 xlim([0 time(end)])
 ylim([min([min(WPI) min(WPF) min(WP2) min(WPB)]) max([max(WPI) max(WPF) max(WP2) max(WPB)])]*1.5)
 

 %%
 % SELECT START/END TIME FOR DATA PROCESSING
  starttime = 25;
  endtime   = 40;
 % counter for start and end time i.e. corresponding to frame no.
 t_start = round(length(time)/time(end)*starttime);
 t_end   = round(length(time)/time(end)*endtime);

 % SMOOTHEN DATA (FROM START TO END TIME)
 s_factor    = 0.01;
 WPI_smooth = smooth(time(t_start:t_end),WPI(t_start:t_end),s_factor,'lowess');
 WPF_smooth = smooth(time(t_start:t_end),WPF(t_start:t_end),s_factor,'lowess');
 WPB_smooth = smooth(time(t_start:t_end),WPB(t_start:t_end),s_factor,'lowess');
 WP2_smooth = smooth(time(t_start:t_end),WP2(t_start:t_end),s_factor,'lowess');
 
 % PLOT SMOOTH DATA
 figure(1)
 subplot(4,1,1)
 hold on
 plot(time(t_start:t_end),WPI_smooth,'k--')
 subplot(4,1,2)
 hold on
 plot(time(t_start:t_end),WPF_smooth,'k--')
 subplot(4,1,3)
 hold on
 plot(time(t_start:t_end),WP2_smooth,'k--')
 subplot(4,1,4)
 hold on
 plot(time(t_start:t_end),WPB_smooth,'k--')

 %%
 % FIND LOCAL MAX AND MIN
 % for WPI
 c2 = 2; % counter for ss range
 d1 = 1; % counter for max surge
 d2 = 1; % counter for min surge
 for time_c = [t_start:t_end-2];
  % find local max
  if WPI_smooth(c2) > WPI_smooth(c2-1) && WPI_smooth(c2) > WPI_smooth(c2+1)
   WPI_max(d1) = WPI_smooth(c2);
   time_WPI_max(d1) = time(c2-1+t_start);
   d1 = d1+1;
  % find local min
  elseif WPI_smooth(c2) < WPI_smooth(c2-1) && WPI_smooth(c2) < WPI_smooth(c2+1)
   WPI_min(d2) = WPI_smooth(c2);
   time_WPI_min(d2) = time(c2-1+t_start);
   d2 = d2+1;
  end
   c2 = c2+1;
 end
 % for WPF
 c2 = 2; % counter for ss range
 d1 = 1; % counter for max surge
 d2 = 1; % counter for min surge
 for time_c = [t_start:t_end-2];
  % find local max
  if WPF_smooth(c2) > WPF_smooth(c2-1) && WPF_smooth(c2) > WPF_smooth(c2+1)
   WPF_max(d1) = WPF_smooth(c2);
   time_WPF_max(d1) = time(c2-1+t_start);
   d1 = d1+1;
  % find local min
  elseif WPF_smooth(c2) < WPF_smooth(c2-1) && WPF_smooth(c2) < WPF_smooth(c2+1)
   WPF_min(d2) = WPF_smooth(c2);
   time_WPF_min(d2) = time(c2-1+t_start);
   d2 = d2+1;
  end
   c2 = c2+1;
 end
 % for WPB
 c2 = 2; % counter for ss range
 d1 = 1; % counter for max surge
 d2 = 1; % counter for min surge
 for time_c = [t_start:t_end-2];
  % find local max
  if WPB_smooth(c2) > WPB_smooth(c2-1) && WPB_smooth(c2) > WPB_smooth(c2+1)
   WPB_max(d1) = WPB_smooth(c2);
   time_WPB_max(d1) = time(c2-1+t_start);
   d1 = d1+1;
  % find local min
  elseif WPB_smooth(c2) < WPB_smooth(c2-1) && WPB_smooth(c2) < WPB_smooth(c2+1)
   WPB_min(d2) = WPB_smooth(c2);
   time_WPB_min(d2) = time(c2-1+t_start);
   d2 = d2+1;
  end
   c2 = c2+1;
 end
 % for WP2
 c2 = 2; % counter for ss range
 d1 = 1; % counter for max surge
 d2 = 1; % counter for min surge
 for time_c = [t_start:t_end-2];
  % find local max
  if WP2_smooth(c2) > WP2_smooth(c2-1) && WP2_smooth(c2) > WP2_smooth(c2+1)
   WP2_max(d1) = WP2_smooth(c2);
   time_WP2_max(d1) = time(c2-1+t_start);
   d1 = d1+1;
  % find local min
  elseif WP2_smooth(c2) < WP2_smooth(c2-1) && WP2_smooth(c2) < WP2_smooth(c2+1)
   WP2_min(d2) = WP2_smooth(c2);
   time_WP2_min(d2) = time(c2-1+t_start);
   d2 = d2+1;
  end
   c2 = c2+1;
 end

 % ARRANGE MAX/MIN PEAK DATA AND CORRESPONDING TIME
 % for WPI
 WPI_max = [WPI_max', time_WPI_max'];
 WPI_min = [WPI_min', time_WPI_min']; 
 % for WPF
 WPF_max = [WPF_max', time_WPF_max'];
 WPF_min = [WPF_min', time_WPF_min']; 
 % for WPB
 WPB_max = [WPB_max', time_WPB_max'];
 WPB_min = [WPB_min', time_WPB_min']; 
 % for WP2
 WP2_max = [WP2_max', time_WP2_max'];
 WP2_min = [WP2_min', time_WP2_min']; 

 % CALCULATE AVERAGE AMPLITUDES
 % for incident probe
 WPI_avg_max   = mean(WPI_max(:,1)); % average max
 WPI_avg_min   = mean(WPI_min(:,1)); % average min
 WPI_avg(c1,1) = mean(WPI_max(:,1)) - mean(WPI_min(:,1));
 % for front probe
 WPF_avg_max   = mean(WPF_max(:,1)); % average max
 WPF_avg_min   = mean(WPF_min(:,1)); % average min
 WPF_avg(c1,1) = mean(WPF_max(:,1)) - mean(WPF_min(:,1));
 % for back probe
 WPB_avg_max   = mean(WPB_max(:,1)); % average max
 WPB_avg_min   = mean(WPB_min(:,1)); % average min
 WPB_avg(c1,1) = mean(WPB_max(:,1)) - mean(WPB_min(:,1));
 % for front probe (2)
 WP2_avg_max   = mean(WP2_max(:,1)); % average max
 WP2_avg_min   = mean(WP2_min(:,1)); % average min
 WP2_avg(c1,1) = mean(WP2_max(:,1)) - mean(WP2_min(:,1));
 
 % CALCULATE AVERAGE WAVE PERIODS
 T_WPI(c1,1) = mean([mean(diff(WPI_max(:,2))) mean(diff(WPI_min(:,2)))]);
 T_WPF(c1,1) = mean([mean(diff(WPF_max(:,2))) mean(diff(WPF_min(:,2)))]);
 T_WPB(c1,1) = mean([mean(diff(WPB_max(:,2))) mean(diff(WPB_min(:,2)))]);
 T_WP2(c1,1) = mean([mean(diff(WP2_max(:,2))) mean(diff(WP2_min(:,2)))]);

 % PLOT LOCAL MAX/MIN
 figure(1)
 subplot(4,1,1)
 hold on
 plot(WPI_max(:,2),WPI_max(:,1),'rv','MarkerSize',6)
 plot(WPI_min(:,2),WPI_min(:,1),'r^','MarkerSize',6)
 plot([time(t_start) time(t_end)],[WPI_avg_max WPI_avg_max],'r--')
 plot([time(t_start) time(t_end)],[WPI_avg_min WPI_avg_min],'r--') 
 xlim([0 time(end)])
 title({['T = ',num2str(T_WPI(c1)),' s, ','H = ',num2str(WPI_avg(c1)),' mm']})
 subplot(4,1,2)
 hold on
 plot(WPF_max(:,2),WPF_max(:,1),'rv','MarkerSize',6)
 plot(WPF_min(:,2),WPF_min(:,1),'r^','MarkerSize',6)
 plot([time(t_start) time(t_end)],[WPF_avg_max WPF_avg_max],'r--')
 plot([time(t_start) time(t_end)],[WPF_avg_min WPF_avg_min],'r--') 
 xlim([0 time(end)])
 title({['T = ',num2str(T_WPF(c1)),' s, ','H = ',num2str(WPF_avg(c1)),' mm']})
 subplot(4,1,3)
 hold on
 plot(WP2_max(:,2),WP2_max(:,1),'rv','MarkerSize',6)
 plot(WP2_min(:,2),WP2_min(:,1),'r^','MarkerSize',6)
 plot([time(t_start) time(t_end)],[WP2_avg_max WP2_avg_max],'r--')
 plot([time(t_start) time(t_end)],[WP2_avg_min WP2_avg_min],'r--') 
 xlim([0 time(end)])
 title({['T = ',num2str(T_WP2(c1)),' s, ','H = ',num2str(WP2_avg(c1)),' mm']})
 subplot(4,1,4)
 hold on
 plot(WPB_max(:,2),WPB_max(:,1),'rv','MarkerSize',6)
 plot(WPB_min(:,2),WPB_min(:,1),'r^','MarkerSize',6)
 plot([time(t_start) time(t_end)],[WPB_avg_max WPB_avg_max],'r--')
 plot([time(t_start) time(t_end)],[WPB_avg_min WPB_avg_min],'r--') 
 xlim([0 time(end)])
 title({['T = ',num2str(T_WPB(c1)),' s, ','H = ',num2str(WPB_avg(c1)),' mm']})
 
 %%
 % USE FFT TO CALCULATE WAVE HEIGHTS AND FREQUENCIES
  
 [out] = FFT(WPI(t_start:t_end),50,1/T_WPI(c1));
 WPI_fft(c1,1) = out(1)*2;
 [out] = FFT(WPF(t_start:t_end),50,1/T_WPF(c1));
 WPF_fft(c1,1) = out(1)*2;
 [out] = FFT(WP2(t_start:t_end),50,1/T_WP2(c1));
 WP2_fft(c1,1) = out(1)*2;
 [out] = FFT(WPB(t_start:t_end),50,1/T_WPB(c1));
 WPB_fft(c1,1) = out(1)*2;
 
 
 %%
 
%  % PLOT SINUSOIDAL WAVE
%  % calculate wavenumber
%  [field] = wavefield('f',1/T_WP2(c1),0.5);
%  k       = cell2mat(field(5,2));
%  Wsin_I = WPI_avg(c1)/2.* cos(2*pi/T_WPI(c1).*time) - (WPI_avg(c1)/2 - WPI_avg_max);
%  Wsin_F = WPF_avg(c1)/2.* cos(2*pi/T_WPF(c1).*time) - (WPF_avg(c1)/2 - WPF_avg_max);
%  Wsin_2 = WP2_avg(c1)/2.* cos(2*pi/T_WP2(c1).*time) - (WP2_avg(c1)/2 - WP2_avg_max);
%  Wsin_B = WPB_avg(c1)/2.* cos(2*pi/T_WPB(c1).*time) - (WPB_avg(c1)/2 - WPB_avg_max);
%  figure(1)
%  subplot(4,1,1)
%  hold on
%  plot(time+time_WPI_max(1),Wsin_I,'c')
%  subplot(4,1,2)
%  hold on
%  plot(time+time_WPF_max(1),Wsin_F,'c')
%  subplot(4,1,3)
%  hold on
%  plot(time+time_WP2_max(1),Wsin_2,'c')
%  subplot(4,1,4)
%  hold on
%  plot(time+time_WPB_max(1),Wsin_B,'c')
 
 
 clearvars -except condition run runnos c1 ss_dat ...
  WPI_avg WPF_avg WPB_avg WP2_avg T_WPI T_WPF T_WPB ...
  WPI_fft WPF_fft WPB_fft WP2_fft
 
 c1 = c1+1;
 
end

