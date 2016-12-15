% function [time,y1,y2,y3,y4] = WaveProbes(runno,doplot)
%
% LJ YIEW
% Created on  Jun 2014
% Last edited Oct 2016
%
% Extracts wave probe data from incident, phase, front and rear probes.
%
% INPUTS:
%  runno  = run number
%  doplot = flag (1/0) for wave profile plots
% 
% OUTPUTS:
%  time        = time data
%  y1,y2,y3,y4 = wave profiles at WPI,WPP,WPF,WPR
%
% FILES NEEDED:
%  AMC_Observations.xlsx
%  wavefield.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time,y1,y2,y3,y4] = WaveProbes(runno,doplot)


if ~exist('runno','var'); runno = 117; end
if ~exist('doplot','var'); doplot = 1; end
 
% read measured data for:
% frequency, wave height, wavelength, wave celerity (f*lambda)
run_dat = xlsread('AMC_Observations','A6:F67');


% USABLE RUNS:
% condition 1: waveheight = 0.02 m, freq = [0.5,2] Hz
% [6;7;8;9;10;11;12;13;14;33;34;35;37;38;67;68;69;70;71]'
% condition 2: waveheight = 0.04 m, freq = [0.5,2] Hz
% [15;16;17;18;19;20;21;22;23;24;39;41;42;43;44;74;77;78;79;80;81;82]'
% condition 3: waveheight = 0.08 m, freq = [0.5,1.5] Hz
% [25;26;27;28;29;30;31;32;83]'
% condition 4: waveheight = [0.01,0.08] m, freq = 1.5 Hz
% [45;46;47;48;49;50]'
% condition 5: waveheight = [0.01,0.1] m, freq = 1.25 Hz
% [51;52;53;54;55;56;57]'
% condition 6: waveheight = [0.005,0.04] m, freq = 1.8 Hz
% [58;59;60;61;62]'
% condition 7: waveheight = [0.01,0.08] m, freq = 1.25 Hz
% [63;64;65]'
% condition 8: waveheight = 0.04 m freq = [0.5,1.5] Hz
% [84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;]'
% condition 9: waveheight = 0.02 m, freq = [0.5,1.5] Hz
% [100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;116]'
% condition 10: waveheight = 0.08 m. freq = [0.5,1.5] Hz
% [117;118;119;120;121;123;125;127;128;129;133;134;135;137;138;140;141;144;145;]'


for run = [runno]' % <== INPUT run number(s)

 clear x_max x_min time_x_max time_x_min

 % adjust format for file name
 if run < 10
  run = ['0',num2str(run)];
 else
  run = [num2str(run)];
 end

 % read wave probe data
 filename  = ['RawData/WaveProbes/R',run,'-02_moving.tsv'];
 waveprobe = fopen(filename,'rt');
 data      = textscan(waveprobe,'%f','HeaderLines',22);
 waveprobe = fopen(filename,'rt');
 corr_data = textscan(waveprobe,'%f','HeaderLines',16);
 data      = data{1};
 corr_data = corr_data{1};
 fclose(waveprobe);
 
 % run parameters
 freq     = run_dat(str2num(run)-83,3);     % frequency
 w_height = run_dat(str2num(run)-83,4);     % wave height
 for j = 1:length(freq)
  [field]   = wavefield('f',freq(j),0.831);
  lambda(j) = cell2mat(field(4,2));
  c(j)      = freq(j)*lambda(j);
 end


 % arrange and extract data
 n = 5; % (corresponds to no. of columns)
 
 cfactor = corr_data(7:10); % correction factor for each WP
 zfactor = corr_data(2:5);  % zero factor for each WP
 
 % extract time
 time = data(1:n:end);
 y1   = (data(2:n:end)-zfactor(1))*cfactor(1);  % incident WP
 y2   = (data(3:n:end)-zfactor(2))*cfactor(2);  % phase WP
 y3   = (data(4:n:end)-zfactor(3))*cfactor(3);  % WP1 (back floe)
 y4   = (data(5:n:end)-zfactor(4))*cfactor(4);  % WP2 (front floe)
 
 
 % calculate steady state period
 phase_dist = 20;  % approx dist between wave machine and floe (MTB length = 35m)
 s_time = phase_dist/c;          % start time
 e_time = 35+(35-phase_dist)/c;  % end time
 
 
 if doplot == 1
  ymax = max(y1);
  figure
  set(gcf,'position',[200 200 1000 800]);
  % plot WP I (incident)
  subplot(4,1,1)
  hold on
  plot([s_time s_time],[-ymax/2*5 ymax/2*5],'r')
  plot([e_time e_time],[-ymax/2*5 ymax/2*5],'r')
  plot(time,y1)
  hold off
  set(gca,'FontSize',12)
  title({['Measured Wave Heights [mm],'],['Run ',num2str(run),...
         ', Frequency = ',num2str(freq),' Hz, Wave Height = ',...
         num2str(w_height),' mm']})
  ylabel('Incident')
  xlim([0 60]),ylim([-ymax/2*1.5 ymax/2*1.5])
  ylim([1.5*min(y1) 1.5*max(y1)])
  box on
  % plot WP P (phase)
  subplot(4,1,2)
  hold on
  plot([s_time s_time],[-ymax/2*5 ymax/2*5],'r')
  plot([e_time e_time],[-ymax/2*5 ymax/2*5],'r')
  plot(time,y2)
  hold off
  set(gca,'FontSize',12)
  ylabel('Phase')
  xlim([0 60]),ylim([-ymax/2*1.5 ymax/2*1.5])
  ylim([1.5*min(y2) 1.5*max(y2)])
  box on
  % plot WP F (front)
  subplot(4,1,3)
  hold on
  plot([s_time s_time],[-ymax/2*5 ymax/2*5],'r')
  plot([e_time e_time],[-ymax/2*5 ymax/2*5],'r')
  plot(time,y4)
  hold off
  set(gca,'FontSize',12)
  ylabel('Front')
  xlabel('Time [s]')
  xlim([0 60]),ylim([-ymax/2*1.5 ymax/2*1.5])
  ylim([1.5*min(y4) 1.5*max(y4)])
  box on
  % plot WP R (rear)
  subplot(4,1,4)
  hold on
  plot([s_time s_time],[-ymax/2*5 ymax/2*5],'r')
  plot([e_time e_time],[-ymax/2*5 ymax/2*5],'r')
  plot(time,y3)
  hold off
  set(gca,'FontSize',12)
  ylabel('Back')
  xlabel('Time [s]')
  xlim([0 60]),ylim([-ymax/2*1.5 ymax/2*1.5])
  ylim([1.5*min(y3) 1.5*max(y3)])
  box on
  
 end
 

end


