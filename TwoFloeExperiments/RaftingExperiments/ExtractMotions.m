% function ExtractMotions(condition)
%
% LJ YIEW
% Created on  Feb 2016
% Last edited Dec 2016
%
% Extracts the motions from the 6-DOF system (DUT rafting experiments)
% Plots the time series of motions in the
%  x axis (surge and drift)
%  z axis (heave)
%  rotations about the y axis (pitch)
% 
% FILES NEEDED:
%  DUT_Runsheet.xlsx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ExtractMotions(condition)

clear all
close all
clc
tic

% READ COLLISION, RAFTING, OVERWASH DATA
cro_dat = xlsread('DUT_Runsheet','E122:T181');
% READ WAVE PARAMETERS DATA
wave_dat = xlsread('DUT_Runsheet','E9:M177');


if ~exist('condition','var')
 condition = 101112;
end


switch condition
 case 2 % mooring tests
  runnos = [156;157;158;159;160]';
 case 10
  runnos = [95;96;97;98;99;100;101;102;103;104]';
  sep = 3;
 case 11
  runnos = [105;106;107;108;109;110;111;112;113;114]';
  sep = 3;
 case 12
  runnos = [115;116;117;118;119;120;121;122]';
  sep = 3;
 case 13
  runnos = [123;124;125;126;127;128;129;130;131;132]';
  sep = 6;
 case 14
  runnos = [133;134;135;136;137;138;139;140;141;142]';
  sep = 6;
 case 15
  runnos = [143;144;145;146;147;148;149;150]';
  sep = 6;
 case 101112 % all tests with separation 30mm
  runnos = 95:122;
  sep = 3;
 case 131415 % all tests with separation 60mm
  runnos = 123:150;
  sep = 6;
end

c1 = 1; % counter for run number

for run =  runnos 
 
 run = run
 
 close all
 clearvars -except condition run runnos c1 ss_dat cro_dat ...
  X_avg Z_avg Yrot_avg ...
  X_err Z_err Yrot_err sep ...
  col_rT col_cF col_cN ...
  plot_raos which_floe

 runno  = num2str(run);

 % READ DATA
 % floe 1
 filename = ['RawData/6-DOF/',runno,'_1_1.txt'];
 dof1      = fopen(filename,'rt');
 data1     = textscan(dof1,'%f');
 data1     = data1{1};
 fclose(dof1);
 % floe 2
 if condition == 1 || condition ==2
 else
  filename = ['RawData/6-DOF/',runno,'_1_2.txt'];
  dof2      = fopen(filename,'rt');
  data2     = textscan(dof2,'%f');
  data2     = data2{1};
  fclose(dof2);
 end

 % ARRANGE AND EXTRACT DATA
 n = 7;   %(corresponds to no. of columns)
 time = linspace(0,40.96,2048);
  X1   = data1(2:n:end)-data1(2); % units: cm 
  Y1   = data1(3:n:end);  %
  Z1   = data1(4:n:end);  % 
  Yrot1 = data1(5:n:end); % units: degrees
  Xrot1 = data1(6:n:end); %
  Zrot1 = data1(7:n:end); % 
 if exist('data2')
  X2   = data2(2:n:end)-data2(2)-sep; % units: cm
  Y2   = data2(3:n:end);  % 
  Z2   = data2(4:n:end);  % 
  Yrot2 = data2(5:n:end); % units: degrees
  Xrot2 = data2(6:n:end); %
  Zrot2 = data2(7:n:end); %
 end
 % set units to mm and degrees
  X1 = X1*10;
  Y1 = Y1*10;
  Z1 = Z1*10;
 if exist('X2')
  X2 = X2*10;
  Y2 = Y2*10; 
  Z2 = Z2*10;
 end
 
%%
 if condition == 2
 else
  % APPLY FILTER TO ADJUST SEPARATION DISTANCE
  for j = 1:length(cro_dat)
   if cro_dat(j,1) == run
    % flags for:
    flag_col = cro_dat(j,14); % collisions
    flag_raf = cro_dat(j,15); % rafting
    flag_ove = cro_dat(j,16); % overwash
   end
  end

  % for runs without collisions
  if flag_col == 0
   loop = 1;
   adjust = 0;
   while loop == 1
    for j = 1:length(X1) 
     if X1(j) - X2(j) < 1.2 % if there is too much overlap
      X2 = X2 - 1;
      loop = 1;
      adjust = 1;
      break
     end
     if j == length(X1)
      if adjust == 0 && min(X1)-min(X2) > sep % too far apart
       X2 = X2 + 1;
       break
      else
       loop = 0; % stop looping
      end
     end
    end
   end
  end

  % for runs with collisions but without rafting
  if flag_col == 1 && flag_raf == 0 
   loop = 1;
   adjust = 0;
   while loop == 1
    for j = 1:length(X1) 
     if X1(j) - X2(j) < -2 % if there is too much overlap
      X2 = X2 - 0.1;
      loop = 1;
      adjust = 1;
      break
     end
     if j == length(X1)
      if adjust == 0 % if they are too far apart
       X2 = X2 + 0.1;
       break
      elseif adjust == 1
       loop = 0;
      end
     end
    end
   end
  end

  % fine tuning
  if run == 99
   X2 = X2 - 0.8;
  elseif run == 96
   X2 = X2 - 3.3;
  elseif run == 107
   X2 = X2 + 10.2;
  elseif run == 108
   X2 = X2 + 0.4;
  elseif run == 116
   X2 = X2 - 64.2;
  elseif run == 123
   X2 = X2 + 10;
  elseif run == 134
   X2 = X2 + 14.5;
  elseif run == 135
   X2 = X2 + 7;
  elseif run == 140
   X2 = X2 + 4;
  elseif run == 145;
   X2 = X2 - 23;
  elseif run == 149
   X2 = X2 + 17;
  end
    
 
%%   

  % DETERMINE RAFTING EVENT
  if flag_col == 1
   c = 1;
   for j = 1:length(X1)
    if X1(j) - X2(j) <= 1 % collision defined when sep <= 1
     col_r(c,:) = [time(j) X1(j) X2(j)];
     c = c+1;
    end
   end
   % DETERMINE COLLISION EVENT
   % filter events, exclude events if they are within 0.1 seconds of another
   c = 1;
   for j = 1:length(col_r)
    if j == 1
     col_c(c,:) = col_r(j,:);
     c = c+1;
    else
     if col_r(j,1) - col_r(j-1,1) > 0.1
      col_c(c,:) = col_r(j,:);
      c = c+1;
     end
    end
   end
  end

  if flag_col == 1
   % CALCULATE RAFTING DURATION
   c = 1;
   for j = 1:length(col_r)
    if col_r(j,1) - col_r(1,1) <= 30 % limit to 30s interval
     dum(c,:) = col_r(j,:);
     c = c+1;
    end
   end
   col_rT(c1,1) = size(dum,1)*0.02;

   % CALCULATE MEAN COLLISION FREQUENCY
   col_cF(c1,1) = 1/mean(diff(col_c(:,1)));
   % NUMBER OF COLLISIONS
   c = 1; clear dum
   for j = 1:size(col_c,1)
    if col_c(j,1) - col_c(1,1) <= 30 % limit to 30s interval
     dum(c,:) = col_c(j,:);
     c = c+1;
    end
   end
   col_cN(c1,1) = size(dum,1);
  else
   col_rT(c1,1) = 0;
   col_cF(c1,1) = 0;
   col_cN(c1,1) = 0;
  end
 
 end

 %%
 
 % PLOT RAW DATA VS TIME
 f1 = figure(1);
 set(gcf,'position',[100 600 500 400])
 hold on
 plot(time,X1,'b')
 if exist('X2')
  plot(time,X2,'k')
 end
 if exist('flag_col')
  if flag_col == 1
   for j = 1:length(col_r)
    plot([col_r(j,1) col_r(j,1)],[col_r(j,2) col_r(j,3)],'r:')
   end
  plot(col_c(:,1),col_c(:,2),'r.','MarkerSize',20)
  end
 end
 ylabel('X [mm]');xlabel('t [s]')
 xlim([0 time(end)])
 %
 f2 = figure(2);
 set(gcf,'position',[1100 600 500 400])
 hold on
 plot(time,Z1,'g')
 if exist('Z2')
  plot(time,Z2,'k')
 end
 ylabel('Z [mm]');xlabel('t [s]')
 xlim([0 time(end)])
 %
 f3 = figure(3);
 set(gcf,'position',[600 100 500 400])
 hold on
 plot(time,Yrot1,'r')
 if exist('Yrot2')
  plot(time,Yrot2,'k')
 end
 ylabel('\beta [deg]');xlabel('t [s]')
 xlim([0 time(end)])
 ylim([-max(abs([min(Yrot1) max(Yrot1)])) max(abs([min(Yrot1) max(Yrot1)]))]*1.2)
 %
 f4 = figure(4);
 set(gcf,'position',[600 600 500 400])
 hold on
 plot(time,Y1,'r')
 if exist('Y2')
  plot(time,Y2,'k')
 end
 ylabel('Y [mm]');xlabel('t [s]')
 xlim([0 time(end)])
 %
 f5 = figure(5);
 set(gcf,'position',[100 100 500 400])
 hold on
 plot(time,Xrot1,'b')
 if exist('Xrot2')
  plot(time,Xrot2,'k')
 end
 ylabel('\alpha [deg]');xlabel('t [s]')
 xlim([0 time(end)])
 ylim([-max(abs([min(Yrot1) max(Yrot1)])) max(abs([min(Yrot1) max(Yrot1)]))]*1.2)
 %
 f6 = figure(6);
 set(gcf,'position',[1100 100 500 400])
 hold on
 plot(time,Zrot1,'g')
 if exist('Zrot2')
  plot(time,Zrot2,'k')
 end
 ylabel('\gamma [deg]');xlabel('t [s]')
 xlim([0 time(end)])
 ylim([-max(abs([min(Yrot1) max(Yrot1)])) max(abs([min(Yrot1) max(Yrot1)]))]*1.2)

 c1 = c1+1;
  
end

toc

