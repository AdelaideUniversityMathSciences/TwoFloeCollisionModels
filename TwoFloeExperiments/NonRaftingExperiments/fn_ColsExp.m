% function [res,out] = fn_ColsExp(what_tests,tt,d_factor,DO_SG,DO_PLT,DO_DSP,xtras,nanvec)
%
% LG BENNETTS & LJ YIEW
% Created on  Jun 2014
% Last edited Oct 2016
%
% Plots the time series motions in the x axis
% Calculates the mean collision frequency and velocity from the two-floe
% non-rafting experiments.
% Algorithm based on Bennetts et al. (2014) AFMC Proceedings.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [res,out] = fn_ColsExp(what_tests,tt,d_factor,DO_SG,DO_PLT,DO_DSP,xtras,nanvec)

%% DEFINE PROBLEM:

if ~exist('what_tests','var')
 what_tests = 'find(run==88)';
 %  what_tests = 'find(run==95)';
 %what_tests = 'find(H==20 & f==0.5)';
end

if ~exist('xtras','var')
 xtras = 'no vel';
end

if ~exist('nanvec','var')
 nanvec = [];
end


%% PRELIMINARIES

if ~exist('DO_SG','var');  DO_SG=0; end
if ~exist('DO_PLT','var'); DO_PLT=1; end
if ~exist('DO_DSP','var'); DO_DSP=1; end

%% PARAMETERS

% Radius  = 200; % [mm]
xsep    = 020; % [mm]
rbm     = 'surge';
epsilon = 1;   % [mm]
eps_h   = 5;   % [mm]
dt      = 400; % [1]
vt      = 10;  % [1]
xh_thr  = 0.5;   % [mm]
%%

tests = fn_findtests(what_tests,'strict');

if isempty(tests)
 cprintf('m','>>> No viable tests\n')
 return
end

%% INITIALISE

f     = zeros(length(tests(:,1)),1);
H     = zeros(length(tests(:,1)),1);
N_col = zeros(length(tests(:,1)),1);
T_col = zeros(length(tests(:,1)),2);
V_col = zeros(length(tests(:,1)),2);

for lp_f = 1:length(tests(:,1))
 
 run = tests(lp_f,1);
 run = sprintf('%03d',run);
 
 if DO_DSP == 1
  cprintf(0.4*[1,1,1],['>>> Run ',run,', Frequency = ',num2str(tests(lp_f,2)),...
   ' [Hz], Wave Height = ',num2str(tests(lp_f,3)),' [mm]\n'])
 end
 
 f(lp_f) = tests(lp_f,2);
 H(lp_f) = tests(lp_f,3);
 [field] = wavefield('f',f(lp_f),0.831);
 lambda(lp_f) = cell2mat(field(4,2));
 
 filename = ['RawData/Qualisys/Run0',run,'_6D.tsv'];
 
 out = fn_getdata(filename,['time ' rbm]);
 
 time = out.time; x = eval(['out.' rbm]); clear out
 
 x(:,2) = x(:,2)-xsep;
 
 % Calculate steady state period
 dist = 13; % approx dist between wave maker and floe (MTB length = 35m)
 c    = f*lambda;
 ss_s = dist/c;          % ss start time
 ss_e = 35+(35-dist)/c;  % ss end time
 
 %% CALCULATE SURGE RAO 
 % especially for runs without collisions
 
  if ~exist('tt','var')
   tt_s = ss_s;
   tt_e = ss_e;
  else
   tt_s = tt.s;
   tt_e = tt.e;
  end

  t_s = round(length(time)/time(end)*(tt_s));
  t_e = round(length(time)/time(end)*(tt_e));

  % smoothen data
  s_factor = 0.005;
  x_s(:,1) = smooth(time(t_s:t_e),x(t_s:t_e,1),s_factor,'lowess');
  x_s(:,2) = smooth(time(t_s:t_e),x(t_s:t_e,2),s_factor,'lowess');

  % remove drift
  if ~exist('d_factor','var')
   d_factor = 0.2;
  end
 %  d_factor = 0.6; % ###
  x_d(:,1) = smooth(time(t_s:t_e),x_s(:,1),d_factor,'lowess');
  x_d(:,2) = smooth(time(t_s:t_e),x_s(:,2),d_factor,'lowess');
  x_s(:,1) = x_s(:,1) - x_d(:,1);
  x_s(:,2) = x_s(:,2) - x_d(:,2) + x(1,2);


  % FLOE F
  % calculate local max and min
  c2 = 2; % counter for ss range
  d1 = 1; % counter for max surge
  d2 = 1; % counter for min surge
  for time_c = [t_s:t_e-2];
   % find local max
   if x_s(c2,1) > x_s(c2-1,1) && x_s(c2,1) > x_s(c2+1,1)
    x_max(d1) = x_s(c2,1);
    t_x_max(d1) = time(c2-1+t_s);
    d1 = d1+1;
   % find local min
   elseif x_s(c2,1) < x_s(c2-1,1) && x_s(c2,1) < x_s(c2+1,1)
    x_min(d2) = x_s(c2,1);
    t_x_min(d2) = time(c2-1+t_s);
    d2 = d2+1;
   end
    c2 = c2+1;
  end
  x_max_F = [x_max', t_x_max'];
  x_min_F = [x_min', t_x_min']; 
  clear x_max x_min
  % calculate average amplitudes
  x_avg_max = mean(x_max_F(:,1)); % average max
  x_avg_min = mean(x_min_F(:,1)); % average min
  x_avg     = mean(x_max_F(:,1)) - mean(x_min_F(:,1));
  x_max     = max(x_max_F(:,1)) - min(x_min_F(:,1));
  x_min     = min(x_max_F(:,1)) - max(x_min_F(:,1));
  x_rao_F   = [x_avg x_max x_min]/H/coth(2*pi/lambda*0.831);


  % FLOE B
  clear x_max t_x_max x_min t_x_min
  % calculate local max and min
  c2 = 2; % counter for ss range
  d1 = 1; % counter for max surge
  d2 = 1; % counter for min surge
  for time_c = [t_s:t_e-2];
   % find local max
   if x_s(c2,2) > x_s(c2-1,2) && x_s(c2,2) > x_s(c2+1,2)
    x_max(d1) = x_s(c2,2);
    t_x_max(d1) = time(c2-1+t_s);
    d1 = d1+1;
   % find local min
   elseif x_s(c2,2) < x_s(c2-1,2) && x_s(c2,2) < x_s(c2+1,2)
    x_min(d2) = x_s(c2,2);
    t_x_min(d2) = time(c2-1+t_s);
    d2 = d2+1;
   end
    c2 = c2+1;
  end
  x_max_B = [x_max', t_x_max'];
  x_min_B = [x_min', t_x_min'];
  clear x_max x_min
  % calculate average amplitudes
  x_avg_max = mean(x_max_B(:,1)); % average max
  x_avg_min = mean(x_min_B(:,1)); % average min
  x_avg     = mean(x_max_B(:,1)) - mean(x_min_B(:,1));
  x_max     = max(x_max_B(:,1)) - min(x_min_B(:,1));
  x_min     = min(x_max_B(:,1)) - max(x_min_B(:,1));
  x_rao_B   = [x_avg x_max x_min]/H/coth(2*pi/lambda*0.831);

  res.F = x_rao_F;
  res.B = x_rao_B;
  
 
 
 %% Plot raw data vs time
 if DO_PLT
  fig = 1;%fn_getfig;
  figure(fig); hold on
  set(gca,'FontSize',14,'box','on')
  if ~isempty(nanvec)
   eval(nanvec)
   inds = find(~isnan(x0(:,1)));
   plot(time(inds),x0(inds,1),'b')
   inds = find(~isnan(x0(:,2)));
   plot(time(inds),x0(inds,2),'g')
   clear x0 inds
  else
   % plot raw and smooth data
   plot(time,x(:,1),'b')
   plot(time,x(:,2),'g')
   if DO_SG == 1
    plot(time(t_s:t_e),x_s(:,1),'k--','LineWidth',1)
    plot(time(t_s:t_e),x_s(:,2),'k--','LineWidth',1)
    plot(time(t_s:t_e),x_d(:,1),'r--','LineWidth',1)
    plot(time(t_s:t_e),x_d(:,2),'r--','LineWidth',1)
    %    % plot local max and min
    plot(x_max_F(:,2),x_max_F(:,1),'bv','MarkerSize',6)
    plot(x_min_F(:,2),x_min_F(:,1),'b^','MarkerSize',6)
    plot(x_max_B(:,2),x_max_B(:,1),'gv','MarkerSize',6)
    plot(x_min_B(:,2),x_min_B(:,1),'g^','MarkerSize',6)
   end
   % Plot steady state window
   plot([ss_s ss_s],[min(x(:,2))*1.2 max(x(:,1))*1.2],'r--')
   plot([ss_e ss_e],[min(x(:,2))*1.2 max(x(:,1))*1.2],'r--')
%    plot([time(t_s) time(t_e)],[x_avg_max x_avg_max],'r--')
%    plot([time(t_s) time(t_e)],[x_avg_min x_avg_min],'r--')
   ylim([min(x(:,2))*1.2 max(x(:,1))*1.2])
   xlim([0 60])
  end
  ylabel('x [mm]'),xlabel('Time [s]')
  title(['Run ',run,', Frequency = ',num2str(tests(lp_f,2)),...
   ' [Hz], Wave Height = ',num2str(tests(lp_f,3)),' [mm]'])
 end

  
 
 %% Identify collisions
 
 % Filter technique: doesn't seem to work
 
 %   fmax     = 1/(time(2)-time(1))/2;
 %   f_filt   = 4*f(lp_f);% [Hz]
 %   [Bh,Ah]  = butter(3,f_filt/fmax,'high');
 %   xh(:,1)  = filter(Bh,Ah,x(:,1));
 %   xh(:,2)  = filter(Bh,Ah,x(:,2));
 %   inds     = find(xh(:,1)>=x_thresh & xh(:,2)>=x_thresh);
 %
 %   inds = find(xh>=x_thresh);
 %
 %   if ~isempty(inds)
 %    inds0=find(diff(inds)>dt);
 %    i0 = 1;
 %    for lp_c=1:length(inds0)
 %     i1=inds0(lp_c);
 %     ic(lp_c) = inds(i0);
 %     i0       = i1+1;
 %    end
 %    i1=length(inds);
 %    ic(length(inds0)+1) = inds(i0);
 %    inds0=ic; clear ic
 %    if DO_PLT
 %     plot(time(inds(inds0)),mean(x(inds(inds0),:),2),'r.','markersize',16)
 %    end
 %   end
 
 % collisions be distance
 
 dx   = x(:,1)-x(:,2);
 inds = find(dx<=epsilon);
 
 % collisions by relaxed distance and high-pass filter
 
 i_b = find(epsilon<dx & dx<=eps_h);
 
 i_b = i_b(find(i_b>1));
 
 fmax     = 1/(time(2)-time(1))/2;
 f_filt   = 4*f(lp_f);% [Hz]
 [Bh,Ah]  = butter(3,f_filt/fmax,'high');
 xh(:,1)  = filter(Bh,Ah,x(:,1));
 xh(:,2)  = filter(Bh,Ah,x(:,2));
 i_b     = intersect(i_b,find(abs(xh(:,1))>=xh_thr & abs(xh(:,2))>=xh_thr));
 
 clear xh fmax f_filt Ah Bh
 
 inds = unique([inds;i_b]); clear i_b
 
 if ~isempty(inds)
  inds0=find(diff(inds)>dt);
  i0 = 1;
  for lp_c=1:length(inds0)
   i1=inds0(lp_c);
   [~,ic(lp_c)]=min(dx(inds(i0:i1)));
   ic(lp_c) = ic(lp_c) + i0-1;
   i0=i1+1;
  end
  i1=length(inds);
  [~,ic(length(inds0)+1)]=min(dx(inds(i0:i1)));
  ic(length(inds0)+1) = ic(length(inds0)+1) + i0-1;
  inds0=ic; clear ic
  if and(DO_PLT,isempty(strfind(xtras,'no col')))
   plot(time(inds(inds0)),mean(x(inds(inds0),:),2),'r.','markersize',16)
  end
  
  N_col(lp_f)   = length(inds0);
  T_col(lp_f,:) = [mean(diff(time(inds(inds0)))), ...
   median(diff(time(inds(inds0))))];
  
 else
  
  N_col(lp_f)   = 0;
  T_col(lp_f,:) = [0,0];
  
 end
 
 if DO_DSP == 1
  cprintf('blue',['>>> numbers cols    = ' int2str(N_col(lp_f)) '\n'])
  cprintf('blue',['>>> mean col freq   = ' num2str(1/T_col(lp_f,1)) ' [Hz]\n'])
  cprintf('blue',['>>> median          = ' num2str(1/T_col(lp_f,2)) ' [Hz]\n'])
 end
 
 
 %% Collision velocity
 
%  % LUCAS' CODE:
%  col_tc = inds(inds0); % collision time count
%  step_t = diff(time);
%  for j = 1:length(col_tc)
%   
%   % pre/post collision velocities
%   col_u1(j) = (x(col_tc(j)-10,1) - x(col_tc(j),1)) / (step_t(1)*10);
%   col_v1(j) = (x(col_tc(j)+10,1) - x(col_tc(j),1)) / (step_t(1)*10);
%   
%   col_u2(j) = (x(col_tc(j)-10,2) - x(col_tc(j),2)) / (step_t(1)*10);
%   col_v2(j) = (x(col_tc(j)+10,2) - x(col_tc(j),2)) / (step_t(1)*10);
%   
%   col_u(j) = (col_u1(j)-col_u2(j))/2;  
%   col_v(j) = (col_v1(j)-col_v2(j))/2;  
%   
%  end
 
  
 % LUKE'S CODE:
 if ~isempty(inds)
  
  for lp_c=1:length(inds0)
   
   inds1  = inds(inds0(lp_c))-vt+1:inds(inds0(lp_c));
   inds11 = inds(inds0(lp_c)):inds(inds0(lp_c))+vt-1;
   
   dum1i = [ones(1,vt); time(inds1).'].'\x(inds1,1);
   dum2i = [ones(1,vt); time(inds1).'].'\x(inds1,2);
   
   dum1f = [ones(1,vt); time(inds11).'].'\x(inds11,1);
   dum2f = [ones(1,vt); time(inds11).'].'\x(inds11,2);
   
   inds2 = inds(inds0(lp_c))-4*vt+1:inds(inds0(lp_c));
   
   if and(DO_PLT,isempty(strfind(xtras,'no vel')))
    plot(time(inds2),dum1i(1)*(1+0*inds2)+dum1i(2)*time(inds2).','m')
    plot(time(inds2),dum2i(1)*(1+0*inds2)+dum2i(2)*time(inds2).','m')
   end
   
   spd_i(lp_c) = (-dum1i(2) + dum2i(2))/2;  % before collision
   spd_f(lp_c) = ( dum1f(2) - dum2f(2))/2;  % after collision
   
   cres_i(lp_c) = spd_f(lp_c)/spd_i(lp_c);
   
   
   clear dum1i dum2i dum1f dum2f
   
  end
  
  V_col(lp_f,1) = mean(spd_i);  V_col(lp_f,2) = median(spd_i);
  
  V0_col(lp_f,1) = spd_i(1);  % first collision
  
  cres(lp_f,1) = mean(cres_i);
  
  if DO_DSP == 1
   cprintf('blue',['>>> mean col vel    = ' num2str(V_col(lp_f,1)) ' [m/s]\n'])
   cprintf('blue',['>>> median          = ' num2str(V_col(lp_f,2)) ' [m/s]\n'])
   cprintf('blue',['>>> mean res coef   = ' num2str(cres(lp_f,1)) ])
  end
  
  clear inds1 inds2 spd_i spd_f
  
 else
  
  V_col(lp_f,1)=nan; V_col(lp_f,2)=nan; V0_col(lp_f,1)=nan; cres(lp_f)=nan;
  
 end
 
 clear inds inds0
 
end

%if DO_PLT; if length(tests(:,1))>1; cascade; end; end

%% OUTPUT

out.f  = f;
out.H  = H;
out.Nc = N_col;
out.Tc = T_col;
out.Vc = V_col;
out.V0 = V0_col;
out.CR = cres;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=fn_getdata(filename,outtyp)

% DO NOT CHANGE n!!! (corresponds to no. of columns if tsv file)
n = 34;

qualysis = fopen(filename,'rt');
data     = textscan(qualysis,'%f','HeaderLines',11);
data     = data{1};
fclose(qualysis);

data = reshape(data,n,length(data)/n).';

if strfind(outtyp,'frame')
 out.frame = data(:,1);
end
if strfind(outtyp,'time')
 out.time  = data(:,2);
end
if strfind(outtyp,'surge')
 out.surge = [data(:,3), data(:,19)];
end
if strfind(outtyp,'sway')
 out.sway  = [data(:,4), data(:,20)];
end
if strfind(outtyp,'heave')
 out.heave = [data(:,5), data(:,21)];
end
if strfind(outtyp,'roll')
 out.roll  = [data(:,6), data(:,22)];
end
if strfind(outtyp,'pitch')
 out.pitch = [data(:,7), data(:,23)];
end
if strfind(outtyp,'yaw')
 out.yaw  = [data(:,8), data(:,24)];
end
if strfind(outtyp,'res')
 out.res   = [data(:,9), data(:,25)];
end
if strfind(outtyp,'rot')
 out.rot   = [data(:,10:18), data(:,26:34)];
end

return

%

function tests = fn_findtests(what_tests,xtra)

if ~exist('xtra','var'); xtra = 'normal'; end

% test #, freq [Hz], Wave height [mm]

tests = ...
 [84, 1.5, 40, 0;
 85, 1.5, 40, 1;
 86, 1.4, 40, 1;
 87, 1.4, 40, 1;
 88, 1.3, 40, 1;
 89, 1.3, 40, 1;
 90, 1.2, 40, 1;
 91, 1.2, 40, 1;
 92, 1.1, 40, 1;
 93, 1.1, 40, 1;
 94, 1.0, 40, 1;
 95, 1.0, 40, 1;
 96, 0.75,40, 1;
 97, 0.75,40, 1;
 98, 0.5, 40, 1;
 99, 0.5, 40, 1;
 100,1.5, 20, 1;
 101,1.5, 20, 1;
 102,1.4, 20, 1;
 103,1.4, 20, 1;
 104,1.3, 20, 1;
 105,1.3, 20, 1;
 106,1.2, 20, 1;
 107,1.2, 20, 1;
 108,1.1, 20, 1;
 109,1.1, 20, 1;
 110,1.0, 20, 1;
 111,1.0, 20, 1;
 112,0.75,20, 1;
 113,0.75,20, 1;
 114,0.5, 20, 1;
 115,0  ,  0, 0;
 116,0.5, 20, 1;
 117,0.5, 80, 1;
 118,0.5, 80, 1;
 119,0.75,80, 1;
 120,0.75,80, 1;
 121,0.75,80, 1;
 122,   0, 0, 0;
 123, 1.0,80, 1;
 124,   0, 0, 0;
 125, 1.0,80, 1;
 126,   0, 0, 0;
 127, 1.1,80, 1;
 128, 1.1,80, 1;
 129, 1.2,80, 0;
 130,   0, 0, 0;
 131, 1.2,80, 0;
 132,   0, 0, 0;
 133, 1.2,80, 1;
 134, 1.0,80, 1;
 135, 1.3,80, 1;
 136, 1.3,80, 0;
 137, 1.3,80, 1;
 138,   0, 0, 0;
 139,   0, 0, 0;
 140, 1.4,80, 1;
 141, 1.4,80, 1;
 142,   0, 0, 0;
 143,   0, 0, 0;
 144,   0, 0, 0;
 145,   0, 0, 0];

run = tests(:,1);
f   = tests(:,2);
H   = tests(:,3);
dud = tests(:,4);

inds = eval(what_tests);

if strcmp(xtra,'normal')
 inds0 = find(H~=0);
 inds = intersect(inds,inds0);
elseif strcmp(xtra,'strict')
 inds0 = find(dud~=0);
 inds = intersect(inds,inds0);
end

tests = tests(inds,:);

return