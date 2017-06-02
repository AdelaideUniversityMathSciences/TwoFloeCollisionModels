%% Parameterise wave scattering and drift forces
% Derive forces from non-rafting collision experiments:
% Runs  : 85,  100, 101, 102, 135, 137, 140, 141
% f[Hz] = 1.5, 1.5, 1.5, 1.4, 1.3, 1.3, 1.4, 1.4
% H[mm] = 40,  20,  20,  20,  80,  80,  80,  80

% run_2Floe_SD(runno,Cd,Fs,Fd,tw) % nb. tw=approx transient start time

% tw = 9.2; ph = 4;
% % run_2Floe_SD(85,0.01,0.03,0.085,tw)
% run_2Floe_SD(85,0.05,0.03,0.085,tw)

% tw = 8.3; ph = 4;
% run_2Floe_SD(86,0.05,0.007,0.055,tw)

% tw = 5; ph = 4;
% run_2Floe_SD(87,0.05,0.007,0.055,tw)

tw = 13; ph = 3; 
% run_2Floe_SD(100,0.001,0.004,0.009,tw)
run_2Floe_SD(100,0.05,0.004,0.009,tw)

% tw = 13; ph = 3; 
% run_2Floe_SD(101,0.001,0.003,0.01,tw)

% tw = 11; ph = 0; 
% run_2Floe_SD(102,0.001,0.003,0.006,tw)

% tw = 20; ph = 2; 
% run_2Floe_SD(135,0.05,0.021,0.14,tw)

% tw = 20.5; ph = 2; 
% run_2Floe_SD(141,0.02,0.04,0.2,tw)

%
hold on
t = linspace(0,60,1000);
% plot(t+ph,-0.2.*1e3./(1+exp(4.*(-t+tw))),'r')
plot(t+ph,-0.2.*1e3./(1+exp(10.*(-t+tw+0))),'k')
ylim('auto')

%% Plot parameterised drift and scattering forces

dat = [85,  100, 101, 102, 135, 141; ... % run number
       1.5, 1.5, 1.5, 1.4, 1.3, 1.4; ... % frequency
       40,  20,  20,  20,  80,  80;  ... % wave height
       0.03, 0.004, 0.003, 0.003, 0.021, 0.04; ... % scattering force
       0.085, 0.009, 0.01, 0.006, 0.14, 0.2;   ... % drift force
       0.0197, 0.0049, 0.0049, 0.00020104, 0.0046, 0.0032]; % mauro's drift
       
       
figure      
hold on      
for j = 1:length(dat(3,:))
 if dat(3,j) == 20
  plot(dat(2,j),dat(4,j),'bv','MarkerSize',14)
 elseif dat(3,j) == 40
  plot(dat(2,j),dat(4,j),'bo','MarkerSize',14)
 elseif dat(3,j) == 80 
  plot(dat(2,j),dat(4,j),'b^','MarkerSize',14)
 end
end
ylabel('Scattering Force [N]')
xlabel('f [Hz]')
xlim([1.2 1.6])

figure      
hold on      
for j = 1:length(dat(3,:))
 if dat(3,j) == 20
  plot(dat(2,j),dat(5,j),'rv','MarkerSize',14)
  plot(dat(2,j),dat(6,j),'kv','MarkerSize',14)
 elseif dat(3,j) == 40
  plot(dat(2,j),dat(5,j),'ro','MarkerSize',14)
  plot(dat(2,j),dat(6,j),'ko','MarkerSize',14)
 elseif dat(3,j) == 80 
  h1=plot(dat(2,j),dat(5,j),'r^','MarkerSize',14);
  h2=plot(dat(2,j),dat(6,j),'k^','MarkerSize',14);
 end
end
ylabel('Drift Force [N]')
xlabel('f [Hz]')
xlim([1.2 1.6])
legend([h1 h2],'Empirical','Mauros Formula')

%% Check drift velocity using drift force *no scattering force, no mooring

% run_2Floe_D(runno,Cd,Fs,Fd,tw)

tw = 9.5; ph = 4;
% run_2Floe_D(85,0.05,0,0.085,tw)
run_2Floe_D(85,0.05,0,0.0197,tw) % using mauros drift force

% tw = 13; ph = 3; % change K = 0.3; C = 0.2;
% % run_2Floe_D(100,0.05,0,0.009,tw)
% run_2Floe_D(100,0.05,0,0.0049,tw) % using mauros drift force

% tw = 13; ph = 3; % change K = 0.3; C = 0.2;
% run_2Floe_D(101,0.05,0,0.01,tw)
% % run_2Floe_D(101,0.001,0,0.0049,tw) % using mauros drift force

% tw = 11; ph = 0; % change K = 0.3; C = 0.2;
% % run_2Floe_D(102,0.05,0,0.006,tw)
% run_2Floe_D(102,0.05,0,0.00020104,tw) % using mauros drift force

% tw = 20; ph = 2; 
% run_2Floe_D(135,0.05,0,0.14,tw)
% run_2Floe_D(135,0.05,0,0.0046,tw) % using mauros drift force

% tw = 20; ph = 2; 
% % run_2Floe_D(141,0.05,0,0.2,tw)
% run_2Floe_D(141,0.05,0,0.0032,tw) % using mauros drift force


%% Plot drift velocity vs steepness, check with single-floe data
clear all
      
dat = [85,  100, 101, 102, 135, 141; ... % run number
       1.5, 1.5, 1.5, 1.4, 1.3, 1.4; ... % frequency
       40,  20,  20,  20,  80,  80;  ... % wave height
       0.03, 0.004, 0.003, 0.003, 0.021, 0.04; ... % scattering force
       0.085, 0.009, 0.01, 0.006, 0.14, 0.2;   ... % drift force
       0.0197, 0.0049, 0.0049, 0.00020104, 0.0046, 0.0032; ... % mauro's drift force
       0.106, 0.0342, 0.1112, 0.0317, 0.1235, 0.1499;   ... % empirical drift velocity
       0.0673, 0.0321, 0.0623, 0.0105, 0.0353, 0.0404]; % mauro's drift velocity      
       % nb. above values generated when drag coeff set to 0.05
      
      
for j = [1,2,4,5,6]%:length(dat(1,:))
 freq = dat(2,j);
 [field] = wavefield('f',freq,0.83);
 lambda = field{4,2};
 k = field{5,2};
 celerity = freq*lambda;
 Vdrift_emp_nd(j) = dat(7,j)/celerity;
 Vdrift_mau_nd(j) = dat(8,j)/celerity;
 steepness(j) = k*dat(3,j)*1e-3/2;
 hold on,plot(steepness(j),Vdrift_emp_nd(j),'k.','MarkerSize',20);
end
%
% hold on
% h1=plot(steepness,Vdrift_emp_nd,'k^','MarkerSize',14);
% h2=plot(steepness,Vdrift_mau_nd,'m*','MarkerSize',14);
% legend([h1 h2],'Empirical','Mauros')
% legend([h1],'Empirical')

%% Check drift force using scattered coeffs from PF model

F = 1.4;%linspace(0.5,1.8,50);
for j = 1:length(F)
rhoF = 1000;
g = 9.81;
% f = 1.4;
f = F(j)
h = 0.83;
[field] = wavefield('f',f,h);
k = field{5,2};
sigma = (2*pi*f)^2/g;
addpath('../../SingleFloeModels/PotentialFlow')
% PF model:
if ~exist('d','var'); d=0.01; end
if ~exist('D','var'); D=0.015; end
if ~exist('L','var'); L=0.2; end
if ~exist('N','var'); N=50; end
if ~exist('M','var'); M=1; end
if ~exist('rho','var'); rho=650; end
if ~exist('A_p0','var'); A_p0=1; end
if ~exist('B_m0','var'); B_m0=0; end
if ~exist('modes','var'); modes=111; end
if ~exist('TS','var'); TS=1; end
[s_s,s_h,s_p] = run_PF_2D(f,h,d,D,L,N,M,rho,A_p0,B_m0,modes,TS);
[a_m_d,b_p_d,~,~,~,~,~,~] = fn_Diffraction(sigma,h,d,L,N,A_p0,B_m0);                                                     
[a_m_s,b_p_s,~,~,~,~,~,~] = fn_Surge(sigma,h,d,L,N,s_s);
[a_m_h,b_p_h,~,~,~,~,~,~] = fn_Heave(sigma,h,d,L,N,s_h);
[a_m_p,b_p_p,~,~,~,~,~,~] = fn_Pitch(sigma,h,d,D,L,N,s_p);
TC(j) = abs(b_p_d(1) + b_p_s(1) + b_p_h(1) + b_p_p(1));
RC(j) = abs(a_m_d(1) + a_m_s(1) + a_m_h(1) + a_m_p(1));
end
% Mauros formula:
AI = 40*1e-3; % wave amplitude
Fd = rhoF*g/2*(RC*AI)^2*(1+2*k*h/sinh(2*k*h))