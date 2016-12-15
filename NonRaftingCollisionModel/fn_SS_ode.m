% function sol = fn_SS_ode(t,X,Xm,WaveParam,FloeParam,Coeff,Mooring,Trans)
%
% LJ YIEW
% Created on  Jul 2013
% Last edited Oct 2016
%
% Sets up the ODEs for the Rumer/Marchenko model (Grotmaack & Meylan, 2006)
%
% INPUTS:
% t         = time
% X         = displacement
% Xm        = equilibrium position of floe
% WaveParam = wave parameters (see code)
% FloeParam = floe parameters (see code)
% Coeff     = drag and added mass coefficients
% Mooring   = spring and damping coefficients
% Trans     = transient wave parameters (see code)
%
% OUTPUTS:
% sol = [displacement,velocity]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sol = fn_SS_ode(t,X,Xm,WaveParam,FloeParam,Coeff,Mooring,Trans)
%% INPUTS:
H     = WaveParam.H;     % wave height (2*wave amplitude)
omega = WaveParam.omega; % angular wave frequency
k     = WaveParam.k;     % wave number
rho   = WaveParam.rho;   % fluid density
h     = WaveParam.h;     % water depth

m = FloeParam.m; % floe mass
A = FloeParam.A; % floe wetted surface area
g = 9.81;

Cd = Coeff.Cd; % drag coefficient
Cm = Coeff.Cm; % added mass coefficient

K = Mooring.K; % spring coefficient
C = Mooring.C; % damping coefficient

f    = Trans.f;  % transient wave forcing coefficients
trans = Trans.t; % flag for transient wave forcing

%% ODE SOLUTION
sol = zeros(size(X));
sol(1) = X(2); % X(2) = FLOE VELOCITY

% X DERIVATIVE OF WAVE PROFILE
% X(1) = DISPLACEMENT
dndx = H.*k./2.*cos(k.*X(1)-omega.*t + f.c);
if trans == 1 % INCLUDE TRANSIENTS
 dndx = 1e-3*f.a.*t.^f.b.*k*cos(k.*X(1) - omega.*t + f.c);
end

% VELOCITY OF FLUID AT FREE SURFACE
Vw = omega.*H./2.*sin(k.*X(1)-omega.*t +f.c).*coth(k.*h);
if trans == 1
 Vw = omega.*1e-3*f.a.*t.^f.b.*sin(k.*X(1)-omega.*t + f.c).*coth(k.*h);
end

% ACCELERATION (dVdt)
sol(2) = 1./(m.*(1+Cm))...
          .*(-m.*g.*dndx ...
             + rho.*Cd.*A.*(tanh(100*(Vw-X(2))).*(Vw-X(2))).*(Vw-X(2)) ...
             - K*(X(1)-Xm) - C*X(2) );          % MOORING FORCES
             


