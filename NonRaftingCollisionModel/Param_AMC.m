function Param = Param_AMC

% THICKNESS [m]
Param.D = 0.015;

% GRAVITATIONAL ACCELERATION [m/s^2]
Param.g = 9.81;            

% FLUID DENSITY [kg/m^3]
Param.rho = 1000;           

% FLOE DENSITY [kg/m^3]
Param.rho_b = 650;

% DRAFT [m]
Param.d = Param.D.*Param.rho_b./Param.rho; 

% FLUID DEPTH [m]
Param.h = 0.831;

% FLOE RADIUS [m]
Param.L = 0.2;

% MASS OF FLOE [kg/m]
Param.M = Param.rho_b*2*Param.L*Param.D; 

% NUMBER OF VERTICAL EIGENFUNCTIONS
Param.N = 20;


return