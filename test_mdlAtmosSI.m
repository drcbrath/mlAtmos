% test mdlAtmosSI

% test altitudes; below, at each, between each, and above breakpoints
H = [  -100; 0.0; 5500; 11000.0; 15500; 20000.0; 26000; 32000.0; 39500; 47000.0; 49000; 51000.0; 61000; 71000.0; 80000; 84852.0; 90000 ];

% test dT offset
dT = 30;

% test custom temperature profile
% (4) custom profile, use profile generated using verified dT offset computations
% construct a custom profile that is equivalent to the std atm, i.e. has linear and isothermal
% layers corresponding to std profile layers with same lapse rates and boundaries 
% at each layer change of the std profile; being equivalent will generate the same
% results for comparison to the std atmosphere while exercising the custom
% temperature profile logic.
% std atmosphere temperature profile breakpoints
Hk     = [   0.0; 11000.0; 20000.0; 32000.0; 47000.0; 51000.0; 71000.0; 84852.0 ];
Tk     = [288.15;  216.65;  216.65;  228.65;  270.65;  270.65;  214.65;  186.95];
% constructed with std atm breakpoints + midpoints for std atm profile ---hence,
% same lapse rates and layer types; more layers, but plotted will show same profile,
% so same results will follow
H_custom = sort([Hk; Hk(1:end-1) + diff(Hk)/2]);
T_custom = interp1(Hk,Tk,H_custom);


% constructors
%
% (1) default
a1 = mdlAtmosSI()

% (2) dT offset
a2 = mdlAtmosSI(dT)

% (3) custom profile

a3 = mdlAtmosSI(H_custom,T_custom)
