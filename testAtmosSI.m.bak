% test AtmosSI
%
% Copyright (c) 2006-2019 Curtis B. Rath, PhD., aka drcbrath
% available under the MIT license from Github, https://github.com/drcbrath/mdlAtmos

% 1. Validate AtmosSI, on std day by comparison to 
%    independent authoritative reference, US Std Atm
% 2. Verify & validate dT modified off std day by comparison to hand calculations 
%    in the two types of layers, linear thermal and isothermal. All other logic is 
%    unchanged, so dT=/=0 does not affect other aspects; hence demonstration at just 
%    two points suffices to validate dT=/=0 computation at any altitude.
% 3. Test custom atmosphere temperature profile. Construct a custom temperature profile 
%    that is equivalent to the standard atmosphere profile in layer types and
%    temperature lapse rates, differing in the breakpoints; hence results at
%    common altitudes must agree. Then validation by comparison to standard day.
% 4. Validate AtmosUS by comparison to verified AtmosSI results multiplied by 
%    conversion factors.
	% 4.1 std day
	% 4.2 dT offset day
	% 4.3 custom profile
%
% compare program oututs to reference results in TestAtmospherePropertiesModel_SI.dat. 
% Ones favorite diff tool might be handy.

% test altitudes; below, at each, between each, and above breakpoints

H = [  -100; 0.0; 5500; 11000.0; 15500; 20000.0; 26000; 32000.0; 39500; 47000.0; 49000; 51000.0; 61000; 71000.0; 80000; 84852.0; 90000 ];
dT = 30;


% (1) standard day, for comparison to US Standard Atmosphere, 1976, reference work
[T_std, rho_std, P_std, a_std, visc_std, theta_std, sigma_std, delta_std, kappa_std] = AtmosSI(H);    % compare to old version


% (2) linear layer and isothermal layers by "hand calcs" with and without dT offset,
%     comparison at reference work table values and for verification of formulation
%     Note: test point for linear thermal layer is selected at top of troposphere
%     so the result will serve as bottom of stratosphere in isothermal layer test calc.

% linear thermal layer
H_L = 11000;       % test point, top of troposphere
[T_2L_dT20, rho_2L_dT20, P_2L_dT20] = AtmosSI(H_L,dT);   % with dT = 20, should match hand calc in linear layer

% isothermal layer
H_I = 15500;   % test point
[T_2I_dT20, rho_2I_dT20, P_2I_dT20] = AtmosSI(H_I,dT);   % with dT = 20, should match hand calc in linear layer


% (3) dT offset, usinsg same dT offset as verified hand calcs in (2), so test points
%  11000 & 15500 should match proving property computation, while full results table
% verifies the rest of the input-copmute-output code
[T_dT, rho_dT, P_dT, a_dT, visc_dT, theta_dT, sigma_dT, delta_dT, kappa_dT] = AtmosSI(H,dT);    % compare to old version


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

[T_cus, rho_cus, P_cus, a_cus, visc_cus, theta_cus, sigma_cus, delta_cus, kappa_cus] = AtmosSI(H,H_custom,T_custom);    % compare to old version


%% Results
fid = fopen('testAtmosSI.dat','wt');

fprintf(fid,'Atmosphere Properties Model Test Data --- Validated against "US Standard Atmosphere, 1976"\n\n');
fprintf(fid,'\n(1) Comparison of: std day, evaluated at test altitudes\n');
fprintf(fid,'     H   T          rho        P          a          theta      sigma      delta      kappa      visc\n');
fprintf(fid,'    (m) (K)        (kg/m^3)   (N/m^2)    (m/s)      (1)        (1)        (1)        (1)        (N*s/m^2)\n');
fprintf(fid,'%7.1f %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n',[H, T_std, rho_std, P_std, a_std, theta_std, sigma_std, delta_std, kappa_std, visc_std]');

fprintf(fid,'\n(2) Comparison to hand calculations, one in linear layer, one in isothermal layer\n');
fprintf(fid,'H, T, P only, since other results follow from these as is verified in (1) above\n');
fprintf(fid,'     H   T          P\n');
fprintf(fid,'    (m) (K)   (N/m^2)\n');
fprintf(fid,'%7.1f %10.4e %10.4e --- linear layer with dT=%4.0f\n',H_L,T_2L_dT20, P_2L_dT20, dT);
fprintf(fid,'%7.1f %10.4e %10.4e --- isothermal layer with dT=%4.0f\n',H_I,T_2I_dT20, P_2I_dT20, dT);

fprintf(fid,'\n(3) Comparison of: dT offset, std+(%f) day, evaluated at test altitudes\n',dT);
fprintf(fid,'Compare at H=11000 & 15500 m, T & P should be as above in (2) "hand calsc". Rest of table results from logic tested in (1) above\n');
fprintf(fid,'     H   T          rho        P          a          theta      sigma      delta      kappa      visc\n');
fprintf(fid,'    (m) (K)        (kg/m^3)   (N/m^2)    (m/s)      (1)        (1)        (1)        (1)        (N*s/m^2)\n');
fprintf(fid,'%7.1f %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n',[H, T_dT, rho_dT, P_dT, a_dT, theta_dT, sigma_dT, delta_dT, kappa_dT, visc_dT]');

fprintf(fid,'\n(4) Comparison of: custom profile, equivalent to std day, evaluated at test altitudes\n');
fprintf(fid,'     H   T          rho        P          a          theta      sigma      delta      kappa      visc\n');
fprintf(fid,'    (m) (K)        (kg/m^3)   (N/m^2)    (m/s)      (1)        (1)        (1)        (1)        (N*s/m^2)\n');
fprintf(fid,'%7.1f %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n',[H, T_cus, rho_cus, P_cus, a_cus, theta_cus, sigma_cus, delta_cus, kappa_cus, visc_cus]');

fclose(fid);
