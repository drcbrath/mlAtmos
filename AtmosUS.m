function [T,rho,P,a,visc,theta,sigma,delta,kappa] = AtmosUS(Hgp,varargin)
%function [T,rho,P,a,visc,theta,sigma,delta,kappa] = AtmosUS(Hgp,varargin)
%
% compute atmosphere properties at given geopotential altitude
%
% AtmosUS(Hgp)
%	Input
%		Hgp   === (ft) geopotential altitude for desired properties
%
% AtmosUS(Hgp,dT)
%	Input
%		Hgp   === (ft) geopotential altitude for desired properties
%		dT    === (degR) scalar, delta temperature to add to temperature profile
% AtmosUS(Hgp,bpHgp,bpT)
%	Input
%		Hgp   === (ft) geopotential altitude for desired properties
%		bpHgp === (ft) vector, geopotential altitude breakpoints for custom profile
%		bpT   === (R) vector, temperature breakpoints for custom profile
%
% Output
%	 T     ===
%   rho   === 
%   P     === 
%   a     === 
%   visc  === 
%   theta === 
%   sigma === 
%   delta === 
%   kappa === 
%
% Copyright (c) 2006-2019 Curtis B. Rath, PhD., aka drcbrath
% available under the MIT license from Github, https://github.com/drcbrath/mdlAtmos

% Note: US Standard Atmosphere model is defined in SI units
% So, here the inputs from US (FSS) to SI, calculations are done is SI, and 
% results converted from SI to US for output
%------- constants -------
Re   = 6369000;       % (m) radius of the earth
GMR  = 0.034163195;   % (degK/m) combined gravity and gas constant of dry air on earth
H0   = 0.0;           % (m) datum, sea level
T0   = 288.15;        % (K), SL std temp
rho0 = 1.225;         % (kg/m^3), SL std density
P0   = 101325;        % (N/m^2), SL std pressure
a0   = 340.2686;      % (m/s), speed of sound at SL std temperature

% std day atmosphere profiles
Hk     = [   0.0; 11000.0; 20000.0; 32000.0; 47000.0; 51000.0; 71000.0; 84852.0 ];                       % (m)
Tk     = [288.15;  216.65;  216.65;  228.65;  270.65;  270.65;  214.65;  186.95];                        % (K)
Tgradk = [-0.0065000; 0.0000000; 0.0010000; 0.0028000; 0.0000000; -0.0028000; -0.0019997; -0.0019997];   % (degK/ft)
  
%------- process alternate input arguments -------
if nargin == 1                                         % std atmosphere
elseif nargin == 2 && isscalar(varargin{1})            % add temperature delta, use std lapse rates
	Tk = Tk + varargin{1}/1.8; 
elseif nargin == 3 && isvector(varargin{1}) && isvector(varargin{2})   % custom temperature profile
	Hk = varargin{1}*0.3048;
	Tk = varargin{2}/1.8;             
	Tgradk = diff(Tk)./diff(Hk);    % compute temperature gradients through atmosphere
	Tgradk = [Tgradk; Tgradk(end)];
	Tsl = Tk(1);                    % note sea level
else
	% error, display usage
end

Hgp = Hgp*0.3048;   % convert from ft to m

%------- compute ratios using established temperature profile -------

% initilize
theta = nan(size(Hgp));
sigma = nan(size(Hgp));
delta = nan(size(Hgp));
kappa = nan(size(Hgp));

% loop through input altitudes
for m = 1:length(Hgp)
    Pks = 1;
    for k = 1:size(Hk,1)-1

        if Hgp(m) > Hk(k+1)                    % h is above layer k

            if Tgradk(k)~=0                    % linear thermal layer
                sr = ( Tk(k+1)./Tk(k) ).^(-GMR./Tgradk(k));
            else                               % isothermal layer
                sr = exp(-GMR*(Hk(k+1)-Hk(k))./Tk(k+1));
            end
            Pks = Pks .* sr;

        else   % have reached layer such that Hk(k) < h < Hk(k+1)

            dh = Hgp(m) - Hk(k);
            T = Tk(k) + Tgradk(k) .* dh;
            theta(m) = T / T0;                 % temperature ratio

            if Tgradk(k)~=0                    % linear thermal layer
                sr = ( T./Tk(k) ).^(-GMR./Tgradk(k));
            else                               % isothermal layer
                sr = exp(-GMR*(Hgp(m)-Hk(k))./Tk(k+1));
            end
            delta(m) = Pks .* sr;              % pressure ratio
            sigma(m) = delta(m) ./ theta(m);   % density ratio
            kappa(m) = sqrt(theta(m));         % sonic speed ratio
            break;

        end
    end
end

% reshape to that of input Hgp
sz    = size(Hgp);
theta = reshape(theta,sz);
sigma = reshape(sigma,sz);
delta = reshape(delta,sz);
kappa = reshape(kappa,sz);

% dimensionalize with US unit sea level std values
T0     = 518.67;        % deg R, SL std temp
rho0   = 0.002377;      % sl/ft^3, SL std density
P0     = 2116.2;        % lbf/ft^2, SL std pressure
a0     = 1116.3;        % ft/s, speed of sound at SL std temperature
T   = T0   * theta;
rho = rho0 * sigma;
P   = P0   * delta;
a   = a0   * kappa;

% viscosity
visc = AirViscUS(T);

end
