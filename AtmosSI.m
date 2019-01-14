function [T,rho,P,a,visc,theta,sigma,delta,kappa] = AtmosSI(Hgp,varargin)
%function [T,rho,P,a,visc,theta,sigma,delta,kappa] = AtmosSI(Hgp,varargin)
%
% compute atmosphere properties at given geopotential altitude
%
% AtmosSI(Hgp)
%	Input
%		Hgp   === (m) geopotential altitude for desired properties
%
% AtmosSI(Hgp,dT)
%	Input
%		Hgp   === (m) geopotential altitude for desired properties
%		dT    === (degK) scalar, delta temperature to add to temperature profile
% AtmosSI(Hgp,bpHgp,bpT)
%	Input
%		Hgp   === (m) geopotential altitude for desired properties
%		bpHgp === (m) vector, geopotential altitude breakpoints for custom profile
%		bpT   === (degK) vector, temperature breakpoints for custom profile
%
% Output
%	T     ===
%   rho   === 
%   P     === 
%   a     === 
%   visc  === 
%   theta === 
%   sigma === 
%   delta === 
%   kappa === 

%------- constants -------
Re   = 6369000;       % (m) radius of the earth
GMR  = 0.034163195;   % (degK/m) combined gravity and gas constant of dry air on earth
H0   = 0.0;           % (m) datum, sea level
T0   = 288.15;        % (K), SL std temp
rho0 = 1.225;         % (kg/m^3), SL std density
P0   = 101325;        % (N/m^2), SL std pressure
a0   = 340.2686;      % (m/s), speed of sound at SL std temperature

% std day atmosphere profiles
Hk     = [   0.0; 11000.0; 20000.0; 32000.0; 47000.0; 51000.0; 71000.0; 84852.0 ];
Tk     = [288.15;  216.65;  216.65;  228.65;  270.65;  270.65;  214.65;  186.95];

% Tgradk = diff(Tk)./diff(Hk);    % compute temperature gradients through atmosphere
% Tgradk = [Tgradk; Tgradk(end)]
Tgradk = [-0.0065000; 0.0000000; 0.0010000; 0.0028000; 0.0000000; -0.0028000; -0.0019997; -0.0019997];
  
%------- process alternate input arguments -------
if nargin == 1                                        % std atmosphere
elseif nargin == 2 && isscalar(varargin{1})            % add temperature delta, use std lapse rates
	Tk = Tk + varargin{1}; 
elseif nargin == 3 && isvector(varargin{1}) && isvector(varargin{2})   % custom temperature profile
	Hk = varargin{1};
	Tk = varargin{2};             
	Tgradk = diff(Tk)./diff(Hk);    % compute temperature gradients through atmosphere
	Tgradk = [Tgradk; Tgradk(end)];
	Tsl = Tk(1);                    % note sea level
else
	% error, display usage
end

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

        if Hgp(m) > Hk(k+1)      % h is above layer k

            if Tkgrad(k)~=0    % linear thermal layer
                sr = ( Tk(k+1)./Tk(k) ).^(-GMR./Tkgrad(k));
            else               % isothermal layer
                sr = exp(-GMR*(Hk(k+1)-Hk(k))./Tk(k+1));
            end
            Pks = Pks .* sr;

        else   % have reached layer such that Hk(k) < h < Hk(k+1)

            dh = Hgp(m) - Hk(k);
            T = Tk(k) + Tgradk(k) .* dh;
            theta(m) = T / T0;              % temperature ratio

            if Tgradk(k)~=0    % linear thermal layer
                sr = ( T./Tk(k) ).^(-GMR./Tgradk(k));
            else               % isothermal layer
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

% dimensionalize
T   = T0 * theta;
rho = rho0 * sigma;
P   = P0 * delta;
a   = a0 * kappa;

% viscosity
visc = AirViscSI(T);

end
