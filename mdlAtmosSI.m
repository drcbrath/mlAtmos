classdef mdlAtmosSI
%
% Copyright (c) 2006-2019 Curtis B. Rath, PhD., aka drcbrath
% available under the MIT license from Github, https://github.com/drcbrath/mdlAtmos

   properties
      T     = nan;   % (K) temperature
      Rho   = nan;   % (kg/m^3) density
      P     = nan;   % (N/m^2) pressure
      a     = nan;   % (m/s) sonic speed
      visc  = nan;   % (N*s/m^2) absolute viscosity
      theta = nan;   % (1) temperature ratio to sea level standard
      sigma = nan;   % (1) density ratio to sea level standard
      delta = nan;   % (1) pressure ratio to sea level standard
      kappa = nan;   % (1) sonic speed ratio to sea level standard
   endproperties
   
   properties (private)
      %------- constants -------
      Re   = 6369000;       % (m) radius of the earth
      GMR  = 0.034163195;   % (degK/m) combined gravity and gas constant of dry air on earth
      H0   = 0.0;           % (m) datum, sea level
      T0   = 288.15;        % (K), SL std temp
      rho0 = 1.225;         % (kg/m^3), SL std density
      P0   = 101325;        % (N/m^2), SL std pressure
      a0   = 340.2686;      % (m/s), speed of sound at SL std temperature

      % std day atmosphere profiles
      stdHk     = [   0.0; 11000.0; 20000.0; 32000.0; 47000.0; 51000.0; 71000.0; 84852.0 ];
      stdTk     = [288.15;  216.65;  216.65;  228.65;  270.65;  270.65;  214.65;  186.95];
      stdTgradk = [-0.0065000; 0.0000000; 0.0010000; 0.0028000; 0.0000000; -0.0028000; -0.0019997; -0.0019997];

      % temperature profile      
      Hk
      Tk
      Tgradk

   endproperties

   methods

      % constructors
      function atmos = mdlAtmosSI()
         % set profiles
         atmos.Hk = std Hk;
         atmos.Tk = stdTk;
         atmos.Tgradk = stdTgradk;
         % initialize atmosphere properties
         atmos.at(0);
      endfunction

      function atmos = mdlAtmosSI(dT)
         if isscalar(dT)
            % set profiles
            atmos.Hk = std Hk;
            atmos.Tk = stdTk+dT;
            atmos.Tgradk = stdTgradk;
            % initialize atmosphere properties
            atmos.at(0);
         else
            dispUsage();
         end
      endfunction
         
      function atmos = mdlAtmosSI(Hk,Tk)
         if isvector(Hk) && isvector(Tk)
            % set profiles
            atmos.Hk = Hk;
            atmos.Tk = Tk;
            atmos.Tgradk = diff(Tk)./diff(Hk);    % compute temperature gradients through atmosphere
            atmos.Tgradk = [atmos.Tgradk atmos.Tgradk(end)];
            % initialize atmosphere properties
            atmos.at(0);
         else
            dispUsage();
         end
      endfunction
   
      % evaluate properties
      function [atmos,T,rho,p,a,visc,theta,sigma,delta,kappa] = at(atmos,Hgp)
         %function [atmos,T,rho,p,a,visc,theta,sigma,delta,kappa] = at(Hgp)
         %
         % compute atmosphere properties at given geopotential altitude
         %
         % atmos.at(Hgp)
         %	Input
         %		Hgp   === (m) geopotential altitude for desired properties
         %
         % Output
         %   atmos === this atmosSI object, with properties set at Hgp
         %	 T     ===
         %   rho   === 
         %   P     === 
         %   a     === 
         %   visc  === 
         %   theta === 
         %   sigma === 
         %   delta === 
         %   kappa === 

         %------- compute ratios using established temperature profile -------

         % initilize
         theta = nan(size(Hgp));
         sigma = nan(size(Hgp));
         delta = nan(size(Hgp));
         kappa = nan(size(Hgp));

         % loop through input altitudes
         for m = 1:length(Hgp)
             Pks = 1;
             for k = 1:size(atmos.Hk,1)-1

                 if Hgp(m) > atmos.Hk(k+1)      % h is above layer k

                     if atmos.Tgradk(k)~=0    % linear thermal layer
                         sr = ( atmos.Tk(k+1)./atmos.Tk(k) ).^(-GMR./atmos.Tgradk(k));
                     else               % isothermal layer
                         sr = exp(-atmos.GMR*(atmos.Hk(k+1)-atmos.Hk(k))./atmos.Tk(k+1));
                     end
                     Pks = Pks .* sr;

                 else   % have reached layer such that Hk(k) < h < Hk(k+1)

                     dh = Hgp(m) - atmos.Hk(k);
                     T = atmos.Tk(k) + atmos.Tgradk(k) .* dh;
                     theta(m) = T / atmos.T0;              % temperature ratio

                     if atmos.Tgradk(k)~=0    % linear thermal layer
                         sr = ( T./atmos.Tk(k) ).^(-GMR./atmos.Tgradk(k));
                     else               % isothermal layer
                         sr = exp(-atmos.GMR*(Hgp(m)-atmos.Hk(k))./atmos.Tk(k+1));
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
         % Note: 2 coeff form of Sutherland's law results differ slightly from 3 coeff form
         % But, it is consistent with & taken from the US Standard Atmosphere, 1976, pg 19

         visc = 1.458e-6 * T.^1.5 ./ (T+110.4);

         atmos.T     = T;
         atmos.rho   = rho;
         atmos.P     = P;
         atmos.a     = a;
         atmos.visc  = visc;
         atmos.theta = theta;
         atmos.sigma = sigma;
         atmos.delta = delta;
         atmos.kappa = kappa;

         end

      endfunction

      function T = T(atmos)
         T = atmos.T;
      end

      function rho = rho(atmos)
         rho = atmos.rho;
      end

      function P = P(atmos)
         P = atmos.P;
      end

      function a = a(atmos)
         a = atmos.a;
      end

      function visc = visc(atmos)
         visc = atmos.visc;
      end

      function theta = theta(atmos)
         theta = atmos.theta;
      end

      function sigma = sigma(atmos)
         sigma = atmos.sigma;
      end

      function delta = delta(atmos)
         delta = atmos.delta;
      end

      function kappa = kappa(atmos)
         kappa = atmos.kappa;
      end

   endmethods

endclassdef
