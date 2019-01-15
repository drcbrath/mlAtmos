classdef mdlAtmosSI

   properties
      T
      Rho
      P
      a
      Visc
      theta
      sigma
      delta
      kappa
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
         % set profiles
         atmos.Hk = std Hk;
         atmos.Tk = stdTk+dT;
         atmos.Tgradk = stdTgradk;
         % initialize atmosphere properties
         atmos.at(0);
      endfunction
         
      function atmos = mdlAtmosSI(Hk,Tk)
         % set profiles
         atmos.Hk = Hk;
         atmos.Tk = Tk;
         atmos.Tgradk = diff(Tk)./diff(Hk);    % compute temperature gradients through atmosphere
         atmos.Tgradk = [atmos.Tgradk atmos.Tgradk(end)];
         % initialize atmosphere properties
         atmos.at(0);
      endfunction
   
      % evaluate properties
      function [T,rho,p,visc,theta,sigma,delta,kappa,visc] = at(Hgp)
         
      endfunction
   
      function [T,rho,p,visc,theta,sigma,delta,kappa,visc] = operator(Hgp)
         [T,rho,p,visc,theta,sigma,delta,kappa,visc] = at(Hgp);
      endfunction

   endmethods

endclassdef
