% test AtmosSI

H = [0; 100; 1000; 10000];

% test std atm
[T,rho,P,a] = AtmosSI(H);
[oT,orho,oP,oa] = AtmSI(H);    % compare to old version

fprintf('[T,oT,T-oT, rho,orho,rho-orho, P,oP,P-oP, a,oa,a-oa]\n');
fprintf('%13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e\n',[T,oT,T-oT, rho,orho,rho-orho, P,oP,P-oP, a,oa,a-oa]')


% test dT atm
dT = 20;   % hot day
[T,rho,P,a] = AtmosSI(H,dT);
[oT,orho,oP,oa] = AtmSI(H,dT);   % compare to old version

fprintf('[T,oT,T-oT, rho,orho,rho-orho, P,oP,P-oP, a,oa,a-oa]\n');
fprintf('%13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e\n',[T,oT,T-oT, rho,orho,rho-orho, P,oP,P-oP, a,oa,a-oa]')

% test custom profile
Hk     = [   0.0; 11000.0; 20000.0; 32000.0; 47000.0; 51000.0; 71000.0; 84852.0 ];
Tk     = [308.15;  236.65;  236.65;  248.65;  290.65;  290.65;  234.65;  206.95];

[T,rho,P,a] = AtmosSI(H,Hk,Tk);
[oT,orho,oP,oa] = AtmSI(H,dT);   % compare to old version

fprintf('[T,oT,T-oT, rho,orho,rho-orho, P,oP,P-oP, a,oa,a-oa]\n');
fprintf('%13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e\n',[T,oT,T-oT, rho,orho,rho-orho, P,oP,P-oP, a,oa,a-oa]')
