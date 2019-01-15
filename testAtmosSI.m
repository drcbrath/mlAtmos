% test AtmosSI

H = [  -100; 0.0; 5500; 11000.0; 15500; 20000.0; 26000; 32000.0; 39500; 47000.0; 49000; 51000.0; 61000; 71000.0; 80000; 84852.0; 90000 ];

% test std atm
[T,rho,P,a] = AtmosSI(H);
[oT,orho,oP,oa] = AtmSI(H);    % compare to old version

fprintf('[H: T,oT,T-oT, rho,orho,rho-orho, P,oP,P-oP, a,oa,a-oa]\n');
fprintf('%10.0f %13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e\n',[H,T,oT,T-oT, rho,orho,rho-orho, P,oP,P-oP, a,oa,a-oa]')


% test dT atm
dT = 20;   % hot day
[T,rho,P,a] = AtmosSI(H,dT);
[oT,orho,oP,oa] = AtmSI(H,dT);   % compare to old version

fprintf('[H: T,oT,T-oT, rho,orho,rho-orho, P,oP,P-oP, a,oa,a-oa]\n');
fprintf('%10.0f %13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e\n',[H,T,oT,T-oT, rho,orho,rho-orho, P,oP,P-oP, a,oa,a-oa]')

% test custom profile
Hk     = [   0.0; 11000.0; 20000.0; 32000.0; 47000.0; 51000.0; 71000.0; 84852.0 ];
Tk     = [308.15;  236.65;  236.65;  248.65;  290.65;  290.65;  234.65;  206.95];

[T,rho,P,a] = AtmosSI(H,Hk,Tk);
[oT,orho,oP,oa] = AtmSI(H,dT);   % compare to old version

fprintf('[H: T,oT,T-oT, rho,orho,rho-orho, P,oP,P-oP, a,oa,a-oa]\n');
fprintf('%10.0f %13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e -- %13.6e %13.6e %13.6e\n',[H,T,oT,T-oT, rho,orho,rho-orho, P,oP,P-oP, a,oa,a-oa]')
