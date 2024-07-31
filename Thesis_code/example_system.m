
% Demonstration of the vector fitting method using simple inverter model
% As is described in section 3.4

close all
clear all
Ns=501;    %Number of frequency samples
Nc=1;      %Size of Y (after reduction)
bigY=zeros(Nc,Nc,Ns);

s=2*pi*1i*logspace(1,5,Ns);

% switching frequencies
fsw = 10e3;
Ts = 1/fsw;

%Inverter side filter ind.
Lf = 0.87e-2;  %Lf = 0.87e-3 to make system unstable;
%Filter capacitor
Cf = 22e-6; %Cf = 22e-6;
%Grid side filter ind.
Lg = 0.22e-3; %Lg = 0.22e-3;
%Grid R
Rg = 0.1;
%Proportional gain of controller
Kp = 5.6;
%Integrator gain of the controller
Kr = 1000;

%Produce impedence frequency response of inverter
for k=1:Ns
  sk=s(k);
  exp_term = (1-(1.5*Ts*sk/2))/(1+(1.5*Ts*sk/2));
  Yc = ((sk.^2*Cf*Lf + 1)./(sk.^3*Cf*Lf*Lg + sk*(Lf+Lg) + Kp*exp_term));
  Yg = 1./(Lg*sk + Rg);
  bigY(:,:,k)=Yc;
end  

%Build inverter analytical transfer function
s_f = tf('s');
%Pade approx. of exponent
exp_term = (1-(1.5*Ts*s_f/2))/(1+(1.5*Ts*s_f/2));
%Transfer function
Yc_sys = (((s_f^2)*Cf*Lf + 1)/((s_f^3)*Cf*Lf*Lg + s_f*(Lf+Lg) + Kp*exp_term));
% Yg_sys = 1/(s_f*Lg + Rg);
sys_tf = Yc_sys;
%================================================
%=           POLE-RESIDUE FITTING               =
%================================================ 
opts.N = 4;              %Order of approximation. (Is used when opts.poles=[]).  
opts.stable = 0; 
opts.poletype='logcmplx';   %Will use logarithmically spaced, complex poles. (Is used when opts.poles=[]).
poles=[]; %[] -->N initial poles are automatically generated as defined by opts.startpoleflag 
opt.weightparam = 3;

opts.screen = 0;
opts.errplot = 0;
opts.cmplx_ss = 1;

[SER,rmserr,bigYfit,opts2]=VFdriver(bigY,s,poles,opts);

A = SER.A;B = SER.B;C = SER.C;D = SER.D;E = SER.E;
sys = ss(A,B,C,D);

% Evaluate HSVs
figure;
s = hsvd(sys);
bar(s/s(1));
set(gca, 'YScale', 'log');
ylabel('Contribution [p.u.]','FontSize', 14)
xlabel('Order (number of states)','FontSize', 14)

% Model reduction using balanced truncation
[sys_b,g] = balreal(sys);
rho = s(4);

elim = find(g<rho);
sys = modred(sys_b,elim,'del');
length(sys.A)
ylim([0,10]);
grid on

figure;
plot(eig(sys_tf),'square','MarkerSize',11)
hold on
plot(eig(sys),'X','MarkerSize',11)
grid on
legend('System TF','Vector fit model','FontSize', 12)
XLabel = xlabel('Imaginary [sec^{-1}]','FontSize', 14);
YLabel = ylabel('Real [sec^{-1}]','FontSize', 14);
% legend('System Model','VF')
% title('System Eigenvalues')
xlim([-14e3,1000]);







