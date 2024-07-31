% Author: Finnur Helgi Malmquist
% Date: 07/30/2024
% This script was created to process and produce a vector fit model of the
% black-box STATCOM and grid systems, using their immitance measurements as
% inputs.
% Note that in order to run this code, the Vector Fitting Toolbox is
% required



% filename = "Zqd_idx1001_14.0MW-222m-50Hz_SDL_6_MVAR_SCR_2_5_PERT_0_25.csv";
filename = "Zqd_idx901_14.0MW-222m-50Hz_SDL_STATCOM_2Hz.csv";
% filename = "Zqd_idx404_14.0MW-222m-50Hz_SDL_7_MVAR_SCR_1_5_WT_PERT_0_25.csv";


%Data processing

data = readtable(filename);

freq = data.Freq_inj_Hz_; s = freq*2*pi*1i;
N = length(freq);
interpolation = N;
f_start = freq(1); f_end = freq(end);

freq_new = linspace(f_start,f_end,interpolation);
s_new = freq_new*2*pi*1i; 
N_new = length(freq_new);


method = 'linear';  %Interpolation method

% Statcom data processing

statcom_response_qq = (cellfun(@str2num,data.Zqq_c));
statcom_response_dq = (cellfun(@str2num,data.Zdq_c));
statcom_response_qd = (cellfun(@str2num,data.Zqd_c));
statcom_response_dd = (cellfun(@str2num,data.Zdd_c));

statcom_resp = ([statcom_response_qq,statcom_response_dq,statcom_response_qd,statcom_response_dd]);

% Data interpolation
statcom_response_qq_int = interp1(freq,statcom_resp(:,1),freq_new,method)';
statcom_response_dq_int = interp1(freq,statcom_resp(:,2),freq_new,method)';
statcom_response_qd_int = interp1(freq,statcom_resp(:,3),freq_new,method)';
statcom_response_dd_int = interp1(freq,statcom_resp(:,4),freq_new,method)';

statcom_response = (reshape(statcom_resp', [2, 2, N]));

statcom_resp_int = ([statcom_response_qq_int,statcom_response_dq_int,statcom_response_qd_int,statcom_response_dd_int]);
statcom_response_int = (reshape(statcom_resp_int', [2, 2, N_new]));

statcom_bigH_orig = pageinv(statcom_response);
statcom_bigH = pageinv(statcom_response_int);

%Signal filtering
for j = 2:1:size((statcom_bigH),3)
        if (abs(freq_new(j)-1*50)<1e-9) || (abs(freq_new(j)-2*50)<1e-9) || (abs(freq_new(j)-3*50)<1e-9)
            statcom_bigH(:,:,j) = 1/2*(statcom_bigH(:,:,j-1)+statcom_bigH(:,:,j+1));
        end 
end

% Grid data processing

grid_response_qq = (cellfun(@str2num,data.Zqq_g));
grid_response_dq = (cellfun(@str2num,data.Zdq_g));
grid_response_qd = (cellfun(@str2num,data.Zqd_g));
grid_response_dd = (cellfun(@str2num,data.Zdd_g));

grid_response_qq_int = interp1(freq,grid_response_qq,freq_new,method)';
grid_response_dq_int = interp1(freq,grid_response_dq,freq_new,method)';
grid_response_qd_int = interp1(freq,grid_response_qd,freq_new,method)';
grid_response_dd_int = interp1(freq,grid_response_dd,freq_new,method)';

grid_resp = ([grid_response_qq,grid_response_dq,grid_response_qd,grid_response_dd]);
grid_response = (reshape(grid_resp', [2, 2, N]));

grid_resp_int = ([grid_response_qq_int,grid_response_dq_int,grid_response_qd_int,grid_response_dd_int]);
grid_response_int = (reshape(grid_resp_int', [2, 2, N_new]));

grid_bigH_orig = grid_response;
grid_bigH = grid_response_int;

for j = 2:1:size((grid_bigH),3)
        if (abs(freq_new(j)-1*50)<1e-9) || (abs(freq_new(j)-2*50)<1e-9) || (abs(freq_new(j)-3*50)<1e-9)
            grid_bigH(:,:,j) = 1/2*(grid_bigH(:,:,j-1)+grid_bigH(:,:,j+1));
        end 
end


figure(1)
n=1;
for i = 1:2
    for j = 1:2
        if n==1
%             title('qq')
%             loglog(freq,squeeze(abs(statcom_bigH_orig(i,j,:))),'-','LineWidth',1.5)
%             loglog(freq_new,squeeze(abs(statcom_bigH(i,j,:))),'-','LineWidth',1.5)
            loglog(freq,squeeze(abs(grid_bigH_orig(i,j,:))),'-','LineWidth',1.5)
%             Legend = legend('qq admittance');
        elseif n==2
%             title('qd')
%             loglog(freq,squeeze(abs(statcom_bigH_orig(i,j,:))),'-','LineWidth',1.5)
%             loglog(freq_new,squeeze(abs(statcom_bigH(i,j,:))),'-','LineWidth',1.5)
            loglog(freq,squeeze(abs(grid_bigH_orig(i,j,:))),'LineWidth',1.5)
%             Legend = legend('qd admittance');
        elseif n==3
%             title('dq')
%             loglog(freq,squeeze(abs(statcom_bigH_orig(i,j,:))),'g--','LineWidth',1.5)
%             loglog(freq_new,squeeze(abs(statcom_bigH(i,j,:))),'--','LineWidth',1.5)
            loglog(freq,squeeze(abs(grid_bigH_orig(i,j,:))),'g--','LineWidth',1.5)
%             Legend = legend('dq admittance');
        else
%             title('dd')
%             loglog(freq,squeeze(abs(statcom_bigH_orig(i,j,:))),'-','LineWidth',1.5)
%             loglog(freq_new,squeeze(abs(statcom_bigH(i,j,:))),'-','LineWidth',1.5)
            loglog(freq,squeeze(abs(grid_bigH_orig(i,j,:))),'--','LineWidth',1.5)
%             Legend = legend('dd admittance');
            
        end
        hold on
        Legend = legend('qq','qd','dq','dd');
        set([Legend, gca], 'FontSize', 12)
        set(gca, 'FontName', 'Helvetica')
        XLabel = xlabel('Frequency [Hz]');
        YLabel = ylabel('Magnitude [p.u.]');
        set([XLabel, YLabel], 'FontName', 'AvantGarde')
        set([XLabel, YLabel], 'FontSize', 14)
        grid on

         n=n+1;
    end
end

%%  Multiport vector fitting of system
% The vector fitting method is applied to both subsystems. The vector
% fitting method can be selected.

opts.poletype='logcmplx';   %Will use logarithmically spaced, complex poles. (Is used when opts.poles=[]).
poles=[];
opts.screen = 0; opt.weightparam = 3;

opts.errplot = 0;
opts.cmplx_ss = 1;
opts.remove_HFpoles = 0; opt.factor_HF = 1;
opts.stable = 1; opts.plot=1;



%Select model order and method for subsystem

opts.N = 22;
method = "MF";
% method = "SIMO";
% method = "PCCF";
% mintol = 1;
% tol = 3;

[system_grid, rmserr_grid] = Vf_driver_driver(grid_bigH,s_new,poles,opts,method,tol,mintol);


%Select model order and method for subsystem
opts.N = 27;

method = "MF";
% method = "SIMO";
% method = "PCCF";
% mintol = 1;
% tol = 3;

[system_statcom,rmserr_statcom] = Vf_driver_driver(statcom_bigH,s_new,poles,opts,method,tol,mintol);


figure(100)
plot(real(eig(system_statcom)),imag(eig(system_statcom))/(2*pi),'X') 
xlabel('Real sec^-1','FontSize', 14); ylabel('Imag Hz','FontSize', 14);
hold on
plot(real(eig(system_grid)),imag(eig(system_grid))/(2*pi),'*')
grid on
xlim([-60, 10]);
ylim([-100,100]);

disp(['STATCOM freq. scan fit model order: ' num2str(length(system_statcom.A))])
disp(['STATCOM freq. scan fit RMSE: ' + string(rmserr_statcom)])
disp(['Grid freq. scan model order: ' num2str(length(system_grid.A))])
disp(['grid freq. scan fit RMSE: ' + string(rmserr_grid)])


%% Order selection using HSV
%HSVs of both subsystems are plotted and evaluated


figure(100)
s = hsvd(system_statcom);
i = find(s==inf);
cumulativeHSV = cumsum(s) / sum(s);

tol = find(cumulativeHSV*100>=99.90 & cumulativeHSV*100<=99.991);
ind = find(cumulativeHSV*100>99.95);
x = ind(1);
x = 35;

yyaxis left
bar(s/s(1));
hold on
set(gca, 'YScale', 'log');
ylabel('Relative Contribution [p.u.]','FontName', 'AvantGarde', 'FontSize', 14)
ylim([1e-7,1.1])
ylimits = ylim;
% xlim([0,62])
xline(x,'--','LineWidth',1.5);
fill([tol(1),tol(1),tol(end),tol(end)],[ylimits(1),ylimits(2),ylimits(2),ylimits(1)],'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
y = 40e-2;
text(x+1, y,'Cut-off point: N = '+string(x), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right','Color', 'black', 'FontWeight', 'bold','Rotation',90);
% text(x+3, y,'Cum. value '+string(cumulativeHSV(x)*100), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right','Color', 'black', 'FontWeight', 'bold','Rotation',90);

yyaxis right
plot(cumulativeHSV*100, '-o')
ylabel('Cumulative contribution [%]','FontName', 'AvantGarde', 'FontSize', 14)
ylim([cumulativeHSV(1)*100,100])
xlabel('Order (number of states)','FontName', 'AvantGarde', 'FontSize', 14)
grid on

figure(101)
s = hsvd(system_grid);
cumulativeHSV = cumsum(s) / sum(s);
x = 29;
tol = find(cumulativeHSV*100>=99.90 & cumulativeHSV*100<=99.991);
ind = find(cumulativeHSV*100>99.95);

yyaxis left
bar(s/s(1));
hold on
set(gca, 'YScale', 'log');
ylabel('Relative Contribution [p.u.]','FontName', 'AvantGarde', 'FontSize', 14)
ylim([1e-5,1.1]);
ylimits = ylim;
% ylim([1e-5,1.1])
xline(x,'--','LineWidth',1.5);

fill([tol(1),tol(1),tol(end),tol(end)],[ylimits(1),ylimits(2),ylimits(2),ylimits(1)],'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
y = 20e-2;
text(x+1, y,'Cut-off point: N = '+string(x), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right','Color', 'black', 'FontWeight', 'bold','Rotation',90);

yyaxis right
plot(cumulativeHSV*100, '-o')
ylabel('Cumulative contribution [%]','FontName', 'AvantGarde', 'FontSize', 14)
ylim([cumulativeHSV(1)*100,100])
grid on
xlabel('Order (number of states)','FontName', 'AvantGarde', 'FontSize', 14)





%% Balanced Truncation
% Desired model order selected, after evaluating HSVs

statcom_rho = 31;
grid_rho = 28;


[sys_b,g] = balreal(system_grid);
p = g(grid_rho);
elim = find(g<p);
system_grid = modred(sys_b,elim,'del');

[sys_b,g] = balreal(system_statcom);
p = g(statcom_rho);
elim = find(g<p);
system_statcom = modred(sys_b,elim,'del');

disp(['New STATCOM freq. scan fit model order: ' num2str(length(system_statcom.A))])
disp(['New Grid freq. scan fit model order: ' num2str(length(system_grid.A))])

figure(1)
plot(real(eig(system_statcom)),imag(eig(system_statcom)/(2*pi)),'X') 
hold on
plot(real(eig(system_grid)),imag(eig(system_grid)/(2*pi)),"diamond") 

xlabel('Real sec^-1');
ylabel('Imaginary [Hz]');
grid on
xlim([-60, 10]);
ylim([-100,100]);
legend('STATCOM poles','Grid poles')


%%  Connect the systems

% close all

% Define the interconnection using specific naming conventions

system_statcom.InputName = {'u1_1', 'u1_2'};
system_statcom.OutputName = {'y1_1', 'y1_2'};

%define state names of systems for later use
n1 = length(system_grid.A);
system_grid.StateName = compose('grid x_%d', 1:n1);

n2 = length(system_statcom.A);
system_statcom.StateName = compose('stat x_%d', 1:n2);

system_grid.InputName = system_statcom.OutputName;
system_grid.OutputName = {'y2_1', 'y2_2'};

% Define the external input and output names for the overall system
externalInputs = {'u1_1', 'u1_2'};
externalOutputs = {'y2_1', 'y2_2'};

% Interconnect the systems using the input and output names
glob_sys = connect(system_statcom,system_grid, externalInputs, externalOutputs);

figure(12)
plot(real(eig(glob_sys)),imag(eig(glob_sys))/(2*pi),"diamond",'MarkerSize',8)
hold on
xlabel('Real sec^{-1}','FontSize', 14);
ylabel('Imaginary [Hz]','FontSize', 14);
grid on
% legend('SCR 2.5','SCR 1.5','FontSize', 12)
legend('Model eigenvalues','FontSize', 12)
title(title_name + ' using ' + method_name + ', relative err: ' + string(rel_err) + ', ' + id_values+'/'+total_values + ' states')












