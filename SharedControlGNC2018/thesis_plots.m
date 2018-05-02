%% Anomalous Actuators 
anom = AnomalyResponseSim('act',0);

SP = anom.simPars;
ts    = SP.ts;
ts_rm = SP.ts_rm;
tfin  = SP.tfin;

PSP = anom.postSimPars;
i_tsrm = PSP.i_tsrm;
tsim   = PSP.tsim;
norm_e = PSP.norm_e;

c2 = [0.466, 0.674, 0.188];
% vfa.pltOpt.legfontsize = 14;
% vfa.pltOpt.fontsize = 16;
% vfa.pltOpt.weight = 'n';
% vfa.pltOpt.fontname = 'Times New Roman';

set(0,'defaultAxesFontSize', 16);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%% Plot 1

figure('Position',[100, 100, 800, 250])
plot(tsim, PSP.r_deg, 'LineWidth', 1.5); hold on; grid on;
plot(tsim, PSP.x_deg(1,:), 'LineWidth', 1.5, 'color', c2);
ylim([-2 25])
line([ts ts],ylim,'Color',[0 0 0],'LineStyle','-.', 'LineWidth', 1);
line([ts_rm ts_rm], ylim,'Color',[0 0 0],'LineStyle','-.', 'LineWidth', 1);
ylabel('$\phi$ (deg)', 'interpreter', 'latex')
xlabel('Time (s)')
title('Bank Angle (deg)', 'interpreter', 'latex')
h = legend('Command ($r$)', 'Output ($\phi$)')
set(h,'Interpreter','Latex'); legend('boxoff')
tightfig()

%% Figure 2: Adaptive Parameters
figure('Position',[100, 100, 800, 250])
plot(tsim, PSP.th1, 'LineWidth', 1.5); hold on; grid on;
plot(tsim, PSP.th2, 'LineWidth', 1.5); plot(tsim(i_tsrm+1:end), PSP.th3, 'LineWidth', 1.5);
plot(tsim, PSP.q, 'LineWidth', 1.5);
ylim([-50 50])
line([ts ts],ylim,'Color',[0 0 0],'LineStyle','-.', 'LineWidth', 1);
line([ts_rm ts_rm], ylim,'Color',[0 0 0],'LineStyle','-.', 'LineWidth', 1);
ylabel('$\theta_i$', 'interpreter', 'latex')
xlabel('Time (s)')
title('Adaptive Parameter Values', 'interpreter', 'latex')
h = legend('$\theta_1$', '$\theta_2$', '$\theta_3$', '$q$')
set(h,'Interpreter','Latex'); legend('boxoff')
tightfig()

%% Figure 3: Errors 
% f2 = figure('Position',[1,1, 800, 240]);
% set(f2,'defaultAxesColorOrder',[c1; c2]);
% yyaxis left
% plot(time, err_norm, 'LineWidth', 1.5); grid on; hold on;
% xlim([0 1500])
% ylim([0 4e-3])
% line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
% line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
% line([SOO2.t_sim(end) SOO2.t_sim(end)],ylim,'Color',[0.6 0.6 0.6],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
% yyaxis right
% plot(time, err_norm2, 'LineWidth', 1.5, 'Color', c2); grid on; hold on;
% ylim([0 15])
% title('Output Error Signals: $\|y(t)-y_m(t)\|_2$ and $\|z(t)-z_{cmd}(t)\|_2$', 'interpreter','latex')
% xlabel('Time (s)')
% h=legend('$\|y(t)-y_m(t)\|_2$', '$\|z(t)-z_{cmd}(t)\|_2$');
% set(h,'fontsize',vfa1.pltOpt.legfontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname,'Interpreter','Latex','Location','SouthEast'); legend('boxoff')
% set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)
% tightfig(f2);
norm_e = zeros(1, length(PSP.tsim));
norm_e(1:i_tsrm) = (squeeze(anom.controlSim.e(1,:,1:i_tsrm)).^2 + squeeze(anom.controlSim.e(2,:,1:i_tsrm)).^2).^(0.5);
norm_e(i_tsrm+1:end) = (squeeze(anom.controlSim.e_3d(1,:,i_tsrm+1:end)).^2 + squeeze(anom.controlSim.e_3d(2,:,i_tsrm+1:end)).^2).^(0.5);

% norm_e_cmd = zeros(1, length(PSP.tsim));
% norm_e_cmd = (squeeze(anom.controlSim.r(1,:,:) - anom.controlSim.x(1,:,:)).^2 + squeeze(anom.controlSim.r(2,:,:) - anom.controlSim.x(2,:,:)).^2).^(0.5);

figure('Position',[100, 100, 800, 250])
% yyaxis left
plot(tsim, norm_e, 'LineWidth', 1.5); hold on; grid on;
ylim([0 0.3])
line([ts ts],ylim,'Color',[0 0 0],'LineStyle','-.', 'LineWidth', 1);
line([ts_rm ts_rm], ylim,'Color',[0 0 0],'LineStyle','-.', 'LineWidth', 1);
ylabel('$\|e\|$ ', 'interpreter', 'latex')
xlabel('Time (s)')
% yyaxis right
% plot(tsim, norm_e_cmd, 'LineWidth', 1.5); hold on; grid on;
title('Magnitude of Error ($e = \left[e_\phi, \quad e_p\right]^T)$', 'interpreter', 'latex')
tightfig()

%% Figure 4: step response

% last 5 seconds of stage 2
i21 = ceil((ts_rm-5)/SP.dT)+1;
i22 = i_tsrm;
theta1_bar_t2 = mean(PSP.th1(i21:i22));
theta2_bar_t2 = mean(PSP.th2(i21:i22));
qbar_t2       = mean(PSP.q(i21:i22));
% last 5 seconds of stage 3
i31 = ceil((tfin-5)/SP.dT)+1;
theta1_bar_t3 = mean(PSP.th1(i31:end));
theta2_bar_t3 = mean(PSP.th2(i31:end));
theta3_bar_t3 = mean(PSP.th3_aug(i31:end));
qbar_t3       = mean(PSP.q(i31:end));

s = tf('s');
rm2 = 8/(s^2+6*s+8);
rm3 = 32/((s^2+6*s+8)*(s+4));
Y_cl_trans_1 = (0.318*qbar_t2)/(s^3 + 2.10*s^2 + (1.10 - 0.318*theta2_bar_t2)*s - 0.318*theta1_bar_t2); % using thetas, q from end of adjustment phase
Y_cl_3 = (0.318*qbar_t3)/(s^3 + (2.10 - 0.318*theta3_bar_t3)*s^2 + (1.10 - 0.318*theta2_bar_t3)*s - 0.318*theta1_bar_t3); % using thetas, q from end of adjustment phase
% step(rm2, Y_cl_trans_1, Y_cl_3, rm3, 8); 
[y1, t1] = step(rm2, 6);
[y2, t2] = step(Y_cl_trans_1, 6); 
[y3, t3] = step(Y_cl_3, 6); 
[y4, t4] = step(rm3, 6); 

figure('Position',[100, 100, 800, 250]);
plot(t1, y1, 'LineWidth', 1.5); hold on; grid on;
plot(t2, y2, 'LineWidth', 1.5, 'LineStyle', '-.');
plot(t3, y3, 'LineWidth', 1.5, 'LineStyle', '--');
plot(t4, y4, 'LineWidth', 1.5);
title('Dynamic Response to Unit Step Command with Anomalous Actuators', 'interpreter', 'latex'); 
xlabel('Time (s)')
ylim([0, 1.2])
h = legend('Second-Order Reference Model', 'Pre-Correction ($t < t_2^*$)', 'Post-Correction ($t > t_2^*$)', 'Third-Order Reference Model', 'location', 'southeast'); 
set(h,'Interpreter','Latex'); legend('boxoff')
tightfig()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Anomalous Delay 
% clear; close all;
anom = AnomalyResponseSim('del',0);

SP = anom.simPars;
ts    = SP.ts;
ts_rm = SP.ts_rm;
tfin  = SP.tfin;

PSP = anom.postSimPars;
i_tsrm = PSP.i_tsrm;
tsim   = PSP.tsim;
norm_e = PSP.norm_e;

c2 = [0.466, 0.674, 0.188];
% vfa.pltOpt.legfontsize = 14;
% vfa.pltOpt.fontsize = 16;
% vfa.pltOpt.weight = 'n';
% vfa.pltOpt.fontname = 'Times New Roman';

set(0,'defaultAxesFontSize', 16);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%% Calcs
bp2 = 1.0;
bm_3 = 32;
ppole_1 = 0;
ppole_2 = 1.10;
ppole_3 = 1.0;

lambda_3d = bp2/bm_3;
q_0_3d = 1/lambda_3d
theta_0_1_3d = ((ppole_1*ppole_2*ppole_3) - 32) / bp2;
theta_0_2_3d = ((ppole_1*ppole_2 + ppole_1*ppole_3 + ppole_2*ppole_3) - 32) / bp2;
theta_0_3_3d = ((ppole_1 + ppole_2 + ppole_3)) / bp2;
Theta_0_3d = [theta_0_1_3d theta_0_2_3d theta_0_3_3d]

%% Plot 1

figure('Position',[100, 100, 800, 250])
% plot(tsim, PSP.r_deg, 'LineWidth', 1.5); hold on; grid on;
plot(t_nom, r_nom, 'LineWidth', 1.5); hold on; grid on;
plot(tsim, PSP.x_deg(1,:), 'LineWidth', 1.5, 'color', c2);
ylim([-2 25])
line([ts ts],ylim,'Color',[0 0 0],'LineStyle','-.', 'LineWidth', 1);
line([tsim(end) tsim(end)], ylim,'Color',[0 0 0],'LineStyle','-.', 'LineWidth', 1);
ylabel('$\phi$ (deg)', 'interpreter', 'latex')
xlabel('Time (s)')
title('Bank Angle (deg)', 'interpreter', 'latex')
h = legend('Command ($r$)', 'Output ($\phi$)')
set(h,'Interpreter','Latex'); legend('boxoff')
tightfig()

%% Figure 2: Adaptive Parameters
figure('Position',[100, 100, 800, 250])
plot(tsim, PSP.th1, 'LineWidth', 1.5); hold on; grid on;
plot(tsim, PSP.th2, 'LineWidth', 1.5); plot(tsim(i_tsrm+1:end), PSP.th3, 'LineWidth', 1.5);
plot(tsim, PSP.q, 'LineWidth', 1.5);
ylim([-40 40])
line([ts ts],ylim,'Color',[0 0 0],'LineStyle','-.', 'LineWidth', 1);
line([ts_rm ts_rm], ylim,'Color',[0 0 0],'LineStyle','-.', 'LineWidth', 1);
ylabel('$\theta_i$', 'interpreter', 'latex')
xlabel('Time (s)')
title('Adaptive Parameter Values', 'interpreter', 'latex')
h = legend('$\theta_1$', '$\theta_2$', '$\theta_3$', '$q$')
set(h,'Interpreter','Latex'); legend('boxoff')
tightfig()

%% Figure 3: Errors 

figure('Position',[100, 100, 800, 250])
% yyaxis left
plot(tsim, norm_e, 'LineWidth', 1); hold on; grid on;
ylim([0 0.3])
line([ts ts],ylim,'Color',[0 0 0],'LineStyle','-.', 'LineWidth', 1);
line([ts_rm ts_rm], ylim,'Color',[0 0 0],'LineStyle','-.', 'LineWidth', 1);
ylabel('$\|e\|$ ', 'interpreter', 'latex')
xlabel('Time (s)')
% yyaxis right
% plot(tsim, norm_e_cmd, 'LineWidth', 1.5); hold on; grid on;
title('Magnitude of Error ($e = \left[e_\phi, \quad e_p\right]^T)$', 'interpreter', 'latex')
tightfig()

%% Figure 4: step response

% last 5 seconds of stage 2
i21 = ceil((ts_rm-5)/SP.dT)+1;
i22 = i_tsrm;
theta1_bar_t2 = mean(PSP.th1(i21:i22));
theta2_bar_t2 = mean(PSP.th2(i21:i22));
qbar_t2       = mean(PSP.q(i21:i22));
% last 5 seconds of stage 3
i31 = ceil((tfin-5)/SP.dT)+1;
theta1_bar_t3 = mean(PSP.th1(i31:end));
theta2_bar_t3 = mean(PSP.th2(i31:end));
theta3_bar_t3 = mean(PSP.th3_aug(i31:end));
qbar_t3       = mean(PSP.q(i31:end));

s = tf('s');
rm2 = 8/(s^2+6*s+8);
rm3 = 32/((s^2+6*s+8)*(s+4));
tau = SP.tau;
Y_cl_trans_1 = (0.318*qbar_t2)/(s^2 + (1.98 - 0.318*theta2_bar_t2 * exp(-tau*s))*s  - 0.318*theta1_bar_t2 * exp(-tau*s)); % using thetas, q from end of adjustment phase
Y_cl_3 = (0.318*qbar_t3)/((1 - 0.318*theta3_bar_t3 * exp(-tau*s))*s^2 + (1.98 - 0.318*theta2_bar_t3 * exp(-tau*s))*s  - 0.318*theta1_bar_t3 * exp(-tau*s)); % using thetas, q from end of adjustment phase
[y1, t1] = step(rm2, 6);
[y2, t2] = step(Y_cl_trans_1, 6); 
[y3, t3] = step(Y_cl_3, 6); 
[y4, t4] = step(rm3, 6); 

figure('Position',[100, 100, 800, 250]);
plot(t1, y1, 'LineWidth', 1.5); hold on; grid on;
plot(t2, y2, 'LineWidth', 1.5, 'LineStyle', '-.');
plot(t3, y3, 'LineWidth', 1.5, 'LineStyle', '--');
plot(t4, y4, 'LineWidth', 1.5);
title('Dynamic Response to Unit Step Command with Sensor Delay', 'interpreter', 'latex'); 
xlabel('Time (s)')
ylim([0, 1.2])
h = legend('Second-Order Reference Model', 'Pre-Correction ($t < t_2^*$)', 'Post-Correction ($t > t_2^*$)', 'Third-Order Reference Model', 'location', 'southeast'); 
set(h,'Interpreter','Latex'); legend('boxoff')
tightfig()
