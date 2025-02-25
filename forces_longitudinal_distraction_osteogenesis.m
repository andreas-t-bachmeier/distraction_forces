%Overall code for calculating Longitudinal Distraction Osteogenesis (LDO) distraction forces:
clc
close all
clear all
l_0_l = 1; % initial gap width for callus formation
l_delta_l = 0.33; % distraction step size for purely lateral strain
DL = 68; % defect length in mm
S = round(DL/l_delta_l); % number of distraction steps
T_s = 8; % time step size in hours
k = 1.5; % scaling factor from Meyers et al.
k_A = 1.1; % area scaling
T_end = S*T_s*60*60; % total duration of distraction in seconds (each step 8 hours)
DT = 1; % time step size
TS = T_s*60*60; % time step size per expansion step in seconds
Td = round(S/(24/T_s)); % distraction days
DS = round(DL/22); % scaling for graphs
% Estimate of the initially distracted area
A_ini = 0.000849; % distracted area in m^2
% Initialization of vectors
t = [0:DT:T_end-1]; % time vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lateral Distraction Region
% Initialization of vectors
e_iL = zeros(1,length(t)); % strain vector
E_tiL = zeros(1,length(t)); % E-modulus vector
s_tiL = zeros(1,length(t)); % stress vector
F_tiL = zeros(1,length(t)); % force vector
a_uiL = zeros(1,length(t)); % equilibrium coefficient
a_fiL = zeros(1,length(t)); % fast coefficient
a_siL = zeros(1,length(t)); % slow coefficient
t_siL = zeros(1,length(t)); % slow time constant
t_fiL = zeros(1,length(t)); % fast time constant
E_0iL = zeros(1,length(t));
A_iL = zeros(1,length(t));
% Initialization of initial values
e_iL(1) = 0; % strain of the callus is 0 at the start before distraction
a_uiL(1) = 0;
tclock_L = 2; % global timer = 2 (initial state at 1, start at 2)
tstep_L = 0; % local timer for each distraction step
% Calculation of strain and E-modulus for each time step
% Calculation of the first step
% step 1
i = 1; % distraction step
while(tclock_L >= 2 && tclock_L <= TS) % global timer for the duration of the current distraction step
    % Calculation of the time constants and coefficients for the distraction step
    a_uiL(tclock_L) = alpha_u(i); % store values in the vector one index after the initial values
    a_fiL(tclock_L) = alpha_f_mod(i,S);
    a_siL(tclock_L) = alpha_s(a_uiL(2),a_fiL(2));
    t_siL(tclock_L) = tau_s(i);
    t_fiL(tclock_L) = 0.97;
    E_0iL(tclock_L) = E_modulus_0(i);
    A_iL(tclock_L) = A_ini;
    % Calculation of the strain
    e_iL(tclock_L) = strain_L(l_delta_l,l_0_l,i,e_iL(1),a_uiL(1)); % substituting the initial values
    % Calculation of the E-modulus
    E_tiL(tclock_L) = E_modulus(E_0iL(2),a_uiL(2),a_fiL(2),a_siL(2),t_fiL(2),t_siL(2),tstep_L,k); % calculation for the first value per step
    % Incrementing the local and global timers
    tstep_L = tstep_L + 1;
    tclock_L = tclock_L + 1;
end
while i < S % Calculation of all subsequent distraction steps
    i = i + 1; % increment distraction step
    tstep_L = 0; % reset the local timer to 0 for the next distraction step
    while (tclock_L > (i-1)*TS && tclock_L <= i*TS)
        a_uiL(tclock_L) = alpha_u(i);
        a_fiL(tclock_L) = alpha_f_mod(i,S);
        a_siL(tclock_L) = alpha_s(a_uiL(((i-1)*TS)+1),a_fiL(((i-1)*TS)+1));
        t_siL(tclock_L) = tau_s(i);
        t_fiL(tclock_L) = 0.97;
        E_0iL(tclock_L) = E_modulus_0(i);
        A_iL(tclock_L) = A_ini;
        e_iL(tclock_L) = strain_L(l_delta_l,l_0_l,i,e_iL((i-1)*TS),a_uiL((i-1)*TS));
        E_tiL(tclock_L) = E_modulus(E_0iL(((i-1)*TS)+1),a_uiL(((i-1)*TS)+1),a_fiL(((i-1)*TS)+1),a_siL(((i-1)*TS)+1),t_fiL(((i-1)*TS)+1),t_siL(((i-1)*TS)+1),tstep_L,k);
    
        tstep_L = tstep_L + 1;
        tclock_L = tclock_L + 1;
    end
end
% Calculation of the stress progression
for x = 1:length(e_iL)
    s_tiL(x) = e_iL(x) * E_tiL(x);
end
% Calculation of the force progression in Newtons (stress in kPa -> multiplied by 1000 to obtain Pascals)
for y = 1:length(s_tiL)
    F_tiL(y) = (s_tiL(y)*1000) * (A_iL(y)*k_A);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization of the results in graphs for strain, E-modulus,
% stress, and force
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lateral Distraction
%figure(1)
subplot(4,1,1)
hold on
grid on
plot(t,E_tiL,'r','LineWidth',2);
set(gca, 'XTick', 0:(length(t)*DS/(Td)):T_end, 'XTickLabel', 0:DS:Td)
%title('E-Modulus')
ax = gca;
ax.GridAlpha = 0.9;
set(gca,'FontSize',15);
set(gca,'FontWeight','bold')
xlabel('t [d]')
ylabel('E [kPa]')
axis([0 S*TS 0 12000])
yticks([0:3000:13000]);
% figure(2)
subplot(4,1,2)
hold on
grid on
plot(t,e_iL,'g','LineWidth',2);
set(gca, 'XTick', 0:(length(t)*DS/(Td)):T_end, 'XTickLabel', 0:DS:Td)
%title('Lateral Strain')
ax = gca;
ax.GridAlpha = 0.9;
set(gca,'FontSize',15);
set(gca,'FontWeight','bold')
xlabel('t [d]')
ylabel('\epsilon []')
axis([0 S*TS 0 0.4])
yticks([0:0.1:0.75]);
subplot(4,1,3)
hold on
grid on
plot(t,s_tiL,'k','LineWidth',2);
set(gca, 'XTick', 0:(length(t)*DS/(Td)):T_end, 'XTickLabel', 0:DS:Td)
%title('Lateral Stress')
ax = gca;
ax.GridAlpha = 0.9;
set(gca,'FontSize',15);
set(gca,'FontWeight','bold')
xlabel('t [d]')
ylabel('\sigma [kPa]')
yticks([0:20:80]);
axis([0 S*TS 0 80])
figure(2)
%subplot(4,1,4)
hold on
grid on
plot(t,F_tiL,'b','LineWidth',2);
%set(gca, 'XTick', 0:(length(t)*DS/(Td)):T_end, 'XTickLabel', 0:DS:Td)
%title('Callus Distraction Force')
ax = gca;
ax.GridAlpha = 0.9;
set(gca,'FontSize',30);
set(gca,'FontWeight','bold')
xlabel('t [d]')
ylabel('CDF [N]')
yticks([0:10:90]);
%axis([0 S*TS 0 90])
set(gcf,'position',[10,10,1900,600])