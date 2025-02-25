% Overall code for calculating Radial Distraction Osteogenesis (RDO) distraction forces in the fibula:
clc
close all
clear all
l_0_l = 1; % initial gap width for callus formation
l_delta_l = 0.33; % distraction step size for purely lateral strain
DL = 10; % distraction length in mm
di = 4; % initial diameter of the expansion unit in mm
S = round(DL/l_delta_l); % number of distraction steps
T_s = 8; % time step size in hours
k = 1.5; % scaling factor from Meyers et al.
k_A = 1.5; % area scaling
T_end = S*T_s*60*60; % total duration of distraction in seconds (each step 8 hours)
DT = 1; % time step size
TS = T_s*60*60; % time step size per expansion step in seconds
Td = round(S/(24/T_s)); % distraction days
DS = round(DL/10); % scaling for graphs
% Estimate of the initially distracted area (see Excel spreadsheet)
A_ini = 0.001051; % distracted area in m^2
LB = 200; % length of the expansion unit (mm)
% Initialization of vectors
t = [0:DT:T_end-1]; % time vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lateral Distraction Region
% Initialization of vectors
e_iC = zeros(1,length(t)); % strain vector
E_tiL = zeros(1,length(t)); % E-modulus vector
s_tiL = zeros(1,length(t)); % stress vector
F_tiL = zeros(1,length(t)); % force vector
p_tiL = zeros(1,length(t)); % pressure vector
A_tiL = zeros(1,length(t)); % area vector
a_uiL = zeros(1,length(t)); % equilibrium coefficient
a_fiL = zeros(1,length(t)); % fast coefficient
a_siL = zeros(1,length(t)); % slow coefficient
t_siL = zeros(1,length(t)); % slow time constant
t_fiL = zeros(1,length(t)); % fast time constant
E_0iL = zeros(1,length(t));
A_iL = zeros(1,length(t));
% Initialization of initial values
e_iC(1) = 0; % callus strain is 0 at the start before distraction
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
    a_fiL(tclock_L) = alpha_f(i);
    a_siL(tclock_L) = alpha_s(a_uiL(2),a_fiL(2));
    t_siL(tclock_L) = tau_s(i);
    t_fiL(tclock_L) = 0.97;
    E_0iL(tclock_L) = E_modulus_0(i);
    A_iL(tclock_L) = A_ini;
    % Calculation of the E-modulus
    E_tiL(tclock_L) = E_modulus(E_0iL(2), a_uiL(2), a_fiL(2), a_siL(2), t_fiL(2), t_siL(2), tstep_L, k); % calculation for the first value per step
    % Incrementing the local and global timers
    tstep_L = tstep_L + 1;
    tclock_L = tclock_L + 1;
end
while i < S % Calculation of all subsequent distraction steps
    i = i + 1; % increment distraction step
    tstep_L = 0; % reset the local timer to 0 for the next distraction step
    while (tclock_L > (i-1)*TS && tclock_L <= i*TS)
        a_uiL(tclock_L) = alpha_u(i);
        a_fiL(tclock_L) = alpha_f(i);
        a_siL(tclock_L) = alpha_s(a_uiL(((i-1)*TS)+1), a_fiL(((i-1)*TS)+1));
        t_siL(tclock_L) = tau_s(i);
        t_fiL(tclock_L) = 0.97;
        E_0iL(tclock_L) = E_modulus_0(i);
        A_iL(tclock_L) = A_ini;
        E_tiL(tclock_L) = E_modulus(E_0iL(((i-1)*TS)+1), a_uiL(((i-1)*TS)+1), a_fiL(((i-1)*TS)+1), a_siL(((i-1)*TS)+1), t_fiL(((i-1)*TS)+1), t_siL(((i-1)*TS)+1), tstep_L, k);
        tstep_L = tstep_L + 1;
        tclock_L = tclock_L + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the distraction vector and the area vector for the circumferential expansion
d_c = zeros(1,S+1); % initialization of distraction vector
A_c = zeros(1,S); % initialization of area vector
delta_d_c = zeros(1,S); % initialization of difference distraction vector
i_d = 0; % loop variable
i_A = 1; % loop variable
wd = 0; % gap width
de = 0; % radius of the expansion unit
% Calculation of the distraction vector
for o = 1:length(d_c)
    wd = i_d * l_delta_l + l_0_l;
    de = i_d * l_delta_l + di;
    d_c(o) = ((asin((wd/2)/(de/2)))/(360*(pi/180))*(2*pi*(de/2)))*2;
    i_d = i_d + 1;
end
% % % % Start with l0 equal to the gap width and not directly curved
d_c(1) = l_0_l;
% Calculation of the difference distraction vector to determine the change in length
for k = 1:length(delta_d_c)
    delta_d_c(k) = (d_c(k+1) - d_c(k));
end
% Calculation of the area vector
for m = 1:length(A_c)
    de = i_A * l_delta_l + di;
    A_c(m) = 2*pi*(de/2*0.001)*(LB*0.001); % lateral surface area of the expansion unit in m^2
    i_A = i_A + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Circumferential Distraction Region
% Initialization of vectors
e_iC = zeros(1,length(t)); % strain vector
s_tiC = zeros(1,length(t)); % stress vector
F_tiC = zeros(1,length(t)); % force vector
a_uiC = zeros(1,length(t)); % equilibrium coefficient
% Initialization of initial values
e_iC(1) = 0; % callus strain is 0 at the start before distraction
a_uiC(1) = 0;
A_tiL(1) = 2*pi*(di/2*0.001)*LB*0.001;
l_0_c = d_c(1);
tclock_C = 2; % global timer = 2 (initial state at 1, start at 2)
tstep_C = 0; % local timer for each distraction step
% Calculation of strain and E-modulus for each time step
% Calculation of the first step
% step 1
i = 1; % distraction step
l_delta_c = delta_d_c(i);
while(tclock_C >= 2 && tclock_C <= TS) % global timer for the duration of the current distraction step
    % Calculation of the time constants and coefficients for the distraction step
    a_uiC(tclock_C) = alpha_u(i); % store values in the vector one index after the initial values
    % Calculation of the strain
    e_iC(tclock_C) = strain_C(l_delta_c, l_0_c, i, e_iC(1), a_uiC(1)); % substituting the initial values
    A_tiL(tclock_C) = A_c(1);
    % Incrementing the local and global timers
    tstep_C = tstep_C + 1;
    tclock_C = tclock_C + 1;
end
while i < S % Calculation of all subsequent distraction steps
    i = i + 1; % increment distraction step
    l_delta_c = delta_d_c(i); % determine distraction distance
    tstep_C = 0; % reset the local timer to 0 for the next distraction step
    while (tclock_C > (i-1)*TS && tclock_C <= i*TS)
        a_uiC(tclock_C) = alpha_u(i);
        e_iC(tclock_C) = strain_C(l_delta_c, l_0_c, i, e_iC((i-1)*TS), a_uiC((i-1)*TS));
        A_tiL(tclock_C) = A_c(i);
        tstep_C = tstep_C + 1;
        tclock_C = tclock_C + 1;
    end
end
% Calculation of the stress progression
for x = 1:length(e_iC)
    s_tiL(x) = e_iC(x) * E_tiL(x);
end
% Calculation of the force progression in Newtons (stress in kPa -> multiplied by 1000 to obtain Pascals)
for y = 1:length(s_tiL)
    F_tiL(y) = (s_tiL(y)*1000) * (A_iL(y)*k_A);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization of the results in graphs for strain, E-modulus,
% stress, and force
figure(1)
subplot(3,1,1)
hold on
grid on
plot(t,E_tiL,'r','LineWidth',3);
set(gca, 'XTick', 0:(length(t)*DS/(Td)):T_end, 'XTickLabel', 0:DS:Td)
%title('E-Modulus')
ax = gca;
ax.GridAlpha = 0.9;
set(gca,'FontSize',30);
set(gca,'FontWeight','bold')
xlabel('t [d]')
ylabel('E [kPa]')
axis([0 S*TS 0 E_tiL(S*TS-(TS-1))*1.1])
yticks([0:500:E_tiL(S*TS-(TS-1))*1.1]);
%figure(2)
subplot(3,1,2)
hold on
grid on
plot(t,e_iC,'g','LineWidth',3);
set(gca, 'XTick', 0:(length(t)*DS/(Td)):T_end, 'XTickLabel', 0:DS:Td)
%title('Lateral Strain')
ax = gca;
ax.GridAlpha = 0.9;
set(gca,'FontSize',30);
set(gca,'FontWeight','bold')
xlabel('t [d]')
ylabel('\epsilon []')
axis([0 S*TS 0 0.4])
yticks([0:0.1:0.4]);
%figure(2)
subplot(3,1,3)
hold on
grid on
plot(t,s_tiL,'k','LineWidth',3);
set(gca, 'XTick', 0:(length(t)*DS/(Td)):T_end, 'XTickLabel', 0:DS:Td)
%title('Lateral Stress')
ax = gca;
ax.GridAlpha = 0.9;
set(gca,'FontSize',30);
set(gca,'FontWeight','bold')
xlabel('t [d]')
ylabel('\sigma [kPa]')
yticks([0:20:s_tiL(S*TS-(TS-1))*1.1]);
axis([0 S*TS 0 s_tiL(S*TS-(TS-1))*1.1])
figure(3)
% subplot(4,1,4)
hold on
grid on
plot(t,F_tiL,'b','LineWidth',2);
set(gca, 'XTick', 0:(length(t)*DS/(Td)):T_end, 'XTickLabel', 0:DS:Td)
%title('Callus Distraction Force')
ax = gca;
ax.GridAlpha = 0.9;
set(gca,'FontSize',30);
set(gca,'FontWeight','bold')
xlabel('t [d]')
ylabel('CDF [N]')
yticks([0:20:F_tiL(S*TS-(TS-1))*1.1]);
axis([0 S*TS 0 F_tiL(S*TS-(TS-1))*1.1])
set(gcf,'position',[10,10,1900,600])