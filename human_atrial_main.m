%% Matlab code of the human atrial myocyte model (2023)
% 
% Reference:
% Enhanced Ca2+-Dependent SK-Channel Gating and Membrane Trafficking in Human Atrial Fibrillation.
% Heijman J, Zhou X, Morotti S, Molina CE, Abu-Taha IH, Tekook M, Jespersen T, Zhang Y, Dobrev S,
% Milting H, Gummert J, Karck M, Kamler M, El-Armouche A, Saljic A, Grandi E, Nattel S, Dobrev D.
% Circ Res. 2023. doi: 10.1161/CIRCRESAHA.122.321858.
 
clear
close all
clc

%% Initial Conditions

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nSR
%load yfin_nSR_0p2Hz
%load yfin_nSR_0p5Hz
load yfin_nSR_1Hz
%load yfin_nSR_2Hz
%load yfin_nSR_3Hz
%load yfin_nSR_4Hz
% cAF
%load yfin_cAF_0p2Hz
%load yfin_cAF_0p5Hz
%load yfin_cAF_1Hz
%load yfin_cAF_2Hz
%load yfin_cAF_3Hz
%load yfin_cAF_4Hz
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SK block
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nSR
%load yfin_nSR_SKblock_0p2Hz
%load yfin_nSR_SKblock_0p5Hz
%load yfin_nSR_SKblock_1Hz
%load yfin_nSR_SKblock_2Hz
%load yfin_nSR_SKblock_3Hz
%load yfin_nSR_SKblock_4Hz
% cAF
%load yfin_cAF_SKblock_0p2Hz
%load yfin_cAF_SKblock_0p5Hz
%load yfin_cAF_SKblock_1Hz
%load yfin_cAF_SKblock_2Hz
%load yfin_cAF_SKblock_3Hz
%load yfin_cAF_SKblock_4Hz
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y0 = yfinal;

%% Simulation parameters
% Stimulation protocol parameters
prot_index = 1; % (-)
% 0) no stimulation 1) normal pacing 2) ERP 3) V-step 4) caffeine
prot_rate = 1; % (Hz)
prot_interval = 100; % (ms) not used
prot_vm = 10; % (mV) not used
prot_par = [prot_index prot_rate prot_interval prot_vm]; % 1 2 3 4

% Cell model parameters
% AF flag
AF_flag = 0; % Set 0 for nSR or 1 for cAF

% ISO flag
ISO_flag = 0; % (boolean)
% CCh concentration
CCh = 0;
% ACh concentration
ACh = 0;
% IK2P parameter
K2P_cond = 1; % with 0 IK2P block (default: 1)
% ISK parameters
SK_cond = 1; % with 0 ISK block (default: 1)
SK_shift = 0; % shift of Ca-sens (default: 0) kd = 10^(shift-3.45) % (mM)

% Ca buffering
Ca_clamp_index = 0; % 0 none, 1 Ca_clamp, 2 BAPTA, 3 EGTA
if Ca_clamp_index == 1 % [Ca]c,sl,j clamped to Cai0
    Cai0 = 500e-6; 
    y0(36) = Cai0; y0(37) = Cai0; y0(38) = Cai0;
end
if Ca_clamp_index == 2 % [Ca]c,sl,j clamped to Cai0_Bapta
    Cai0_Bapta = 5*100e-6;
    y0(36) = Cai0_Bapta; y0(37) = Cai0_Bapta; y0(38) = Cai0_Bapta;
end
if Ca_clamp_index == 3 % [Ca]c clamped to Cai0_Egta
    Cai0_Egta = 500e-6;
    y0(38) = Cai0_Egta;
end

% Na buffering
Na_clamp_index = 0; % 0 none, 1 Na_clamp
if Na_clamp_index == 1 % [Na]c clamped to Nai0
    Nai0 = 8; 
    y0(34) = Nai0;
end
    
cell_par = [AF_flag ISO_flag CCh ACh K2P_cond SK_cond SK_shift Ca_clamp_index Na_clamp_index]; % 5-13

% Sensitivity analysis parameters
% SA_par: 1) GNa 2) GNaL 3) GNaB 4) vNaK 5) Gtof 6) GKr 7) GKs
%         8) GKur 9) GK1 10) GKp 11) GK2P 12) GKACh 13) GSK 14) KdSK
%         15) GClCa 16) GClB 17) GCaL 18) GCaB 19) vPMCA 20) vNCX
%         21) vSERCA 22) vRel 23) vLeak
% Additional scaling factors: 24) GK2P-cAF 25) GSK-cAF
% 26) GCaB-cAF 27) vNCX-cAF 28) koCa-cAF 29) vSERCA-cAF 30) kmf-cAF
% 31) kmf 32) kmnai

SA_par = ones(1,33);

% Parameter array for passing nondefault conditions
p = [prot_par cell_par SA_par];

% Simulation duration
duration = 10e3; % (ms)

%% Single Run Simulation
tic
tspan = [0; duration];
options = odeset('RelTol',1e-6,'MaxStep',1,'Stats','off'); 
[t,y] = ode15s(@human_atrial_model,tspan,y0,options,p);
toc

%% Final Conditions
yfinal = y(end,:);

%save yfin_XXX yfinal

%% Figures
figure(1);
subplot(4,1,1); hold on,plot(t,y(:,39)); ylabel('Em (mV)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,2); hold on,plot(t,y(:,38)); ylabel('[Ca]i (mM)');
%plot(t,y(:,37),'r',t,y(:,36),'m')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,3); hold on,plot(t,y(:,34)); ylabel('[Na]i (mM)'); %
%plot(t,y(:,33),'r')%,t,y(:,32),'m')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,4); hold on,plot(t,y(:,31)); ylabel('[Ca]SR (mM)');
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel('Time (ms)');
set(gcf,'color','w')

%% Calc and Plot Currents
currents = human_atrial_calcCurrents(t,y,p);

%% Plot all currents
%currents = human_atrial_calcCurrents(t,y,p);
I_Na = currents(:,5);
I_NaL = currents(:,6);
I_NaK = currents(:,7);
I_to = currents(:,8);
I_Kr = currents(:,9);
I_Ks = currents(:,10);
I_Kur = currents(:,11);
I_K1 = currents(:,12);
I_K2P = currents(:,13);
I_SK = currents(:,14);
I_Kp = currents(:,15);
I_ClCa = currents(:,16);
I_CaL = currents(:,17);
I_PMCA = currents(:,18);
I_NCX = currents(:,19);
I_CaBKG = currents(:,20);
I_NaBKG = currents(:,21);
I_ClBKG = currents(:,22);

figure(2),set(gcf,'color','w')
subplot(3,1,1); hold on,plot(t,y(:,39)); ylabel('Em (mV)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(3,1,2); hold on,plot(t,I_to,t,I_Kr,t,I_Ks,t,I_Kur,t,I_K1,t,I_K2P,t,I_SK,t,I_Kp,t,I_NaK,t,I_NaBKG);
ylabel('(A/F)');
set(gca,'box','off','tickdir','out','fontsize',12)
legend('Ito','IKr','IKs','IKur','IK1','IK2P','ISK','IKp','INaK','INaBKG')
subplot(3,1,3); hold on,plot(t,I_Na,t,I_NaL,t,I_CaL,t,I_PMCA,t,I_NCX,t,I_ClCa,t,I_CaBKG,t,I_ClBKG);
ylabel('(A/F)');
set(gca,'box','off','tickdir','out','fontsize',12)
legend('INa','INaL','ICaL','IPMCA','INCX','IClCa','ICaBKG','IClBKG')
xlabel('Time (ms)');

figure(3),set(gcf,'color','w')
subplot(4,1,1); hold on,plot(t,y(:,39)); ylabel('Em (mV)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,2); hold on,plot(t,I_CaL,t,I_to);
ylabel('(A/F)');
set(gca,'box','off','tickdir','out','fontsize',12)
legend('ICaL','Ito')
subplot(4,1,3); hold on,plot(t,I_Kur,t,I_SK,t,I_K2P,t,I_K1);
ylabel('(A/F)');
set(gca,'box','off','tickdir','out','fontsize',12)
legend('IKur','ISK','IK2P','IK1')
subplot(4,1,4); hold on,plot(t,I_NCX,t,I_NaK,t,I_NaL);
ylabel('(A/F)');
set(gca,'box','off','tickdir','out','fontsize',12)
legend('INCX','INaK','INaL')
xlabel('Time (ms)');

figure(4),set(gcf,'color','w')
subplot(3,1,1); hold on,plot(t,y(:,39)); ylabel('Em (mV)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(3,1,2); hold on,plot(t,y(:,38)*1e6);
ylabel('[Ca]i (nM)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(3,1,3); hold on,plot(t,I_SK);
ylabel('ISK (A/F)');
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel('Time (ms)');
