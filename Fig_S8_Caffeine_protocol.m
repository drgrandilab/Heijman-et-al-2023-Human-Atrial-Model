clear
close all
clc

%% Initial Conditions

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nSR
%load yfin_nSR_0p2Hz
load yfin_nSR_0p5Hz
%load yfin_nSR_1Hz
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
prot_index = 4; % (-)
% 0) no stimulation 1) normal pacing 2) ERP 3) V-step 4) caffeine
prot_rate = 0.5; % (Hz)
prot_interval = 100; % (ms) not used
prot_vm = 10; % (mV) not used
prot_par = [prot_index prot_rate prot_interval prot_vm]; % 1 2 3 4

% Cell model parameters
% AF flag
AF_flag = 0; % 0: nSR; 1: cAF
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

%     % BAPTA
%     Ca_clamp_index = 2; % [Ca]c,sl,j clamped to Cai0_Bapta
%     Cai0_Bapta = 5*100e-6;
%     y0(36) = Cai0_Bapta; y0(37) = Cai0_Bapta; y0(38) = Cai0_Bapta;
    
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
duration = 20e3; % (ms)

%% Single Run Simulation
tic
tspan = [0; duration];
options = odeset('RelTol',1e-6,'MaxStep',1,'Stats','off'); 
[t,y] = ode15s(@human_atrial_model,tspan,y0,options,p);
toc

% Final Conditions
yfinal = y(end,:);
%save yfin_XXX yfinal

% Calc and Plot Currents
currents = human_atrial_calcCurrents(t,y,p);
I_NCX = currents(:,19);

figure(1); set(gcf,'color','w')
subplot(4,1,1); hold on,plot(t,y(:,39),'k'); ylabel('Em (mV)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,2); hold on,plot(t,y(:,38),'k'); ylabel('[Ca]i (mM)');
%plot(t,y(:,37),'r',t,y(:,36),'m')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,3); hold on,plot(t,I_NCX,'k'); ylabel('INCX (A/F)'); %
%plot(t,y(:,33),'r')%,t,y(:,32),'m')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,4); hold on,plot(t,y(:,31),'k'); ylabel('[Ca]SR (mM)');
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel('Time (ms)');

%% cAF
% AF flag
AF_flag = 1; % 0: nSR; 1: cAF

load yfin_cAF_0p5Hz
%load yfin_cAF_1Hz
y0 = yfinal;

cell_par = [AF_flag ISO_flag CCh ACh K2P_cond SK_cond SK_shift Ca_clamp_index Na_clamp_index]; % 5-13
p = [prot_par cell_par SA_par];

% Single Run Simulation
tic
tspan = [0; duration];
options = odeset('RelTol',1e-6,'MaxStep',1,'Stats','off'); 
[t,y] = ode15s(@human_atrial_model,tspan,y0,options,p);
toc

% Final Conditions
yfinal = y(end,:);
%save yfin_XXX yfinal

% Calc and Plot Currents
currents = human_atrial_calcCurrents(t,y,p);
I_NCX = currents(:,19);

figure(1);set(gcf,'color','w')
subplot(4,1,1); hold on,plot(t,y(:,39),'r'); ylabel('Em (mV)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,2); hold on,plot(t,y(:,38),'r'); ylabel('[Ca]i (mM)');
%plot(t,y(:,37),'r',t,y(:,36),'m')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,3); hold on,plot(t,I_NCX,'r'); ylabel('INCX (A/F)'); %
%plot(t,y(:,33),'r')%,t,y(:,32),'m')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,4); hold on,plot(t,y(:,31),'r'); ylabel('[Ca]SR (mM)');
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel('Time (ms)');

subplot(4,1,1);xlim([6e3 16e3]); title('Caffeine administration at t = 10 s')
subplot(4,1,2);xlim([6e3 16e3]); legend('nSR','cAF')
subplot(4,1,3);xlim([6e3 16e3]);
subplot(4,1,4);xlim([6e3 16e3]);
