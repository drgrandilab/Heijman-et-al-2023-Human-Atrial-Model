% This file plots the results of linear regression analysis in nSR vs cAF.

close all
clear
clc

%% Select outputs (plot)
% 1) dVm_max 2) Vm_max 3) -Vm_min 4) AP_amp 5) APD90
% 6) APD70 7) APD50 8) APD30 9) Ca_max 10) Ca_min
% 11) CaT_amp 12) CaT_rise 13) CaT_decay_50 14) CaT_decay_63 15) Na_min
% 16) CaSR_max 17) CaSR_min 18) B3_serca 19) B3_ncx 20) B3_pmca
% 21) ERP

%output_6 = [5 7 3 11 15 21]; % APD90 APD50 RMP CaTamp Na ERP
%output_6 = [5 7 11 21]; % APD90 APD50 CaTamp ERP

output_6 = [5 7 3 11]; % APD90 APD50 RMP CaTamp

%output_6 = [5 21]; % APD90 ERP
%output_6 = [3 7]; % RMP APD50
%output_6 = [11 15]; % CaTamp Na

if length(output_6)==1 || length(output_6)==2
    sp1 = 2; sp2 = 1;
end
if length(output_6)==3 || length(output_6)==4
    sp1 = 2; sp2 = 2;
end
if length(output_6)>4
    sp1 = 2; sp2 = 3;
end

%% Load parameters
% load matrix all_parameters (columns: N parameters, rows: N trials)
load perturbation_matrix_1000_0p1 % sigma 0.1
% all_parameters
% parameter_names

% SA_par: 1) GNa 2) GNaL 3) GNaB 4) INaK 5) Gtof 6) GKr 7) GKs
%         8) GKur 9) GK1 10) GKp 11) GK2P 12) GKACh 13) GSK 14) KdSK
%         15) GClCa 16) GClB 17) GCaL (PCa, PK, PNa) 18) GCaB 19) IPMCA
%         20) INCX 21) SERCA 22) SR_release 23) SR_leak

% Remove GK,ACh from the analysis
all_parameters = all_parameters(:,[1:11 13:end]);
parameter_names = parameter_names([1:11 13:end]);

[N_trials N_pars] = size(all_parameters);

%% Load outputs
freq = 1;
SK_index = 1;

% Baseline
if freq == 1 && SK_index == 1
    % 1 Hz
    load output_matrix_nSR_1Hz_500s
    load output_ERP_nSR_1Hz_500s
    %all_outputs_nSR_0 = [all_outputs,output_ERP];
    all_outputs_nSR_1 = [all_outputs,output_ERP];

    load output_matrix_cAF_1Hz_500s
    load output_ERP_cAF_1Hz_500s
    %all_outputs_cAF_0 = [all_outputs,output_ERP];
    all_outputs_cAF_1 = [all_outputs,output_ERP];
end

%% Output properties
output_names = [output_names,output_name];
output_units = [output_units,output_unit];
N_outputs = length(output_names);

%% nSR
all_outputs = all_outputs_nSR_1;

color = [0 0 0]; % BLACK

[N_trials N_pars] = size(all_parameters);

% Check basic properties
% % AP amp > 10 mV
% crit_1 = find(all_outputs(:,4)>10);
% AP amp > 10 mV & ERP > 25 ms
crit_1 = find( (all_outputs(:,4)>10).*(all_outputs(:,21)>25) > 0 );

good_count = length(crit_1);

all_good_parameters = all_parameters(crit_1,:);
all_good_outputs = all_outputs(crit_1,:);

% Define X and Y
X = log(all_good_parameters);
Y = log(all_good_outputs);

% Call the PLS routine - nipals algorithm (2003)
[T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,R2x,R2y] = PLS_nipals(X,Y,rank(X));

% figure
% set(gcf,'color','w','Position',[50,100,1500,750])
% for dex4 = 1:length(output_6),
%     output_dex4 = output_6(dex4);
%     subplot(2,3,dex4)
%     bar(Bpls(:,output_dex4),'FaceColor',color)
%     title(output_names(output_dex4))
%     set(gca,'box','off','tickdir','out','fontsize',10)
%     set(gca,'XTick',1:N_pars)
%     set(gca,'XTickLabel',parameter_names)
%     set(gca,'XLim',[0 N_pars+1])
%     rotateXLabels( gca(), 90)
% end

Bpls_nSR_1 = Bpls;

%% cAF
all_outputs = all_outputs_cAF_1;

color = [1 0 0]; % RED

[N_trials N_pars] = size(all_parameters);

% Check basic properties
% % AP amp > 10 mV
% crit_1 = find(all_outputs(:,4)>10);
% AP amp > 10 mV & ERP > 25 ms
crit_1 = find( (all_outputs(:,4)>10).*(all_outputs(:,21)>25) > 0 );

good_count = length(crit_1);

all_good_parameters = all_parameters(crit_1,:);
all_good_outputs = all_outputs(crit_1,:);

% Define X and Y
X = log(all_good_parameters);
Y = log(all_good_outputs);

% Call the PLS routine - nipals algorithm (2003)
[T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,R2x,R2y] = PLS_nipals(X,Y,rank(X));

% figure
% set(gcf,'color','w','Position',[50,100,1500,750])
% for dex4 = 1:length(output_6),
%     output_dex4 = output_6(dex4);
%     subplot(2,3,dex4)
%     bar(Bpls(:,output_dex4),'FaceColor',color)
%     title(output_names(output_dex4))
%     set(gca,'box','off','tickdir','out','fontsize',10)
%     set(gca,'XTick',1:N_pars)
%     set(gca,'XTickLabel',parameter_names)
%     set(gca,'XLim',[0 N_pars+1])
%     rotateXLabels( gca(), 90)
% end

Bpls_cAF_1 = Bpls;

%% Multiple plot
color_nSR_1 = [0 0 0]; % BLACK
color_nSR_1_e = [0 0 0]; % BLACK

color_cAF_1 = [1 0 0]; % RED
color_cAF_1_e = [1 0 0]; % RED

%% normal SK, nSR vs. cAF 
figure
set(gcf,'color','w','Position',[50,100,1500,750])
for dex4 = 1:length(output_6)
    output_dex4 = output_6(dex4);
    subplot(sp1,sp2,dex4)
    
    Bpls = [Bpls_nSR_1(:,output_dex4) Bpls_cAF_1(:,output_dex4)];
    hb = bar(Bpls);%,'FaceColor',color)
    set(hb(1), 'FaceColor',color_nSR_1, 'EdgeColor',color_nSR_1_e)
    set(hb(2), 'FaceColor',color_cAF_1, 'EdgeColor',color_cAF_1_e)

    title(output_names(output_dex4))
    set(gca,'box','off','tickdir','out','fontsize',10)
    set(gca,'XTick',1:N_pars)
    set(gca,'XTickLabel',parameter_names)
    set(gca,'XLim',[0 N_pars+1])
    rotateXLabels( gca(), 90)
end
