% AF flag
AF_flag = 1; % 0: nSR; 1: cAF
colorarray = [1 0 0];

if SK_cond == 0
    if AF_flag == 0
        % nSR
        load yfin_nSR_SKblock_0p2Hz
        ICs_matrix(1,:) = yfinal;
        load yfin_nSR_SKblock_0p5Hz
        ICs_matrix(2,:) = yfinal;
        load yfin_nSR_SKblock_1Hz
        ICs_matrix(3,:) = yfinal;
        load yfin_nSR_SKblock_2Hz
        ICs_matrix(4,:) = yfinal;
        load yfin_nSR_SKblock_3Hz
        ICs_matrix(5,:) = yfinal;
        load yfin_nSR_SKblock_4Hz
        ICs_matrix(6,:) = yfinal;
        load yfin_nSR_SKblock_5Hz
        ICs_matrix(7,:) = yfinal;
        load yfin_nSR_SKblock_6Hz
        ICs_matrix(8,:) = yfinal;
    else
        % cAF
        load yfin_cAF_SKblock_0p2Hz
        ICs_matrix(1,:) = yfinal;
        load yfin_cAF_SKblock_0p5Hz
        ICs_matrix(2,:) = yfinal;
        load yfin_cAF_SKblock_1Hz
        ICs_matrix(3,:) = yfinal;
        load yfin_cAF_SKblock_2Hz
        ICs_matrix(4,:) = yfinal;
        load yfin_cAF_SKblock_3Hz
        ICs_matrix(5,:) = yfinal;
        load yfin_cAF_SKblock_4Hz
        ICs_matrix(6,:) = yfinal;
        load yfin_cAF_SKblock_5Hz
        ICs_matrix(7,:) = yfinal;
        load yfin_cAF_SKblock_6Hz
        ICs_matrix(8,:) = yfinal;
    end    
end

y0_matrix = ICs_matrix;

cell_par = [AF_flag ISO_flag CCh ACh K2P_cond SK_cond SK_shift Ca_clamp_index Na_clamp_index]; % 5-13

% Parameter array for passing nondefault conditions
p = [prot_par cell_par SA_par];

HR = FREQ; % Hz
outputs = zeros(length(HR),20);

tic
for ii = 1:length(HR)%-2
    prot_rate = HR(ii) %#ok<NOPTS> % (Hz)
    prot_par = [prot_index prot_rate prot_interval prot_vm]; % 1 2 3 4
    p = [prot_par cell_par SA_par];
    
    tspan = [0; simDuration(ii)];
    
    y0 = y0_matrix(ii,:);
    [t,y] = ode15s(@human_atrial_model,tspan,y0,options,p);
    currents = human_atrial_calcCurrents(t,y,p);
    
    figure(11),subplot(2,4,ii),plot(t,y(:,39),'Color',colorarray)
    figure(12),subplot(2,4,ii),plot(t,y(:,38)*1e6,'Color',colorarray)
    
    if prot_rate == 0.5
        figure(10),hold on,subplot(2,4,1),plot(t,y(:,39),'Color',colorarray)
        figure(10),hold on,subplot(2,4,5),plot(t,y(:,38)*1e6,'Color',colorarray)
    elseif prot_rate == 1
        figure(10),hold on,subplot(2,4,2),plot(t,y(:,39),'Color',colorarray)
        figure(10),hold on,subplot(2,4,6),plot(t,y(:,38)*1e6,'Color',colorarray)
    elseif prot_rate == 2
        figure(10),hold on,subplot(2,4,3),plot(t,y(:,39),'Color',colorarray)
        figure(10),hold on,subplot(2,4,7),plot(t,y(:,38)*1e6,'Color',colorarray)
    elseif prot_rate == 3
        figure(10),hold on,subplot(2,4,4),plot(t,y(:,39),'Color',colorarray)
        figure(10),hold on,subplot(2,4,8),plot(t,y(:,38)*1e6,'Color',colorarray)
    end
    
    time = t; % (ms)
    Vm = y(:,39); % (mV)
    Ca = y(:,38); % (mM)
    CaSR = y(:,31); % mM
    Na = y(:,34); % (mM)
    dVm = currents(:,1); % (mV/ms)
    %dVm_delta = (Vm(2:end)-Vm(1:end-1))./(time(2:end)-time(1:end-1));
    %dVm = [dVm_delta; dVm_delta(end)];
    Jserca  = currents(:,2);
    Jncx = currents(:,3);
    Jpmca = currents(:,4);
    period = 1000/prot_rate; % ms
    AP_index = 2; % w/ 1 first AP, otherwise last AP

    % 1) dVm_max 2) Vm_max 3) -Vm_min 4) AP_amp 5) APD90 6) APD70 7) APD50 8) APD30
    % 9) Ca_max 10) Ca_min 11) CaT_amp 12) CaT_rise 13) CaT_decay_50 14) CaT_decay_63
    % 15) Na_min 16) CaSR_max 17) CaSR_min 18) B3_serca 19) B3_ncx 20) B3_pmca
    outputs(ii,:) = function_beat_analysis_B3_alt_2019(time,Vm,Ca,CaSR,Na,dVm,Jserca,Jncx,Jpmca,period,AP_index,Ca_clamp_index);
end
toc

figure(11),set(gcf,'color','w')

subplot(2,4,1),set(gca,'box','off','tickdir','out','fontsize',12)
title('0.2 Hz'),ylabel('Em (mV)'),xlabel('Time (ms)')
subplot(2,4,2),set(gca,'box','off','tickdir','out','fontsize',12)
title('0.5 Hz'),ylabel('Em (mV)'),xlabel('Time (ms)')
subplot(2,4,3),set(gca,'box','off','tickdir','out','fontsize',12)
title('1 Hz'),ylabel('Em (mV)'),xlabel('Time (ms)')
subplot(2,4,4),set(gca,'box','off','tickdir','out','fontsize',12)
title('2 Hz'),ylabel('Em (mV)'),xlabel('Time (ms)')
subplot(2,4,5),set(gca,'box','off','tickdir','out','fontsize',12)
title('3 Hz'),ylabel('Em (mV)'),xlabel('Time (ms)')
subplot(2,4,6),set(gca,'box','off','tickdir','out','fontsize',12)
title('4 Hz'),ylabel('Em (mV)'),xlabel('Time (ms)')
subplot(2,4,7),set(gca,'box','off','tickdir','out','fontsize',12)
title('5 Hz'),ylabel('Em (mV)'),xlabel('Time (ms)')
subplot(2,4,8),set(gca,'box','off','tickdir','out','fontsize',12)
title('6 Hz'),ylabel('Em (mV)'),xlabel('Time (ms)')

figure(12),set(gcf,'color','w')

subplot(2,4,1),set(gca,'box','off','tickdir','out','fontsize',12)
title('0.2 Hz'),ylabel('[Ca]i (nM)'),xlabel('Time (ms)')
subplot(2,4,2),set(gca,'box','off','tickdir','out','fontsize',12)
title('0.5 Hz'),ylabel('[Ca]i (nM)'),xlabel('Time (ms)')
subplot(2,4,3),set(gca,'box','off','tickdir','out','fontsize',12)
title('1 Hz'),ylabel('[Ca]i (nM)'),xlabel('Time (ms)')
subplot(2,4,4),set(gca,'box','off','tickdir','out','fontsize',12)
title('2 Hz'),ylabel('[Ca]i (nM)'),xlabel('Time (ms)')
subplot(2,4,5),set(gca,'box','off','tickdir','out','fontsize',12)
title('3 Hz'),ylabel('[Ca]i (nM)'),xlabel('Time (ms)')
subplot(2,4,6),set(gca,'box','off','tickdir','out','fontsize',12)
title('4 Hz'),ylabel('[Ca]i (nM)'),xlabel('Time (ms)')
subplot(2,4,7),set(gca,'box','off','tickdir','out','fontsize',12)
title('5 Hz'),ylabel('[Ca]i (nM)'),xlabel('Time (ms)')
subplot(2,4,8),set(gca,'box','off','tickdir','out','fontsize',12)
title('6 Hz'),ylabel('[Ca]i (nM)'),xlabel('Time (ms)')

%% Figure 3
figure(3),set(gcf,'color','w')

subplot(4,5,1),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,1),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('UV (mV/ms)')%, xlabel('HR (Hz)')
subplot(4,5,2),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,2),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('AP peak (mV)')%, xlabel('HR (Hz)')
subplot(4,5,3),hold on,plot(HR(HR_low:HR_up),-outputs(HR_low:HR_up,3),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('Em rest (mV)')%, xlabel('HR (Hz)')
subplot(4,5,4),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,4),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('AP amp (mV)')%, xlabel('HR (Hz)')
subplot(4,5,5),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,5),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('APD90 (ms)')%, xlabel('HR (Hz)')

subplot(4,5,6),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,6),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('APD70 (ms)')%, xlabel('HR (Hz)')
subplot(4,5,7),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,7),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('APD50 (ms)')%, xlabel('HR (Hz)')
subplot(4,5,8),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,8),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('APD30 (ms)')%, xlabel('HR (Hz)')
subplot(4,5,9),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,9)*1e6,'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('syst [Ca] (nM)')%, xlabel('HR (Hz)')
subplot(4,5,10),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,10)*1e6,'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('diast [Ca] (nM)')%, xlabel('HR (Hz)')

subplot(4,5,11),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,11)*1e6,'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('CaT amp (nM)')%, xlabel('HR (Hz)')
subplot(4,5,12),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,12),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('CaT ttp (ms)')%, xlabel('HR (Hz)')
subplot(4,5,13),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,13),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('CaT t50 (ms)')%, xlabel('HR (Hz)')
subplot(4,5,14),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,14),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('CaT tau (ms)')%, xlabel('HR (Hz)')
subplot(4,5,15),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,15),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('diast [Na] (mM)')%, xlabel('HR (Hz)')

subplot(4,5,16),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,16),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('CaSR max (mM)'), xlabel('HR (Hz)')
subplot(4,5,17),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,17),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('CaSR min (mM)'), xlabel('HR (Hz)')
subplot(4,5,18),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,18),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('B3 serca (%)'), xlabel('HR (Hz)')
subplot(4,5,19),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,19),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('B3 ncx (%)'), xlabel('HR (Hz)')
subplot(4,5,20),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,20),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('B3 pmca (%)'), xlabel('HR (Hz)')

%% Figure 4
figure(4),set(gcf,'color','w')

subplot(2,3,1),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,5),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('APD90 (ms)')%, xlabel('HR (Hz)')
subplot(2,3,2),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,7),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('APD50 (ms)')%, xlabel('HR (Hz)')
subplot(2,3,3),hold on,plot(HR(HR_low:HR_up),-outputs(HR_low:HR_up,3),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('Em rest (mV)')%, xlabel('HR (Hz)')
subplot(2,3,4),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,11)*1e6,'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('CaT amp (nM)'), xlabel('HR (Hz)')
subplot(2,3,5),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,10)*1e6,'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('CaT min (nM)'), xlabel('HR (Hz)')
subplot(2,3,6),hold on,plot(HR(HR_low:HR_up),outputs(HR_low:HR_up,15),'Color',colorarray)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('diast [Na] (mM)'), xlabel('HR (Hz)')
