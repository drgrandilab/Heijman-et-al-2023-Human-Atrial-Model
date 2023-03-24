function output = human_atrial_model(t,y,p,runType)

%% State variables
% 1       2       3       4       5       6       7       8      9       10  
% m       h       j       d       f       fcaBj   fcaBsl  xkur   ykur    xtof  
% 11      12      13      14      15      16      17      18     19      20   
% ytof    xkr     xks     RyRr    RyRo    RyRi    NaBj    NaBsl  TnCL    TnCHc
% 21      22      23      24      25      26      27      28     29      30
% TnCHm   CaM     Myoc    Myom    SRB     SLLj    SLLsl   SLHj   SLHsl   Csqnb
% 31      32      33      34      35      36      37      38     39      40
% Ca_sr   Naj     Nasl    Nai     Ki      Caj     Casl    Cai    Vm      m(late)
% 41      42
% h(late) xkp2

ydot = zeros(size(y));

%% Cell Parameters
%epi = 1; % EPI or ENDO?
RA = 0; % Right ATRIUM
AF = p(5); % AF
ISO = p(6); % ISO (boolean)
CCh = p(7); % [uM]
ACh = p(8); % [uM]
K2P_cond = p(9);
SK_cond = p(10);
SK_shift = p(11);
Ca_clamp = p(12);
Na_clamp = p(13);

% Scaling factors for sensitivity analysis
SA_par = p(14:end);
if length(SA_par) < 33
    %SA_par = [SA_par 1 1];
    SA_par = [SA_par 1 1 1 1 1 1 1 1 1 1];
end

% Optimization
%parameters_min = ones(1,33);
parameters_min = [0.999248086266895,0.990815925018498,0.996228305841765,0.900939968026087,1.93500438838323,0.822474143493678,0.401661576953456,0.559778787899607,1.08447907955945,0.492465534668169,1.72973882387353,1,0.459661014874920,0.794287503008617,0.505257758179443,1.78778706282580,1.47482899087547,0.998031797040116,0.571769611427024,1.21545276029896,0.894199394767172,0.967002482576837,0.725942111546702,1.89170366530244,1.05312785823022,2.91948962912352,0.986167083837026,0.290250017745365,0.194539767775644,1.48233908592898,0.906532786812203,0.669969161314351,1.46825825179351];
SA_par = SA_par.*parameters_min;

% Caffeine protocol
CaffeineFlag = 0;
if p(1) == 4 && t > 10e3
    CaffeineFlag = 1;
end

%% Model Parameters
% Constants
R = 8314;       % [J/kmol*K]  
Frdy = 96485;   % [C/mol]  
Temp = 310;     % [K]
FoRT = Frdy/R/Temp;
Qpow = (Temp-310)/10;

%% Geometry
% Capacitance
Acell = 11e3; % [um^2]
Cmem = Acell*1e-14; % [F] 110 pF membrane capacitance
%Cmem = 1.1e-10;   % [F] membrane capacitance 1.3810e-10 in ventricles

% Fractional currents in compartments
Fjunc = 0.11; Fsl = 1-Fjunc;
Fjunc_CaL = 0.9; Fsl_CaL = 1-Fjunc_CaL;

% Cell dimensions and volume
cellLength = 100;     % cell length [um]113;%100
cellRadius = 10.25;   % cell radius [um]12;%10.25
Vcell = pi*cellRadius^2*cellLength*1e-15; % [L]
Vmyo = 0.65*Vcell;
Vsr = 0.035*Vcell;
Vsl = 0.02*Vcell;
Vjunc = 0.0539*0.01*Vcell; 

% Diffusion rates
J_ca_juncsl = 1/1.2134e12; % [L/msec]           8.2413e-13
J_ca_slmyo = 1/2.68510e11; % [L/msec]           3.2743e-12
J_na_juncsl = 1/(1.6382e12/3*100); % [L/msec]   1.8313e-14
J_na_slmyo = 1/(1.8308e10/3*100);  % [L/msec]   1.6386e-12

%% Ion concentrations
% Fixed concentrations     
Cli = 15;   % Intracellular Cl  [mM]
Clo = 150;  % Extracellular Cl  [mM]
Ko = 5.4;   % Extracellular K   [mM]
Nao = 140;  % Extracellular Na  [mM]
Cao = 1.8;  % Extracellular Ca  [mM]
    %Ko = 4;   % Extracellular K   [mM]
    %Nao = 147;  % Extracellular Na  [mM]
    %Cao = 2;  % Extracellular Ca  [mM]
Mgi = 1;    % Intracellular Mg  [mM]

% Nernst Potentials
ena_junc = (1/FoRT)*log(Nao/y(32));     % [mV]
ena_sl = (1/FoRT)*log(Nao/y(33));       % [mV]
ek = (1/FoRT)*log(Ko/y(35));	        % [mV]
eca_junc = (1/FoRT/2)*log(Cao/y(36));   % [mV]
eca_sl = (1/FoRT/2)*log(Cao/y(37));     % [mV]
ecl = (1/FoRT)*log(Cli/Clo);            % [mV]

%% I_Na
% % Grandi et al, 2011
% GNa = SA_par(1)*23*(1-0.1*AF);  % [mS/uF] Grandi
% 
% mss = 1 / ((1 + exp( -(56.86 + y(39)) / 9.03 ))^2);
% taum = 0.1292 * exp(-((y(39)+45.79)/15.54)^2) + 0.06487 * exp(-((y(39)-4.823)/51.12)^2);                 
%  
% ah = (y(39) >= -40) * (0)... 
%    + (y(39) < -40) * (0.057 * exp( -(y(39) + 80) / 6.8 )); 
% bh = (y(39) >= -40) * (0.77 / (0.13*(1 + exp( -(y(39) + 10.66) / 11.1 )))) ...
%    + (y(39) < -40) * ((2.7 * exp( 0.079 * y(39)) + 3.1*10^5 * exp(0.3485 * y(39)))); 
% tauh = 1 / (ah + bh); 
% hss = 1 / ((1 + exp( (y(39) + 71.55)/7.43 ))^2);
%  
% aj = (y(39) >= -40) * (0) ...
%     +(y(39) < -40) * (((-2.5428 * 10^4*exp(0.2444*y(39)) - 6.948*10^-6 * exp(-0.04391*y(39))) * (y(39) + 37.78)) / ...
%                      (1 + exp( 0.311 * (y(39) + 79.23) )));
% bj = (y(39) >= -40) * ((0.6 * exp( 0.057 * y(39))) / (1 + exp( -0.1 * (y(39) + 32) ))) ...
%    + (y(39) < -40) * ((0.02424 * exp( -0.01052 * y(39) )) / (1 + exp( -0.1378 * (y(39) + 40.14) ))); 
% tauj = 1 / (aj + bj);
% jss = 1 / ((1 + exp( (y(39) + 71.55)/7.43 ))^2);         
% 
% ydot(1) = (mss - y(1)) / taum;
% ydot(2) = (hss - y(2)) / tauh;
% ydot(3) = (jss - y(3)) / tauj;

% Courtemanche et al, 1998
%GNa = SA_par(1)*7.8*(1-0.1*AF);  % [mS/uF] Courtemanche
GNa = SA_par(1)*9*(1-0.1*AF);  % [mS/uF] Courtemanche (adjusted for dV/dt max) MOD1

am = (y(39) == -47.13) * (3.2)... 
   + (y(39) ~= -47.13) * (0.32 * (y(39)+47.13) / (1 - exp( -0.1*(y(39)+47.13)))); 
bm = 0.08*exp(-y(39)/11);

ah = (y(39) >= -40) * (0)... 
   + (y(39) < -40) * (0.135 * exp( -(y(39) + 80) / 6.8 )); 
bh = (y(39) >= -40) * (1 / (0.13*(1 + exp( -(y(39) + 10.66) / 11.1 )))) ...
   + (y(39) < -40) * ((3.56 * exp( 0.079 * y(39)) + 3.1*10^5 * exp(0.35 * y(39)))); 

aj = (y(39) >= -40) * (0) ...
    +(y(39) < -40) * (( (-127140*exp(0.2444*y(39)) - 3.474*10^-5*exp(-0.04391*y(39))) * (y(39) + 37.78)) / ...
                     (1 + exp( 0.311 * (y(39) + 79.23) ) ));
bj = (y(39) >= -40) * ((0.3 * exp(-2.535*10^-7*y(39))) / (1 + exp( -0.1 * (y(39) + 32) ))) ...
   + (y(39) < -40) * ((0.1212 * exp( -0.01052 * y(39) )) / (1 + exp( -0.1378 * (y(39) + 40.14) )));

ydot(1) = am*(1-y(1))-bm*y(1);
ydot(2) = ah*(1-y(2))-bh*y(2);
ydot(3) = aj*(1-y(3))-bj*y(3);
    
I_Na_junc = Fjunc*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_junc);
I_Na_sl = Fsl*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_sl);
I_Na = I_Na_junc + I_Na_sl; % #ok<NASGU>

%% Late I_Na
%GNaL = SA_par(2)*0.0025*AF; % [mS/uF]
GNaL = SA_par(2)*(0.0005+0.0020*AF); % [mS/uF] (adjusted) MOD1

aml = 0.32*(y(39)+47.13)/(1-exp(-0.1*(y(39)+47.13)));
bml = 0.08*exp(-y(39)/11);
hlinf = 1/(1+exp((y(39)+91)/6.1));
tauhl = 600;
ydot(40) = aml*(1-y(40))-bml*y(40);
ydot(41) = (hlinf-y(41))/tauhl;

I_NaL_junc = Fjunc*GNaL*y(40)^3*y(41)*(y(39)-ena_junc);
I_NaL_sl = Fsl*GNaL*y(40)^3*y(41)*(y(39)-ena_sl);
I_NaL = I_NaL_junc + I_NaL_sl; 

%% I_nabk: Na Background Current
GNaB = SA_par(3)*0.597e-3; % [mS/uF] 

I_nabk_junc = Fjunc*GNaB*(y(39)-ena_junc);
I_nabk_sl = Fsl*GNaB*(y(39)-ena_sl);
I_nabk = I_nabk_junc + I_nabk_sl; 

%% I_nak: Na/K Pump Current
IbarNaK = SA_par(4)*1.26; % [uA/uF]
KmNaip = 11*(1-0.25*ISO); % [mM]
KmKo = 1.5; % [mM]
%Q10NaK = 1.63;  
%Q10KmNai = 1.39;

sigma = (exp(Nao/67.3)-1)/7;
fnak = 1/(1+0.1245*exp(-0.1*y(39)*FoRT)+0.0365*sigma*exp(-y(39)*FoRT));

I_nak_junc = Fjunc*IbarNaK*fnak*Ko / (1+(KmNaip/y(32))^4) /(Ko+KmKo);
I_nak_sl = Fsl*IbarNaK*fnak*Ko / (1+(KmNaip/y(33))^4) /(Ko+KmKo);
I_nak = I_nak_junc + I_nak_sl;

%% I_to: Transient Outward K Current
%GtoFast = SA_par(5)*(1-0.7*AF)*0.165;
GtoFast = SA_par(5)*(1-0.8*AF*(1-RA)-0.45*AF*RA)*0.165; % Updated from Caballero 2010 (-45% RA, -80% LA)

xtoss = 1/( 1 + exp( -(y(39)+1)/11 ) );
tauxtof = 3.5*exp(-((y(39)/30)^2))+1.5;
ytoss = 1/( 1 + exp( (y(39)+40.5)/11.5) ) ;
tauytof = 25.635*exp(-(((y(39)+52.45)/15.8827)^2))+24.14; %14.14

ydot(10) = (xtoss-y(10))/tauxtof;
ydot(11) = (ytoss-y(11))/tauytof;

I_tof = GtoFast*y(10)*y(11)*(y(39)-ek);
I_to = I_tof;

%% I_kr: Rapidly Activating K Current
gkr = SA_par(6)*0.035*sqrt(Ko/5.4);

xrss = 1/(1+exp(-(y(39)+10)/5));
tauxr = 550/(1+exp((-22-y(39))/9))*6/(1+exp((y(39)-(-11))/9))+230/(1+exp((y(39)-(-40))/20));
ydot(12) = (xrss-y(12))/tauxr;
rkr = 1/(1+exp((y(39)+74)/24));
I_kr = gkr*y(12)*rkr*(y(39)-ek);

%% I_ks: Slowly Activating K Current
gks_junc = SA_par(7)*(1+1*AF+2*ISO)*0.0035;
gks_sl = SA_par(7)*(1+1*AF+2*ISO)*0.0035;
pNaK = 0.01833; 
eks = (1/FoRT)*log((Ko+pNaK*Nao)/(y(35)+pNaK*y(34)));

xsss = 1/(1+exp(-(y(39)+40*ISO+3.8)/14.25));
tauxs = 990.1/(1+exp(-(y(39)+40*ISO+2.436)/14.12));
ydot(13) = (xsss-y(13))/tauxs;

I_ks_junc = Fjunc*gks_junc*y(13)^2*(y(39)-eks);
I_ks_sl = Fsl*gks_sl*y(13)^2*(y(39)-eks);                                                                                                                                   
I_ks = I_ks_junc + I_ks_sl;

%% I_kur: Ultra-Rapid Delayed Rectifier Outward K Current
%Gkur = SA_par(8)*(1-0.5*AF)*(1+2*ISO)*0.045*(1+0.2*RA);
Gkur = SA_par(8)*(1+2*ISO)*0.045*(1+0.2*RA*(1-AF)-0.5*AF); % Updated from Caballero 2010

xkurss = 1/( 1 + exp( (y(39)+6)/-8.6) );
tauxkur = 9/(1+exp((y(39)+5)/12.0))+0.5;
ykurss = 1/( 1 + exp( (y(39)+7.5)/10) );
tauykur = 590/(1+exp((y(39)+60)/10))+3050;

ydot(8) = (xkurss-y(8))/tauxkur;
ydot(9) = (ykurss-y(9))/tauykur;
I_kur = Gkur*y(8)*y(9)*(y(39)-ek);

%% I_ki: Time-Independent K Current
%Gki = SA_par(9)*(1+1*AF)*(2.1*0.0525)*sqrt(Ko/5.4);
Gki = SA_par(9)*(1+SA_par(33)*1*AF)*(2.1*0.0525)*sqrt(Ko/5.4); % MOD2

% Na-dependence (Voigt-Heijman et al 2013)
aki_sl = (0.1+0.9/(1+(y(33)/7)^2))*1/(1+exp(0.2385*(y(39)-ek-59.215)));
aki_j = (0.1+0.9/(1+(y(32)/7)^2))*1/(1+exp(0.2385*(y(39)-ek-59.215)));
% Na-dependence (Schmidt et al 2015)
% aki_sl = 1/(1+(y(33)/7)^2);
% aki_j = 1/(1+(y(32)/7)^2);

bki = (0.49124*exp(0.08032*(y(39)+5.476-ek)) + exp(0.06175*(y(39)-ek-594.31))) /(1 + exp(-0.5143*(y(39)-ek+4.753)));
kiss_sl = aki_sl/(aki_sl+bki);
kiss_j = aki_j/(aki_j+bki);

I_ki_sl = Fsl*Gki*kiss_sl*(y(39)-ek);
I_ki_j = Fjunc*Gki*kiss_j*(y(39)-ek);
I_ki = I_ki_j + I_ki_sl;

%% I_kp: Plateau K current
gkp = SA_par(10)*0.002;

kp_kp = 1/(1+exp(7.488-y(39)/5.98));
I_kp = gkp*kp_kp*(y(39)-ek);

%% I_k2p: K2P3.1 K current (from Schmidt et al 2015) 
%gk2p = SA_par(11)*K2P_cond*0.0050*(1+(0.0145/0.005-1)*AF); % 2.9-fold increase in AF
gk2p = SA_par(11)*K2P_cond*0.0050*(1+(SA_par(24)*0.0145/0.005-1)*AF);

xk2pss = 0.2+0.8/(1+exp(-(y(39)-10+15*AF)/14));
tauxk2p = 2+40/(1+exp((y(39)+25)^2/80));
ydot(42) = (xk2pss-y(42))/tauxk2p;

I_k2p = gk2p*y(42)*(y(39)-ek);

%% I_k,ach: Muscarinic-Receptor-Activated K Current
% Grandi (ACh)
% I_KAch = 1/(1+(0.03/Ach)^2.1)*(0.08+0.4./(1+exp((y(39)+91)/12))).*(y(39)-ek);

% Schmidt (CCh)
% fkach = 1.5/(1+(9/y(33))^4); % Na-dependence (NaSL)
% rkach = 0.055+0.4/(1+exp((y(39)-ek+9.53)/17.18));
% if AF == 0
%     Gkach = 0.10*(1+fkach)/0.9411765;
% else
%     Gkach = 0.05/0.9411765;
% end
% I_kach = Gkach*sqrt(Ko/5.4)*(CCh/(CCh+0.125))*rkach*(y(39)-ek);

if CCh == 0 && ACh == 0 % none
    Gkach_sl = SA_par(12)*0;
    Gkach_j = SA_par(12)*0;
    d_kach = 0;
    r_kach = 0;
elseif CCh > 0 && ACh == 0 % CCh
    % Na-dependence
    fkach_sl = 1.5/(1+(9/y(33))^4); % Na-dependence sl
    fkach_j = 1.5/(1+(9/y(32))^4); % Na-dependence j
    if AF == 0 
        Gkach_sl = SA_par(12)*(RA*1+(1-RA)*0.3)*0.10*(1+fkach_sl)/0.9411765*sqrt(Ko/5.4);
        Gkach_j = SA_par(12)*(RA*1+(1-RA)*0.3)*0.10*(1+fkach_j)/0.9411765*sqrt(Ko/5.4);
    else
        Gkach_sl = SA_par(12)*(RA*0.5+(1-RA)*0.3)*0.10/0.9411765*sqrt(Ko/5.4);
        Gkach_j = SA_par(12)*(RA*0.5+(1-RA)*0.3)*0.10/0.9411765*sqrt(Ko/5.4);
    end
    % V-dependence
    r_kach = 0.055+0.4/(1+exp((y(39)-ek+9.53)/17.18));
    % Dose-dependence
    d_kach = CCh/(CCh+0.125);
elseif CCh == 0 && ACh > 0 % ACh
	% Na-dependence
    fkach_sl = 1.5/(1+(9/y(33))^4); % Na-dependence sl
    fkach_j = 1.5/(1+(9/y(32))^4); % Na-dependence j
    if AF == 0
        Gkach_sl = SA_par(12)*(RA*1+(1-RA)*0.3)*5*0.10*(1+fkach_sl)/0.9411765*sqrt(Ko/5.4);
        Gkach_j = SA_par(12)*(RA*1+(1-RA)*0.3)*5*0.10*(1+fkach_j)/0.9411765*sqrt(Ko/5.4);
    else
        Gkach_sl = SA_par(12)*(RA*0.5+(1-RA)*0.3)*5*0.10/0.9411765*sqrt(Ko/5.4);
        Gkach_j = SA_par(12)*(RA*0.5+(1-RA)*0.3)*5*0.10/0.9411765*sqrt(Ko/5.4);
    end
    % V-dependence
    r_kach = 0.08+0.4/(1+exp((y(39)+91)/12));
    % Dose-dependence
    d_kach = 1/(1+(0.03/ACh)^2.1);
end

I_kach_sl = Fsl*Gkach_sl*d_kach*r_kach*(y(39)-ek);
I_kach_j = Fjunc*Gkach_j*d_kach*r_kach*(y(39)-ek);

I_kach = I_kach_sl + I_kach_j;

%% I_sk: Small-Conductance Ca-Activated K Current
par_sk = [0.0506381114404388,0.273335569451572,2.96381060498817,0.199981221802789,0.279328126521496,-86.9289059836381,0.00636311816933264,5.22915055145375];

%gsk = SA_par(13)*SK_cond*par_sk(1)*(1+(par_sk(8)-1)*AF);
gsk = SA_par(13)*SK_cond*par_sk(1)*(1+(SA_par(25)*par_sk(8)-1)*AF);

kdsk = SA_par(14)*(10^(SK_shift-3.45)); % kdsk = (10^(ISK_shift-3.3)); % (mM)
gsk_ca_junc = 1/(1+ exp((log10(kdsk)-log10(y(36)))/0.3));
gsk_ca_sl = 1/(1+ exp((log10(kdsk)-log10(y(37)))/0.3));

%     % Ca clamp for SK channels
%     gsk_ca = 500e-6;
%     gsk_ca_junc = 1/(1+ exp((log10(kdsk)-log10(gsk_ca))/0.3));
%     gsk_ca_sl = 1/(1+ exp((log10(kdsk)-log10(gsk_ca))/0.3));

%     % Min Ca seen by SKs is 350 nM
%     if y(36) < 350e-6
%         gsk_ca_junc = 1/(1+ exp((log10(kdsk)-log10(350e-6))/0.3));
%     end
%     if y(37) < 350e-6
%         gsk_ca_sl = 1/(1+ exp((log10(kdsk)-log10(350e-6))/0.3));
%     end

gsk_vm = par_sk(2)/(1+exp((y(39)-ek+par_sk(3))*par_sk(4))) + par_sk(5)/(1+exp((-(y(39)-ek+par_sk(6))*par_sk(7))));

I_sk_junc = Fjunc*gsk*gsk_ca_junc*gsk_vm*(y(39)-ek);
I_sk_sl = Fsl*gsk*gsk_ca_sl*gsk_vm*(y(39)-ek);                                                                                                                                
I_sk = I_sk_junc + I_sk_sl;

%% I_ClCa and I_Clbk: Ca-activated and Background Cl Currents
GClCa = SA_par(15)*0.0548;     % [mS/uF]
KdClCa = 100e-3;               % [mM]
%GClB = SA_par(16)*9e-3;        % [mS/uF]
%GClB = SA_par(16)*0.5*9e-3;        % [mS/uF] MOD1
GClB = SA_par(16)*0.75*9e-3;        % [mS/uF] MOD2
GClCFTR = 0; % 4.9e-3*ISO;     % [mS/uF]

I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/y(36))*(y(39)-ecl);
I_ClCa_sl = Fsl*GClCa/(1+KdClCa/y(37))*(y(39)-ecl);
I_ClCa = I_ClCa_junc+I_ClCa_sl;

I_Clbk = GClB*(y(39)-ecl);

I_ClCFTR = GClCFTR*(y(39)-ecl);

%% I_Ca: L-type Ca Current
pNa = SA_par(17)*(1+0.5*ISO)*(1-0.5*AF)*0.75e-8; % [cm/sec]
pCa = SA_par(17)*(1+0.5*ISO)*(1-0.5*AF)*2.7e-4; % [cm/sec]
pK = SA_par(17)*(1+0.5*ISO)*(1-0.5*AF)*1.35e-7; % [cm/sec]
Q10CaL = 1.8;

dss = 1/(1+exp(-(y(39)+3*ISO+9)/6));
taud = dss*(1-exp(-(y(39)+3*ISO+9)/6))/(0.035*(y(39)+3*ISO+9)); 
fss = 1/(1+exp((y(39)+3*ISO+30)/7))+0.2/(1+exp((50-y(39)-3*ISO)/20));
tauf = 1/(0.0197*exp( -(0.0337*(y(39)+3*ISO+25))^2 )+0.02);
ydot(4) = (dss-y(4))/taud;
ydot(5) = (fss-y(5))/tauf;
ydot(6) = 1.7*y(36)*(1-y(6))-1*11.9e-3*y(6); % fCa_junc
ydot(7) = 1.7*y(37)*(1-y(7))-1*11.9e-3*y(7); % fCa_sl
%fcaCaMSL = 0.1/(1+(0.01/y(37)));
%fcaCaj = 0.1/(1+(0.01/y(36)));
fcaCaMSL = 0;
fcaCaj = 0;
ibarca_j = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(36)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibarca_sl = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(37)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibark = pK*(y(39)*Frdy*FoRT)*(0.75*y(35)*exp(y(39)*FoRT)-0.75*Ko) /(exp(y(39)*FoRT)-1);
ibarna_j = pNa*(y(39)*Frdy*FoRT) *(0.75*y(32)*exp(y(39)*FoRT)-0.75*Nao) /(exp(y(39)*FoRT)-1);
ibarna_sl = pNa*(y(39)*Frdy*FoRT) *(0.75*y(33)*exp(y(39)*FoRT)-0.75*Nao) /(exp(y(39)*FoRT)-1);

I_Ca_junc = (Fjunc_CaL*ibarca_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45;
I_Ca_sl = (Fsl_CaL*ibarca_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*0.45;
I_Ca = I_Ca_junc+I_Ca_sl;
I_CaK = (ibark*y(4)*y(5)*(Fjunc_CaL*(fcaCaj+(1-y(6)))+Fsl_CaL*(fcaCaMSL+(1-y(7))))*Q10CaL^Qpow)*0.45;
I_CaNa_junc = (Fjunc_CaL*ibarna_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45;
I_CaNa_sl = (Fsl_CaL*ibarna_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*.45;
I_CaNa = I_CaNa_junc+I_CaNa_sl;
I_Catot = I_Ca + I_CaK + I_CaNa; 

%% I_cabkg: Ca Background Current
%GCaB = SA_par(18)*6.0643e-4;   % [uA/uF]
%GCaB = SA_par(18)*6.0643e-4*(1+SA_par(26)*0*AF)*(1-CaffeineFlag); % [uA/uF]
GCaB = SA_par(18)*6.0643e-4*(1+SA_par(26)*0.5*AF)*(1-CaffeineFlag); % [uA/uF] MOD2

I_cabk_junc = Fjunc*GCaB*(y(39)-eca_junc);
I_cabk_sl = Fsl*GCaB*(y(39)-eca_sl);
I_cabk = I_cabk_junc + I_cabk_sl;

%% I_pca: Sarcolemmal Ca Pump Current
IbarSLCaP = SA_par(19)*0.0471; % [uA/uF]
KmPCa = 0.5e-3;     % [mM] 
Q10SLCaP = 2.35;    % [none]

I_pca_junc = Fjunc*Q10SLCaP^Qpow*IbarSLCaP*y(36)^1.6/(KmPCa^1.6+y(36)^1.6);
I_pca_sl = Fsl*Q10SLCaP^Qpow*IbarSLCaP*y(37)^1.6/(KmPCa^1.6+y(37)^1.6);
I_pca = I_pca_junc + I_pca_sl;

%% I_ncx: Na/Ca Exchanger flux
%IbarNCX = SA_par(20)*(1+0.4*AF)*3.15; % [uA/uF]
IbarNCX = SA_par(20)*(1+SA_par(27)*0.4*AF)*3.15; % [uA/uF]
KmCai = 3.59e-3;    % [mM]
KmCao = 1.3;        % [mM]
KmNai = 12.29;      % [mM]
KmNao = 87.5;       % [mM]
ksat = 0.27;        % [none]  
nu = 0.35;          % [none]
Kdact = 0.384e-3;   % [mM]
Q10NCX = 1.57;      % [none]

Ka_junc = 1/(1+(Kdact/y(36))^2);
Ka_sl = 1/(1+(Kdact/y(37))^2);
s1_junc = exp(nu*y(39)*FoRT)*y(32)^3*Cao;
s1_sl = exp(nu*y(39)*FoRT)*y(33)^3*Cao;
s2_junc = exp((nu-1)*y(39)*FoRT)*Nao^3*y(36);
s3_junc = KmCai*Nao^3*(1+(y(32)/KmNai)^3) + KmNao^3*y(36)*(1+y(36)/KmCai)+KmCao*y(32)^3+y(32)^3*Cao+Nao^3*y(36);
s2_sl = exp((nu-1)*y(39)*FoRT)*Nao^3*y(37);
s3_sl = KmCai*Nao^3*(1+(y(33)/KmNai)^3) + KmNao^3*y(37)*(1+y(37)/KmCai)+KmCao*y(33)^3+y(33)^3*Cao+Nao^3*y(37);

I_ncx_junc = Fjunc*IbarNCX*Q10NCX^Qpow*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1+ksat*exp((nu-1)*y(39)*FoRT));
I_ncx_sl = Fsl*IbarNCX*Q10NCX^Qpow*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1+ksat*exp((nu-1)*y(39)*FoRT));
I_ncx = I_ncx_junc + I_ncx_sl;

%% SR fluxes: Calcium Uptake, Release, and Leak
%Vmax_SRCaP = SA_par(21)*5.3114e-3*(1-SA_par(29)*0*AF)*(1-CaffeineFlag);  % [mM/msec] (286 umol/L cytosol/sec)
%Kmf = SA_par(31)*(2.5-1.25*ISO)*0.246e-3*(1+SA_par(30)*0*AF); % [mM]
Vmax_SRCaP = SA_par(21)*5.3114e-3*(1-SA_par(29)*0.25*AF)*(1-CaffeineFlag); % [mM/msec] MOD2
Kmf = SA_par(31)*(2.5-1.25*ISO)*0.246e-3*(1+SA_par(30)*0.25*AF); % [mM] MOD2
Kmr = 1.7;               % [mM]L cytosol
Q10SRCaP = 2.6;          % [none]
hillSRCaP = 1.787;       % [mM]
ks = SA_par(22)*25;      % [1/ms]      
%koCa = 10+20*AF+10*ISO*(1-AF); % [mM^-2 1/ms]
%koCa = 10*(1+2*AF)*(1+1*ISO*(1-AF)); % [mM^-2 1/ms]
koCa = 10*(1+SA_par(28)*2*AF)*(1+1*ISO*(1-AF))*(1+CaffeineFlag*6.5); % [mM^-2 1/ms]
kom = 0.06;              % [1/ms]
kiCa = 0.5;              % [1/mM/ms]
kim = 0.005;             % [1/ms]
%ec50SR = 0.45*(1+SA_par(32)*0*AF);           % [mM]
ec50SR = 0.45*(1+SA_par(32)*0.25*AF);           % [mM] MOD2
MaxSR = 15;              % [mM]
MinSR = 1;               % [mM]
kleak = SA_par(23)*(1+0.25*AF)*5.348e-6;

kCaSR = MaxSR - (MaxSR-MinSR)/(1+(ec50SR/y(31))^2.5);
koSRCa = koCa/kCaSR;
kiSRCa = kiCa*kCaSR;
RI = 1-y(14)-y(15)-y(16);
ydot(14) = (kim*RI-kiSRCa*y(36)*y(14))-(koSRCa*y(36)^2*y(14)-kom*y(15));   % R
ydot(15) = (koSRCa*y(36)^2*y(14)-kom*y(15))-(kiSRCa*y(36)*y(15)-kim*y(16));% O
ydot(16) = (kiSRCa*y(36)*y(15)-kim*y(16))-(kom*y(16)-koSRCa*y(36)^2*RI);   % I
J_SRCarel = ks*y(15)*(y(31)-y(36)); % [mM/ms]

J_serca = Q10SRCaP^Qpow*Vmax_SRCaP*((y(38)/Kmf)^hillSRCaP-(y(31)/Kmr)^hillSRCaP)...
    /(1+(y(38)/Kmf)^hillSRCaP+(y(31)/Kmr)^hillSRCaP); % [mM/ms]

J_SRleak = kleak*(y(31)-y(36)); % [mM/ms]

%% Sodium and Calcium Buffering
% koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
Bmax_Naj = 7.561;       % [mM] % Na buffering
Bmax_Nasl = 1.65;       % [mM]
koff_na = 1e-3;         % [1/ms]
kon_na = 0.1e-3;        % [1/mM/ms]
Bmax_TnClow = 70e-3;    % [mM]                      % TnC low affinity
koff_tncl = (1+0.5*ISO)*19.6e-3;    % [1/ms] 
kon_tncl = 32.7;        % [1/mM/ms]
Bmax_TnChigh = 140e-3;  % [mM]                      % TnC high affinity 
koff_tnchca = 0.032e-3; % [1/ms] 
kon_tnchca = 2.37;      % [1/mM/ms]
koff_tnchmg = 3.33e-3;  % [1/ms] 
kon_tnchmg = 3e-3;      % [1/mM/ms]
Bmax_CaM = 24e-3;       % [mM] **? about setting to 0 in c-code**   % CaM buffering
koff_cam = 238e-3;      % [1/ms] 
kon_cam = 34;           % [1/mM/ms]
Bmax_myosin = 140e-3;   % [mM]                      % Myosin buffering
koff_myoca = 0.46e-3;   % [1/ms]
kon_myoca = 13.8;       % [1/mM/ms]
koff_myomg = 0.057e-3;  % [1/ms]
kon_myomg = 0.0157;     % [1/mM/ms]
Bmax_SR = 19*.9e-3;     % [mM] (Bers text says 47e-3) 19e-3
koff_sr = 60e-3;        % [1/ms]
kon_sr = 100;           % [1/mM/ms]
Bmax_SLlowsl = 37.4e-3*Vmyo/Vsl;        % [mM]   	% SL buffering
Bmax_SLlowj = 4.6e-3*Vmyo/Vjunc*0.1;    % [mM]
koff_sll = 1300e-3;     % [1/ms]
kon_sll = 100;          % [1/mM/ms]
Bmax_SLhighsl = 13.4e-3*Vmyo/Vsl;       % [mM] 
Bmax_SLhighj = 1.65e-3*Vmyo/Vjunc*0.1;  % [mM]
koff_slh = 30e-3;       % [1/ms]
kon_slh = 100;          % [1/mM/ms]
Bmax_Csqn = 140e-3*Vmyo/Vsr;            % [mM]      % Csqn buffering
koff_csqn = 65;         % [1/ms] 
kon_csqn = 100;         % [1/mM/ms] 

% Junctional and SL Na Buffers
ydot(17) = kon_na*y(32)*(Bmax_Naj-y(17))-koff_na*y(17);        % NaBj      [mM/ms]
ydot(18) = kon_na*y(33)*(Bmax_Nasl-y(18))-koff_na*y(18);       % NaBsl     [mM/ms]

% Cytosolic Ca Buffers
ydot(19) = kon_tncl*y(38)*(Bmax_TnClow-y(19))-koff_tncl*y(19);            % TnCL      [mM/ms]
ydot(20) = kon_tnchca*y(38)*(Bmax_TnChigh-y(20)-y(21))-koff_tnchca*y(20); % TnCHc     [mM/ms]
ydot(21) = kon_tnchmg*Mgi*(Bmax_TnChigh-y(20)-y(21))-koff_tnchmg*y(21);   % TnCHm     [mM/ms]
ydot(22) = kon_cam*y(38)*(Bmax_CaM-y(22))-koff_cam*y(22);                 % CaM       [mM/ms]
ydot(23) = kon_myoca*y(38)*(Bmax_myosin-y(23)-y(24))-koff_myoca*y(23);    % Myosin_ca [mM/ms]
ydot(24) = kon_myomg*Mgi*(Bmax_myosin-y(23)-y(24))-koff_myomg*y(24);      % Myosin_mg [mM/ms]
ydot(25) = kon_sr*y(38)*(Bmax_SR-y(25))-koff_sr*y(25);                    % SRB       [mM/ms]
J_CaB_cytosol = ydot(19)+ydot(20)+ydot(22)+ydot(23)+ydot(25);

% Junctional and SL Ca Buffers
ydot(26) = kon_sll*y(36)*(Bmax_SLlowj-y(26))-koff_sll*y(26);       % SLLj      [mM/ms]
ydot(27) = kon_sll*y(37)*(Bmax_SLlowsl-y(27))-koff_sll*y(27);      % SLLsl     [mM/ms]
ydot(28) = kon_slh*y(36)*(Bmax_SLhighj-y(28))-koff_slh*y(28);      % SLHj      [mM/ms]
ydot(29) = kon_slh*y(37)*(Bmax_SLhighsl-y(29))-koff_slh*y(29);     % SLHsl     [mM/ms]
J_CaB_junction = ydot(26)+ydot(28);
J_CaB_sl = ydot(27)+ydot(29);

%% Ion concentrations
% SR Ca Concentrations
ydot(30) = kon_csqn*y(31)*(Bmax_Csqn-y(30))-koff_csqn*y(30);       % Csqn      [mM/ms]
ydot(31) = J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel)-ydot(30);         % Ca_sr     [mM/ms]

% Sodium Concentrations
I_Na_tot_junc = I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc+I_NaL_junc;   % [uA/uF]
I_Na_tot_sl = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl+I_NaL_sl;   % [uA/uF]

ydot(32) = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(y(33)-y(32))-ydot(17);
ydot(33) = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(y(32)-y(33))...
   +J_na_slmyo/Vsl*(y(34)-y(33))-ydot(18);
ydot(34) = (J_na_slmyo/Vmyo*(y(33)-y(34)))*(1-Na_clamp);   % [mM/msec] 

% Potassium Concentration
I_K_tot = I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp+I_kur+I_kach+I_k2p+I_sk; % [uA/uF] %SVP: added IKur
ydot(35) = 0; % -I_K_tot*Cmem/(Vmyo*Frdy);         % [mM/msec]

% Calcium Concentrations
I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc;                   % [uA/uF]
I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl;            % [uA/uF]
ydot(36) = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(y(37)-y(36))...
    -J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc;  % Ca_j
ydot(37) = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(y(36)-y(37))...
    +J_ca_slmyo/Vsl*(y(38)-y(37))-J_CaB_sl;   % Ca_sl
ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol+J_ca_slmyo/Vmyo*(y(37)-y(38));             % [mM/msec]
if Ca_clamp == 1 || Ca_clamp == 2 % Ca_clamp, BAPTA
    ydot(36) = 0; ydot(37) = 0; ydot(38) = 0; 
elseif Ca_clamp == 3 % EGTA
    ydot(38) = 0; 
end

%% Simulation type
switch p(1)
    case 0          % no stimulation
        I_app = 0;
    case 1          % pace w/ current injection at rate 'rate' (Hz)
		rate = p(2); % Hz
        period = 1000/rate; % ms
        if mod(t,period) <= 5
            I_app = 12.5;
        else
            I_app = 0.0;
        end  
    case 2      % ERP
        rate = p(2); % Hz
        rec_interval = p(3);
        if t <= 5
            I_app = 12.5;
        elseif t > 5 && t <= rec_interval
            I_app = 0.0;
        elseif t > rec_interval && t <= rec_interval+5
            %if rate == 0.5 && AF == 0
            %    DTE = 12.5*0.125;
            %else
            %    DTE = 12.5*0.2;
            %end
                % new DTE (opt17)
                DTE = 12.5*0.3;
            I_app = 2*DTE;
        else
            I_app = 0.0;
        end  
	case 3      % Voltage step
        rate = p(2); % Hz
        period = 1000/rate; % ms
        step_duration = p(3);
        V_test = p(4);
        V_hold1 = -80; T_hold1 = 5;
        V_hold2 = V_test; T_hold2 = step_duration;
	    if mod(t,period) <= T_hold1 %#ok<ALIGN>
            V_clamp = V_hold1;
        elseif mod(t,period) > T_hold1 && mod(t,period) <= T_hold1+T_hold2
            V_clamp = V_hold2;
        else
            V_clamp = V_hold1;
        end
		R_clamp = 0.02;
		I_app = (V_clamp-y(39))/R_clamp; 
    case 4 % caffeine
		rate = p(2); % Hz
        period = 1000/rate; % ms
        if mod(t,period) <= 5 && t < 10e3 
            I_app = 12.5;
        else
            I_app = 0.0;
        end
%         if t > 10e3
%             V_clamp = -75;
%             R_clamp = 0.02;
%             I_app = (V_clamp-y(39))/R_clamp; 
%         end
end  

%% Membrane Potential
I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;
I_Cl_tot = I_ClCa + I_Clbk + I_ClCFTR;
I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl;
I_tot = I_Na_tot + I_Cl_tot + I_Ca_tot + I_K_tot; % [uA/uF]
%ydot(39) = -(I_tot-I_app);
ydot(39) = -(I_tot-I_app)*(1-CaffeineFlag);
vmax = ydot(39);

%% Output adjustment depending on the function call
if (nargin == 3)
    output = ydot;
elseif (nargin == 4) && strcmp(runType,'ydot')
    output = ydot;
elseif (nargin == 4) && strcmp(runType,'rates')
    output = r;
elseif (nargin == 4) && strcmp(runType,'currents')
    %currents = [vmax J_serca*Vsr/Vmyo -I_ncx*2*Cmem/(2*Frdy*Vmyo) I_pca*Cmem/(2*Frdy*Vmyo)];
    currents = [vmax J_serca*Vsr/Vmyo -I_ncx*2*Cmem/(2*Frdy*Vmyo) I_pca*Cmem/(2*Frdy*Vmyo) I_Na I_NaL I_nak I_tof I_kr I_ks I_kur I_ki I_k2p I_sk I_kp I_ClCa I_Catot I_pca I_ncx I_cabk I_nabk I_Clbk];
    output = currents;
end
