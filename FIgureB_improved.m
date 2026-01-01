%2270069 전희선
%% ====== Figure B imporved =======
clear; clc;
close all;

gas = IdealGasMix('gri30.xml');

N  = nSpecies(gas);
R  = gasconstant;       

% Species index
iNH3 = speciesIndex(gas, 'NH3');
iH2  = speciesIndex(gas, 'H2');
iN2  = speciesIndex(gas, 'N2');
iO2  = speciesIndex(gas, 'O2');
iH2O = speciesIndex(gas, 'H2O');

% Molar masses [kg/kmol]
M = molarMasses(gas);
M_NH3 = M(iNH3);
M_H2  = M(iH2);
M_N2  = M(iN2);
M_O2  = M(iO2);
M_H2O = M(iH2O);

% 2. 문제 조건 (PCEC mode)
T_amb   = 298.15;     % ambient [K]
P_amb   = 1e5;        % 1 bar [Pa]
P_sys   = 50e5;       % 50 bar [Pa]

T_sys  = 500 + 273.15;  % [K]
eps_HX  = 0.85;          % HEX effectiveness

steam_util   = 0.80;

% 3. Stream 배열 준비
nDot = zeros(29,1);     % kmol/s
mDot = zeros(29,1);
T    = zeros(29,1);     
P    = zeros(29,1);
X    = zeros(N,29);

%% Stream 14 : Steam feed
nDot(14) = 1.0;         % kmol/s
mDot(14) = nDot(14)*M_H2O;
T(14) = T_amb;
P(14) = P_sys;

X14=zeros(N,1); 
X14(iH2O)=1;
X(:,14)=X14;

set(gas, 'T', T(14), 'P', P(14), 'X', X(:,14));
h14=enthalpy_mass(gas);
%% Stream 15 : After HEX1 (Steam preheating) without T15
P(15) = P_sys;
nDot(15)=nDot(14);
mDot(15)=mDot(14);
X(:,15)=X(:,14);

%% Stream 16 : After Electric Heater (up to cell T)
T(16) = T_sys;
P(16) = P_sys;
nDot(16)=nDot(15);
mDot(16)=mDot(15);
X(:,16)=X(:,15);

set(gas, 'T', T(16), 'P', P(16), 'X', X(:,16));
h16=enthalpy_mass(gas);
%% Stream 17 : Anode outlet (Steam electrolysis)
n_H2O_in = nDot(16);
n_H2O_react = steam_util * n_H2O_in;  % 80% steam utilization

nO2_17  = 0.5*n_H2O_react;
nH2O_17 = n_H2O_in - n_H2O_react;

nDot(17) = nO2_17 + nH2O_17;
mDot(17) = nO2_17*M_O2 + nH2O_17*M_H2O;
T(17) = T_sys;
P(17) = P_sys;

X17 = zeros(N,1);
X17(iO2) = nO2_17  / nDot(17);
X17(iH2O) = nH2O_17 / nDot(17);
X(:,17)=X17;

set(gas, 'T', T(17), 'P', P(17), 'X', X(:,17));
h17=enthalpy_mass(gas);
%% Stream 18
P(18) = P_sys;
nDot(18)=nDot(17);
mDot(18)=mDot(17);
X(:,18)=X17;

%% Initial Guess of T15_in
T15_in = (25+500)/2 +273.15;
while (1)
    T15_in_guess = T15_in;
    % HEX1_2
    % Hot side: 17 -> 17_out
    % Cold side: 15_in -> 15
    [T(15), T17_out, Q_dot1B] = HEX(gas, ...
                T15_in_guess, P_sys, mDot(14), X(:,14), ...   % Cold side
                T(17), P(17), mDot(17), X(:,17), ...   % Hot side
                eps_HX, P_sys);
    % HEX1_1
    % Hot side: 17_out -> 18
    % Cold side: 14 -> 15_in
    [T15_in, T(18), Q_dot1A] = HEX(gas, ...
                T(14), P_sys, mDot(14), X(:,14), ...   % Cold side
                T17_out, P_sys, mDot(17), X(:,17), ...   % Hot side
                eps_HX, P_sys);
    if T15_in_guess - T15_in < 1
        break
    end
end
%% Stream 19 & 20
% nO2_17  = 0.5*n_H2O_react = 0.4
% nH2O_17 = n_H2O_in - n_H2O_react = 0.2
nDot(19) = nO2_17;
mDot(19) = nO2_17*M_O2;

nDot(20) = nH2O_17;
mDot(20) = nH2O_17*M_H2O;

T(19)=T_amb; T(20)=T_amb;
P(19)=P_sys; P(20)=P_sys;

X19 = zeros(N,1); 
X19(iO2) = 1; 
X(:,19)=X19;
set(gas, 'T', T(19), 'P', P(19), 'X', X(:,19));
h19=enthalpy_mass(gas);

X20 = zeros(N,1);
X20(iH2O) = 1; 
X(:,20)=X20;
set(gas, 'T', T(20), 'P', P(20), 'X', X(:,20));
h20=enthalpy_mass(gas);

%% Stream 21 : Fresh N2 feed
T(21)=T_amb;
P(21)=P_sys;

nN2_21 = n_H2O_react/3;
nDot(21) = nN2_21;
mDot(21) = nDot(21)*M_N2;
X21=zeros(N,1); 
X21(iN2)=1;
X(:,21)=X21;
set(gas, 'T', T(21), 'P', P(21), 'X', X(:,21));
h21=enthalpy_mass(gas);

%% Stream 28 : Recycle gas (from condenser + HEX3)
nN2_28=0; nH2_28=0;
nDot(28)=0.0; mDot(28)=0.0;
T(28)=T_amb;
P(28)=P_sys;
X28 = zeros(N,1);
X(:,28)=zeros(N,1);

iteration = 1;
%%
% ================= Iteration history storage =================
iter_hist      = [];
T28_hist       = [];
nDot23_H2_hist = [];
nDot23_N2_hist = [];

while(1)
    % initial guess of stream 28
    T28_old = T(28); 
    nDot28_old = nDot(28);
    %% Stream 22 : Mixer (21 + 28)
    % mass flow
    % 수정한 부분 별 
    nN2_22 = nN2_28 + nN2_21;
    nH2_22 = nH2_28;
    nDot(22) = nN2_22 + nH2_22;
    mDot(22) = nN2_22*M_N2 + nH2_22*M_H2;
    P(22)=P_sys;
    
    % composition
    X22 = zeros(N,1); 
    X22(iN2) = nN2_22 / nDot(22);
    X22(iH2) = nH2_22 / nDot(22);
    X(:,22)=X22;
    
    % temperature by energy balance
    set(gas,'T',T(21),'P',P_sys,'X',X(:,21)); 
    h21 = enthalpy_mass(gas);
    if nDot(28)>0
        set(gas,'T',T28_old,'P',P_sys,'X',X(:,28)); 
        h28 = enthalpy_mass(gas);
    else
        h28 = 0;
    end
    
    H_in_B = mDot(21)*h21 + mDot(28)*h28; 
    h22 = H_in_B/mDot(22);
    
    set(gas,'P',P(22),'X',X(:,22),'H',h22 );
    % T22 solving (adiabatic mixer)
    T(22)=temperature(gas);
    %% Stream 23 without T23
    X(:,23)=X22;
    nDot(23)=nDot(22);
    mDot(23)=mDot(22);
    P(23)=P_sys;
    
   % ================= Store iteration history =================
    iter_hist(end+1) = iteration;
    
    % T28 history
    T28_hist(end+1) = T(28);
    
    % Stream 23: mole flow rates [kmol/s]
    nDot23_H2 = nDot(23) * X(iH2,23);
    nDot23_N2 = nDot(23) * X(iN2,23);
    
    nDot23_H2_hist(end+1) = nDot23_H2;
    nDot23_N2_hist(end+1) = nDot23_N2;
    %% Stream 24
    T(24) = T_sys;
    P(24) = P_sys;
   
    % Stream 24 (NH3 cracking equilibrium)
    T(24) = T_sys;
    P(24) = P_sys;

    nNH3_gen = (2/3) * n_H2O_react;   % electrochemical NH3
    nNH3_24_old = nNH3_gen;
    nN2_24_old = max(nN2_22 - 0.5*nNH3_gen, 0);  % 소비된 N2
    nH2_24_old  = nH2_22;  % 넘어온 H2
    
    nDot_24_old = nNH3_24_old + nN2_24_old + nH2_24_old;
    mDot_24_old = nNH3_24_old*M_NH3 + nN2_24_old*M_N2 + nH2_24_old*M_H2;
    
    X24 = zeros(N,1);
    X24(iNH3)=nNH3_24_old;
    X24(iN2)=nN2_24_old;
    X24(iH2)=nH2_24_old;
    X24 = X24/sum(X24);
    set(gas,'T',T(24),'P',P(24),'X',X24);  % cathode inlet
    equilibrate(gas,'TP'); % NH3 cracking equilibrium
    % 초기 NH3 100% → T_cell, P_sys에서 평형
    XeqB = moleFractions(gas);

    xNH3_B = XeqB(iNH3);
    xN2_B  = XeqB(iN2);
    xH2_B  = XeqB(iH2);
    
    % Conversion x 계산
    % X_NH3 = (A - x)/(B + x)
    % xB = (A-B*X_NH3)/(1+X_NH3)
    A = nNH3_24_old;
    B = nDot_24_old;    % nNH3_24_old+nN2_24_old+nH2_24_old
    xB = (A - B*xNH3_B) / (1 + xNH3_B);

    % 반응물 nNH3_1 기준 scaling
    nNH3_24 = nNH3_24_old * (1 - xB);
    nN2_24  = nNH3_24_old * 0.5 * xB + nN2_24_old;
    nH2_24  = nNH3_24_old * 1.5 * xB + nH2_24_old;

    nDot(24) = nNH3_24 + nN2_24 + nH2_24 ;
    mDot(24) = nNH3_24*M_NH3 + nN2_24*M_N2 + nH2_24*M_H2;
    
    X24 = zeros(N,1);
    X24(iNH3) = nNH3_24 / nDot(24);
    X24(iN2) = nN2_24 / nDot(24);
    X24(iH2) = nH2_24 / nDot(24);

    X(:,24)=X24;
    %% Stream 25 withour T25
    X(:,25)=X(:,24);
    nDot(25)=nDot(24);
    mDot(25)=mDot(24);
    P(25)=P_sys;
    
    %%  ==== HX2 =====
    % Hot side: 24 -> 25
    % Cold side: 22 -> 23
    [T(23), T(25), Q_dot2B] = HEX(gas, ...
                    T(22), P(22), mDot(22), X(:,22), ...   % Cold side
                    T(24), P(24), mDot(24), X(:,24), ...   % Hot side
                    eps_HX, P_sys);
    %% Stream 26 without T26
    X(:,26)=X(:,25);
    nDot(26)=nDot(25);
    mDot(26)=mDot(25);
    nNH3_26 = nDot(26) * X(iNH3,26);
    nN2_26  = nDot(26) * X(iN2,26);
    nH2_26  = nDot(26) * X(iH2,26);
    P(26)=P_sys;
    
    %% Stream 27 without T27
    T(27)=T_amb;
    P(27)=P_sys;
    nN2_27 = nN2_26;
    nH2_27 = nH2_26;
    nDot(27) = nN2_27 + nH2_27;
    mDot(27) = nN2_27*M_N2 + nH2_27*M_H2;
    
    X27 = zeros(N,1);
    X27(iN2) = nN2_27  / nDot(27);
    X27(iH2) = nH2_27  / nDot(27);
    X(:,27) = X27;
    
    %% Stream 29
    T(29)=T_amb;
    P(29)=P_sys;
    nNH3_29 = nNH3_26;
    nDot(29) = nNH3_29;
    mDot(29) = nNH3_29 * M_NH3;
    
    X29 = zeros(N,1);
    X29(iNH3) = nNH3_29/nDot(29); 
    X(:,29) = X29;
    
    %% === HEX3 ===
    % Hot side: 25 -> 26
    % Cold side: 27 -> 28
    [T(28), T(26), Q_dot3B] = HEX(gas, ...
                    T(27), P(27), mDot(27), X(:,27), ...   % Cold side
                    T(25), P(25), mDot(25), X(:,25), ...   % Hot side
                    eps_HX, P_sys);

    % stream 28 update
    mDot(28) = mDot(27);
    nDot(28) = nDot(27);
    X(:,28) = X(:,27);
    nN2_28 = nDot(28) * X(iN2,28);
    nH2_28 = nDot(28) * X(iH2,28);
    iteration = iteration +1;

    if abs(T(28)-T28_old) < 0.1...
            && abs(nDot(28)-nDot28_old)/nDot28_old < 1e-6
        break;
    else
        T(28) = (T(28)+T28_old)/2;
    end
    
end

%% ================= Convergence Plots =================

% --- T28 ---
figure;
plot(iter_hist, T28_hist, '-o','LineWidth',1.5);
xlabel('Iteration');
ylabel('T_{28} [K]');
title('Convergence of Recycle Temperature T_{28}');
grid on;

% --- H2 mole flow rate ---
figure;
plot(iter_hist, nDot23_H2_hist, '-s','LineWidth',1.5);
xlabel('Iteration');
ylabel('ṅ_{H2,23}  [kmol/s]');
title('Convergence of H2 Mole Flow Rate at Stream 23');
grid on;

% --- N2 mole flow rate ---
figure;
plot(iter_hist, nDot23_N2_hist, '-d','LineWidth',1.5);
xlabel('Iteration');
ylabel('ṅ_{N2,23}  [kmol/s]');
title('Convergence of N2 Mole Flow Rate at Stream 23');
grid on;

%% ===============================
% FINAL POST-CONVERGENCE UPDATE
% (수렴된 T28 기준으로 전체 stream 재계산)
% ===============================
% stream 28 
nN2_28 = nDot(28) * X(iN2,28);
nH2_28 = nDot(28) * X(iH2,28);

% Stream 22 : Mixer (21 + 28)
nN2_22 = nN2_28 + nN2_21;
nH2_22 = nH2_28;
nDot(22) = nN2_22 + nH2_22;
mDot(22) = nN2_22*M_N2 + nH2_22*M_H2;
P(22)=P_sys;

X22 = zeros(N,1); 
X22(iN2) = nN2_22 / nDot(22);
X22(iH2) = nH2_22 / nDot(22);
X(:,22)=X22;

set(gas,'T',T(21),'P',P_sys,'X',X(:,21)); h21 = enthalpy_mass(gas);
set(gas,'T',T(28),'P',P_sys,'X',X(:,28)); h28 = enthalpy_mass(gas);

H_in_B = mDot(21)*h21 + mDot(28)*h28;
h22 = H_in_B/mDot(22);

set(gas,'P',P(22),'X',X(:,22),'H',h22 );
T(22)=temperature(gas);
h22=enthalpy_mass(gas);

% Stream 23 without T23
X(:,23)=X22;
nDot(23)=nDot(22);
mDot(23)=mDot(22);
P(23)=P_sys;

% Stream 24 (NH3 cracking equilibrium)
T(24) = T_sys;
P(24) = P_sys;

nNH3_gen = (2/3) * n_H2O_react;   % electrochemical NH3
nNH3_24_old = nNH3_gen;
nN2_24_old  = nN2_22 - 0.5*nNH3_gen;  % 소비된 N2
nH2_24_old  = nH2_22;                 % 넘어온 H2

nDot_24_old = nNH3_24_old+nN2_24_old+nH2_24_old;
mDot_24_old = nNH3_24_old*M_NH3 + nN2_24_old*M_N2 + nH2_24_old*M_H2;

X24 = zeros(N,1);
X24(iNH3)=nNH3_24_old;
X24(iN2)=nN2_24_old;
X24(iH2)=nH2_24_old;
X24 = X24/sum(X24);

set(gas,'T',T_sys,'P',P_sys,'X',X24);
equilibrate(gas,'TP');

% 초기 NH3 100% → T_cell, P_sys에서 평형
XeqB = moleFractions(gas);

xNH3_B = XeqB(iNH3);
xN2_B  = XeqB(iN2);
xH2_B  = XeqB(iH2);

% Conversion x 계산
% X_NH3 = (A - x)/(B + x)
% xB = (A-B*X_NH3)/(1+X_NH3)
A = nNH3_24_old;
B = nDot_24_old;    % nNH3_24_old+nN2_24_old+nH2_24_old
xB = (A - B*xNH3_B) / (1 + xNH3_B);

% 반응물 nNH3_1 기준 scaling
nNH3_24 = nNH3_24_old * (1 - xB);
nN2_24  = nNH3_24_old * 0.5 * xB + nN2_24_old;
nH2_24  = nNH3_24_old * 1.5 * xB + nH2_24_old;

nDot(24) = nNH3_24 + nN2_24 + nH2_24 ;
mDot(24) = nNH3_24*M_NH3 + nN2_24*M_N2 + nH2_24*M_H2;

X24 = zeros(N,1);
X24(iNH3) = nNH3_24 / nDot(24);
X24(iN2) = nN2_24 / nDot(24);
X24(iH2) = nH2_24 / nDot(24);

X(:,24) = X24;
set(gas, 'T', T(24), 'P', P(24), 'X', X(:,24));
h24=enthalpy_mass(gas);

% Stream 25 without T25
X(:,25)=X(:,24);
nDot(25)=nDot(24);
mDot(25)=mDot(24);
P(25)=P_sys;

% HX2 (최종)
[T(23), T(25), Q_dot2B] = HEX(gas,...
    T(22),P(22),mDot(22),X(:,22),...
    T(24),P(24),mDot(24),X(:,24),...
    eps_HX,P_sys);

set(gas, 'T', T(25), 'P', P(25), 'X', X(:,25));
h25=enthalpy_mass(gas);
set(gas, 'T', T(23), 'P', P(23), 'X', X(:,23));
h23=enthalpy_mass(gas);

% Stream 26 without T26
X(:,26)=X(:,25);
nDot(26)=nDot(25);
mDot(26)=mDot(25);
nNH3_26 = nDot(26) * X(iNH3,26);
nN2_26  = nDot(26) * X(iN2,26);
nH2_26  = nDot(26) * X(iH2,26);
P(26)=P_sys;

% Stream 27
T(27) = T_amb;
P(27) = P_sys;
nN2_27 = nN2_26;
nH2_27 = nH2_26;
nDot(27) = nN2_27 + nH2_27;
mDot(27) = nN2_27*M_N2 + nH2_27*M_H2;

X27=zeros(N,1);
X27(iN2)=nN2_27/nDot(27);
X27(iH2)=nH2_27/nDot(27);

X(:,27)=X27;
set(gas, 'T', T(27), 'P', P(27), 'X', X(:,27));
h27=enthalpy_mass(gas);

% HX3 (최종 T26, T28 재확정)
[T(28), T(26), Q_dot3B] = HEX(gas,...
    T(27),P(27),mDot(27),X(:,27),...
    T(25),P(25),mDot(25),X(:,25),...
    eps_HX,P_sys);
set(gas, 'T', T(28), 'P', P(28), 'X', X(:,27));
h28=enthalpy_mass(gas);
set(gas, 'T', T(26), 'P', P(26), 'X', X(:,26));
h26=enthalpy_mass(gas);

% Stream 29
T(29)=T_amb;
P(29)=P_sys;
nNH3_29 = nNH3_26;
nDot(29) = nNH3_29;
mDot(29) = nNH3_29 * M_NH3;

X29 = zeros(N,1);
X29(iNH3) = nNH3_29/nDot(29); 
X(:,29) = X29;
set(gas, 'T', T(29), 'P', P(29), 'X', X(:,29));
h29=enthalpy_mass(gas);

% stream 28 update
mDot(28) = mDot(27);
nDot(28) = nDot(27);
X(:,28) = X(:,27);
nN2_28 = nDot(28) * X(iN2,28);
nH2_28 = nDot(28) * X(iH2,28);
%% 출력
fprintf("\n========= PCEC Mode : Stream Properties =========\n");

for j=14:29
    fprintf("Stream %2d: nDot = %.4f kmol/s, mDot = %.4f kg/s, T = %.2f K,  P = %.2f bar\n", ...
        j, nDot(j), mDot(j), T(j), P(j)/1e5);
    fprintf("   X(NH3)=%.3f  X(H2)=%.3f  X(N2)=%.3f  X(O2)=%.3f  X(H2O)=%.3f\n", ...
        X(iNH3,j), X(iH2,j), X(iN2,j), X(iO2,j), X(iH2O,j));
end

fprintf("\n========= mass conservation =========\n");
fprintf("전체: mDot(14) + mDot(21):%.2f, mDot(19) + mDot(20) + mDot(29):%.2f\n",mDot(14)+mDot(21), mDot(19)+mDot(20)+mDot(29));
fprintf("HEX1: mDot(14) + mDot(17):%.2f, mDot(18) + mDot(15):%.2f\n",mDot(14)+mDot(17), mDot(18)+mDot(15));
fprintf("Water Condensor: mDot(18):%.2f, mDot(19) + mDot(20):%.2f\n",mDot(18), mDot(19)+mDot(20));
fprintf("mDot(15):%.2f, mDot(16):%.2f\n",mDot(15), mDot(16));
fprintf("PCEC: mDot(16) + mDot(23):%.2f, mDot(17) + mDot(24):%.2f\n",mDot(16)+mDot(23), mDot(17)+mDot(24));
fprintf("HEX2: mDot(22) + mDot(24):%.2f, mDot(23) + mDot(25):%.2f\n",mDot(22)+mDot(24), mDot(23)+mDot(25));
fprintf("Mixer: mDot(22):%.2f, mDot(21) + mDot(28):%.2f\n",mDot(22), mDot(21)+mDot(28));
fprintf("HEX3: mDot(25) + mDot(27):%.2f, mDot(26) + mDot(28):%.2f\n",mDot(25)+mDot(27), mDot(26)+mDot(28));
fprintf("Ammonia Condensor: mDot(26):%.2f, mDot(27) + mDot(29):%.2f\n",mDot(26), mDot(27)+mDot(29));
%% ======================== PCEC ELECTRICITY INPUT & SYSTEM EFFICIENCY =========================
fprintf("\n========= PCEC SYSTEM EFFICIENCY =========\n");
% Constants / LHV
LHV_NH3 = 18.6e6;    % [J/kg]    (18.6 MJ/kg) 
F       = 96485e3;        % [C/kmol]

% Reaction Gibbs Energy at cell conditions
% H2O -> H2 + 0.5 O2
Xp = zeros(N,1);

% H2O
Xp(iH2O)=1;
set(gas,'T',T_sys,'P',1e5,'X',Xp);
g_H2O = gibbs_mole(gas);   % [J/kmol]

% N2
Xp(:)=0; Xp(iN2)=1;
set(gas,'T',T_sys,'P',1e5,'X',Xp);
g_N2 = gibbs_mole(gas);

% O2
Xp(:)=0; Xp(iO2)=1;
set(gas,'T',T_sys,'P',1e5,'X',Xp);
g_O2 = gibbs_mole(gas);

% NH3
Xp(:)=0; Xp(iNH3)=1;
set(gas,'T',T_sys,'P',1e5,'X',Xp);
g_NH3 = gibbs_mole(gas);

% ΔG0 for steam electrolysis
dG0 = 3*g_H2O + g_N2 - 2*g_NH3 - 1.5*g_O2;     % [J/kmol], endothermic
E0  = -dG0 / (6*F);               % [V] reversible voltage

% Nernst Voltage
% Using stream 17 (anode outlet) & steam feed stream 16

% NH3
PCEC_NH3_1 = X(iNH3,23); % in
PCEC_NH3_2 = X(iNH3,24); % out
PCEC_NH3 = 0.5*PCEC_NH3_1 + 0.5*PCEC_NH3_2; 
% O2
PCEC_O2_1 = X(iO2, 16); % in
PCEC_O2_2 = X(iO2, 17); % out
PCEC_O2 = 0.5*PCEC_O2_1 + 0.5*PCEC_O2_2;
% H2O
PCEC_H2O_1 = X(iH2O, 16); % in
PCEC_H2O_2 = X(iH2O, 17); % out
PCEC_H2O = 0.5*PCEC_H2O_1 + 0.5*PCEC_H2O_2;
% N2
PCEC_N2_1 = X(iN2,23); % in
PCEC_N2_2 = X(iN2,24); % out
PCEC_N2 = 0.5*PCEC_N2_1 + 0.5*PCEC_N2_2;

% partial pressure
pNH3  = max(PCEC_NH3*P_sys, 1e-12);
pO2  = max(PCEC_O2*P_sys,  1e-12);
pH2O = max(PCEC_H2O*P_sys, 1e-12);
pN2 = max(PCEC_N2*P_sys, 1e-12);

Q = (pNH3^2 * pO2^1.5) / (pH2O^3 * pN2);    % Nernst for 3H2O + N2 -> 2NH3 + 1.5O2
Erev = E0 - (R*T_sys/(6*F))*log(Q);

% Electrical current & power
I_PCEC = 2 * F * n_H2O_react;      % [A]
W_PCEC = Erev * I_PCEC;            % [W]

% Electric heater power (16 <- 15)
set(gas,'T',T(15),'P',P(15),'X',X(:,15)); h15 = enthalpy_mass(gas);
set(gas,'T',T(16),'P',P(16),'X',X(:,16)); h16 = enthalpy_mass(gas);

W_heater = mDot(15)*(h16 - h15);   % [W] = kg/s * (J/kg)

% Hydrogen production power
LHV_out = mDot(29)*LHV_NH3;   % Stream 29

% Efficiency definitions
eta_PCEC  = LHV_out / (W_PCEC + W_heater);        % including heating

% --------------------
% Print Results
% --------------------
fprintf("Reversible voltage Erev         : %.3f V\n", Erev);
fprintf("PCEC current I                  : %.3e A\n", I_PCEC);
fprintf("PCEC electric power W_PCEC      : %.3e W\n", W_PCEC);
fprintf("Electric heater power Q_heater  : %.3e W\n", W_heater);
fprintf("Hydrogen LHV power out          : %.3e W\n", LHV_out);
fprintf("PCEC system efficiency η        : %.4f (%.2f %%)\n", eta_PCEC, eta_PCEC*100);

%% ==== PCEC HEAT TRANSFER ====

% Anode inlet: stream 16
set(gas,'T',T(16),'P',P(16),'X',X(:,16));
h16 = enthalpy_mass(gas);
% Cathode inlet : stream 23
set(gas,'T',T(23),'P',P(23),'X',X(:,23));
h23 = enthalpy_mass(gas);

% Anode outlet: stream 17
set(gas,'T',T(17),'P',P(17),'X',X(:,17));
h17 = enthalpy_mass(gas);
% Cathode outlet: stream 24
set(gas,'T',T(24),'P',P(24),'X',X(:,24));
h24 = enthalpy_mass(gas);

% Total PCFC heat transfer
Q_dot_PCEC = -(mDot(16)*h16 + mDot(23)*h23)...
             +(mDot(17)*h17 + mDot(24)*h24)...
             -W_PCEC; % Q_in - Q_out

fprintf("\n========= PCFC HEAT TRANSFER =========\n");
fprintf("Work done by PCFC                                 : %.3e W\n", W_PCEC);
fprintf("TOTAL PCEC heat transfer(in) (PCEC 흡수열 (Q_in)) : %.3e W\n", Q_dot_PCEC);

fprintf("\n========= PCFC Partial Pressure =========\n");
fprintf("NH3, Partial pressure : %.3e Pa\n", pNH3);
fprintf("O2, Partial pressure : %.3e Pa\n", pO2);
fprintf("H2O, Partial pressure : %.3e Pa\n", pH2O);
fprintf("N2, Partial pressure : %.3e Pa\n", pN2);
fprintf("P_sys : %.3e Pa\n", P_sys);

fprintf("\n========= 비교=========\n");
fprintf("HEX1 Heat Transfer             : %.3e W\n", Q_dot1B);
fprintf("HEX2 Heat Transfer             : %.3e W\n", Q_dot2B);
fprintf("PCEC electric power W_PCFC     : %.3e W\n", W_PCEC);
fprintf("Electric heater power Q_heater : %.3e W\n", W_heater);
fprintf("Total PCEC heat transfer       : %.3e W\n", Q_dot_PCEC);

%% function
function [T_c_out, T_h_out, Q_dot] = HEX(gas, ...
    T_c_in, P_c, mDot_c, X_c, ...
    T_h_in, P_h, mDot_h, X_h, ...
    eps_HX, P_sys)

% =============================================================
% Heat Exchanger with effectiveness (mass basis)
% Cold side  : in (c_in)  -> out (c_out)
% Hot side   : in (h_in)  -> out (h_out)
%
% All enthalpies are mass-based [J/kg]
% =============================================================

% --- Cold side inlet ---
set(gas, 'T', T_c_in, 'P', P_c, 'X', X_c);
h_c_in = enthalpy_mass(gas);
% --- Cold side outlet (max heating assumption: to hot inlet T) ---
set(gas, 'T', T_h_in, 'P', P_c, 'X', X_c);
h_c_out_assump = enthalpy_mass(gas);

% --- Hot side inlet ---
set(gas, 'T', T_h_in, 'P', P_h, 'X', X_h);
h_h_in = enthalpy_mass(gas);

% --- Hot side outlet (max cooling assumption: to cold inlet T) ---
set(gas, 'T', T_c_in, 'P', P_h, 'X', X_h);
h_h_out_assump = enthalpy_mass(gas);
% --- Maximum possible heat transfer ---
Q_dot_C_max = mDot_c * (h_c_out_assump - h_c_in);   % Cold side capacity
Q_dot_H_max = mDot_h * (h_h_in - h_h_out_assump);   % Hot side capacity

Q_dot_max = min(Q_dot_C_max, Q_dot_H_max);
Q_dot     = eps_HX * Q_dot_max;

% --- Outlet enthalpies ---
h_c_out = h_c_in + Q_dot / mDot_c;
h_h_out = h_h_in - Q_dot / mDot_h;

% enthalpy → temperature
set(gas, 'P', P_sys, 'H', h_c_out, 'X', X_c);
T_c_out = temperature(gas);

set(gas, 'P', P_sys, 'H', h_h_out, 'X', X_h);
T_h_out = temperature(gas);

end