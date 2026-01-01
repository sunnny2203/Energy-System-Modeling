%2270069 전희선
%% ====== Figure A =======
clear; clc;

% 1. Cantera 초기화
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

% 2. 문제 조건 (PCFC mode)
T_amb   = 298.15;     % ambient [K]
P_amb   = 1e5;        % 1 bar [Pa]
P_sys   = 50e5;       % 50 bar [Pa]

T_sys  = 500 + 273.15;  % [K]
eps_HX  = 0.85;          % HEX effectiveness
fuel_util   = 0.80;

% 3. Stream 배열 준비
nDot = zeros(13,1);     % kmol/s
mDot = zeros(13,1);
T    = zeros(13,1);     
P    = zeros(13,1);
X    = zeros(N,13);

%% Stream 1 : NH3 feed (pure NH3)
nNH3_1 = 1.0;      % kmol/s
nDot(1) = nNH3_1;
mDot(1) = nNH3_1*M_NH3;

T(1) = T_amb;
P(1) = P_sys;

X1 = zeros(N,1);
X1(iNH3) = 1;
X(:,1) = X1;
set(gas, 'T',T(1),'P',P(1),'X',X(:,1));
h1=enthalpy_mass(gas);
%% Stream 2 without T2
P(2) = P_sys;
X(:,2)   = X1;
nDot(2) = nDot(1);
mDot(2) = mDot(1);
%% Stream 3 : PCFC anode outlet (fuel utilization 80%)
T(3) = T_sys;
P(3) = P_sys;

% NH3 cracking
set(gas,'T',T_sys,'P',P_sys,'X',X(:,2));
equilibrate(gas,'TP');
XeqA = moleFractions(gas);

xNH3_A = XeqA(iNH3);
xN2_A  = XeqA(iN2);
xH2_A  = XeqA(iH2);

% Conversion x 계산
% X_NH3 = (1 - x)/(1 + x)
xA = (1 - xNH3_A) / (1 + xNH3_A);
nNH3_2p = nNH3_1 * (1 - xA);
nN2_2p  = nNH3_1 * 0.5 * xA;
nH2_2p  = nNH3_1 * 1.5 * xA;

nDot_2p = nNH3_2p + nN2_2p + nH2_2p;
mDot_2p = nNH3_2p*M_NH3 + nN2_2p*M_N2 + nH2_2p*M_H2;

X2p = zeros(N,1);
X2p(iNH3,1) = nNH3_2p/nDot_2p;
X2p(iN2,1) = nN2_2p/nDot_2p;
X2p(iH2,1) = nH2_2p/nDot_2p;

nH2_3  = nH2_2p * (1 - fuel_util);
nN2_3  = nN2_2p;
nNH3_3 = nNH3_2p;

nDot(3) = nH2_3 + nN2_3 + nNH3_3;
mDot(3) = nH2_3*M_H2 + nN2_3*M_N2 + nNH3_3*M_NH3;

X3 = zeros(N,1);
X3(iH2)  = nH2_3  / nDot(3);
X3(iN2)  = nN2_3  / nDot(3);
X3(iNH3) = nNH3_3 / nDot(3);
X(:,3)   = X3;

% 반응된 H2
nH2_reacted = nH2_2p * fuel_util;

%% Stream 4 : Pure O2 for combustor (stoichiometric)
nO2_4 = 0.5 * nH2_3 + 0.75*nNH3_3; 
nDot(4) = nO2_4;
mDot(4) = nO2_4*M_O2;

T(4) = T_amb;
P(4) = P_sys;

X4 = zeros(N,1);
X4(iO2) = 1;
X(:,4) = X4;

%% Stream 5 : Combustor adiabatic temperature

% Stream 3 enthalpies (pure species)
set(gas,'T',T(3),'P',P(3),'X',X(:,3));
h3= enthalpy_mass(gas);

% Stream 4
set(gas,'T',T(4),'P',P(4),'X', X(:,4));
h4=enthalpy_mass(gas);

% 총 유입 엔탈피
H_in_A = mDot(3)*h3 + mDot(4)*h4;
mDot(5) =mDot(3) +mDot(4);

%
nH2O_5 = nH2_3 + 1.5*nNH3_3;
nN2_5  = nN2_3 + 0.5*nNH3_3;

nDot(5) = nH2O_5 + nN2_5;
mDot_5 = nH2O_5*M_H2O + nN2_5*M_N2;

X5 = zeros(N,1);
X5(iH2O) = nH2O_5 / nDot(5);
X5(iN2)  = nN2_5  / nDot(5);
X(:,5)= X5;

P(5) = P_sys;

set(gas,'P',P(5),'X',X(:,5),'H',H_in_A/mDot(5));
% T5 solving (adiabatic combustor)
T(5)=temperature(gas);
h5=enthalpy_mass(gas);
%% stream 6 without T6
nH2O_6 = nH2O_5;
nN2_6 = nN2_5 ; 
nDot(6)= nDot(5);
mDot(6)= mDot(5);

X(:,6) = X(:,5);      % 조성 변화 없음 (단순 sensible HX)
P(6)=P_sys;
%% HEX1

set(gas, 'T', T(1), 'P', P(1), 'X', X(:,1));
h_c_in = enthalpy_mass(gas);
% --- Cold side outlet (max heating assumption: to hot inlet T) ---
set(gas, 'T', T(5), 'P', P(2), 'X', X(:,1));
equilibrate(gas,'TP');
h_c_out_assump = enthalpy_mass(gas);

% --- Hot side inlet ---
set(gas, 'T', T(5), 'P', P(5), 'X', X(:,5));
h_h_in = enthalpy_mass(gas);

% --- Hot side outlet (max cooling assumption: to cold inlet T) ---
set(gas, 'T', T(1), 'P', P(6), 'X', X(:,5));
h_h_out_assump = enthalpy_mass(gas);
% --- Maximum possible heat transfer ---
Q_dot_C_max = mDot(1) * (h_c_out_assump - h_c_in);   % Cold side capacity
Q_dot_H_max = mDot(5) * (h_h_in - h_h_out_assump);   % Hot side capacity

Q_dot_max = min(Q_dot_C_max, Q_dot_H_max);
Q_dot1A     = eps_HX * Q_dot_max;

% --- Outlet enthalpies ---
h_c_out = h_c_in + Q_dot1A / mDot(1);
h_h_out = h_h_in - Q_dot1A / mDot(5);

% enthalpy → temperature
set(gas, 'P', P(2), 'H', h_c_out, 'X', X(:,1));
equilibrate(gas,'HP');
T(2) = temperature(gas);
h2 = enthalpy_mass(gas);

set(gas, 'P', P(6), 'H', h_h_out, 'X', X(:,5));
T(6) = temperature(gas);
h6 = enthalpy_mass(gas);
%% Streams 7,8 : Condenser
T(7)=T_amb; P(7)=P_sys;
T(8)=T_amb; P(8)=P_sys;

nH2O_7 = nH2O_5;
nN2_8  = nN2_5;

nDot(7)=nH2O_7;   X7=zeros(N,1); X7(iH2O)=1; X(:,7)=X7;
mDot(7)=nDot(7)*M_H2O;
set(gas, 'T', T(7), 'P', P(7), 'X', X(:,7));
h7=enthalpy_mass(gas);

nDot(8)=nN2_8;    X8=zeros(N,1); X8(iN2)=1;  X(:,8)=X8;
mDot(8)=nDot(8)*M_N2;
set(gas, 'T', T(8), 'P', P(8), 'X', X(:,8));
h8=enthalpy_mass(gas);

%% Stream 9 : Air inlet (ambient)
nO2_cons = 0.5*nH2_reacted;
nO2_9 = 3*nO2_cons;     % 300% excess air → 3배 O2

nDot(9) = nO2_9 / 0.21;
mDot(9) = nDot(9)*(0.21*M_O2+0.79*M_N2);

T(9)=T_amb; P(9)=P_amb;
X9=zeros(N,1); X9(iO2)=0.21; X9(iN2)=0.79; X(:,9)=X9;
set(gas, 'T', T(9), 'P', P(9), 'X', X(:,9));
h9=enthalpy_mass(gas);

%% Stream 10: compressed air
nDot(10)=nDot(9);
mDot(10)=mDot(9);
X(:,10)=X9;
P(10)=P_sys;
T(10)=T_amb; % isothermal
set(gas, 'T', T(10), 'P', P(10), 'X', X(:,10));
h10=enthalpy_mass(gas);

%% Stream 12: cathode outlet (cell 온도 = T_cell)
T(12)=T_sys; 
P(12)=P_sys;

% Use stream 10 composition and flowrate for cathode inlet
nO2_10 = 0.21*nDot(10);
nN2_10 = 0.79*nDot(10);
nH2O_10 = 0;

nO2_12 = nO2_10 - nO2_cons;
nN2_12 = nN2_10;
nH2O_12 = nH2_reacted;

nDot(12)=nO2_12 + nN2_12 + nH2O_12;
mDot(12)=nO2_12*M_O2 + nN2_12*M_N2 + nH2O_12*M_H2O;

X12=zeros(N,1);
X12(iO2)=nO2_12/nDot(12);
X12(iN2)=nN2_12/nDot(12);
X12(iH2O)=nH2O_12/nDot(12);
X(:,12)=X12;
set(gas, 'T', T(12), 'P', P(12), 'X', X(:,12));
h12=enthalpy_mass(gas);

%% Stream 11: cold outlet, without T11
P(11)=P_sys;
nDot(11)=nDot(10);
mDot(11)=mDot(10);
X(:,11)=X(:,10);

%% Stream 13: hot outlet, without T13
P(13)=P_sys;
nDot(13)=nDot(12);
mDot(13)=mDot(12);
X(:,13)=X(:,12);

%% ==== HX2 =====
% Hot side: 12 -> 13
% Cold side: 10 -> 11

[T(11), T(13), Q_dot2A] = HEX(gas, ...
    T(10), P(10), mDot(10), X(:,10), ...   % Cold side
    T(12), P(12), mDot(12), X(:,12), ...   % Hot side
    eps_HX, P_sys);

set(gas, 'T', T(13), 'P', P(13), 'X', X(:,13));
h13=enthalpy_mass(gas);

set(gas, 'T', T(11), 'P', P(11), 'X', X(:,11));
h11=enthalpy_mass(gas);
%% 출력
fprintf("\n========= PCFC Mode : Stream Properties =========\n");

for j=1:13
    fprintf("Stream %2d: nDot = %.4f kmol/s, mDot = %.4f kg/s, T = %.2f K,  P = %.2f bar\n", ...
        j, nDot(j), mDot(j), T(j), P(j)/1e5);
    fprintf("   X(NH3)=%.3f  X(H2)=%.3f  X(N2)=%.3f  X(O2)=%.3f  X(H2O)=%.3f\n", ...
        X(iNH3,j), X(iH2,j), X(iN2,j), X(iO2,j), X(iH2O,j));
end

fprintf("\n========= mass conservation =========\n");
fprintf("total sys: mDot(1) + mDot(4) + mDot(9):%.2f, mDot(7) + mDot(8) + mDot(13):%.2f\n",mDot(1)+mDot(4)+mDot(9), mDot(7)+mDot(8)+mDot(13));
fprintf("HEX2: mDot(1) + mDot(6):%.2f, mDot(2) + mDot(5):%.2f\n",mDot(1)+mDot(6), mDot(2)+mDot(5));
fprintf("combustor: mDot(3) + mDot(4):%.2f, mDot(5):%.2f\n",mDot(3)+mDot(4), mDot(5));
fprintf("water condensor: mDot(6):%.2f, mDot(7) + mDot(8):%.2f\n",mDot(6), mDot(7)+mDot(8));
fprintf("PCFC: mDot(2) + mDot(11):%.2f, mDot(3) + mDot(12):%.2f\n",mDot(2)+mDot(11), mDot(3)+mDot(12));
fprintf("HEX2: mDot(12) + mDot(10):%.2f, mDot(13) + mDot(11):%.2f\n",mDot(12)+mDot(10), mDot(13)+mDot(11));

%% ====== PCFC ELECTRICAL POWER & SYSTEM EFFICIENCY ======

% LHV 정의
LHV_NH3 = 18.6e6;    % [J/kg]    (18.6 MJ/kg) 

% Faraday constant: kmol 기준으로 맞추기
F = 96485e3;         % [C/kmol]
% (Cantera의 R, gibbs_mole이 kmol 기준이므로 F도 kmol 기준이 제일 깔끔함)

% Gibbs energies at T_cell, 1 bar
Xp = zeros(N,1);

Xp(iH2)=1;
set(gas,'T',T_sys,'P',1e5,'X',Xp);   
g_H2 = gibbs_mole(gas);   % [J/kmol]

Xp(:)=0; Xp(iO2)=1;
set(gas,'T',T_sys,'P',1e5,'X',Xp);
g_O2 = gibbs_mole(gas);   % [J/kmol]

Xp(:)=0; Xp(iH2O)=1;
set(gas,'T',T_sys,'P',1e5,'X',Xp);
g_H2O = gibbs_mole(gas);  % [J/kmol]

% ΔG0 and reversible voltage
dG0 = g_H2O - g_H2 - 0.5*g_O2;   % [J/kmol]
E0  = -dG0 / (2*F);              % [V], J/kmol / C/kmol

% Nernst 전압
% H2: PCFC anode의 stream 2p 사용 (과제 조건)

PCFC_H2_1   = X2p(iH2,1);   % mole fraction at 2p
PCFC_H2_2   = X(iH2, 3);
PCFC_H2     = 0.5* PCFC_H2_1 + 0.5* PCFC_H2_2;

% O2: cathode side air (stream 11 or 10; 11 사용)
PCFC_O2_1   = X(iO2, 11); % in
PCFC_O2_2   = X(iO2, 12); % out
PCFC_O2     = 0.5*PCFC_O2_1 + 0.5*PCFC_O2_2;

% H2O: cathode outlet (stream 12)
PCFC_H2O_1  = X(iH2O, 11); % in 
PCFC_H2O_2  = X(iH2O, 12); % out
PCFC_H2O    = 0.5*PCFC_H2O_1 + 0.5*PCFC_H2O_2;

% partial pressure [Pa]
pH2  = PCFC_H2  * P_sys;
pO2  = PCFC_O2  * P_sys;
pH2O = PCFC_H2O * P_sys;

% to prevent the pressure to be zero
pH2  = max(pH2,  1e-12);
pO2  = max(pO2,  1e-12);
pH2O = max(pH2O, 1e-12);

% H2 + 1/2 O2 -> H2O 에 대한 Nernst 식
Q    = pH2O / (pH2 * sqrt(pO2));
Erev = E0 - (R*T_sys/(2*F))*log(Q);

% PCFC current & power
% nH2_reacted : [kmol/s], F: [C/kmol] 이므로 곱하면 A가 됨
I_PCFC = 2 * F * nH2_reacted;   % [A]
W_PCFC = Erev * I_PCFC;         % [W]

% Compressor power (isothermal, reversible)
W_comp = nDot(9) * R * T_amb * log(P(10)/P(9));  % [W]

% --- 연료 입력 (NH3 LHV, 질량 기준) ---
Fuel_LHV_in =  xA* mDot(1) * LHV_NH3;  % [kg/s] * [J/kg] = [J/s]

% 효율
eta_PCFC = (W_PCFC - W_comp) / Fuel_LHV_in;

fprintf("\n========= PCFC CELL & SYSTEM EFFICIENCY =========\n");
fprintf("Reversible voltage Erev        : %.3f V\n", Erev);
fprintf("PCFC current I                 : %.3e A\n", I_PCFC);
fprintf("PCFC electric power W_PCFC     : %.3e W\n", W_PCFC);
fprintf("Compressor power W_comp        : %.3e W\n", W_comp);
fprintf("연료 입력                      : %.3e W\n", Fuel_LHV_in);
fprintf("PCFC system efficiency η       : %.4f (%.2f %%)\n", eta_PCFC, eta_PCFC*100);

%% ==== PCFC HEAT TRANSFER ====

% Anode inlet: stream 2p
set(gas,'T',T_sys,'P',P_sys,'X',X2p);
h2p = enthalpy_mass(gas);
% Cathode inlet : stream 11
set(gas,'T',T(11),'P',P(11),'X',X(:,11));
h11 = enthalpy_mass(gas);

% Anode outlet: stream3
set(gas,'T',T(3),'P',P(3),'X',X(:,3));
h3=enthalpy_mass(gas);
% Cathode outlet: stream 12
set(gas,'T',T(12),'P',P(12),'X',X(:,12));
h12 = enthalpy_mass(gas);

% Total PCFC heat transfer
Q_dot_PCFC = -(mDot_2p*h2p + mDot(11)*h11)...
             +(mDot(3)*h3 + mDot(12)*h12)...
             +W_PCFC; % Q_in

fprintf("\n========= PCFC HEAT TRANSFER =========\n");
fprintf("Work done by PCFC                                 : %.3e W\n", W_PCFC);
fprintf("TOTAL PCFC heat transfer(in) (PCFC 흡수열 (Q_in)) : %.3e W\n", Q_dot_PCFC);

fprintf("\n========= PCFC Partial Pressure =========\n");
fprintf("H2, Partial pressure : %.3e Pa\n", pH2);
fprintf("O2, Partial pressure : %.3e Pa\n", pO2);
fprintf("H2O, Partial pressure : %.3e Pa\n", pH2O);
fprintf("P_sys : %.3e Pa\n", P_sys);

fprintf("\n========= 비교=========\n");
fprintf("HEX1 Heat Transfer             : %.3e W\n", Q_dot1A);
fprintf("HEX2 Heat Transfer             : %.3e W\n", Q_dot2A);
fprintf("PCFC electric power W_PCFC     : %.3e W\n", W_PCFC);
fprintf("Compressor power W_comp        : %.3e W\n", W_comp);
fprintf("Total PCFC heat transfer       : %.3e W\n", Q_dot_PCFC);

%% Check Energy Consesrvation
fprintf("\n================ ENERGY CONSERVATION CHECK =================\n");

% Water Condenser
H_in  = mDot(6)*h6;
H_out = mDot(7)*h7 + mDot(8)*h8;

res_wc = H_out - H_in;

fprintf("\n[Water Condenser]\n");
fprintf("Hin  = %.6e W\n", H_in);
fprintf("Hout = %.6e W\n", H_out);
fprintf("Residual = %.6e W\n", res_wc);

% External Cracker: HEX1
H_in  = mDot(1)*h1;
H_out = mDot_2p*h2p;

res_cracker = H_out - H_in;

fprintf("\n[External Cracker: HEX1]\n");
fprintf("Hin  = %.6e W\n", H_in);
fprintf("Hout = %.6e W\n", H_out);
fprintf("Residual = %.6e W\n", res_cracker);

% Catalytic Combustor
H_in  = mDot(3)*h3 + mDot(4)*h4;
H_out = mDot(5)*h5;

res_comb = H_out - H_in;

fprintf("\n[Catalytic Combustor]\n");
fprintf("Hin  = %.6e W\n", H_in);
fprintf("Hout = %.6e W\n", H_out);
fprintf("Residual = %.6e W\n", res_comb);

% PCFC Cell
H_in  = mDot_2p*h2p + mDot(11)*h11;
H_out = mDot(3)*h3 + mDot(12)*h12;

Q_in = Q_dot_PCFC;
W_out = W_PCFC;

res_pcfc = H_out - (H_in + Q_in - W_out);

fprintf("\n[PCFC]\n");
fprintf("Hin  = %.6e W\n", H_in);
fprintf("Hout = %.6e W\n", H_out);
fprintf("Qin  = %.6e W\n", Q_in);
fprintf("Wout = %.6e W\n", W_out);
fprintf("Residual = %.6e W\n", res_pcfc);

% HEX2
H_in  = mDot(10)*h10 + mDot(12)*h12;
H_out = mDot(11)*h11 + mDot(13)*h13;

res_hex2 = H_out - H_in;

fprintf("\n[HEX2]\n");
fprintf("Hin  = %.6e W\n", H_in);
fprintf("Hout = %.6e W\n", H_out);
fprintf("Residual = %.6e W\n", res_hex2);

% Compressor
H_in  = mDot(9)*h9;
H_out = mDot(10)*h10;
W_in  = W_comp;

res_comp = H_out - (H_in - W_in);

fprintf("\n[Compressor]\n");
fprintf("Hin  = %.6e W\n", H_in);
fprintf("Hout = %.6e W\n", H_out);
fprintf("Win  = %.6e W\n", W_in);
fprintf("Residual = %.6e W\n", res_comp);

% OVERALL SYSTEM
H_in  = mDot(1)*h1 + mDot(9)*h9;
H_out = mDot(7)*h7 + mDot(8)*h8 + mDot(13)*h13;

Q_in = Q_dot_PCFC;    % PCFC에서 흡수한 열
W_out = W_PCFC - W_comp;

res_sys = H_out - (H_in + Q_in - W_out);

fprintf("\n[Energy Check] OVERALL SYSTEM\n");
fprintf("Hin  = %.6e W\n", H_in);
fprintf("Hout = %.6e W\n", H_out);
fprintf("Qin  = %.6e W\n", Q_in);
fprintf("Wnet = %.6e W\n", W_out);
fprintf("Residual = %.6e W\n", res_sys);
fprintf("Relative error = %.3e %%\n", abs(res_sys)/max(abs(H_out),1)*100);

%%
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