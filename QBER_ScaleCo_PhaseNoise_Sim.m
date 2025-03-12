clear;
clc;


%==========================================================================
%Simulation Parameters
global P_dBm;       %Peak transmitted power
global Omega_z_P; 
global Omega_z_G;
global B            %Bit rate
global G_A_dB;      %The gain of the optical amplifier
global G_TX_S_dB    %TX Telescope Gain at Satellite
global G_RX_G_dB    %RX Telescope Gain at Ground Station
global G_TX_P_dB    %TX Telescope Gain at HAP
global G_RX_P_dB    %RX Telescope Gain at HAP
global ModDepth;    %Modulation depth
global lengthBit;   %Bit length for Monte Carlo simulation
global delta_fIF;   %Received signal frequency

P_dBm=29;

Omega_z_P=1475;
Omega_z_G=50;
B=10^9;

G_TX_S_dB=141;              
G_RX_G_dB=132;          

G_TX_P_dB=98;               
G_RX_P_dB=100;              

G_A_dB=0;                      
ModDepth=0.45;

lengthBit=4*10^7;              
delta_fIF=1.7*10^7;         

C2n_Weak=5*10^-17;  %Refractive index structure coefficients
C2n_Strong=7*10^-10;

ScaleCo=0:0.5:3;    %D-T scale coefficient


% Note: G_A_dB=1, G_TX_P_dB=98 -> Strong condition needs to be fitted
% G_A_dB=0, G_TX_P_dB=98 -> Better

%==========================================================================
%Calculate QBER via Gamma-Gamma Channels with weak and strong turbulence
%Weak condition
QBER_Direct_Weak=zeros(1,length(ScaleCo));
P_sift_Direct_Weak=zeros(1,length(ScaleCo));

QBER_Direct_Phase_Weak=zeros(1,length(ScaleCo));
P_sift_Direct_Phase_Weak=zeros(1,length(ScaleCo));

% QBER_HAP_Phase_Weak=zeros(1,length(ScaleCo));
% P_sift_HAP_Phase_Weak=zeros(1,length(ScaleCo));

for i=1:length(ScaleCo)
    [QBER_Direct_Weak(i),P_sift_Direct_Weak(i)]=calculateQBER_Direct(ScaleCo(i),C2n_Weak);
    [QBER_Direct_Phase_Weak(i),P_sift_Direct_Phase_Weak(i)]=calculateQBER_Direct_Phase(ScaleCo(i),C2n_Weak);
    % [QBER_HAP_Phase_Weak(i),P_sift_HAP_Phase_Weak(i)]=calculateQBER_HAP_Phase(ScaleCo(i),C2n_Weak);
end


%==========================================================================
%Strong condition
QBER_Direct_Strong=zeros(1,length(ScaleCo));
P_sift_Direct_Strong=zeros(1,length(ScaleCo));

QBER_Direct_Phase_Strong=zeros(1,length(ScaleCo));
P_sift_Direct_Phase_Strong=zeros(1,length(ScaleCo));

% QBER_HAP_Phase_Strong=zeros(1,length(ScaleCo));
% P_sift_HAP_Phase_Strong=zeros(1,length(ScaleCo));

for i=1:length(ScaleCo)
    [QBER_Direct_Strong(i),P_sift_Direct_Strong(i)]=calculateQBER_Direct(ScaleCo(i),C2n_Strong);
    [QBER_Direct_Phase_Strong(i),P_sift_Direct_Phase_Strong(i)]=calculateQBER_Direct_Phase(ScaleCo(i),C2n_Strong);
    % [QBER_HAP_Phase_Strong(i),P_sift_HAP_Phase_Strong(i)]=calculateQBER_HAP_Phase(ScaleCo(i),C2n_Strong);
end


%==========================================================================
%Simulation
[ScaleCo_Direct_Sim,QuantumBER_Direct_Sim_Weak,Psift_Direct_Sim_Weak,QuantumBER_Direct_Sim_Strong,Psift_Direct_Sim_Strong]=QBER_Direct_Sim();
[ScaleCo_Direct_Phase_Sim,QuantumBER_Direct_Phase_Sim_Weak,Psift_Direct_Phase_Sim_Weak,QuantumBER_Direct_Phase_Sim_Strong,Psift_Direct_Phase_Sim_Strong]=QBER_Direct_Phase_Sim();
% [ScaleCo_HAP_Phase_Sim,QuantumBER_HAP_Phase_Sim_Weak,Psift_HAP_Phase_Sim_Weak,QuantumBER_HAP_Phase_Sim_Strong,Psift_HAP_Phase_Sim_Strong]=QBER_HAP_Phase_Sim();


%==========================================================================
%Plot funciton of Gamma-Gamma channels with weak and strong turbulence
figure(1)
subplot(1,2,1)
semilogy(ScaleCo,QBER_Direct_Weak,'r-+',ScaleCo,P_sift_Direct_Weak,'b-+','LineWidth',2);
hold on 
semilogy(ScaleCo,QBER_Direct_Phase_Weak,'r-o',ScaleCo,P_sift_Direct_Phase_Weak,'b-o','LineWidth',2);
hold off
% hold on
% semilogy(ScaleCo,QBER_HAP_Phase_Weak,'r-*',ScaleCo,P_sift_HAP_Phase_Weak,'b-*','LineWidth',2); 
% hold off

hold on
semilogy(ScaleCo_Direct_Sim,QuantumBER_Direct_Sim_Weak,'ko',ScaleCo_Direct_Sim,Psift_Direct_Sim_Weak,'ko','LineWidth',2);
hold off
hold on
semilogy(ScaleCo_Direct_Phase_Sim,QuantumBER_Direct_Phase_Sim_Weak,'ko',ScaleCo_Direct_Phase_Sim,Psift_Direct_Phase_Sim_Weak,'ko','LineWidth',2);
hold off
% hold on
% semilogy(ScaleCo_HAP_Phase_Sim,QuantumBER_HAP_Phase_Sim_Weak,'ko',ScaleCo_HAP_Phase_Sim,Psift_HAP_Phase_Sim_Weak,'ko','LineWidth',2);
% hold off

grid on
xlabel('D-T scale coefficient, \xi');
ylabel('Probability');
% legend('QBER','P_{sift}','QBER_{(phase fluctuation)}','P_{sift (phase fluctuation)}','QBER_{(phase fluctuation + HAP)}','P_{sift (phase fluctuation + HAP)}','M-C Simulation','Location','southwest');
legend('QBER','P_{sift}','QBER_{(phase fluctuation)}','P_{sift (phase fluctuation)}','M-C Simulation','Location','southwest');
axis([0,3,1.e-5,1.e-0]);


subplot(1,2,2)
semilogy(ScaleCo,QBER_Direct_Strong,'r-+',ScaleCo,P_sift_Direct_Strong,'b-+','LineWidth',2);
hold on
semilogy(ScaleCo,QBER_Direct_Phase_Strong,'r-o',ScaleCo,P_sift_Direct_Phase_Strong,'b-o','LineWidth',2);
hold off
% hold on 
% semilogy(ScaleCo,QBER_HAP_Phase_Strong,'r-*',ScaleCo,P_sift_HAP_Phase_Strong,'b-*','LineWidth',2);
% hold off

hold on
semilogy(ScaleCo_Direct_Sim,QuantumBER_Direct_Sim_Strong,'ko',ScaleCo_Direct_Sim,Psift_Direct_Sim_Strong,'ko','LineWidth',2);
hold off
hold on
semilogy(ScaleCo_Direct_Phase_Sim,QuantumBER_Direct_Phase_Sim_Strong,'ko',ScaleCo_Direct_Phase_Sim,Psift_Direct_Phase_Sim_Strong,'ko','LineWidth',2);
hold off
% hold on
% semilogy(ScaleCo_HAP_Phase_Sim,QuantumBER_HAP_Phase_Sim_Strong,'ko',ScaleCo_HAP_Phase_Sim,Psift_HAP_Phase_Sim_Strong,'ko','LineWidth',2);
% hold off

grid on
xlabel('D-T scale coefficient, \xi');
ylabel('Probability');
% legend('QBER','P_{sift}','QBER_{(phase fluctuation)}','P_{sift (phase fluctuation)}','QBER_{(phase fluctuation + HAP)}','P_{sift (phase fluctuation + HAP)}','M-C Simulation','Location','southwest');
legend('QBER','P_{sift}','QBER_{(phase fluctuation)}','P_{sift (phase fluctuation)}','M-C Simulation','Location','southwest');
axis([0,3,1.e-5,1.e-0]);