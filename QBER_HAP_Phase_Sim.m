function [ScaleCo_HAP_Phase_Sim,QuantumBER_HAP_Phase_Sim_Weak,Psift_HAP_Phase_Sim_Weak,QuantumBER_HAP_Phase_Sim_Strong,Psift_HAP_Phase_Sim_Strong]=QBER_HAP_Phase_Sim()
    %Simulation Parameters
    C2n_Weak=5*10^-17;  %Refractive index structure coefficient
    C2n_Strong=7*10^-9;

    ScaleCo_HAP_Phase_Sim=0:0.5:2.5;

    %Weak and strong turbulence
    [QuantumBER_HAP_Phase_Sim_Weak,Psift_HAP_Phase_Sim_Weak,QuantumBER_HAP_Phase_Sim_Strong,Psift_HAP_Phase_Sim_Strong]=Simulation_HAP_Phase(C2n_Weak,C2n_Strong,ScaleCo_HAP_Phase_Sim);
end