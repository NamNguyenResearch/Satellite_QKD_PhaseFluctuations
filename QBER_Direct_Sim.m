function [ScaleCo_Direct_Sim,QuantumBER_Direct_Sim_Weak,Psift_Direct_Sim_Weak,QuantumBER_Direct_Sim_Strong,Psift_Direct_Sim_Strong]=QBER_Direct_Sim()
    %Simulation Parameters
    C2n_Weak=5*10^-17;  %Refractive index structure coefficient
    C2n_Strong=7*10^-10;

    ScaleCo_Direct_Sim=0:0.5:2.5;

    %Weak and strong turbulence
    [QuantumBER_Direct_Sim_Weak,Psift_Direct_Sim_Weak,QuantumBER_Direct_Sim_Strong,Psift_Direct_Sim_Strong]=Simulation_Direct(C2n_Weak,C2n_Strong,ScaleCo_Direct_Sim);
end

