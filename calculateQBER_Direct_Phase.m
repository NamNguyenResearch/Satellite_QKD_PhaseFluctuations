function [QBER_Gamma_C2n_Phase,P_sift_Gamma_C2n_Phase]=calculateQBER_Direct_Phase(ScaleCo,C2n)
    P_sift_Gamma_C2n_Phase=calculateJointProbabilite_Phase(0,0,ScaleCo,C2n)+calculateJointProbabilite_Phase(0,1,ScaleCo,C2n)...
                     +calculateJointProbabilite_Phase(1,0,ScaleCo,C2n)+calculateJointProbabilite_Phase(1,1,ScaleCo,C2n);
    P_error_Gamma_C2n_Phase=calculateJointProbabilite_Phase(0,1,ScaleCo,C2n)+calculateJointProbabilite_Phase(1,0,ScaleCo,C2n);
    QBER_Gamma_C2n_Phase=P_error_Gamma_C2n_Phase/P_sift_Gamma_C2n_Phase;
end