function [QBER,P_sift]=calculateQBER_HAP_Phase(ScaleCo,C2n)
    P_sift=calculateJointProbabilite_HAP_Phase(0,0,ScaleCo,C2n)+calculateJointProbabilite_HAP_Phase(0,1,ScaleCo,C2n)...
          +calculateJointProbabilite_HAP_Phase(1,0,ScaleCo,C2n)+calculateJointProbabilite_HAP_Phase(1,1,ScaleCo,C2n);
    P_error=calculateJointProbabilite_HAP_Phase(0,1,ScaleCo,C2n)+calculateJointProbabilite_HAP_Phase(1,0,ScaleCo,C2n);
    
    QBER=P_error/P_sift;
end