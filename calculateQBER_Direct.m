function [QBER,P_sift]=calculateQBER_Direct(ScaleCo,C2n)
    P_sift=jointProbabiliteGamma_Direct(0,0,ScaleCo,C2n)+jointProbabiliteGamma_Direct(0,1,ScaleCo,C2n)...
          +jointProbabiliteGamma_Direct(1,0,ScaleCo,C2n)+jointProbabiliteGamma_Direct(1,1,ScaleCo,C2n);
    P_error=jointProbabiliteGamma_Direct(0,1,ScaleCo,C2n)+jointProbabiliteGamma_Direct(1,0,ScaleCo,C2n);
    
    QBER=P_error/P_sift;
end