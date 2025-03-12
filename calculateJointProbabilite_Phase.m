function [jointProbabiliteGamma_Phase]=calculateJointProbabilite_Phase(bit_Alice,bit_Bob,ScaleCo,C2n) 
    global bit_Alice_Glo;
    global bit_Bob_Glo;
    global ScaleCo_Glo;
    global C2n_Glo;
    
    bit_Alice_Glo=bit_Alice;
    bit_Bob_Glo=bit_Bob;
    ScaleCo_Glo=ScaleCo;
    C2n_Glo=C2n;

    tmp1=0;
    tmp2=0;
    
    for n=-10:1:10
        a=2*n*pi-0.5*pi;
        b=2*n*pi+0.5*pi;
        
        tmp1=tmp1+quadl('jointProbabiliteGamma_Direct_Phase',a,b);    
    end
    
    for n=-10:1:10
        a=2*n*pi-1.5*pi;
        b=2*n*pi-0.5*pi;
        tmp2=tmp2+quadl('jointProbabiliteGamma_Direct_Phase',a,b);
    end
    
    jointProbabiliteGamma_Phase=tmp1+tmp2;
end