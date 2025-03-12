function [QuantumBER_Direct_Weak,Psift_Direct_Weak,QuantumBER_Direct_Strong,Psift_Direct_Strong]=Simulation_HAP_Phase(C2n_Weak,C2n_Strong,ScaleCo_HAP_Phase)
    %Parameter simulation
    global P_dBm;                  %Peak transmitted power
    global Omega_z_P;  
    global Omega_z_G;
    global ModDepth;
    global B;
    global C2n_0;                  %The value of C2n at the ground level   
    global H_G;                    %The height of the ground node (m)
    global Omega;                  %Windspeed
    global G_TX_S_dB
    global G_RX_G_dB
    global G_TX_P_dB
    global G_RX_P_dB

    global G_A_dB;
    
    H_G=5; 
    Omega=21;
    
    global delta_fIF;
    global lengthBit;              %Bit length for Monte Carlo simulation
    
    
    %======================================================================
    %FSO channel 
    %Variance of background noise
    sigmab_P_2=4.435*10^(-28);  %Satellite-to-HAP
    sigmab_G_2=1.445*10^(-25);  %HAP-to-ground
    
    sigma1=0.43;                %Attenuation coefficient (dB/km)
    sigma=sigma1/4.343;  
    
    B0=125*10^9;                %Optical bandwidth

    G_A=10.^(G_A_dB./10); 
    
    
    %======================================================================
    %LEO Satellite (Alice)
    lamda=1550*10^-9;           %Wavelength
    H_S=610*10^3;               %LEO satellite altitude
    ZenithAngleDegree_S=50;     %Zenith angle
    ZenithAngle_S=ZenithAngleDegree_S./180.*pi;
    G_TX_S=10.^(G_TX_S_dB./10);
    deltaf=B;                   %Bandwidth
    
    
    %======================================================================
    %Relay (HAP)
    H_P=20*10^3;               %HAP altitude
    ZenithAngleDegree_P=50;    %Zenith angle
    ZenithAngle_P=ZenithAngleDegree_P./180.*pi;
    a_P=0.05;                  %The radius of the detection qperture
    n_sp=5;                    %ASE parameter

    G_TX_P=10.^(G_TX_P_dB./10);
    G_RX_P=10.^(G_RX_P_dB./10);
    
    
    %======================================================================
    %Bob (Vehicle or UAV)
    a_G=0.31;                  %The radius of the dection aperture
    G_RX_G=10.^(G_RX_G_dB./10);
    nguy=0.62;                 %Quantum efficiency
    kA=0.7;                    %Ionization factor
    M=10;                      %Avalanche Multiplication Factor
    FA=kA.*M+(2-1./M).*(1-kA); %Excess noise factor
    Fn=2;                      %Amplifier noise figure
    RL=1000;                   %Load resistance
    T=300;                     %Temperature
    
    P=(10.^(P_dBm./10)).*10.^-3;
    q=1.6*10^-19;              %Electron charge
    kB=1.38*10^-23;            %Boltzmann's constant
    h=6.626*10^-34;            %Planck's constant 
    c=3*10^8;                  %Speed of Light
    R=(nguy.*q)./(h.*c./lamda);%Detector responsivity 
    
    
    %======================================================================
    dataTransmitted_Weak=randi([0,1],1,lengthBit);

    QuantumBER_Direct_Weak=zeros(1,length(ScaleCo_HAP_Phase));
    Psift_Direct_Weak=zeros(1,length(ScaleCo_HAP_Phase));
    
    dataTransmitted_Strong=randi([0,1],1,lengthBit);

    QuantumBER_Direct_Strong=zeros(1,length(ScaleCo_HAP_Phase));
    Psift_Direct_Strong=zeros(1,length(ScaleCo_HAP_Phase));

    for i=1:length(ScaleCo_HAP_Phase)
        %Calculate QBER via Gamma-Gamma Channels 
        iTransmitted_Weak=zeros(1,lengthBit);
        iReceive_Weak=zeros(1,lengthBit);
        dataReceive_Weak=zeros(1,lengthBit);
    
        lengthSiftKey_Weak=0;
        ErrorSiftKey_Weak=0;
        
        iTransmitted_Strong=zeros(1,lengthBit);
        iReceive_Strong=zeros(1,lengthBit);
        dataReceive_Strong=zeros(1,lengthBit);
    
        lengthSiftKey_Strong=0;
        ErrorSiftKey_Strong=0;

        for j=1:lengthBit 
            %==============================================================
            %Calculate QBER
            %Free-space loss
            L_S=(H_S-H_P)./cos(ZenithAngle_S);
            FSL=(4.*pi.*L_S./lamda).^2;


            %==============================================================
            %Path loss
            L_P=H_P./cos(ZenithAngle_P);
            hl=exp(-sigma.*L_P./1000);  


            %==============================================================
            %The fraction of the power collected by the detector at HAP
            r=0;
            v_P=sqrt(pi).*a_P./(sqrt(2).*Omega_z_P);
            A_0_P=(erf(v_P)).^2;
            Omega_zeq_P_2=Omega_z_P.^2.*(sqrt(pi).*erf(v_P))./(2.*v_P.*exp(-v_P.^2));
            h_p_P=A_0_P.*exp(-2.*(r.^2)./Omega_zeq_P_2);
           


            %==============================================================
            %The fraction of the power collected by the detector at ground station
            v_G=sqrt(pi).*a_G./(sqrt(2).*Omega_z_G);
            A_0_G=(erf(v_G)).^2;
            Omega_zeq_G_2=Omega_z_G.^2.*(sqrt(pi).*erf(v_G))./(2.*v_G.*exp(-v_G.^2));
            h_p_G=A_0_G.*exp(-2.*(r.^2)./Omega_zeq_G_2);


            %==============================================================
            %Refrective index structure coefficient
            k=2.*pi./lamda;
            
            tmp=0;
            
            %Weak condition
            if tmp==0 
                C2n_0=C2n_Weak;
                integral_Weak=quad(@calculateSigma_R_2,H_G,H_P);
                sigma_R_2_Weak=2.25.*k.^(7/6).*(sec(ZenithAngle_P)).^(11./6).*integral_Weak;
                tmp=1;
            end
            
            %Strong condition
            if tmp==1 
                C2n_0=C2n_Strong;
                integral_Strong=quad(@calculateSigma_R_2,H_G,H_P);
                sigma_R_2_Strong=2.25.*k.^(7/6).*(sec(ZenithAngle_P)).^(11./6).*integral_Strong;
            end
            

            %==============================================================
            %Received background light power at HAP and GS
            N_0_P=2.*sigmab_P_2;
            N_0_G=2.*sigmab_G_2;

            P_b_P=N_0_P.*B0;
            P_b_G=N_0_G.*B0;


            %==============================================================
            %ASE noise
            P_a=2.*n_sp.*h.*c./lamda.*B0;
            
            
            %==============================================================
            %Peak received power at Bob
            P_r_G=1./FSL.*P.*G_TX_S   .*h_p_P.*G_RX_P.*G_A.*G_TX_P   .*h_p_G.*G_RX_G.*hl;

         
            %==============================================================
            %Phase noise
            f_s=B; %The statistical standard deviation of the received signal frequency
            sigma_teta_e_2=2.*pi.*delta_fIF./f_s;
            teta_e=normrnd(0,sqrt(sigma_teta_e_2));
            
            
            %==============================================================
            %Calculate gamma-gamma variable
           
            %=============================
            %WEAK CONDITION
            
            alpha_Weak=(exp(0.49*sigma_R_2_Weak/(1+1.11*(sqrt(sigma_R_2_Weak))^(12/5))^(7/6))-1)^(-1);
            beta_Weak=(exp(0.51*sigma_R_2_Weak/(1+0.69*(sqrt(sigma_R_2_Weak))^(12/5))^(5/6))-1)^(-1);

            g1_Weak=gamrnd(alpha_Weak,1/alpha_Weak);
            g2_Weak=gamrnd(beta_Weak,1/beta_Weak);
            h_t_Weak=g1_Weak*g2_Weak;
            
            
            %==============================================================
            %Variance  
            sigmaN_2_Weak=2*q*FA*(M^2)*R*(1/4*R*M*P_r_G*ModDepth*h_t_Weak+P_b_P*G_RX_P*...
            G_A*G_TX_P*h_p_G*G_RX_G*hl*h_t_Weak+P_a*G_RX_P*G_A*G_TX_P*h_p_G*G_RX_G*hl*h_t_Weak+...
            P_b_G*G_RX_G)*deltaf+4*kB*T*Fn*deltaf/RL;

            
            %Calculate n(t) (received noise)
            n_t_Weak=normrnd(0,sqrt(sigmaN_2_Weak));
         
            if dataTransmitted_Weak(j)==0
                iTransmitted_Weak(j)=-1./4.*R.*M.*P_r_G.*ModDepth;
                
            else
                iTransmitted_Weak(j)=1./4.*R.*M.*P_r_G.*ModDepth;
            end
            
            
            %==============================================================
            %Threshold
            d0_Weak=-1./4.*R.*M.*P_r_G.*ModDepth-ScaleCo_HAP_Phase(i)*sqrt(sigmaN_2_Weak);
            d1_Weak=1./4.*R.*M.*P_r_G.*ModDepth+ScaleCo_HAP_Phase(i)*sqrt(sigmaN_2_Weak);
            
            
            %==============================================================
            %Calculate QBER
            iReceive_Weak(j)=iTransmitted_Weak(j)*cos(teta_e)*h_t_Weak+n_t_Weak;

            if iReceive_Weak(j)<=d0_Weak 
                dataReceive_Weak(j)=0;
                lengthSiftKey_Weak=lengthSiftKey_Weak+1;
            elseif iReceive_Weak(j)>=d1_Weak
                dataReceive_Weak(j)=1;
                lengthSiftKey_Weak=lengthSiftKey_Weak+1;
            else
                dataReceive_Weak(j)=2;
            end
            
            if dataReceive_Weak(j)~=2 
                if dataReceive_Weak(j)~=dataTransmitted_Weak(j)
                    ErrorSiftKey_Weak=ErrorSiftKey_Weak+1;
                end
            end 
            
            %=============================
            %STRONG CONDITION
            alpha_Strong=(exp(0.49*sigma_R_2_Strong/(1+1.11*(sqrt(sigma_R_2_Strong))^(12/5))^(7/6))-1)^(-1);
            beta_Strong=(exp(0.51*sigma_R_2_Strong/(1+0.69*(sqrt(sigma_R_2_Strong))^(12/5))^(5/6))-1)^(-1);

            g1_Strong=gamrnd(alpha_Strong,1/alpha_Strong);
            g2_Strong=gamrnd(beta_Strong,1/beta_Strong);
            h_t_Strong=g1_Strong*g2_Strong;
            
            
            %==============================================================
            %Variance  
            sigmaN_2_Strong=2*q*FA*(M^2)*R*(1/4*R*M*P_r_G*ModDepth*h_t_Strong+P_b_P*G_RX_P*...
            G_A*G_TX_P*h_p_G*G_RX_G*hl*h_t_Strong+P_a*G_RX_P*G_A*G_TX_P*h_p_G*G_RX_G*hl*h_t_Strong+...
            P_b_G*G_RX_G)*deltaf+4*kB*T*Fn*deltaf/RL;
                 
            
            %Calculate n(t) (received noise)
            n_t_Strong=normrnd(0,sqrt(sigmaN_2_Strong));
         
            if dataTransmitted_Strong(j)==0
                iTransmitted_Strong(j)=-1./4.*R.*M.*P_r_G.*ModDepth;
                
            else
                iTransmitted_Strong(j)=1./4.*R.*M.*P_r_G.*ModDepth;
            end
            
            
            %==============================================================
            %Threshold
            d0_Strong=-1./4.*R.*M.*P_r_G.*ModDepth-ScaleCo_HAP_Phase(i)*sqrt(sigmaN_2_Strong);
            d1_Strong=1./4.*R.*M.*P_r_G.*ModDepth+ScaleCo_HAP_Phase(i)*sqrt(sigmaN_2_Strong);
            
            
            %==============================================================
            %Calculate QBER
            iReceive_Strong(j)=iTransmitted_Strong(j)*cos(teta_e)*h_t_Strong+n_t_Strong;

            if iReceive_Strong(j)<=d0_Strong 
                dataReceive_Strong(j)=0;
                lengthSiftKey_Strong=lengthSiftKey_Strong+1;
            elseif iReceive_Strong(j)>=d1_Strong
                dataReceive_Strong(j)=1;
                lengthSiftKey_Strong=lengthSiftKey_Strong+1;
            else
                dataReceive_Strong(j)=2;
            end
            
            if dataReceive_Strong(j)~=2 
                if dataReceive_Strong(j)~=dataTransmitted_Strong(j)
                    ErrorSiftKey_Strong=ErrorSiftKey_Strong+1;
                end
            end 
        end
    
        
        %=============================
        %WEAK CONDITION
        QuantumBER_Direct_Weak(i)=ErrorSiftKey_Weak/lengthSiftKey_Weak;
        Psift_Direct_Weak(i)=lengthSiftKey_Weak/lengthBit;
       

        %=============================
        %STRONG CONDITION
        QuantumBER_Direct_Strong(i)=ErrorSiftKey_Strong/lengthSiftKey_Strong;
        Psift_Direct_Strong(i)=lengthSiftKey_Strong/lengthBit;
    end
end