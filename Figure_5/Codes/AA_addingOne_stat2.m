

% differentiate

%dr = find(Stim_order3>= 41);
%Stim_order3(dr,:) = [];
%VectoS3(dr,:) =[];
%Vecto_FR3(dr,:)=[];
%%
%winp = 10;
%Alp_FR = D_Alph_FR(winp,:)';
%P_time = Alp_FR;
%save('P_time_.mat','P_time')

%%
%winp = 10;
%Del_FR =  D_Delta_FR(winp,:)';
%P_time = Del_FR;
%save('P_time_.mat','P_time')

%%
%
VectoS3w = VectoS3;
Vecto_FR3w = Vecto_FR3;
%
% differentiate
pttFR = size(Vecto_FR3w,1);

% Plotting Spectrogram for each 
Fg = 120;
NumS = 10
Wind = 1000

params.Fs = 1000;
params.fpass = [0 Fg]; 
params.tapers = [10 19];
params.trialave = 0;
params.pad=1;
params.err = 1;

ptt = size(Vecto_FR3w,1);
pt = size(Vecto_FR3w,2);
%frk = 30001;

for p = 2:pt
if (VectoS3w(1,p-1) <= 1) && (VectoS3w(1,p)>=2);
    frk = p 
end
end


%

for con = 5:5:NumS
    
    
    % cutting required parts of vector
    
    
    for s = 1:pttFR
                Vecto_FR3t(s,:) = Vecto_FR3w(s,(frk-con*Wind):(frk+con*Wind)-1);
    end

                
         % Collecting data for power matrix       
                
               tsec = size(Vecto_FR3t,2);
                tbef = 1:tsec/2;
                taft = tsec/2+1:tsec-1;  
        
               %detecting spectrogram
                

                
        for sx = 1:pttFR
                [FRb,f] = mtspectrumc( Vecto_FR3t(sx,tbef), params );
                [FRa,f] = mtspectrumc( Vecto_FR3t(sx,taft), params );
                V_FRb(sx,:) = 10*log10(FRb);
                V_FRa(sx,:) = 10*log10(FRa);
        end

                % f / X Hz
                points_f = round(length(f)/Fg);
                Delta_p = (1*points_f:4*points_f);
                Theta_p = (5*points_f:8*points_f);
                Spindle_p = (9*points_f:16*points_f);
                Beta_p = (16*points_f:29*points_f);
                Gamma_p = (30*points_f:50*points_f);
                GammaX_p = (50*points_f:80*points_f);
                GammaXX_p = (110*points_f:119*points_f);
                
                % Calculate Ch
                
              
              for fr = 1:pttFR
                Delt_FRb(con,fr) = mean(V_FRb(fr,Delta_p));
                Delt_FRa(con,fr) = mean(V_FRa(fr,Delta_p));
                Thet_FRb(con,fr) = mean(V_FRb(fr,Theta_p));
                Thet_FRa(con,fr) = mean(V_FRa(fr,Theta_p));
                Alph_FRb(con,fr) = mean(V_FRb(fr,Spindle_p));
                Alph_FRa(con,fr) = mean(V_FRa(fr,Spindle_p));
                Beta_FRb(con,fr) = mean(V_FRb(fr,Beta_p));
                Beta_FRa(con,fr) = mean(V_FRa(fr,Beta_p));
                Gamma_FRb(con,fr) = mean(V_FRb(fr,Gamma_p));
                Gamma_FRa(con,fr) = mean(V_FRa(fr,Gamma_p));
                GammaX_FRb(con,fr) = mean(V_FRb(fr,GammaX_p));
                GammaX_FRa(con,fr) = mean(V_FRa(fr,GammaX_p));
                GammaXX_FRb(con,fr) = mean(V_FRb(fr,GammaXX_p));
                GammaXX_FRa(con,fr) = mean(V_FRa(fr,GammaXX_p));
              end
              
               % Calcuate diff
                
                  D_Delta_FR(con,:) = Delt_FRa(con,:)-Delt_FRb(con,:);
                  D_Theta_FR(con,:) = Thet_FRa(con,:)-Thet_FRb(con,:);
                  D_Alph_FR(con,:) = Alph_FRa(con,:)-Alph_FRb(con,:);
                  D_Beta_FR(con,:) = Beta_FRa(con,:)-Beta_FRb(con,:);
                  D_Gamma_FR(con,:) = Gamma_FRa(con,:)-Gamma_FRb(con,:);
                  D_GammaX_FR(con,:) = GammaX_FRa(con,:)-GammaX_FRb(con,:);
                  D_GammaXX_FR(con,:) = GammaXX_FRa(con,:)-GammaXX_FRb(con,:);
   
                % zerowing stacks
                
  
                Vecto_FR3t = [];
                VectoS3t = [];

                V_FRb = [];

                V_FRa = [];
    
end
%
piq = 5;

FR_D1 = D_Delta_FR(piq,:);
FR_T1 = D_Theta_FR(piq,:);
FR_A1 = D_Alph_FR(piq,:);
FR_B1 = D_Beta_FR(piq,:);
FR_G1 = D_Gamma_FR(piq,:);
FR_GX1 = D_GammaX_FR(piq,:);

%
%
ds =(5:5:NumS)/5;

FR_D = mean(D_Delta_FR,2);
FR_T = mean(D_Theta_FR,2);
FR_A = mean(D_Alph_FR,2);
FR_B = mean(D_Beta_FR,2);
FR_G = mean(D_Gamma_FR,2);
FR_GX = mean(D_GammaX_FR,2);
FR_GXX = mean(D_GammaXX_FR,2);


%FR
d_FR = FR_D(5:5:NumS);
t_FR = FR_T(5:5:NumS);
a_FR = FR_A(5:5:NumS);
b_FR = FR_B(5:5:NumS);
g_FR = FR_G(5:5:NumS);
gx_FR = FR_GX(5:5:NumS);
gxx_FR = FR_GXX(5:5:NumS);



%FR

AC_BFR_Pfx = [];
AC_BFR_Pfx(1,ds) = d_FR;
AC_BFR_Pfx(2,ds) = t_FR;
AC_BFR_Pfx(3,ds) = a_FR;
AC_BFR_Pfx(4,ds) = b_FR;
AC_BFR_Pfx(5,ds) = g_FR;
AC_BFR_Pfx(6,ds) = gx_FR;
AC_BFR_Pfx(7,ds) = gxx_FR;


printmat(AC_BFR_Pfx , 'BandsFR', 'Delta Theta Alpha Beta Gamma GammaX GammaXX', '5sec 10sec 15sec 20sec 25sec 30sec')
%title('Difference before and after stimulation')
%
%

%AR_FR

winp = 10;

%

%Del_ARb - 1
%Del_ARa - 2


Del_FRf(:,1) =  Delt_FRb(winp,:)';
Del_FRf(:,2) =  Delt_FRa(winp,:)';


The_FRf(:,1) = Thet_FRb(winp,:)';
The_FRf(:,2) = Thet_FRa(winp,:)';


Alp_FRf(:,1) = Alph_FRb(winp,:)';
Alp_FRf(:,2) = Alph_FRa(winp,:)';


Bet_FRf(:,1) = Beta_FRb(winp,:)';
Bet_FRf(:,2) = Beta_FRa(winp,:)';

Gam_FRf(:,1) = Gamma_FRb(winp,:)';
Gam_FRf(:,2) = Gamma_FRa(winp,:)';

GamX_FRf(:,1) = GammaX_FRb(winp,:);
GamX_FRf(:,2) = GammaX_FRa(winp,:);

%

sec = 5;
Nit = sqrt(size(FR_D1,2));

FR_D1_sem = std(D_Delta_FR(sec,:))/Nit;
FR_T1_sem = std(D_Theta_FR(sec,:))/Nit;
FR_A1_sem = std(D_Alph_FR(sec,:))/Nit;
FR_B1_sem = std(D_Beta_FR(sec,:))/Nit;
FR_G1_sem = std(D_Gamma_FR(sec,:))/Nit;
FR_GX1_sem = std(D_GammaX_FR(sec,:))/Nit;

%FR_D1_sem
%FR_A1_sem
%FR_B1_sem



%% All 
der = FR_A1'

%% Arch caudal 
Dc_ArcA = Del_AR;
Tc_ArcA = The_AR;
Ac_ArcA = Alp_AR;
Bc_ArcA = Bet_AR;
Gc_ArcA = Gam_AR;
GXc_ArcA = GamX_AR;
GXXc_ArcA = GamXX_AR;

Dc_ArcF = Del_FR;
Tc_ArcF = The_FR;
Ac_ArcF = Alp_FR;
Bc_ArcF = Bet_FR;
Gc_ArcF = Gam_FR;
GXc_ArcF = GamX_FR;
GXXc_ArcF = GamXX_FR;

%% Arch rostral
Dr_ArcA = Del_AR;
Tr_ArcA = The_AR;
Ar_ArcA = Alp_AR;
Br_ArcA = Bet_AR;
Gr_ArcA = Gam_AR;
GXr_ArcA = GamX_AR;
GXXr_ArcA = GamXX_AR;

Dr_ArcF = Del_FR;
Tr_ArcF = The_FR;
Ar_ArcF = Alp_FR;
Br_ArcF = Bet_FR;
Gr_ArcF = Gam_FR;
GXr_ArcF = GamX_FR;
GXXr_ArcF = GamXX_FR;

%% Halo caudal
Dc_HalA = Del_AR;
Tc_HalA = The_AR;
Ac_HalA = Alp_AR;
Bc_HalA = Bet_AR;
Gc_HalA = Gam_AR;
GXc_HalA = GamX_AR;
GXXc_HalA = GamXX_AR;

Dc_HalF = Del_FR;
Tc_HalF = The_FR;
Ac_HalF = Alp_FR;
Bc_HalF = Bet_FR;
Gc_HalF = Gam_FR;
GXc_HalF = GamX_FR;
GXXc_HalF = GamXX_FR;

%% Halo rostral
Dr_HalA = Del_AR;
Tr_HalA = The_AR;
Ar_HalA = Alp_AR;
Br_HalA = Bet_AR;
Gr_HalA = Gam_AR;
GXr_HalA = GamX_AR;
GXXr_HalA = GamXX_AR;

Dr_HalF = Del_FR;
Tr_HalF = The_FR;
Ar_HalF = Alp_FR;
Br_HalF = Bet_FR;
Gr_HalF = Gam_FR;
GXr_HalF = GamX_FR;
GXXr_HalF = GamXX_FR;

%% Control caudal
Dc_ConA = Del_AR;
Tc_ConA = The_AR;
Ac_ConA = Alp_AR;
Bc_ConA = Bet_AR;
Gc_ConA = Gam_AR;
GXc_ConA = GamX_AR;
GXXc_ConA = GamXX_AR;

Dc_ConF = Del_FR;
Tc_ConF = The_FR;
Ac_ConF = Alp_FR;
Bc_ConF = Bet_FR;
Gc_ConF = Gam_FR;
GXc_ConF = GamX_FR;
GXXc_ConF = GamXX_FR;

%% Halo rostral
Dr_ConA = Del_AR;
Tr_ConA = The_AR;
Ar_ConA = Alp_AR;
Br_ConA = Bet_AR;
Gr_ConA = Gam_AR;
GXr_ConA = GamX_AR;
GXXr_ConA = GamXX_AR;

Dr_ConF = Del_FR;
Tr_ConF = The_FR;
Ar_ConF = Alp_FR;
Br_ConF = Bet_FR;
Gr_ConF = Gam_FR;
GXr_ConF = GamX_FR;
GXXr_ConF = GamXX_FR;

%% Lets do anova
%%
% Data


% Delta

%Big1 = [Dc_ArcA Dc_ArcF Dr_ArcA Dr_ArcF Dc_HalA Dc_HalF Dr_HalA Dr_HalF ...
%    Dc_ConA Dc_ConF Dr_ConA Dr_ConF];

% Theta

%Big1 = [Tc_ArcA Tc_ArcF Tr_ArcA Tr_ArcF Tc_HalA Tc_HalF Tr_HalA Tr_HalF ...
%    Tc_ConA Tc_ConF Tr_ConA Tr_ConF];

% Alpha

%Big1 = [Ac_ArcA Ac_ArcF Ar_ArcA Ar_ArcF Ac_HalA Ac_HalF Ar_HalA Ar_HalF ...
  %  Ac_ConA Ac_ConF Ar_ConA Ar_ConF];

% Beta

%Big1 = [Bc_ArcA Bc_ArcF Br_ArcA Br_ArcF Bc_HalA Bc_HalF Br_HalA Br_HalF ...
 %   Bc_ConA Bc_ConF Br_ConA Br_ConF];

% Gamma

Big1 = [Gc_ArcA Gc_ArcF Gr_ArcA Gr_ArcF Gc_HalA Gc_HalF Gr_HalA Gr_HalF ...
    Gc_ConA Gc_ConF Gr_ConA Gr_ConF];

% GammaX

%Big1 = [GXc_ArcA GXc_ArcF GXr_ArcA GXr_ArcF GXc_HalA GXc_HalF GXr_HalA GXr_HalF ...
 %  GXc_ConA GXc_ConF GXr_ConA GXr_ConF];

% GammaXX

%Big1 = [GXXc_ArcA GXXc_ArcF GXXr_ArcA GXXr_ArcF GXXc_HalA GXXc_HalF GXXr_HalA GXXr_HalF ...
 %   GXXc_ConA GXXc_ConF GXXr_ConA GXXr_ConF];

%% Virus

% Delta

Dc_ArcA_n = ones(1,size(Dc_ArcA,2));
Dc_ArcF_n = ones(1,size(Dc_ArcF,2));
Dr_ArcA_n = ones(1,size(Dr_ArcA,2));
Dr_ArcF_n = ones(1,size(Dr_ArcF,2));

Dc_HalA_n = ones(1,size(Dc_HalA,2))*2;
Dc_HalF_n = ones(1,size(Dc_HalF,2))*2;
Dr_HalA_n = ones(1,size(Dr_HalA,2))*2;
Dr_HalF_n = ones(1,size(Dr_HalF,2))*2;

Dc_ConA_n = ones(1,size(Dc_ConA,2))*3;
Dc_ConF_n = ones(1,size(Dc_ConF,2))*3;
Dr_ConA_n = ones(1,size(Dr_ConA,2))*3;
Dr_ConF_n = ones(1,size(Dr_ConF,2))*3;


Big2 = [Dc_ArcA_n Dc_ArcF_n Dr_ArcA_n Dr_ArcF_n Dc_HalA_n Dc_HalF_n Dr_HalA_n Dr_HalF_n ...
    Dc_ConA_n Dc_ConF_n Dr_ConA_n Dr_ConF_n];

%%
% Theta

Tc_ArcA_n = ones(1,size(Tc_ArcA,2));
Tc_ArcF_n = ones(1,size(Tc_ArcF,2));
Tr_ArcA_n = ones(1,size(Tr_ArcA,2));
Tr_ArcF_n = ones(1,size(Tr_ArcF,2));

Tc_HalA_n = ones(1,size(Tc_HalA,2))*2;
Tc_HalF_n = ones(1,size(Tc_HalF,2))*2;
Tr_HalA_n = ones(1,size(Tr_HalA,2))*2;
Tr_HalF_n = ones(1,size(Tr_HalF,2))*2;

Tc_ConA_n = ones(1,size(Tc_ConA,2))*3;
Tc_ConF_n = ones(1,size(Tc_ConF,2))*3;
Tr_ConA_n = ones(1,size(Tr_ConA,2))*3;
Tr_ConF_n = ones(1,size(Tr_ConF,2))*3;


Big2 = [Tc_ArcA_n Tc_ArcF_n Tr_ArcA_n Tr_ArcF_n Tc_HalA_n Tc_HalF_n Tr_HalA_n Tr_HalF_n ...
    Tc_ConA_n Tc_ConF_n Tr_ConA_n Tr_ConF_n];

%%
% Alpha

Ac_ArcA_n = ones(1,size(Ac_ArcA,2));
Ac_ArcF_n = ones(1,size(Ac_ArcF,2));
Ar_ArcA_n = ones(1,size(Ar_ArcA,2));
Ar_ArcF_n = ones(1,size(Ar_ArcF,2));

Ac_HalA_n = ones(1,size(Ac_HalA,2))*2;
Ac_HalF_n = ones(1,size(Ac_HalF,2))*2;
Ar_HalA_n = ones(1,size(Ar_HalA,2))*2;
Ar_HalF_n = ones(1,size(Ar_HalF,2))*2;

Ac_ConA_n = ones(1,size(Ac_ConA,2))*3;
Ac_ConF_n = ones(1,size(Ac_ConF,2))*3;
Ar_ConA_n = ones(1,size(Ar_ConA,2))*3;
Ar_ConF_n = ones(1,size(Ar_ConF,2))*3;


Big2 = [Ac_ArcA_n Ac_ArcF_n Ar_ArcA_n Ar_ArcF_n Ac_HalA_n Ac_HalF_n Ar_HalA_n Ar_HalF_n ...
    Ac_ConA_n Ac_ConF_n Ar_ConA_n Ar_ConF_n];
%%
% Beta

Bc_ArcA_n = ones(1,size(Bc_ArcA,2));
Bc_ArcF_n = ones(1,size(Bc_ArcF,2));
Br_ArcA_n = ones(1,size(Br_ArcA,2));
Br_ArcF_n = ones(1,size(Br_ArcF,2));

Bc_HalA_n = ones(1,size(Bc_HalA,2))*2;
Bc_HalF_n = ones(1,size(Bc_HalF,2))*2;
Br_HalA_n = ones(1,size(Br_HalA,2))*2;
Br_HalF_n = ones(1,size(Br_HalF,2))*2;

Bc_ConA_n = ones(1,size(Bc_ConA,2))*3;
Bc_ConF_n = ones(1,size(Bc_ConF,2))*3;
Br_ConA_n = ones(1,size(Br_ConA,2))*3;
Br_ConF_n = ones(1,size(Br_ConF,2))*3;


Big2 = [Bc_ArcA_n Bc_ArcF_n Br_ArcA_n Br_ArcF_n Bc_HalA_n Bc_HalF_n Br_HalA_n Br_HalF_n ...
    Bc_ConA_n Bc_ConF_n Br_ConA_n Br_ConF_n];
%%
% Gamma

Gc_ArcA_n = ones(1,size(Gc_ArcA,2));
Gc_ArcF_n = ones(1,size(Gc_ArcF,2));
Gr_ArcA_n = ones(1,size(Gr_ArcA,2));
Gr_ArcF_n = ones(1,size(Gr_ArcF,2));

Gc_HalA_n = ones(1,size(Gc_HalA,2))*2;
Gc_HalF_n = ones(1,size(Gc_HalF,2))*2;
Gr_HalA_n = ones(1,size(Gr_HalA,2))*2;
Gr_HalF_n = ones(1,size(Gr_HalF,2))*2;

Gc_ConA_n = ones(1,size(Gc_ConA,2))*3;
Gc_ConF_n = ones(1,size(Gc_ConF,2))*3;
Gr_ConA_n = ones(1,size(Gr_ConA,2))*3;
Gr_ConF_n = ones(1,size(Gr_ConF,2))*3;


Big2 = [Gc_ArcA_n Gc_ArcF_n Gr_ArcA_n Gr_ArcF_n Gc_HalA_n Gc_HalF_n Gr_HalA_n Gr_HalF_n ...
    Gc_ConA_n Gc_ConF_n Gr_ConA_n Gr_ConF_n];
%%
% GammaX

GXc_ArcA_n = ones(1,size(GXc_ArcA,2));
GXc_ArcF_n = ones(1,size(GXc_ArcF,2));
GXr_ArcA_n = ones(1,size(GXr_ArcA,2));
GXr_ArcF_n = ones(1,size(GXr_ArcF,2));

GXc_HalA_n = ones(1,size(GXc_HalA,2))*2;
GXc_HalF_n = ones(1,size(GXc_HalF,2))*2;
GXr_HalA_n = ones(1,size(GXr_HalA,2))*2;
GXr_HalF_n = ones(1,size(GXr_HalF,2))*2;

GXc_ConA_n = ones(1,size(GXc_ConA,2))*3;
GXc_ConF_n = ones(1,size(GXc_ConF,2))*3;
GXr_ConA_n = ones(1,size(GXr_ConA,2))*3;
GXr_ConF_n = ones(1,size(GXr_ConF,2))*3;


Big2 = [GXc_ArcA_n GXc_ArcF_n GXr_ArcA_n GXr_ArcF_n GXc_HalA_n GXc_HalF_n GXr_HalA_n GXr_HalF_n ...
    GXc_ConA_n GXc_ConF_n GXr_ConA_n GXr_ConF_n];
%%
% GammaXX

GXXc_ArcA_n = ones(1,size(GXXc_ArcA,2));
GXXc_ArcF_n = ones(1,size(GXXc_ArcF,2));
GXXr_ArcA_n = ones(1,size(GXXr_ArcA,2));
GXXr_ArcF_n = ones(1,size(GXXr_ArcF,2));

GXXc_HalA_n = ones(1,size(GXXc_HalA,2))*2;
GXXc_HalF_n = ones(1,size(GXXc_HalF,2))*2;
GXXr_HalA_n = ones(1,size(GXXr_HalA,2))*2;
GXXr_HalF_n = ones(1,size(GXXr_HalF,2))*2;

GXXc_ConA_n = ones(1,size(GXXc_ConA,2))*3;
GXXc_ConF_n = ones(1,size(GXXc_ConF,2))*3;
GXXr_ConA_n = ones(1,size(GXXr_ConA,2))*3;
GXXr_ConF_n = ones(1,size(GXXr_ConF,2))*3;


Big2 = [GXXc_ArcA_n GXXc_ArcF_n GXXr_ArcA_n GXXr_ArcF_n GXXc_HalA_n GXXc_HalF_n GXXr_HalA_n GXXr_HalF_n ...
    GXXc_ConA_n GXXc_ConF_n GXXr_ConA_n GXXr_ConF_n];

%% Caudal/Rostral

% % Delta

Dc_ArcA_cr = ones(1,size(Dc_ArcA,2));
Dc_ArcF_cr = ones(1,size(Dc_ArcF,2));
Dr_ArcA_cr = ones(1,size(Dr_ArcA,2))*4;
Dr_ArcF_cr = ones(1,size(Dr_ArcF,2))*4;

Dc_HalA_cr = ones(1,size(Dc_HalA,2));
Dc_HalF_cr = ones(1,size(Dc_HalF,2));
Dr_HalA_cr = ones(1,size(Dr_HalA,2))*4;
Dr_HalF_cr = ones(1,size(Dr_HalF,2))*4;

Dc_ConA_cr = ones(1,size(Dc_ConA,2));
Dc_ConF_cr = ones(1,size(Dc_ConF,2));
Dr_ConA_cr = ones(1,size(Dr_ConA,2))*4;
Dr_ConF_cr = ones(1,size(Dr_ConF,2))*4;


Big3 = [Dc_ArcA_cr Dc_ArcF_cr Dr_ArcA_cr Dr_ArcF_cr Dc_HalA_cr Dc_HalF_cr Dr_HalA_cr Dr_HalF_cr ...
    Dc_ConA_cr Dc_ConF_cr Dr_ConA_cr Dr_ConF_cr];

%%
% % Theta

Tc_ArcA_cr = ones(1,size(Tc_ArcA,2));
Tc_ArcF_cr = ones(1,size(Tc_ArcF,2));
Tr_ArcA_cr = ones(1,size(Tr_ArcA,2))*4;
Tr_ArcF_cr = ones(1,size(Tr_ArcF,2))*4;

Tc_HalA_cr = ones(1,size(Tc_HalA,2));
Tc_HalF_cr = ones(1,size(Tc_HalF,2));
Tr_HalA_cr = ones(1,size(Tr_HalA,2))*4;
Tr_HalF_cr = ones(1,size(Tr_HalF,2))*4;

Tc_ConA_cr = ones(1,size(Tc_ConA,2));
Tc_ConF_cr = ones(1,size(Tc_ConF,2));
Tr_ConA_cr = ones(1,size(Tr_ConA,2))*4;
Tr_ConF_cr = ones(1,size(Tr_ConF,2))*4;


Big3 = [Tc_ArcA_cr Tc_ArcF_cr Tr_ArcA_cr Tr_ArcF_cr Tc_HalA_cr Tc_HalF_cr Tr_HalA_cr Tr_HalF_cr ...
    Tc_ConA_cr Tc_ConF_cr Tr_ConA_cr Tr_ConF_cr];

%%
% % Alpha

Ac_ArcA_cr = ones(1,size(Ac_ArcA,2));
Ac_ArcF_cr = ones(1,size(Ac_ArcF,2));
Ar_ArcA_cr = ones(1,size(Ar_ArcA,2))*4;
Ar_ArcF_cr = ones(1,size(Ar_ArcF,2))*4;

Ac_HalA_cr = ones(1,size(Ac_HalA,2));
Ac_HalF_cr = ones(1,size(Ac_HalF,2));
Ar_HalA_cr = ones(1,size(Ar_HalA,2))*4;
Ar_HalF_cr = ones(1,size(Ar_HalF,2))*4;

Ac_ConA_cr = ones(1,size(Ac_ConA,2));
Ac_ConF_cr = ones(1,size(Ac_ConF,2));
Ar_ConA_cr = ones(1,size(Ar_ConA,2))*4;
Ar_ConF_cr = ones(1,size(Ar_ConF,2))*4;


Big3 = [Ac_ArcA_cr Ac_ArcF_cr Ar_ArcA_cr Ar_ArcF_cr Ac_HalA_cr Ac_HalF_cr Ar_HalA_cr Ar_HalF_cr ...
    Ac_ConA_cr Ac_ConF_cr Ar_ConA_cr Ar_ConF_cr];
%%
% % Beta

Bc_ArcA_cr = ones(1,size(Bc_ArcA,2));
Bc_ArcF_cr = ones(1,size(Bc_ArcF,2));
Br_ArcA_cr = ones(1,size(Br_ArcA,2))*4;
Br_ArcF_cr = ones(1,size(Br_ArcF,2))*4;

Bc_HalA_cr = ones(1,size(Bc_HalA,2));
Bc_HalF_cr = ones(1,size(Bc_HalF,2));
Br_HalA_cr = ones(1,size(Br_HalA,2))*4;
Br_HalF_cr = ones(1,size(Br_HalF,2))*4;

Bc_ConA_cr = ones(1,size(Bc_ConA,2));
Bc_ConF_cr = ones(1,size(Bc_ConF,2));
Br_ConA_cr = ones(1,size(Br_ConA,2))*4;
Br_ConF_cr = ones(1,size(Br_ConF,2))*4;


Big3 = [Bc_ArcA_cr Bc_ArcF_cr Br_ArcA_cr Br_ArcF_cr Bc_HalA_cr Bc_HalF_cr Br_HalA_cr Br_HalF_cr ...
    Bc_ConA_cr Bc_ConF_cr Br_ConA_cr Br_ConF_cr];

%%
% % Gamma

Gc_ArcA_cr = ones(1,size(Gc_ArcA,2));
Gc_ArcF_cr = ones(1,size(Gc_ArcF,2));
Gr_ArcA_cr = ones(1,size(Gr_ArcA,2))*4;
Gr_ArcF_cr = ones(1,size(Gr_ArcF,2))*4;

Gc_HalA_cr = ones(1,size(Gc_HalA,2));
Gc_HalF_cr = ones(1,size(Gc_HalF,2));
Gr_HalA_cr = ones(1,size(Gr_HalA,2))*4;
Gr_HalF_cr = ones(1,size(Gr_HalF,2))*4;

Gc_ConA_cr = ones(1,size(Gc_ConA,2));
Gc_ConF_cr = ones(1,size(Gc_ConF,2));
Gr_ConA_cr = ones(1,size(Gr_ConA,2))*4;
Gr_ConF_cr = ones(1,size(Gr_ConF,2))*4;


Big3 = [Gc_ArcA_cr Gc_ArcF_cr Gr_ArcA_cr Gr_ArcF_cr Gc_HalA_cr Gc_HalF_cr Gr_HalA_cr Gr_HalF_cr ...
    Gc_ConA_cr Gc_ConF_cr Gr_ConA_cr Gr_ConF_cr];

%%
% % GammaX

GXc_ArcA_cr = ones(1,size(GXc_ArcA,2));
GXc_ArcF_cr = ones(1,size(GXc_ArcF,2));
GXr_ArcA_cr = ones(1,size(GXr_ArcA,2))*4;
GXr_ArcF_cr = ones(1,size(GXr_ArcF,2))*4;

GXc_HalA_cr = ones(1,size(GXc_HalA,2));
GXc_HalF_cr = ones(1,size(GXc_HalF,2));
GXr_HalA_cr = ones(1,size(GXr_HalA,2))*4;
GXr_HalF_cr = ones(1,size(GXr_HalF,2))*4;

GXc_ConA_cr = ones(1,size(GXc_ConA,2));
GXc_ConF_cr = ones(1,size(GXc_ConF,2));
GXr_ConA_cr = ones(1,size(GXr_ConA,2))*4;
GXr_ConF_cr = ones(1,size(GXr_ConF,2))*4;


Big3 = [GXc_ArcA_cr GXc_ArcF_cr GXr_ArcA_cr GXr_ArcF_cr GXc_HalA_cr GXc_HalF_cr GXr_HalA_cr GXr_HalF_cr ...
    GXc_ConA_cr GXc_ConF_cr GXr_ConA_cr GXr_ConF_cr];
%%

% % GammaXX

GXXc_ArcA_cr = ones(1,size(GXXc_ArcA,2));
GXXc_ArcF_cr = ones(1,size(GXXc_ArcF,2));
GXXr_ArcA_cr = ones(1,size(GXXr_ArcA,2))*4;
GXXr_ArcF_cr = ones(1,size(GXXr_ArcF,2))*4;

GXXc_HalA_cr = ones(1,size(GXXc_HalA,2));
GXXc_HalF_cr = ones(1,size(GXXc_HalF,2));
GXXr_HalA_cr = ones(1,size(GXXr_HalA,2))*4;
GXXr_HalF_cr = ones(1,size(GXXr_HalF,2))*4;

GXXc_ConA_cr = ones(1,size(GXXc_ConA,2));
GXXc_ConF_cr = ones(1,size(GXXc_ConF,2));
GXXr_ConA_cr = ones(1,size(GXXr_ConA,2))*4;
GXXr_ConF_cr = ones(1,size(GXXr_ConF,2))*4;


Big3 = [GXXc_ArcA_cr GXXc_ArcF_cr GXXr_ArcA_cr GXXr_ArcF_cr GXXc_HalA_cr GXXc_HalF_cr GXXr_HalA_cr GXXr_HalF_cr ...
    GXXc_ConA_cr GXXc_ConF_cr GXXr_ConA_cr GXXr_ConF_cr];

%% Auditory/Frontal

% % Delta

Dc_ArcA_af = ones(1,size(Dc_ArcA,2));
Dc_ArcF_af = ones(1,size(Dc_ArcF,2))*5;
Dr_ArcA_af = ones(1,size(Dr_ArcA,2));
Dr_ArcF_af = ones(1,size(Dr_ArcF,2))*5;

Dc_HalA_af = ones(1,size(Dc_HalA,2));
Dc_HalF_af = ones(1,size(Dc_HalF,2))*5;
Dr_HalA_af = ones(1,size(Dr_HalA,2));
Dr_HalF_af = ones(1,size(Dr_HalF,2))*5;

Dc_ConA_af = ones(1,size(Dc_ConA,2));
Dc_ConF_af = ones(1,size(Dc_ConF,2))*5;
Dr_ConA_af = ones(1,size(Dr_ConA,2));
Dr_ConF_af = ones(1,size(Dr_ConF,2))*5;


Big4 = [Dc_ArcA_af Dc_ArcF_af Dr_ArcA_af Dr_ArcF_af Dc_HalA_af Dc_HalF_af Dr_HalA_af Dr_HalF_af ...
    Dc_ConA_af Dc_ConF_af Dr_ConA_af Dr_ConF_af];
%%
% % Theta

Tc_ArcA_af = ones(1,size(Tc_ArcA,2));
Tc_ArcF_af = ones(1,size(Tc_ArcF,2))*5;
Tr_ArcA_af = ones(1,size(Tr_ArcA,2));
Tr_ArcF_af = ones(1,size(Tr_ArcF,2))*5;

Tc_HalA_af = ones(1,size(Tc_HalA,2));
Tc_HalF_af = ones(1,size(Tc_HalF,2))*5;
Tr_HalA_af = ones(1,size(Tr_HalA,2));
Tr_HalF_af = ones(1,size(Tr_HalF,2))*5;

Tc_ConA_af = ones(1,size(Tc_ConA,2));
Tc_ConF_af = ones(1,size(Tc_ConF,2))*5;
Tr_ConA_af = ones(1,size(Tr_ConA,2));
Tr_ConF_af = ones(1,size(Tr_ConF,2))*5;


Big4 = [Tc_ArcA_af Tc_ArcF_af Tr_ArcA_af Tr_ArcF_af Tc_HalA_af Tc_HalF_af Tr_HalA_af Tr_HalF_af ...
    Tc_ConA_af Tc_ConF_af Tr_ConA_af Tr_ConF_af];

%%
% % Alpha

Ac_ArcA_af = ones(1,size(Ac_ArcA,2));
Ac_ArcF_af = ones(1,size(Ac_ArcF,2))*5;
Ar_ArcA_af = ones(1,size(Ar_ArcA,2));
Ar_ArcF_af = ones(1,size(Ar_ArcF,2))*5;

Ac_HalA_af = ones(1,size(Ac_HalA,2));
Ac_HalF_af = ones(1,size(Ac_HalF,2))*5;
Ar_HalA_af = ones(1,size(Ar_HalA,2));
Ar_HalF_af = ones(1,size(Ar_HalF,2))*5;

Ac_ConA_af = ones(1,size(Ac_ConA,2));
Ac_ConF_af = ones(1,size(Ac_ConF,2))*5;
Ar_ConA_af = ones(1,size(Ar_ConA,2));
Ar_ConF_af = ones(1,size(Ar_ConF,2))*5;


Big4 = [Ac_ArcA_af Ac_ArcF_af Ar_ArcA_af Ar_ArcF_af Ac_HalA_af Ac_HalF_af Ar_HalA_af Ar_HalF_af ...
    Ac_ConA_af Ac_ConF_af Ar_ConA_af Ar_ConF_af];
%%
% % Beta

Bc_ArcA_af = ones(1,size(Bc_ArcA,2));
Bc_ArcF_af = ones(1,size(Bc_ArcF,2))*5;
Br_ArcA_af = ones(1,size(Br_ArcA,2));
Br_ArcF_af = ones(1,size(Br_ArcF,2))*5;

Bc_HalA_af = ones(1,size(Bc_HalA,2));
Bc_HalF_af = ones(1,size(Bc_HalF,2))*5;
Br_HalA_af = ones(1,size(Br_HalA,2));
Br_HalF_af = ones(1,size(Br_HalF,2))*5;

Bc_ConA_af = ones(1,size(Bc_ConA,2));
Bc_ConF_af = ones(1,size(Bc_ConF,2))*5;
Br_ConA_af = ones(1,size(Br_ConA,2));
Br_ConF_af = ones(1,size(Br_ConF,2))*5;


Big4 = [Bc_ArcA_af Bc_ArcF_af Br_ArcA_af Br_ArcF_af Bc_HalA_af Bc_HalF_af Br_HalA_af Br_HalF_af ...
    Bc_ConA_af Bc_ConF_af Br_ConA_af Br_ConF_af];
%%
% % Gamma

Gc_ArcA_af = ones(1,size(Gc_ArcA,2));
Gc_ArcF_af = ones(1,size(Gc_ArcF,2))*5;
Gr_ArcA_af = ones(1,size(Gr_ArcA,2));
Gr_ArcF_af = ones(1,size(Gr_ArcF,2))*5;

Gc_HalA_af = ones(1,size(Gc_HalA,2));
Gc_HalF_af = ones(1,size(Gc_HalF,2))*5;
Gr_HalA_af = ones(1,size(Gr_HalA,2));
Gr_HalF_af = ones(1,size(Gr_HalF,2))*5;

Gc_ConA_af = ones(1,size(Gc_ConA,2));
Gc_ConF_af = ones(1,size(Gc_ConF,2))*5;
Gr_ConA_af = ones(1,size(Gr_ConA,2));
Gr_ConF_af = ones(1,size(Gr_ConF,2))*5;


Big4 = [Gc_ArcA_af Gc_ArcF_af Gr_ArcA_af Gr_ArcF_af Gc_HalA_af Gc_HalF_af Gr_HalA_af Gr_HalF_af ...
    Gc_ConA_af Gc_ConF_af Gr_ConA_af Gr_ConF_af];
%%
% % GammaX

GXc_ArcA_af = ones(1,size(GXc_ArcA,2));
GXc_ArcF_af = ones(1,size(GXc_ArcF,2))*5;
GXr_ArcA_af = ones(1,size(GXr_ArcA,2));
GXr_ArcF_af = ones(1,size(GXr_ArcF,2))*5;

GXc_HalA_af = ones(1,size(GXc_HalA,2));
GXc_HalF_af = ones(1,size(GXc_HalF,2))*5;
GXr_HalA_af = ones(1,size(GXr_HalA,2));
GXr_HalF_af = ones(1,size(GXr_HalF,2))*5;

GXc_ConA_af = ones(1,size(GXc_ConA,2));
GXc_ConF_af = ones(1,size(GXc_ConF,2))*5;
GXr_ConA_af = ones(1,size(GXr_ConA,2));
GXr_ConF_af = ones(1,size(GXr_ConF,2))*5;


Big4 = [GXc_ArcA_af GXc_ArcF_af GXr_ArcA_af GXr_ArcF_af GXc_HalA_af GXc_HalF_af GXr_HalA_af GXr_HalF_af ...
    GXc_ConA_af GXc_ConF_af GXr_ConA_af GXr_ConF_af];
%%

% % GammaXX

GXXc_ArcA_af = ones(1,size(GXXc_ArcA,2));
GXXc_ArcF_af = ones(1,size(GXXc_ArcF,2))*5;
GXXr_ArcA_af = ones(1,size(GXXr_ArcA,2));
GXXr_ArcF_af = ones(1,size(GXXr_ArcF,2))*5;

GXXc_HalA_af = ones(1,size(GXXc_HalA,2));
GXXc_HalF_af = ones(1,size(GXXc_HalF,2))*5;
GXXr_HalA_af = ones(1,size(GXXr_HalA,2));
GXXr_HalF_af = ones(1,size(GXXr_HalF,2))*5;

GXXc_ConA_af = ones(1,size(GXXc_ConA,2));
GXXc_ConF_af = ones(1,size(GXXc_ConF,2))*5;
GXXr_ConA_af = ones(1,size(GXXr_ConA,2));
GXXr_ConF_af = ones(1,size(GXXr_ConF,2))*5;


Big4 = [GXXc_ArcA_af GXXc_ArcF_af GXXr_ArcA_af GXXr_ArcF_af GXXc_HalA_af GXXc_HalF_af GXXr_HalA_af GXXr_HalF_af ...
    GXXc_ConA_af GXXc_ConF_af GXXr_ConA_af GXXr_ConF_af];


%%

Data = Big1';
Grp1 = Big2';
Grp2 = Big3';
Grp3 = Big4';

%%

%% ANOVA
[~, ANOVAtbl, stats] = anovan(Data, {Grp1, Grp2}, 'model','interaction','varnames',{'Opsin','Time'});


%%


c = multcompare(stats, 'Ctype', 'hsd', 'Dimension', [1 2], 'Display', 'on');
c


%%
A = table2array(A14)
Data = A(:,1);
Grp1 = A(:,2);
Grp2 = A(:,3);

