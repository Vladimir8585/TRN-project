VectoS3 = VectoS3f;
Vecto_FR3 = Vecto_FR3f;

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
con = 10; % period of seconds to find
    
    
    % cutting required parts of vector
    
    
    for s = 1:pttFR
                Vecto_FR3t(s,:) = Vecto_FR3w(s,(frk-con*Wind):(frk+con*Wind)-1);
    end

                
         % Collecting data for power matrix       
                
               tsec = size(Vecto_FR3t,2);
                tbef = 1:(tsec/2) ;
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
                Beta_p = (16*points_f:30*points_f);
                Gamma_p = (31*points_f:50*points_f);
                GammaX_p = (50*points_f:80*points_f);
                Average_p = (100*points_f:110*points_f);
                
                % Calculate Ch
                
              
              for fr = 1:pttFR
                Delt_FRb(fr,:) = mean(V_FRb(fr,Delta_p));
                Delt_FRa(fr,:) = mean(V_FRa(fr,Delta_p));
                Thet_FRb(fr) = mean(V_FRb(fr,Theta_p));
                Thet_FRa(fr) = mean(V_FRa(fr,Theta_p));
                Alph_FRb(fr) = mean(V_FRb(fr,Spindle_p));
                Alph_FRa(fr) = mean(V_FRa(fr,Spindle_p));
                Beta_FRb(fr) = mean(V_FRb(fr,Beta_p));
                Beta_FRa(fr) = mean(V_FRa(fr,Beta_p));
                Gamma_FRb(fr) = mean(V_FRb(fr,Gamma_p));
                Gamma_FRa(fr) = mean(V_FRa(fr,Gamma_p));
                GammaX_FRb(fr) = mean(V_FRb(fr,GammaX_p));
                GammaX_FRa(fr) = mean(V_FRa(fr,GammaX_p));
                Average_FRb(fr) = sum(V_FRb(fr,Average_p));
                Average_FRa(fr) = sum(V_FRa(fr,Average_p));
              end

              
              
               % Calcuate diffeference

               for cc = 1:pttFR
                
                  D_Delta_FR(cc) = Delt_FRa(cc)-Delt_FRb(cc);
                  D_Theta_FR(cc) = Thet_FRa(cc)-Thet_FRb(cc);
                  D_Alph_FR(cc) = Alph_FRa(cc)-Alph_FRb(cc);
                  D_Beta_FR(cc) = Beta_FRa(cc)-Beta_FRb(cc);
                  D_Gamma_FR(cc) = Gamma_FRa(cc)-Gamma_FRb(cc);
                  D_GammaX_FR(cc) = GammaX_FRa(cc)-GammaX_FRb(cc);
                  D_Average_FR(cc) = Average_FRa(cc)-Average_FRb(cc);

               end

                
 %Calculate mean difference

d_FR = mean(D_Delta_FR);
t_FR = mean(D_Theta_FR);
a_FR = mean(D_Alph_FR);
b_FR = mean(D_Beta_FR);
g_FR = mean(D_Gamma_FR);
gx_FR = mean(D_GammaX_FR);
Aver_FR = mean(D_Average_FR);

% Calculate % normalized difference nomalized by baseline (previous 10 seconds - con = 10)

                 for kk = 1:pttFR
                
                  D_Delta_FR(kk) = (Delt_FRa(kk)-Delt_FRb(kk))/Delt_FRb(kk)*100;
                  D_Theta_FR(kk) = (Thet_FRa(kk)-Thet_FRb(kk))/Thet_FRb(kk)*100;
                  D_Alph_FR(kk) = (Alph_FRa(kk)-Alph_FRb(kk))/Alph_FRb(kk)*100;
                  D_Beta_FR(kk) = (Beta_FRa(kk)-Beta_FRb(kk))/Beta_FRb(kk)*100;
                  D_Gamma_FR(kk) = (Gamma_FRa(kk)-Gamma_FRb(kk))/Gamma_FRb(kk)*100;
                  D_GammaX_FR(kk) = (GammaX_FRa(kk)-GammaX_FRb(kk))/GammaX_FRb(kk)*100;
                  D_Average_FR(kk) = Average_FRa(kk)-Average_FRb(kk);

                 end


% % Calculate raw data before and after stimulation 

                D_b = mean(Delt_FRb);
                D_a = mean(Delt_FRa);
                T_b = mean(Thet_FRb);
                T_a = mean(Thet_FRa);
                A_b = mean(Alph_FRb);
                A_a = mean(Alph_FRa);
                B_b = mean(Beta_FRb);
                B_a = mean(Beta_FRa);
                G_b = mean(Gamma_FRb);
                G_a = mean(Gamma_FRa);
                G_xb = mean(GammaX_FRb);
                G_xa = mean(GammaX_FRa);
                Av_b = mean(Average_FRb);
                AV_a = mean(Average_FRa);



% Deliver mean difference

AC_BFR_Pfx = [];
AC_BFR_Pfx(1) = d_FR;
%AC_BFR_Pfx(2) = t_FR;
AC_BFR_Pfx(2) = a_FR;
AC_BFR_Pfx(3) = b_FR;
%AC_BFR_Pfx(5) = g_FR;
%AC_BFR_Pfx(6) = gx_FR;
%AC_BFR_Pfx(7) = Aver_FR;


printmat(AC_BFR_Pfx , 'Delta Alpha Beta')
%title('Difference before and after stimulation')

% Deliver % normalized difference

Pfx = [];
Pfx(1) = mean(D_Delta_FR);
Pfx(2) = mean(D_Alph_FR);
Pfx(3) = mean(D_Beta_FR);

printmat(Pfx , 'Delta% Alpha% Beta%')
%title('Normalized difference before and after stimulation')


% Deliver real difference

BFR = [];
                BFR(1,1) = D_b;
                BFR(1,2) = D_a;
                %BFR(2,1) = T_b;
                %BFR(2,2) = T_a;
                BFR(2,1) = A_b;
                BFR(2,2) = A_a;
                BFR(3,1) = B_b;
                BFR(3,2) = B_a;
                %BFR(5,1) = G_b;
                %BFR(5,2) = G_a;
                %BFR(6,1) = G_xb;
                %BFR(6,2) = G_xa;
                %BFR(7,1) = Av_b;
                %BFR(7,2) = AV_a;


printmat(BFR, 'BandsFR', 'Delta Alpha Beta ', 'Before After')
%title('Real data  before and after stimulation')

%%
