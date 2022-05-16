%% Plotting Spectrogram for each 
%Fg = 4;


params.Fs = 1000;
params.fpass = [16 30]; 
params.tapers = [3 5];
params.trialave = 0;
params.pad=1;
params.err = 1;

trl = size(Vecto_FR3,1);
tsec = size(Vecto_FR3,2);
tbef = 1:tsec/2;
taft = tsec/2+1:tsec;

[S,f] = mtspectrumc( Vecto_FR3(1,tbef), params );
Pw = 10*log10(S);
%Pw_sm = smooth(Pw,0.02,'lowess');
m_cel = size(Pw,1);

%% 

%V_ARb = zeros(trl,m_cel);
%V_ALb = zeros(trl,m_cel);
%V_FLb = zeros(trl,m_cel);
V_FRb = zeros(trl,m_cel);
%V_EMGb = zeros(trl,m_cel);

%V_ARa = zeros(trl,m_cel);
%V_ALa = zeros(trl,m_cel);
%V_FLa = zeros(trl,m_cel);
V_FRa = zeros(trl,m_cel);
%V_EMGa = zeros(trl,m_cel);

for kl = 1:trl
    % detecting spectrogram
    %[ARb,f] = mtspectrumc( Vecto_AR3(kl,tbef), params );
    %[ALb,f] = mtspectrumc( Vecto_AL3(kl,tbef), params );
    [FRb,f] = mtspectrumc( Vecto_FR3(kl,tbef), params );
    %[FLb,f] = mtspectrumc( Vecto_FL3(kl,tbef), params );
    %[EMGb,f] = mtspectrumc( Vecto_EMG3(kl,tbef), params );
    
    %[ARa,f] = mtspectrumc( Vecto_AR3(kl,taft), params );
    %[ALa,f] = mtspectrumc( Vecto_AL3(kl,taft), params );
    [FRa,f] = mtspectrumc( Vecto_FR3(kl,taft), params );
    %[FLa,f] = mtspectrumc( Vecto_FL3(kl,taft), params );
    %[EMGa,f] = mtspectrumc( Vecto_EMG3(kl,taft), params );
    
    %log and  non smooth
    %placing into matrixes
    
   %V_ARb(kl,:) = 10*log10(ARb);
   % V_ALb(kl,:) = 10*log10(ALb);
    %V_FLb(kl,:) = 10*log10(FLb);
    V_FRb(kl,:) = 10*log10(FRb);
    %V_EMGb(kl,:) = 10*log10(EMGb);

    %V_ARa(kl,:) = 10*log10(ARa);
    %V_ALa(kl,:) = 10*log10(ALa);
    %V_FLa(kl,:) = 10*log10(FLa);
    V_FRa(kl,:) = 10*log10(FRa);
    %V_EMGa(kl,:) = 10*log10(EMGa);
    
end

%% get mean 

%Vec_ARb = mean(V_ARb);
%Vec_ALb = mean(V_ALb);
%Vec_FLb = mean(V_FLb,2);
Vec_FRb = mean(V_FRb,2);
%Vec_EMGb = mean(V_EMGb,2);

%Vec_ARa = mean(V_ARa);
%Vec_ALa = mean(V_ALa);
%Vec_FLa = mean(V_FLa,2);
Vec_FRa = mean(V_FRa,2);
%Vec_EMGa = mean(V_EMGa,2);



%%


%x = Vec_ARb;
%y = Vec_ARa;
 
 
%[h,p] = ttest(x,y,'Alpha',0.01)


x = Vec_FRb;
y = Vec_FRa;

 
[h,p] = ttest(x,y,'Alpha',0.05)

%x = Vec_FLb;
%y = Vec_FLa;


%[h,p] = ttest(x,y,'Alpha',0.01)


%%


%save('Vec_FRb_142st__S.mat','Vec_FRb')
%save('Vec_FRa_142st__AL.mat','Vec_FRa')
%save('Vec_FLb_312st__AR.mat','Vec_FLb')
%save('Vec_FLa_312st__FL.mat','Vec_FLa')
%%



%Arch_AR = [-0.33	-0.35 0.39	-0.64	-1.49	0.6	-1.65	-0.29];
%Arch_AL = [-0.16	0.4	-0.4	0.15	-0.69	-1.28	0.96	0.13	-2.24	-0.32];
%Arch_FL = [-0.29	-0.25	-0.12	0.94	-0.83	0.13	-1.24	-1.04	-1.65	-0.29];
%Arch_FR = [-0.3	-0.04	-0.38	0.72	-0.62	0.21	-1.42	-1.49	0.6	-1.65	-0.29];



%Halo_AR = [-1.03	-0.67	-0.34	-0.14];
%Halo_AL = [-0.07	-1.18	-0.61	-0.35];
%Halo_FL = [-1.21	0.03	-1	-1.1	-0.56];
%Halo_FR = [-1.03	0.09	-1.07	-1.06	-0.64	-0.31];

%Con_AR = [0.47	0.35	0.58	-0.81];
%Con_FL = [0.43	0.63	0.07	0.39];
% Con_FR = [0.34	0.64	0.07	0.62];









