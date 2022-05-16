

%Vecto_FR3 = Vecto_FR3f; 
%VectoS3 = VectoS3f;
%% Plotting Spectrogram for each 
Fg = 90;


params.Fs = 1000;
params.fpass = [0 50]; 
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

% 

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
    %V_ALb(kl,:) = 10*log10(ALb);
    %V_FLb(kl,:) = 10*log10(FLb);
    V_FRb(kl,:) = 10*log10(FRb);
    %V_EMGb(kl,:) = 10*log10(EMGb);

    %V_ARa(kl,:) = 10*log10(ARa);
    %V_ALa(kl,:) = 10*log10(ALa);
    %V_FLa(kl,:) = 10*log10(FLa);
    V_FRa(kl,:) = 10*log10(FRa);
    %V_EMGa(kl,:) = 10*log10(EMGa);
    
end

% get mean 

%Vec_ARb = mean(V_ARb);
%Vec_ALb = mean(V_ALb);
%Vec_FLb = mean(V_FLb);
Vec_FRb = mean(V_FRb);
%Vec_EMGb = mean(V_EMGb);

%Vec_ARa = mean(V_ARa);
%Vec_ALa = mean(V_ALa);
%Vec_FLa = mean(V_FLa);
Vec_FRa = mean(V_FRa);
%Vec_EMGa = mean(V_EMGa);

% Plotting figures

figure
plot(f,Vec_FRb,'r','LineWidth',2)
hold on
plot(f,Vec_FRa,'g','LineWidth',2)
ylim([5 15])
xlim([16 30])
set(gca,'fontsize',14)
% subplot(2,2,4)
% plot(f,Vec_FLb,'r','LineWidth',2)
% hold on
% plot(f,Vec_FLa,'g','LineWidth',2)
% FLstr = sprintf('Frontal Left, N = %d',trl);
% title(FLstr)
% axis xy; xlabel('Hz'); ylabel('Power')
% legend('Stim Off','Stim On')


%%
% figure
% plot(f,Vec_FLb,'r','LineWidth',2)
% hold on
% plot(f,Vec_FLa,'g','LineWidth',2)
% FLstr = sprintf('Frontal Left, N = %d',trl);
% title(FLstr)
% axis xy; xlabel('Hz'); ylabel('Power')
% legend('Stim Off','Stim On')
