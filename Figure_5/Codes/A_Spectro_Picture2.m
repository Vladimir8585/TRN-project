%% Spectrogram picture

movingwin = [1 0.5];
params.Fs = 1000;
params.fpass = [16 30];
params.tapers = [3 5];
params.trialave = 1;
params.err = 1;
maxDbG = 5; 


%% Right vector length
rt = size(Vecto_FR3,2);
[ARV,Gt,Gf] = mtspecgramc(Vecto_FR3(1,:),movingwin,params);

%% for loop  construction and summation
ktk = size(Vecto_FR3,1)
rg = size(ARV,1);
gr = size(ARV,2);
%ALvecto = zeros(ktk,rg,gr);
%ARvecto = zeros(ktk,rg,gr);
%FLvecto = zeros(ktk,rg,gr);
FRvecto = zeros(ktk,rg,gr);

for tj = 1:ktk

%[ALvecto(tj,:,:),Gt,Gf] = mtspecgramc(Vecto_AL3(tj,:),movingwin,params);
%[ARvecto(tj,:,:),Gt,Gf] = mtspecgramc(Vecto_AR3(tj,:),movingwin,params);
%[FLvecto(tj,:,:),Gt,Gf] = mtspecgramc(Vecto_FL3(tj,:),movingwin,params);
[FRvecto(tj,:,:),Gt,Gf] = mtspecgramc(Vecto_FR3(tj,:),movingwin,params);

end



%% try 2 with loop


%sumAL = squeeze(mean(ALvecto,1));
%sumAR = squeeze(mean(ARvecto,1));
%sumFL = squeeze(mean(FLvecto,1));
sumFR = squeeze(mean(FRvecto,1));
%sumEMG = mean(Vecto_EMG3,1);


Tvec = 1:length(VectoS3);
trl = size(VectoS3,1);


%%

% subplot(6,1,1)
% for ls = 1:ktk
%     plot(Tvec,VectoS3(ls,:));
%     hold on
% end
% Sstr = sprintf('Averaged spectrogram, N = %d',trl);
% title(Sstr)
% axis([0000, rt, 0, 3.5])
% subplot(6,1,2)
%     plot(Tvec,sumEMG);
% axis([0, rt, -200, 200])

% subplot(6,1,3);
% %imagesc(Gt,Gf,10*log10(sumAL'),[0 maxDbG])
% %axis xy; xlabel('Time(s)'); ylabel('Freq(Hz)')
% %title('Auditory left')
% subplot(6,1,4);
% imagesc(Gt,Gf,10*log10(sumAR'),[22 maxDbG])
% axis xy; xlabel('Time(s)'); ylabel('Freq(Hz)')
% title('Auditory right')
% subplot(6,1,5);
% imagesc(Gt,Gf,10*log10(sumFL'),[22 maxDbG])
% axis xy; xlabel('Time(s)'); ylabel('Freq(Hz)')
% title('Frontal left')
% subplot(6,1,6);
% imagesc(Gt,Gf,10*log10(sumFR'),[22 maxDbG])
% axis xy; xlabel('Time(s)'); ylabel('Freq(Hz)')
% title('Frontal right')

figure
subplot(3,1,1)
imagesc(Gt,Gf,10*log10(sumFR'),[8 17])
axis xy; xlabel('Time(s)'); ylabel('Freq(Hz)')
%title('Frontal right')
%colorbar
