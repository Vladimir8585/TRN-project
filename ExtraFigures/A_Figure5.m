%VectoS3 = VectoS3f;
%Vecto_FR3 = Vecto_FR3f;
%Vecto_EMG3 = Vecto_FR3f;

%% Trimming code

ptt = size(Vecto_EMG3,1)
pt = size(Vecto_EMG3,2)

for s = 1:ptt
 for p = 2: pt
    if (VectoS3(s,p-1) <= 1) && (VectoS3(s,p)>=2);
        
        %Vecto_AL3t(s,:) = Vecto_AL3(s,(p-9900):(p+9999));
        %Vecto_AR3t(s,:) = Vecto_AR3(s,(p-10000):(p+9999));
        %Vecto_FL3t(s,:) = Vecto_FL3(s,(p-10000):(p+9999));
        Vecto_FR3t(s,:) = Vecto_FR3(s,(p-9900):(p+9999));
        Vecto_EMG3t(s,:) = Vecto_EMG3(s,(p-9900):(p+9999));
        VectoS3t(s,:) = VectoS3(s,(p-9900):(p+9999));
    end
 end
end



%Vecto_AL3 = Vecto_AL3t;
%Vecto_AR3 = Vecto_AR3t;
%Vecto_FL3 = Vecto_FL3t;
Vecto_FR3 = Vecto_FR3t;
Vecto_EMG3 = Vecto_EMG3t;
VectoS3 = VectoS3t;


%%
% Choose number of vector to use

con = 14;

EEG_FR = Vecto_FR3(con,:);
EMG_FR = Vecto_EMG3(con,:);
tx = size(EEG_FR,2);
time = 1:tx;

% Bandpass FIR filter for Delta (1-4 Hz)
%AlphaFilt = designfilt('bandpassfir', 'FilterOrder', 20,'CutoffFrequency1',9,'CutoffFrequency2',16, 'SampleRate',1500 );
%fvtool(h)
%DeltFilt = designfilt('lowpassiir','FilterOrder',20, 'PassbandFrequency',4,'PassbandRipple',0.2, 'SampleRate',500);
%fvtool(DeltFilt)

Lp_d = 1;
Hp_d = 4;
Ny= 400;
n = 200; % order
hd = fir1(n,[Lp_d/Ny Hp_d/Ny],'bandpass');

Delta_f = filtfilt(hd,1,EEG_FR);
Delta_h_f= hilbert(Delta_f);



% Bandpass FIR filter for Alpha (9-16 Hz)

Lp_a = 9;
Hp_a = 16;
Ny= 400;
n = 200; % order
h = fir1(n,[Lp_a/Ny Hp_a/Ny],'bandpass');

Alpha_f = filtfilt(h,1,EEG_FR);
Alph_h_f = abs(hilbert(Alpha_f));


% Spectrogram picture

movingwin = [1 0.01];
params.Fs = 1000;
params.fpass = [0 40];
params.tapers = [0.2 1];
params.trialave = 1;
params.err = 1;
maxDbG = 5; 


% Right vector length
rt = size(Vecto_FR3,2);
[ARV,Gt,Gf] = mtspecgramc(Vecto_FR3(con,:),movingwin,params);

%figure
%subplot(3,1,1)
%imagesc(Gt,Gf,10*log10(ARV'),[18 28])
%axis xy; xlabel('Time(s)'); ylabel('Freq(Hz)')
%colorbar

%
figure
subplot(7,1,1);
plot(time,EMG_FR,'k') 
axis([0,tx,-500,500]);
subplot(7,1,2);
plot(time,EEG_FR,'k') 
axis([0,tx,-700,700]);
figure
subplot(7,1,[3,4]);
imagesc(Gt,Gf,10*log10(ARV'),[18 27])
axis xy; xlabel('Time(s)'); ylabel('Freq(Hz)')
%colorbar
figure
subplot(7,1,5);
plot(time,Delta_f);
axis([0,tx,-300,300]);
subplot(7,1,6);
plot(time,Alpha_f);
hold on
plot(time,Alph_h_f);
axis([0,tx,-150,150]);
%subplot(7,1,7);
%plot(time,VectoS3);
%axis([0,tx,-3,3]);


