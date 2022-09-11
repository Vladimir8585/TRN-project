%% define sequences
y = board_adc_data;
x = amplifier_data;
%x = DD2;
%y = Y2;

%% allocate each channel to each sequence
thr = 3 ; % threshold 

EMG = x(1,:);
%AL = x(5,:);
%AR = x(2,:);
%FL = x(3,:);
FR = x(4,:);
tr = 1:length(EMG);
tk = length(EMG);

%EMG = x(6,:);
%AL = x(10,:);
%AR = x(7,:);
%FL = x(8,:);
%FR = x(9,:);
%tr = 1:length(EMG);
%tk = length(EMG);

%EMG = x(5,:);
%AL = x(9,:);
%AR = x(6,:);
%FL = x(7,:);
%FR = x(8,:);
%tr = 1:length(AL);
%tk = length(AL);


%% Calculate Spectograms

% calculate spectogram 30-150 

movingwin = [5 1];
params.Fs = 1000;
params.fpass = [30 110];
params.tapers = [3 5];
params.trialave = 1;
params.err = 1;
maxDbG = 15;

[FR38,Gt,Gf] = mtspecgramc(FR,movingwin,params);

% calculate spectogramms 0-30

movingwin = [5 1];
params.Fs = 1000;
params.fpass = [0 35];
params.tapers = [3 5];
params.trialave = 1;
params.err = 1;
maxDb = 28;

[FR03,t,f] = mtspecgramc(FR,movingwin,params);

%% EMG scorring

% Making band filter 60-200 Hz

params.Fs = 1000;
params.fpass = [60 200]; % bandpass 
params.tapers = [10 19];
params.trialave = 0;
params.pad=1;
params.err = 1; 

% Introducing 5 sec windows

hx = EMG;
Win_width = 5000; % window windth
Nr = round(length(hx)/Win_width); % making trace equal to N * window number
Nr = Nr-1; % trimming 5 last seconds
tx = Nr*Win_width; % lenght of the trace
Xt = hx(1,1:tx); % EMG trace equal to N * window number
NyLimit = 500; % limit for Power detection
vt = 1:tx; % time of the trace

drills = zeros(Nr,Win_width);
EMG_Powers = ones(1,tx);

for i = 1:Nr
    drills(i,:) = Xt(1,((i-1)*Win_width+1) :(i*Win_width));
    
a = drills(i,:); 
[S,f] = mtspectrumc( a, params );
Ss = 10*log10(S);

EMG_Powers(1,((i-1)*Win_width+1):(i*Win_width)) = max(Ss);

end

%% Delta and Theta - window power

ht = FR;
N_win = 5000; % window windth
Nrt = round(length(ht)/N_win); % making trace equal to N * window number
Nrt = Nrt-1;
txt = Nrt*N_win; % lenght of the trace
Xtt = ht(1,1:txt); %  trace equal to N * window number
vtt = 1:txt;

% Delta 

params.Fs = 1000;
params.fpass = [1 4]; 
params.tapers = [3 5];
params.trialave = 0;
params.pad=1;
params.err = 1; 


mills = zeros(Nr,N_win);
D_Powers = ones(1,tx);

for i = 1:Nr
    mills(i,:) = FR(1,((i-1)*N_win+1) :(i*N_win)); %% Frontal Right!!
    
b = mills(i,:); 
[D,f] = mtspectrumc( b, params );
Ds = 10*log10(D);
D_Powers(1,((i-1)*N_win+1):(i*N_win)) = max(Ds);

end

% Theta

params.Fs = 1000;
params.fpass = [5 8]; % bandpass
params.tapers = [3 5];
params.trialave = 0;
params.pad=1;
params.err = 1; 


vills = zeros(Nr,N_win);
T_Powers = ones(1,tx);

for i = 1:Nr
    vills(i,:) = FR(1,((i-1)*N_win+1) :(i*N_win)); %% Frontal right
    
c = vills(i,:); 
[T,f] = mtspectrumc( c, params );
Ts = 10*log10(T);
T_Powers(1,((i-1)*N_win+1):(i*N_win)) = max(Ts);

end

%% Ratio

RatioTD = T_Powers./D_Powers;

%% Start filtering Sleep from Awake by EMG
% filter 15 seconds periods

% wake traces > 15 seconds

Filt_EMG = zeros(1,Nr);
for yy = 3:Nr
    if (EMG_Powers((yy-2)*Win_width)>= thr) && (EMG_Powers((yy-1)*Win_width) >= thr) && (EMG_Powers((yy)*Win_width) >= thr)
    Filt_EMG(yy) = 1;
    end
end



for i = 1:Nr
    
EMG_PowersF(1,((i-1)*Win_width+1):(i*Win_width)) = EMG_Powers(i*Win_width)*Filt_EMG(i);

end

%% Signal trace trimming

Sig = y(1:tx);

% Plot wake trace

figure
subplot(7,1,1);
plot(vt,EMG(1:tx),'k') 
axis([0,tx,-2000,2000]);
subplot(7,1,2);
plot(vt,EMG_PowersF,'k') 
axis([0,tx,-3,30]);
subplot(7,1,3);
plot(vt,Sig,'k')
axis([0,tx,-1,4]);
subplot(7,1,4)
plot(vtt,T_Powers,'r')
hold on
plot(vtt,D_Powers,'g')
axis([0,tx,15,40]);
subplot(7,1,5);
plot(vtt,RatioTD);
axis([0,tx,0.5,1.5]);
subplot(7,1,6);
imagesc(Gt,Gf,10*log10(FR38'),[0 maxDbG])
axis xy; xlabel('Time(s)'); ylabel('Freq(Hz)')
title('Frontal EEG')


%% Differentiation of SWS and REM


Thr = 23 ; % threshold Delta
Rthr = 1.1; % threshold Theta

Filt_D = zeros(1,Nr);
Filt_T = zeros(1,Nr);
for yy = 3:Nr
    if (D_Powers((yy-2)*Win_width)>= Thr) && (D_Powers((yy-1)*Win_width) >= Thr) && (D_Powers((yy)*Win_width) >= Thr) && ...
            (EMG_Powers((yy)*Win_width)<thr) && (RatioTD((yy)*Win_width) < Rthr)
           
    Filt_D(yy) = 1;
    elseif (RatioTD((yy)*Win_width) >= Rthr) && (EMG_Powers((yy)*Win_width)<thr) && (D_Powers((yy-1)*Win_width) >= Thr) || ...
            (RatioTD((yy)*Win_width) >= Rthr) && (EMG_Powers((yy)*Win_width)<thr) && (RatioTD((yy-1)*Win_width) >= Rthr)
    Filt_T(yy) = 1;
    end
end

UpTh = 39; % Threshold for noise

for ty = 1:Nr
    if (D_Powers((ty)*Win_width) > UpTh)
        Filt_D(ty) = 0;
    end
end


% making vectors


for i = 1:Nr
    
D_PowersF(1,((i-1)*Win_width+1):(i*Win_width)) = D_Powers(i*Win_width)*Filt_D(i);
T_PowersF(1,((i-1)*Win_width+1):(i*Win_width)) = T_Powers(i*Win_width)*Filt_T(i);

end


%% Signam differentiation between Sleep and Wake
Sig = y(1:tx);

% Wake trace

Wake_ind = find(EMG_PowersF > 0);
Wake_ind_On = find(EMG_PowersF >0 & Sig >= 2);
Wake_ind_Off = find(EMG_PowersF >0 & Sig < 2);

Wake_tr = FR(Wake_ind);
Wake_tr_On = FR(Wake_ind_On);
Wake_tr_Off = FR(Wake_ind_Off);

Wake_tim = length(Wake_tr);
Wake_tim_On = length(Wake_tr_On);
Wake_tim_Off = length(Wake_tr_Off);

Twake = 1:(Wake_tim);
Twake_On = 1:(Wake_tim_On);
Twake_Off = 1:(Wake_tim_Off);

%% Sleep trace

Sleep_ind = find(EMG_PowersF == 0);
Sleep_ind_On = find(EMG_PowersF == 0 & Sig >= 2);
Sleep_ind_Off = find(EMG_PowersF == 0 & Sig < 2);

Sleep_tr = FR(Sleep_ind);
Sleep_tr_On = FR(Sleep_ind_On);
Sleep_tr_Off = FR(Sleep_ind_Off);

Sleep_tim = length(Sleep_tr);
Sleep_tim_On = length(Sleep_tr_On);
Sleep_tim_Off = length(Sleep_tr_Off);

Tsleep = 1:(Sleep_tim);
Tsleep_On = 1:(Sleep_tim_On);
Tsleep_Off = 1:(Sleep_tim_Off);


%% REM traces

REM_ind = find(T_PowersF > 0);
REM_ind_On = find(T_PowersF > 0 & Sig >=2);
REM_ind_Off = find(T_PowersF > 0 & Sig <2);

REM_tr = FR(REM_ind);
REM_tr_On = FR(REM_ind_On);
REM_tr_Off = FR(REM_ind_Off);

REM_tim = length(REM_tr);
REM_tim_On = length(REM_tr_On);
REM_tim_Off = length(REM_tr_Off);

Rsleep = 1:(REM_tim);
Rsleep_On = 1:(REM_tim_On);
Rsleep_Off = 1:(REM_tim_Off);

%% SWS traces

SWS_ind = find(D_PowersF > 0);
SWS_ind_On = find(D_PowersF > 0 & Sig >=2);
SWS_ind_Off = find(D_PowersF > 0 & Sig <2);

SWS_tr = FR(SWS_ind);
SWS_tr_On = FR(SWS_ind_On);
SWS_tr_Off = FR(SWS_ind_Off);

SWS_tim = length(SWS_tr);
SWS_tim_On = length(SWS_tr_On);
SWS_tim_Off = length(SWS_tr_Off);

Ssleep = 1:(SWS_tim);
Ssleep_On = 1:(SWS_tim_On);
Ssleep_Off = 1:(SWS_tim_Off);

%% making automatic scorage  Scorr_vector = 0 - Sleep

Scorr_vector = zeros(1,Nr);
for f = 1: Nr
    if ((Filt_EMG(f)) == 0 ) && ((Filt_T(f)) == 1)
        Scorr_vector(1,f) = 0;
    elseif ((Filt_EMG(f)) == 0 )
        Scorr_vector(1,f) = 0;
    else 
        Scorr_vector(1,f) = 1;
    end
end

Scorr = zeros(1,tx);
for i = 1:Nr
    
Scorr(1,((i-1)*Win_width+1):(i*Win_width)) = Scorr_vector(i);

end


%% Making State matrix
Scorr_vector2 = zeros(1,Nr);
for f = 1: Nr
    if ((Filt_EMG(f)) == 0 ) && ((Filt_T(f)) == 1)
        Scorr_vector2(1,f) = 2;
    elseif ((Filt_EMG(f)) == 0 )
        Scorr_vector2(1,f) = 1;
    else 
        Scorr_vector2(1,f) = 0;
    end
end

% expanding score vector to proper size

Scorr2 = zeros(1,tx);
for i = 1:Nr
    
Scorr2(1,((i-1)*Win_width+1):(i*Win_width)) = Scorr_vector2(i);

end

N_row = 0;
for g = 3:Nr
    if Scorr2((g-1)*Win_width) ~= Scorr2((g)*Win_width)
        N_row = N_row +1;
    end
end

%% Plotting signals

params.Fs = 1000;
params.fpass = [0 50]; 
params.tapers = [3 5];
params.trialave = 0;
params.pad=1;
params.err = 1;


[S1,fw] = mtspectrumc(Wake_tr, params );
W = 10*log10(S1);
W_sm = smooth(W,0.02,'lowess');


[S2,fs] = mtspectrumc(SWS_tr, params );
Sl = 10*log10(S2);
Sl_sm = smooth(Sl,0.02,'lowess');


[S3,fr] = mtspectrumc(REM_tr, params );
R = 10*log10(S3);
R_sm = smooth(R,0.02,'lowess');


figure
plot(fw,W_sm,'b','LineWidth',2)
hold on
plot(fs,Sl_sm,'g','LineWidth',2)
plot(fr,R_sm,'r','LineWidth',2)
ylim([20 28])
xlim([16 30])
set(gca,'fontsize',14)

%%

% All duration
Sl_dur = (Sleep_tim)/60000; % min
Wake_dur = (Wake_tim)/60000; % min
SWS_dur = (SWS_tim)/60000; % sec
Rem_dur = (REM_tim)/1000; % sec

% On duration

Sl_dur_On = (Sleep_tim_On)/60000; % min
Wake_dur_On = (Wake_tim_On)/60000; % min
SWS_dur_On = (SWS_tim_On)/60000; % sec
Rem_dur_On = (REM_tim_On)/1000; % sec




%% plot Nr 2

leg = size(vtt,2)

figure
subplot(7,1,1);
plot(vt,Sig,'k')
axis([0,tx,-1,4]);
subplot(7,1,3);
plot(vtt,FR(1:leg),'k')
axis([0,tx,-1200, 1200])
subplot(7,1,2)
plot(vtt,EMG(1:leg),'k')
axis([0,tx,-1200, 1200])
subplot(7,1,4);
plot(vt,EMG_PowersF) 
axis([0,tx,-3,30]);
subplot(7,1,5)
plot(vtt,EMG_PowersF) 
hold on
plot(vtt,D_PowersF,'g')
plot(vtt,T_PowersF,'r')
axis([0,tx,-5,40]);
figure
subplot(3,1,1)
plot(vtt, Scorr)
axis([0,tx,-1,3]);
subplot(3,1,3)
imagesc(t,f,10*log10(FR03'),[0 maxDb])
axis xy; ylabel('Freq(Hz)')


%% Finding episodes for sleep Ptime

%

% finding starting and ending points

Sl_start = [];
Sl_end = [];

Nz = 1;


for z = 2: tx
    if ((Scorr(z-1)) == 1 ) && ((Scorr(z)) == 0)
        Sl_start(Nz) = z;
        Nz = Nz +1;
    elseif ((Scorr(z-1)) == 0 ) && ((Scorr(z)) == 1) || ((Scorr(z-1)) == 0 ) && ((Scorr(z)) == 2)
        Sl_end(Nz) = z;

    end
end
%

% trimming unfinished points


if   (Sl_end(1))<(Sl_start(1)) 
          Sl_end(1) = [];
end

if     (Sl_end(end))<(Sl_start(end)) 
           Sl_start(end) = [];
end

size(Sl_end)
size(Sl_start)
kr = [Sl_end',Sl_start'];


% Making matrix of the starting time/ending time/difference(in sec)

P_dif = (Sl_end - Sl_start);
aP_time = (P_dif/1000)';
%Sc_mat = [Sl_start' Sl_end' P_dif'];

%% Finding epizodes for wake times


W_start = [];
W_end = [];

Pr = 1;


for z = 2: tx
    if ((Scorr(z-1)) == 0 ) && ((Scorr(z)) == 1)
        W_start(Pr) = z;
        Pr = Pr +1;
    elseif ((Scorr(z-1)) == 1 ) && ((Scorr(z)) == 0) || ((Scorr(z-1)) == 2 ) && ((Scorr(z)) == 0)
        W_end(Pr) = z;

    end
end
%

% trimming unfinished points


if   (W_end(1))<(W_start(1)) 
          W_end(1) = [];
end

if     (W_end(end))<(W_start(end)) 
           W_start (end) = [];
end

size(W_end);
size(W_start);
kz = [W_end',W_start'];


% Making matrix of the starting time/ending time/difference(in sec)

W_dif = (W_end - W_start);
Awa_time = (W_dif/1000)';
%Sc_mat = [Sl_start' Sl_end' P_dif'];

%% TO collect epizodes -  P_time and Awa_time

%FR

AC_BFR_Pfx = [];
AC_BFR_Pfx(1,1) = Sl_dur;
AC_BFR_Pfx(1,2) = Sl_dur_On;
AC_BFR_Pfx(2,1) = Wake_dur;
AC_BFR_Pfx(2,2) = Wake_dur_On;
AC_BFR_Pfx(3,1) = SWS_dur;
AC_BFR_Pfx(3,2) = SWS_dur_On;
AC_BFR_Pfx(4,1) = Rem_dur;
AC_BFR_Pfx(4,2) = Rem_dur_On;

printmat(AC_BFR_Pfx , 'States', 'Sleep Wake SWS REM', 'Overall On_Only ')
%title('Difference before and after stimulation')


