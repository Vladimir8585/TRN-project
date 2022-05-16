%% define sequences

thr = 3 ; % threshold 

y = board_adc_data;
x = amplifier_data;
%x = DD2;
%y = Y2;

%% allocate each channel to each sequence

%EMG = x(1,:);
AL = x(4,:);
AR = x(1,:);
FL = x(2,:);
FR = x(3,:);
tr = 1:length(AL);
tk = length(AL);

%EMG = x(6,:);
%AL = x(10,:);
%AR = x(7,:);
%FL = x(8,:);
%FR = x(9,:);
%tr = 1:length(EMG);
%tk = length(EMG);

%EMG = x(6,:);
%AL = x(8,:);
%AR = x(5,:);
%FL = x(6,:);
%FR = x(7,:);
%tr = 1:length(AL);
%tk = length(AL);







%% calculate spectogram 30-150 

movingwin = [5 1];
params.Fs = 1000;
params.fpass = [30 110];
params.tapers = [3 5];
params.trialave = 1;
params.err = 1;
maxDbG = 15;

[AL38,Gt,Gf] = mtspecgramc(AL,movingwin,params);
[AR38,Gt,Gf] = mtspecgramc(AR,movingwin,params);
[FL38,Gt,Gf] = mtspecgramc(FL,movingwin,params);
[FR38,Gt,Gf] = mtspecgramc(FR,movingwin,params);

%% calculate spectogramms 0-30

movingwin = [5 1];
params.Fs = 1000;
params.fpass = [0 35];
params.tapers = [3 5];
params.trialave = 1;
params.err = 1;
maxDb = 30;

[AL03,~,f] = mtspecgramc(AL,movingwin,params);
[AR03,t,f] = mtspecgramc(AR,movingwin,params);
[FL03,t,f] = mtspecgramc(FL,movingwin,params);
[FR03,t,f] = mtspecgramc(FR,movingwin,params);

%% EMG scorring

%% Making band filter 60-200 Hz

params.Fs = 1000;
params.fpass = [60 200]; % bandpass 
params.tapers = [10 19];
params.trialave = 0;
params.pad=1;
params.err = 1; 

%%
%Sleep scoring code for auditory traces

hx = FR;
Win_width = 5000; % window windth
Nr = round(length(hx)/Win_width); % making trace equal to N * window number
Nr = Nr-1
tx = Nr*Win_width; % lenght of the trace
Xt = hx(1,1:tx); %  trace equal to N * window number
NyLimit = 500; % limit for Power detection
vt = 1:tx; % time of the trace

%%

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
Nrt = Nrt-1
txt = Nrt*N_win; % lenght of the trace
Xtt = ht(1,1:txt); %  trace equal to N * window number
vtt = 1:txt;

%% Delta 

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

%% Theta

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


Filt_EMG = zeros(1,Nr);
for yy = 3:Nr
    if (EMG_Powers((yy-2)*Win_width)>= thr) && (EMG_Powers((yy-1)*Win_width) >= thr) && (EMG_Powers((yy)*Win_width) >= thr)
    Filt_EMG(yy) = 1;
    end
end

% making On/Off EMG trace


for i = 1:Nr
    
EMG_PowersF(1,((i-1)*Win_width+1):(i*Win_width)) = EMG_Powers(i*Win_width)*Filt_EMG(i);

end

%% Signal trace trimming

Sig = y(1:tx);

% Wake trace

Wake_ind = find(EMG_PowersF > 0);
Wake_ind_On = find(EMG_PowersF >0 & Sig >2);
Wake_ind_Off = find(EMG_PowersF >0 & Sig <2);

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
Sleep_ind_On = find(EMG_PowersF == 0 & Sig >2);
Sleep_ind_Off = find(EMG_PowersF == 0 & Sig <2);

Sleep_tr = FR(Sleep_ind);
Sleep_tr_On = FR(Sleep_ind_On);
Sleep_tr_Off = FR(Sleep_ind_Off);

Sleep_tim = length(Sleep_tr);
Sleep_tim_On = length(Sleep_tr_On);
Sleep_tim_Off = length(Sleep_tr_Off);

Tsleep = 1:(Sleep_tim);
Tsleep_On = 1:(Sleep_tim_On);
Tsleep_Off = 1:(Sleep_tim_Off);

% Find Ratio and Delta trace during sleep

SleepD = D_Powers(Sleep_ind);
SleepT = T_Powers(Sleep_ind);
Sleep_RatioTD = RatioTD(Sleep_ind);

%%
%% Signal trace trimming

Sig = y(1:tx);

% Plot wake trace

figure
subplot(7,1,1);
plot(vt,AL(1:tx)) 
axis([0,tx,-700,700]);
subplot(7,1,2);
plot(vt,EMG_PowersF) 
axis([0,tx,-3,30]);
subplot(7,1,3);
plot(vt,Sig,'k')
axis([0,tx,-1,4]);
subplot(7,1,4)
plot(vtt,T_Powers,'r')
hold on
plot(vtt,D_Powers)
axis([0,tx,15,40]);
subplot(7,1,5);
plot(vtt,RatioTD);
axis([0,tx,0.5,1.5]);
subplot(7,1,6);
imagesc(Gt,Gf,10*log10(FR38'),[0 maxDbG])
axis xy; xlabel('Time(s)'); ylabel('Freq(Hz)')
title('Frontal EEG')
%colorbar;
subplot(7,1,7);
imagesc(t,f,10*log10(FR03'),[0 maxDb])
axis xy; ylabel('Freq(Hz)')
%colorbar;


%% Differentiation of SWS and REM


Thr = 25 ; % threshold 
Rthr = 1.1
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

UpTh = 39;

for ty = 1:Nr
    if (D_Powers((ty)*Win_width) > UpTh);
        Filt_D(ty) = 0;
    end
end


% making vectors


for i = 1:Nr
    
D_PowersF(1,((i-1)*Win_width+1):(i*Win_width)) = D_Powers(i*Win_width)*Filt_D(i);
T_PowersF(1,((i-1)*Win_width+1):(i*Win_width)) = T_Powers(i*Win_width)*Filt_T(i);

end


%% REM traces

REM_ind = find(T_PowersF > 0);
REM_ind_On = find(T_PowersF > 0 & Sig >2);
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
SWS_ind_On = find(D_PowersF > 0 & Sig >2);
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

%% making automatic scorage 

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

%  Making scorring matrix

Scorr_time = zeros(N_row,3);
N_rows = 0;
for p = 3:Nr
    if Scorr2((p-1)*Win_width) ~= Scorr2((p)*Win_width)
        N_rows = N_rows +1;
        Scorr_time(N_rows,1) = (p)*Win_width;
        Scorr_time(N_rows,2) = Scorr2((p-1)*Win_width);
        Scorr_time(N_rows,3) = Scorr2((p)*Win_width);
    end
end


%%  Detect staring ponts

rl = length(Sig);

Scor_sig = zeros(1,rl);

for p = 2: rl
    if (Sig(1,p-1) <= 1) && (Sig(1,p)>=2);
        Scor_sig(p) = 1;
    elseif (Sig(1,p-1) >= 2) && (Sig(1,p)<=2);
        Scor_sig(p) = 2;
    end
end

Up = find(Scor_sig ==1);
Down = find(Scor_sig ==2);

%% trimming unfinished signal parts

if   (Down(1))<(Up(1)) 
          Scor_sig(Down(1)) = 0;
elseif     (Down(end))<(Up(end))
           Scor_sig(Up(end)) = 0;
end

%% 
N2_rows = 2*length(Up);
Scorr_sig2 = zeros(N2_rows,2);
N_rt = 0;
for b = 1: rl
    if Scor_sig(b) ~= 0
        N_rt = N_rt +1;
        Scorr_sig2(N_rt,1) = b;
        Scorr_sig2(N_rt,2) = Scor_sig(b);
    end
end
    

%% Differentiate between time
time_r = zeros(N2_rows,1);
pf = 0;
for j = 1: N2_rows 
    if (Scorr_sig2(j,2))==2
        time_r(j,1) = Scorr_sig2(j,1)- Scorr_sig2(j-1,1);
    elseif (Scorr_sig2(j,2))==1
        pf = pf+1;
        time_r(j,1) = pf;
        
    end
end

%% combine matrixes

Scorr_sig3 = [Scorr_sig2(:,1) Scorr_sig2(:,2) time_r];


%% making shure that we have enough space between

if Scorr_sig3(1,1) < 60000
    Scorr_sig3(2,:) = [];
    Scorr_sig3(1,:) = [];
elseif Scorr_sig3(end,1) > (txt- 60000)
    Scorr_sig3(end-1,:) = [];
    Scorr_sig3(end,:) = [];
end

N3_rows = size(Scorr_sig3);
%% Automatic sleep and awake signal detection

Vecto_AL6 = [];
Vecto_AR6 = [];
Vecto_FL6 = [];
Vecto_FR6 = [];
Vecto_EMG6 = [];
VectoS6 = [];


Vecto_AL3 = [];
Vecto_AR3 = [];
Vecto_FL3 = [];
Vecto_FR3 = [];
Vecto_EMG3 = [];
VectoS3 = [];


Vecto_AL6w = [];
Vecto_AR6w = [];
Vecto_FL6w = [];
Vecto_FR6w = [];
Vecto_EMG6w = [];
VectoS6w = [];


Vecto_AL3w = [];
Vecto_AR3w = [];
Vecto_FL3w = [];
Vecto_FR3w = [];
Vecto_EMG3w = [];
VectoS3w = [];


Nk60 = 0;
Nk30 = 0;
Nk60w = 0;
Nk30w = 0;

Order60 = [];
Order30 = [];
Order60w = [];
Order30w = [];

Stim_order6 = [];
Stim_order3 = [];
Stim_order6w = [];
Stim_order3w = [];



for a = 1:N3_rows-2
    if Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=59999 && ...
            sum(Scorr((Scorr_sig3(a,1)-60000):(Scorr_sig3(a,1)+120000))) == 0;
        Nk60 = Nk60+1;
        
        Vecto_AL6(Nk60,:) = AL(1,(Scorr_sig3(a,1)-60000):(Scorr_sig3(a,1)+119999));
        Vecto_AR6(Nk60,:) = AR(1,(Scorr_sig3(a,1)-60000):(Scorr_sig3(a,1)+119999));
        Vecto_FL6(Nk60,:) = FL(1,(Scorr_sig3(a,1)-60000):(Scorr_sig3(a,1)+119999));
        Vecto_FR6(Nk60,:) = FR(1,(Scorr_sig3(a,1)-60000):(Scorr_sig3(a,1)+119999));
        Vecto_EMG6(Nk60,:) = AL(1,(Scorr_sig3(a,1)-60000):(Scorr_sig3(a,1)+119999));
        VectoS6(Nk60,:) = Sig(1,(Scorr_sig3(a,1)-60000):(Scorr_sig3(a,1)+119999));
        
        Stim_order6(Nk60,1) = (a-1)/2 +1;
        Order60(Nk60,1) = Scorr_sig3(a,1);
        
    elseif Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&...
            sum(Scorr((Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)+60000))) == 0;
        Nk30 = Nk30+1;
        
        Vecto_AL3(Nk30,:) = AL(1,(Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)+59999));
        Vecto_AR3(Nk30,:) = AR(1,(Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)+59999));
        Vecto_FL3(Nk30,:) = FL(1,(Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)+59999));
        Vecto_FR3(Nk30,:) = FR(1,(Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)+59999));
        Vecto_EMG3(Nk30,:) = AL(1,(Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)+59999));
        VectoS3(Nk30,:) = Sig(1,(Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)+59999));
        
        Stim_order3(Nk30,1) = (a-1)/2 +1;
        Order30(Nk30,1) = Scorr_sig3(a,1);
        
    elseif Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=59999 && ...
            sum(Scorr((Scorr_sig3(a,1)-60000):(Scorr_sig3(a,1)+120000))) >= 179999;
        Nk60w = Nk60w+1;
        
        Vecto_AL6w(Nk60w,:) = AL(1,(Scorr_sig3(a,1)-60000):(Scorr_sig3(a,1)+119999));
        Vecto_AR6w(Nk60w,:) = AR(1,(Scorr_sig3(a,1)-60000):(Scorr_sig3(a,1)+119999));
        Vecto_FL6w(Nk60w,:) = FL(1,(Scorr_sig3(a,1)-60000):(Scorr_sig3(a,1)+119999));
        Vecto_FR6w(Nk60w,:) = FR(1,(Scorr_sig3(a,1)-60000):(Scorr_sig3(a,1)+119999));
        Vecto_EMG6w(Nk60w,:) = AL(1,(Scorr_sig3(a,1)-60000):(Scorr_sig3(a,1)+119999));
        VectoS6w(Nk60w,:) = Sig(1,(Scorr_sig3(a,1)-60000):(Scorr_sig3(a,1)+119999));
        
        Stim_order6w(Nk60w,1) = (a-1)/2 +1;
        Order60w(Nk60w,1) = Scorr_sig3(a,1);
        
    elseif Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&...
            sum(Scorr((Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)+60000))) >= 89999;
        Nk30w = Nk30w+1;
        
        Vecto_AL3w(Nk30w,:) = AL(1,(Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)+59999));
        Vecto_AR3w(Nk30w,:) = AR(1,(Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)+59999));
        Vecto_FL3w(Nk30w,:) = FL(1,(Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)+59999));
        Vecto_FR3w(Nk30w,:) = FR(1,(Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)+59999));
        Vecto_EMG3w(Nk30w,:) = AL(1,(Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)+59999));
        VectoS3w(Nk30w,:) = Sig(1,(Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)+59999));
        
        Stim_order3w(Nk30w,1) = (a-1)/2 +1;
        Order30w(Nk30w,1) = Scorr_sig3(a,1);

    end
end

%% making scatter to work

if size(Order30,1) == 0;
    Z30 = 0;
else
    Z30 = zeros(size(Order30,1),1);
end

if size(Order60,1) == 0;
    Z60 = 0;
else
    Z60 = zeros(size(Order60,1),1);
end

if size(Order30w,1) == 0;
    Z30w = 0;
else
    Z30w = zeros(size(Order30w,1),1);
end

if size(Order60w,1) == 0;
    Z60w = 0;
else
    Z60w = zeros(size(Order60w,1),1);
end

if isempty(Order60w)
    Order60w = 0;
end
if isempty(Order30w)
    Order30w = 0;
end
    
if isempty(Order30)
    Order30  = 0;
end

if isempty(Order60)
    Order60  = 0;
end


%% plot Nr 2

figure
subplot(4,1,1);
plot(vt,EMG_PowersF) 
hold on
plot(vtt,D_PowersF,'g')
plot(vtt,T_PowersF,'r')
axis([0,tx,-5,40]);
subplot(4,1,2)
plot(vtt,Sig);
hold on
scatter(Order30,Z30,'g')
scatter(Order60,Z60,'r')
scatter(Order30w,Z30w,'k')
scatter(Order60w,Z60w,'m') 
axis([0,tx,-1,4]);
subplot(4,1,3)
plot(vtt, Scorr)
axis([0,tx,-1,3]);
subplot(4,1,4)
imagesc(t,f,10*log10(FR03'),[0 maxDb])
axis xy; ylabel('Freq(Hz)')


%% Matrix which says how many traces we got

Scorr_time

Stim_order3
Stim_order3w 
Stim_order6 
Stim_order6w

Scorr_sig3


M_Numbers = []; %zeros(30,4);
M_Numbers(1,1) = size(Order30,1);
M_Numbers(1,2) = size(Order30w,1);
M_Numbers(1,3) = size(Order60,1);
M_Numbers(1,4) = size(Order60w,1);

printmat(M_Numbers, 'Number of detected traces', 'Number ', 'Sleep_30 Wake_30 Sleep_60 Wake_60 ')


