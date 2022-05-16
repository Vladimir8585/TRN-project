% I did not change properly this code for a while, it got parts which are
% not related to interested outcome. I used to use it for different
% outcomes brakets will let you know that this part mostly is not required
% []



%% define sequences
%y = board_adc_data;
%x = amplifier_data;
x = DD2;
y = Y2;

%% allocate each channel to each sequence
thr = 3 ; % threshold for EMG (sleep/wake detection)

EMG = x(1,:);
AL = x(4,:);
%AR = x(2,:);
FL = x(2,:);
FR = x(3,:);
tr = 1:length(EMG);
tk = length(EMG);

%EMG = x(6,:);
%AL = x(10,:);
%AR = x(7,:);
%FL = x(8,:);
%FR = x(9,:);
%tr = 1:length(EMG);
%tk = length(EMG);

%EMG = x(6,:);
%AL = x(9,:);
%AR = x(7,:);
%FL = x(8,:);
%FR = x(9,:);
%tr = 1:length(AR);
%tk = length(AR);


%% calculate spectogram 30-110  [not required  for all channels]

movingwin = [5 1];
params.Fs = 1000;
params.fpass = [30 110];
params.tapers = [3 5];
params.trialave = 1;
params.err = 1;
maxDbG = 15;

[AL38,Gt,Gf] = mtspecgramc(AL,movingwin,params);
%[AR38,Gt,Gf] = mtspecgramc(AR,movingwin,params);
[FL38,~,~] = mtspecgramc(FL,movingwin,params);
[FR38,~,~] = mtspecgramc(FR,movingwin,params);

%% calculate spectogramms 0-30 [not required  for all channels]

movingwin = [5 1];
params.Fs = 1000;
params.fpass = [0 35];
params.tapers = [3 5];
params.trialave = 1;
params.err = 1;
maxDb = 30;

[AL03,t,f] = mtspecgramc(AL,movingwin,params);
%[AR03,t,f] = mtspecgramc(AR,movingwin,params);
[FL03,~,~] = mtspecgramc(FL,movingwin,params);
[FR03,~,~] = mtspecgramc(FR,movingwin,params);

%% EMG scorring

%% Making band filter 60-200 Hz

params.Fs = 1000;
params.fpass = [60 200]; % bandpass 
params.tapers = [10 19];
params.trialave = 0;
params.pad=1;
params.err = 1; 

%%
%Sleep scoring (EMG windows and powers for each windows)

hx = EMG;
Win_width = 5000; % window windth
Nr = round(length(hx)/Win_width); % making trace equal to N * window number
Nr = Nr-1 
tx = Nr*Win_width; % lenght of the trace
Xt = hx(1,1:tx); %  trace equal to N * window number
NyLimit = 500; % limit for Power detection
vt = 1:tx; % time of the trace

%% Detect power of each window and expand it to the size of natural trace

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
% if EMG power > threshold for 15 seconds, I count it as a sleep


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

% Wake trace    [ not required, were used in different circumstances]

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

%% Sleep trace [ not required, were used in different circumstances]

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

% Plot wake trace

figure
subplot(7,1,1);
plot(vt,EMG(1:tx)) 
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


Thr = 25 ; % threshold for the Delta power
Rthr = 1.1 % threshold for Delta/Theta ratio
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


%% REM traces [not required in the interested outcome]

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

%% SWS traces [not required in the interested outcome]

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

%% making automatic scorage for sleep state change detection (combining REM and SWS to sleep vector)

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
%% making automatic scorage for Awake state chnage detection

Scorr_vectoR = zeros(1,Nr);
for f = 1: Nr
    if ((Filt_EMG(f)) == 0 ) && ((Filt_T(f)) == 1)
        Scorr_vectoR(1,f) = 1;
    elseif ((Filt_EMG(f)) == 0 )
        Scorr_vectoR(1,f) = 1;
    else 
        Scorr_vectoR(1,f) = 0;
    end
end

ScorR = zeros(1,tx);
for i = 1:Nr
    
ScorR(1,((i-1)*Win_width+1):(i*Win_width)) = Scorr_vectoR(i);

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

% expanding score vector to proper size of the trace [not used in this
% outcome I guess, but produce score matrix]

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


%%  Detect staring ponts and ending points of the light output

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

%% making score matrix with signal incorporated
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
    

%% Differentiate between time (how long was light signal was applied for)
time_s = zeros(N2_rows/2,1);
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


%% combine matrixes [when the is the mistake , code stops here...and I go over Scorr_sig2
%

Scorr_sig3 = [Scorr_sig2(:,1) Scorr_sig2(:,2) time_r];
%find(Scorr_sig2(:,1) ==   627647)
%Scorr_sig2(73,:) =[]

%% making shure that we have enough space at the begining and at the end

if Scorr_sig3(1,1) < 60000
    Scorr_sig3(2,:) = [];
    Scorr_sig3(1,:) = [];
elseif Scorr_sig3(end,1) > (txt- 60000)
    Scorr_sig3(end-1,:) = [];
    Scorr_sig3(end,:) = [];
end

N3_rows = size(Scorr_sig3);
%% Automatic sleep and awake signal detection (collecting vectors 
% which agree with inclusion criteria

Vecto_AL3 = [];
Vecto_AR3 = [];
Vecto_FL3 = [];
Vecto_FR3 = [];
Vecto_EMG3 = [];
VectoS3 = [];

Vecto_AL3f = [];
Vecto_AR3f = [];
Vecto_FL3f = [];
Vecto_FR3f = [];
Vecto_EMG3f = [];
VectoS3f = [];

Vecto_AL3n = [];
Vecto_AR3n = [];
Vecto_FL3n = [];
Vecto_FR3n = [];
Vecto_EMG3n = [];
VectoS3n = [];


Nk30 = 0;
Nk30f = 0;
Nk30n = 0;

Zk30 = 0;
Zk30f = 0;
Zk30n = 0;


Order30 = [];
Order30f = [];
Order30n = [];


Stim_order3 = [];
Stim_order3f = [];
Stim_order3n = [];


K30 = [];
K30f = [];
K30n = [];

Kp30 = [];
Kp30f = [];
Kp30n = [];

for a = 1:N3_rows-2
     
        
    if Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&... %% for 30 sec stimulation (same stuff)
            sum(Scorr((Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)))) == 0 && ...
            sum(Scorr((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000))) ~= 0;
        Nk30 = Nk30+1;
        
        K30(Nk30,1) = find(Scorr((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000)),1,'first');
        
        Vecto_AL3(Nk30,:) = AL(1,(Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)+59999));
        %Vecto_AR3(Nk30,:) = AR(1,(Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)+59999));
        Vecto_FL3(Nk30,:) = FL(1,(Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)+59999));
        Vecto_FR3(Nk30,:) = FR(1,(Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)+59999));
        Vecto_EMG3(Nk30,:) = EMG(1,(Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)+59999));
        VectoS3(Nk30,:) = Sig(1,(Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)+59999));
        
        Stim_order3(Nk30,1) = (a-1)/2 +1;
        K30(Nk30,2)= (a-1)/2 +1;
        Order30(Nk30,1) = Scorr_sig3(a,1);
        
    elseif Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&... %% 
            sum(Scorr((Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)))) == 0 && ...
            sum(Scorr((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000))) == 0;
        Zk30 = Zk30+1;
        Kp30(Zk30,1)=61000; 
        Kp30(Zk30,2)= (a-1)/2 +1;
    end
end

for a = 1:N3_rows-2
        
    if Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&... %% 
            sum(Scorr((Scorr_sig3(a,1)-20000):(Scorr_sig3(a,1)))) == 0 && ...
            sum(Scorr((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000))) ~= 0;
        Nk30f = Nk30f+1;
        
        K30f(Nk30f,1) = find(Scorr((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000)),1,'first');
        
        Vecto_AL3f(Nk30f,:) = AL(1,(Scorr_sig3(a,1)-20000):(Scorr_sig3(a,1)+59999));
        %Vecto_AR3f(Nk30f,:) = AR(1,(Scorr_sig3(a,1)-20000):(Scorr_sig3(a,1)+59999));
        Vecto_FL3f(Nk30f,:) = FL(1,(Scorr_sig3(a,1)-20000):(Scorr_sig3(a,1)+59999));
        Vecto_FR3f(Nk30f,:) = FR(1,(Scorr_sig3(a,1)-20000):(Scorr_sig3(a,1)+59999));
        Vecto_EMG3f(Nk30f,:) = EMG(1,(Scorr_sig3(a,1)-20000):(Scorr_sig3(a,1)+59999));
        VectoS3f(Nk30f,:) = Sig(1,(Scorr_sig3(a,1)-20000):(Scorr_sig3(a,1)+59999));
        
        Stim_order3f(Nk30f,1) = (a-1)/2 +1;
        K30f(Nk30f,2) = (a-1)/2 +1;
        Order30f(Nk30f,1) = Scorr_sig3(a,1);
        
    elseif Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&... %% 
            sum(Scorr((Scorr_sig3(a,1)-20000):(Scorr_sig3(a,1)))) == 0 && ...
            sum(Scorr((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000))) == 0;
        Zk30f = Zk30f+1;
        Kp30f(Zk30f,1) = 61000;
        Kp30f(Zk30f,2)= (a-1)/2 +1;
    end
end

for a = 1:N3_rows-2

        
    if Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&... %% 
            sum(ScorR((Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)))) == 0 && ...
            sum(ScorR((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000))) ~= 0;
        Nk30n = Nk30n+1;
        
        K30n(Nk30n,1) = find(ScorR((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000)),1,'first');
        
        Vecto_AL3n(Nk30n,:) = AL(1,(Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)+59999));
        %Vecto_AR3n(Nk30n,:) = AR(1,(Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)+59999));
        Vecto_FL3n(Nk30n,:) = FL(1,(Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)+59999));
        Vecto_FR3n(Nk30n,:) = FR(1,(Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)+59999));
        Vecto_EMG3n(Nk30n,:) = EMG(1,(Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)+59999));
        VectoS3n(Nk30n,:) = Sig(1,(Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)+59999));
        
        Stim_order3n(Nk30n,1) = (a-1)/2 +1;
        K30n(Nk30n,2) = (a-1)/2 +1;
        Order30n(Nk30n,1) = Scorr_sig3(a,1);
        
    elseif Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&... %% 
            sum(ScorR((Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)))) == 0 && ...
            sum(ScorR((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000))) == 0;
        Zk30n = Zk30n+1;
        Kp30n(Zk30n,1) = 61000;
        Kp30n(Zk30n,2)= (a-1)/2 +1;
          
    end
end

%% making scatter (circles representing signal start) to work 

if size(Order30,1) == 0;
    Z30 = 0;
else
    Z30 = zeros(size(Order30,1),1);
end

if size(Order30f,1) == 0;
    Z30f = 0;
else
    Z30f = zeros(size(Order30f,1),1);
end

if size(Order30n,1) == 0;
    Z30n = 0;
else
    Z30n = zeros(size(Order30n,1),1);
end


if isempty(Order30f)
    Order30f = 0;
end
    
if isempty(Order30)
    Order30  = 0;
end

if isempty(Order30n)
    Order30n = 0;
end


%% Ratio - Calculating sleeping /wake time (ratio) during stimulation/after stimulation,
Pef = find(Scorr_sig3(:,2) == 2);  
Polo = Scorr_sig3(Pef,3);

fif = size(Up,2);
Pk = 0;
Ps = [];


for kk = 1:fif-2
    
    Up1(kk) = numel(find(Scorr(Up(kk):(Up(kk)+(Polo(kk)/2)))==0));
    Up2(kk) = numel(find(Scorr(Up(kk)+(Polo(kk)/2):(Up(kk)+(Polo(kk))))==0));
    Do1(kk) = numel(find(Scorr(Down(kk):(Down(kk) + (Polo(kk)/2)))==0));
    Do2(kk) = numel(find(Scorr((Down(kk)+(Polo(kk)/2)):(Down(kk) + Polo(kk)))==0));
    
    Us(kk) = numel(find(Scorr(Up(kk):(Up(kk)+Polo(kk)))==0));
    Uw(kk) = numel(find(Scorr(Up(kk):(Up(kk)+Polo(kk)))==1));
    D1s(kk) = numel(find(Scorr(Down(kk):(Down(kk) + Polo(kk)))==0));
    D1w(kk) = numel(find(Scorr(Down(kk):(Down(kk) + Polo(kk)))==1));
    D2s(kk) = numel(find(Scorr((Down(kk)+Polo(kk)):(Down(kk)+2*Polo(kk)))==0));
    D2w(kk) = numel(find(Scorr((Down(kk)+Polo(kk)):(Down(kk)+2*Polo(kk)))==1));
    
    %ratio
    
    
   
    rUp1(kk) = Up1(kk)/(Polo(kk)/2);
    rUp2(kk) = Up2(kk)/(Polo(kk)/2);
    rDo1(kk) = Do1(kk)/(Polo(kk)/2);
    rDo2(kk) = Do2(kk)/(Polo(kk)/2);
    
    
    rUs(kk) = Us(kk)/Polo(kk);
    rUw(kk) = Uw(kk)/Polo(kk);
    rD1s(kk) = D1s(kk)/Polo(kk);
    rD1w(kk) = D1w(kk)/Polo(kk);
    rD2s(kk) = D2s(kk)/Polo(kk);
    rD2w(kk) = D2w(kk)/Polo(kk);
    
    Pk = Pk +1;
    Ps(kk) = Pk;
    
end

%% mean ration  [Us Uw ; D1s D1w; D2s D2w]

MeanRatio = [mean(rUs) mean(rUw) mean(rD1s) mean(rD1w) mean(rD2s) mean(rD2w) mean(rUp1) mean(rUp2) mean(rDo1) mean(rDo2)];
MeanTime = [mean(Us) mean(Uw) mean(D1s) mean(D1w) mean(D2s) mean(D2w) ];
NumberAll = [Us' Uw' D1s' D1w' D2s' D2w' Up1' Up2' Do1' Do2' Ps'];
% final ration meanUs/meanD1s 
Fration = [mean(Us)/mean(D1s) mean(Uw)/mean(D1w) mean(Us)/mean(D2s) mean(Uw)/mean(D2w)];


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
scatter(Order30f,Z30f,'r')
scatter(Order30n,Z30n,'m')
axis([0,tx,-1,4]);
subplot(4,1,3)
plot(vtt, Scorr)
axis([0,tx,-1,3]);
subplot(4,1,4)
imagesc(t,f,10*log10(FR03'),[0 maxDb])
axis xy; ylabel('Freq(Hz)')


%% Matrixes which portray how many traces we got


Ratios = [];

Ratios(1,1) = MeanRatio(1);
Ratios(1,2) = MeanRatio(2);
Ratios(1,3) = MeanRatio(3);
Ratios(1,4) = MeanRatio(4);
Ratios(1,5) = MeanRatio(5);
Ratios(1,6) = MeanRatio(6);

Ratios(2,1) = MeanTime(1);
Ratios(2,2) = MeanTime(2);
Ratios(2,3) = MeanTime(3);
Ratios(2,4) = MeanTime(4);
Ratios(2,5) = MeanTime(5);
Ratios(2,6) = MeanTime(6);

printmat(Ratios, 'Sleep/wake ratio', 'MeanRatio MeanTime ', 'UnderS UnderW Off1_S Off1_W Off2_S Off2_W ')

%if size(K60,1) == 0;
%    K60(:,1) = 0;
%end
%if size(K60f,1) == 0;
%    K60f(:,1) = 0;
%end
%if size(K60n,1) == 0;
%    K60n(:,1) = 0;
%end

if size(K30,1) == 0;
    K30(:,1) = 0;
end
if size(K30f,1) == 0;
    K30f(:,1) = 0;
end
if size(K30n,1) == 0;
    K30n(:,1) = 0;
end


if size(Kp30,1) == 0;
    Kp30(:,1) = 0;
end
if size(Kp30f,1) == 0;
    Kp30f(:,1) = 0;
end
if size(Kp30n,1) == 0;
    Kp30n(:,1) = 0;
end


Pik = [];

Pik(1,1) = numel(K30(:,1));
Pik(1,2) = numel(Kp30(:,1));
Pik(1,3) = numel(K30f(:,1));
Pik(1,4) = numel(Kp30f(:,1));
Pik(1,5) = numel(K30n(:,1));
Pik(1,6) = numel(Kp30n(:,1));


printmat(Pik, 'Scorred', 'N',  'K30 Kp30 K30f Kp30f K30n Kp30n ')
