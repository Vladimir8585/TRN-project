
%VectoS3 = VectoS;
%Vecto_FR3 = VectoFR;
%Vecto_FR3 = Vecto_FL3;
%% Spindle scorring with Halassa

%% Signal trimming

trim_min = 29999;
trim_plus = 14999;
dend = trim_min + trim_plus;

ptt = size(VectoS3,1);
pt = size(VectoS3,2);

for s = 1:ptt
 for p = 2: pt
    if (VectoS3(s,p-1) <= 1) && (VectoS3(s,p)>=2);
        
        Vecto_FR3t(s,:) = Vecto_FR3(s,(p-trim_min):(p+trim_plus));
        Vecto_EMG3t(s,:) = Vecto_EMG3(s,(p-trim_min):(p+trim_plus));
        VectoS3t(s,:) = VectoS3(s,(p-trim_min):(p+trim_plus));
    end
 end
end

Vecto_FR3 = Vecto_FR3t;
Vecto_EMG3 = Vecto_EMG3t;
VectoS3 = VectoS3t;

%% Halassa spindle detection

ktk = size(Vecto_FR3,1);
gr = size(Vecto_FR3,2);
TTg = 1:gr;

% Firr filter

Lp = 9;
Hp = 15;
Ny= 500;
n = 400; % order
h = fir1(n,[Lp/Ny Hp/Ny],'bandpass');

fFR = zeros(ktk,gr);

for df = 1:ktk 

    fFR(df,:) = filtfilt(h,1,Vecto_FR3(df,:));
end

% Hilbert and envelope

en_fFR = abs(hilbert(fFR'));

% Smoothing
Nv = 1000; % smoothing for Nv seconds
coef = ones(1,Nv)/Nv;
delay = (length(coef)-1)/2;

cl_fFR = filter(coef,1,en_fFR);

% Mean and SD for threshold
trim_raw = 15000; % Mean + SD collection

% trial by trial mean + SD
fR_m = mean(cl_fFR(5000:trim_raw,:));  
fR_sd = std(cl_fFR(5000:trim_raw,:));
meansd_fR = fR_m+fR_sd;

% Overall mean + SD

O_fR_m = mean(cl_fFR(1:trim_raw)); 
O_fR_sd = std(cl_fFR(1:trim_raw));
O_mean_sd_fR = O_fR_m + O_fR_sd;


fR_thr = zeros(gr,ktk);
fR_thr_c = zeros(gr,ktk);

for pol = 1:ktk

fR_thr(1:gr,pol) = meansd_fR(pol);
fR_thr_c(1:gr,pol) = O_mean_sd_fR;

end

%%  Find areas higher than threshold for changing threshold

area_ch = zeros(gr,ktk);

for krb = 1:ktk
    for zd = 1: gr
        if   cl_fFR(zd,krb) > meansd_fR(krb) 
        area_ch(zd,krb) = 1;
        else 
        area_ch(zd,krb) = 0;
        end
    end
end

%%   Detect staring and ending points


Setpoints = zeros(gr,ktk);
Ups = zeros(ktk);
Downs = zeros(ktk);
up_count = 0;
down_count =0;

for trz = 1:ktk
    for len = 2: gr
            if (area_ch(len-1,trz) == 0) && (area_ch(len,trz) == 1)
                up_count = up_count + 1;
                Setpoints(len,trz) = 1;
                Ups(up_count,trz) = len;
            elseif (area_ch(len-1,trz) == 1) && (area_ch(len,trz) == 0)
                down_count = down_count + 1;
                Setpoints(len,trz) = 2;
                Downs(down_count,trz) = len;
            end
    end
end

%% %% trimming unfinished signal parts and cleaning it from too short (500) or too long signals (3000)

for gen = 1:ktk
    up = Ups(:,gen);
    down = Downs(:,gen);
    up(up==0) = [];
    down(down==0) = [];
        if (down(1) < up(1))  
            Setpoints(down(1),gen) = 0;
            down(1) =[];
        elseif (down(end) < up(end))
            Setpoints(up(end),gen) = 0;
            up(end) =[];
        end

    difr =  down - up;

        for ip = 1:length(difr)

            if ((difr(ip)>=500) && (difr(ip)<=3000))
            Setpoints(down(ip),gen) = 0;
            else 
            Setpoints(up(ip),gen) = 0;
            Setpoints(down(ip),gen) = 0;
            end  

        end
end

%%

con = 13

figure
subplot(7,1,1)
plot(TTg,Vecto_EMG3(con,:),'k');
xlim([0 gr]);
%title('Raw EMG')
subplot(7,1,2)
plot(TTg,Vecto_FR3(con,:),'k');
xlim([0 gr]);
%title('Raw EEG')
subplot(7,1,3)
plot(TTg,fFR(con,:));
hold on
plot(TTg,en_fFR(:,con))
xlim([0 gr]);
%title('Fir filter/envelope')
figure
subplot(7,1,4)
plot((TTg-delay/Nv),cl_fFR(:,con));
hold on
pd = plot(TTg,fR_thr(:,con),'k');
xlim([0 gr]);
%legend([pd],'1*SD','Location','NorthWest');
%title('Smoothed envelope and + 1SD threshold')
subplot(7,1,5)
plot((TTg-delay/Nv),area_ch(:,con));
axis([0,gr,-1,2]);
figure
%title('Areas higher than threshold')
subplot(7,1,6)
plot(TTg,Setpoints(:,con));
axis([0,gr,-1,2]);
%title('Sleep spindles filtered (> 500 ms and < 3000ms)')
subplot(7,1,7)
plot(TTg,VectoS3);
axis([0,gr,-1,4]);
%title('Light signal onset')

%% Raster plot


Spindles = Setpoints';
xline = ones(ktk+1)*trim_min;
yline = 0:size(Spindles,1);

figure
hold all;

for gg = 1:size(Spindles,1)
    spikePos = find(Spindles(gg,:) == 1);
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [gg-0.2 gg+0.2], 'b');
    end
    line(xline,yline,'Color','green')
end

%% Cut spindles in each 12 sec bin. 12 before stim and 12 sec after.

cut = 12; % second to cut

timer = cut*1000;

n_spindles = Spindles(:,(trim_min-timer):(trim_min+timer)-1);

%% Sliced spindles raster plot


n_xline = ones(ktk+1)*timer;
n_yline = 0:size(n_spindles,1);

figure
hold all;

for gg = 1:size(n_spindles,1)
    spikePos = find(n_spindles(gg,:) == 1);
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [gg-0.2 gg+0.2], 'b');
    end
    line(n_xline,n_yline,'Color','green')
end
%%  Count spindles in each bin

Pre_sp = [];
Post_sp = [];

for vv = 1:size(n_spindles,1)
    Pre_sp(vv) = sum(n_spindles(vv,1:timer));
    Post_sp(vv) = sum(n_spindles(vv,(timer+1):end));
end
%% Calculate means (overall number/number of trials)

Pre_sum = sum(Pre_sp)
Post_sum = sum(Post_sp)

Pre_prob = Pre_sum/size(n_spindles,1)
Post_prob = Post_sum/size(n_spindles,1)



