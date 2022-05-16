%%% Plotting A_Barchanging, saving vectors and clearing from artificial noise


%% 30 sec

rt60 = 100000;
Tvec = 1:length(VectoS3);
Snt = size(VectoS3,1);
figure
str = sprintf('30 second awake ,N = %d',Snt);
title(str)
subplot(6,1,1)
for ls = 1:size(VectoS3,1)
    plot(Tvec,VectoS3(ls,:));
    hold on
end
title(str);
axis([0000, rt60, 0, 3.5])
subplot(6,1,2)
for lar = 1:size(Vecto_AR3,1)
    plot(Tvec,Vecto_AR3(lar,:));
    hold on
end
title('Auditory Right');
axis([0000, rt60, -700, 700])
subplot(6,1,3)
 for lal = 1:size(Vecto_AL3,1)
    plot(Tvec,Vecto_AL3(lal,:));
    hold on
end
title('Auditory Left');
axis([0000, rt60, -700, 700])
subplot(6,1,4)
for lfl = 1:size(Vecto_FL3,1)
    plot(Tvec,Vecto_FL3(lfl,:));
    hold on
end
title('Frontal Left');
axis([0, rt60, -700, 700])
subplot(6,1,5)
for lfr = 1:size(Vecto_FR3,1)
    plot(Tvec,Vecto_FR3(lfr,:));
    hold on
end
title('Frontal Right');
axis([0, rt60, -700, 700])
subplot(6,1,6)
for lfr = 1:size(Vecto_EMG3,1)
    plot(Tvec,Vecto_EMG3(lfr,:));
    hold on
end
title('EMG');
axis([0, rt60, -700, 700])




%% Filtering Overall noise

for k = size(Vecto_FL3,1):-1:1
    for l = size(Vecto_FL3,2):-1:2
    if (Vecto_FL3(k,l)) >= 700;
        Vecto_AL3(k,:)=[];
        Vecto_FR3(k,:)=[];
        Vecto_FL3(k,:)=[];
        Vecto_AR3(k,:)=[];
        VectoS3(k,:)=[];
        Vecto_EMG3(k,:)=[];
        Stim_order3(k) = [];
        K3(k,:) = [];
    end
    end
end


size(Vecto_FL3)
size(Vecto_EMG3)
size(VectoS3)

%% Filtering 10 sec noise

ptt = size(Vecto_FL3,1);
pt = size(Vecto_FL3,2);
win = 9900;

%for s = 1:ptt
 for p = 2: pt
    if (VectoS3(1,p-1) <= 1) && (VectoS3(1,p)>=2);        
        pq = p-1;
    end
 end



for k = size(Vecto_FL3,1):-1:1
    for l = pq-win:pq+win-1;
    if (Vecto_FL3(k,l)) >= 800;
        Vecto_AL3(k,:)=[];
        Vecto_FR3(k,:)=[];
        Vecto_FL3(k,:)=[];
        Vecto_AR3(k,:)=[];
        VectoS3(k,:)=[];
        Vecto_EMG3(k,:)=[];
        Stim_order3(k) = [];
        K3(k,:) = [];
    end
    end
end


size(Vecto_FL3)
size(Vecto_EMG3)
size(VectoS3)

%% Filtering EMG noise

ptt = size(Vecto_FL3,1);
pt = size(Vecto_FL3,2);

%for s = 1:ptt
 for p = 2: pt
    if (VectoS3(1,p-1) <= 1) && (VectoS3(1,p)>=2);        
        pq = p-1;
    end
 end
%end


for k = size(Vecto_FL3,1):-1:1
    for l = 1:pq
    if (Vecto_EMG3(k,l)) >= 700;
        Vecto_AL3(k,:)=[];
        Vecto_FR3(k,:)=[];
        Vecto_FL3(k,:)=[];
        Vecto_AR3(k,:)=[];
        VectoS3(k,:)=[];
        Vecto_EMG3(k,:)=[];
        Stim_order3(k) = [];
        K3(k,:) = [];
    end
    end
end


size(Vecto_FL3)
size(Vecto_EMG3)
size(VectoS3)

%% Window filtering

ptt = size(Vecto_AL3,1);
pt = size(Vecto_AL3,2);

%for s = 1:ptt
 for p = 2: pt
    if (VectoS3(1,p-1) <= 1) && (VectoS3(1,p)>=2);        
        pq = p-1;
    end
 end
%end

win = 10000;

for k = size(Vecto_AL3,1):-1:1
    for l = pq-win:pq+win-1
    if (Vecto_AL3(k,l)) >= 700;
        Vecto_AL3(k,:)=[];
        Vecto_FR3(k,:)=[];
        Vecto_FL3(k,:)=[];
        Vecto_AR3(k,:)=[];
        VectoS3(k,:)=[];
        Vecto_EMG3(k,:)=[];
        Stim_order3(k) = [];
        K3(k,:) = [];
    end
    end
end


size(Vecto_AL3)
size(Vecto_EMG3)
size(VectoS3)


%%

ptt = size(Vecto_FL3,1);
pt = size(Vecto_FL3,2);

for s = 1:ptt
for p = 2: pt
    if (VectoS3(s,p-1) <= 1) && (VectoS3(s,p)>=2);        
        pq(s) = p-1;
    end
end
end

for s = ptt:-1:1
    if pq(s) ~= 20000
        Vecto_AL3(s,:)=[];
        Vecto_FR3(s,:)=[];
        Vecto_FL3(s,:)=[];
        Vecto_AR3(s,:)=[];
        VectoS3(s,:)=[];
        Vecto_EMG3(s,:)=[];
        Stim_order3(s) = [];
        K3(s,:) = [];
    end
end


%%
%% 30 w

rt60 = 100000;
Tvec = 1:length(VectoS3w);
Snt = size(VectoS3w,1);
figure
str = sprintf('30 second awake ,N = %d',Snt);
title(str)
subplot(6,1,1)
for ls = 1:size(VectoS3w,1)
    plot(Tvec,VectoS3w(ls,:));
    hold on
end
title(str);
axis([0000, rt60, 0, 3.5])
subplot(6,1,2)
for lar = 1:size(Vecto_AR3w,1)
    plot(Tvec,Vecto_AR3w(lar,:));
    hold on
end
title('Auditory Right');
axis([0000, rt60, -700, 700])
subplot(6,1,3)
 for lal = 1:size(Vecto_AL3w,1)
    plot(Tvec,Vecto_AL3w(lal,:));
    hold on
end
title('Auditory Left');
axis([0000, rt60, -700, 700])
subplot(6,1,4)
for lfl = 1:size(Vecto_FL3w,1)
    plot(Tvec,Vecto_FL3w(lfl,:));
    hold on
end
title('Frontal Left');
axis([0, rt60, -700, 700])
subplot(6,1,5)
for lfr = 1:size(Vecto_FR3w,1)
    plot(Tvec,Vecto_FR3w(lfr,:));
    hold on
end
title('Frontal Right');
axis([0, rt60, -700, 700])
subplot(6,1,6)
for lfr = 1:size(Vecto_EMG3w,1)
    plot(Tvec,Vecto_EMG3w(lfr,:));
    hold on
end
title('EMG');
axis([0, rt60, -700, 700])

%%
%%
%% 30n 

rt60 = 100000;
Tvec = 1:length(VectoS3n);
Snt = size(VectoS3n,1);
figure
str = sprintf('30 second awake ,N = %d',Snt);
title(str)
subplot(6,1,1)
for ls = 1:size(VectoS3n,1)
    plot(Tvec,VectoS3n(ls,:));
    hold on
end
title(str);
axis([0000, rt60, 0, 3.5])
subplot(6,1,2)
for lar = 1:size(Vecto_AR3n,1)
    plot(Tvec,Vecto_AR3n(lar,:));
    hold on
end
title('Auditory Right');
axis([0000, rt60, -700, 700])
subplot(6,1,3)
 for lal = 1:size(Vecto_AL3n,1)
    plot(Tvec,Vecto_AL3n(lal,:));
    hold on
end
title('Auditory Left');
axis([0000, rt60, -700, 700])
subplot(6,1,4)
for lfl = 1:size(Vecto_FL3n,1)
    plot(Tvec,Vecto_FL3n(lfl,:));
    hold on
end
title('Frontal Left');
axis([0, rt60, -700, 700])
subplot(6,1,5)
for lfr = 1:size(Vecto_FR3n,1)
    plot(Tvec,Vecto_FR3n(lfr,:));
    hold on
end
title('Frontal Right');
axis([0, rt60, -700, 700])
subplot(6,1,6)
for lfr = 1:size(Vecto_EMG3n,1)
    plot(Tvec,Vecto_EMG3n(lfr,:));
    hold on
end
title('EMG');
axis([0, rt60, -700, 700])

%%
ptt = size(Vecto_FL3n,1);
pt = size(Vecto_FL3n,2);
win = 9999;

%for s = 1:ptt
 for p = 2: pt
    if (VectoS3n(1,p-1) <= 1) && (VectoS3n(1,p)>=2);        
        pq = p-1;
    end
 end



for k = size(Vecto_FL3n,1):-1:1
    for l = pq-win:pq+win-1;
    if (Vecto_FL3n(k,l)) >= 800;
        Vecto_AL3n(k,:)=[];
        Vecto_FR3n(k,:)=[];
        Vecto_FL3n(k,:)=[];
        Vecto_AR3n(k,:)=[];
        VectoS3n(k,:)=[];
        Vecto_EMG3n(k,:)=[];
        Stim_order3n(k) = [];
        K3n(k,:) = [];
    end
    end
end


size(Vecto_AL3n)
size(Vecto_EMG3n)
size(VectoS3n)

