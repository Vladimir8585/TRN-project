%% Plotting A_Barchanging, saving vectors and clearing from artificial noise
%

%% Plot 30_60 seconds sleep

rt30 = 90000;
Tvec = 1:length(VectoS3);
Sn = size(VectoS3,1);
figure
str = sprintf('30 second sleep ,N = %d',Sn);
title(str)
subplot(6,1,1)
for ls = 1:size(VectoS3,1)
    plot(Tvec,VectoS3(ls,:));
    hold on
end
title(str);
axis([0000, rt30, 0, 3.5])
subplot(6,1,2)
for lar = 1:size(Vecto_AR3,1)
    plot(Tvec,Vecto_AR3(lar,:));
    hold on
end
title('Auditory Right');
axis([0000, rt30, -500, 500])
subplot(6,1,3)
%for lal = 1:size(Vecto_AL3,1)
%    plot(Tvec,Vecto_AL3(lal,:));
%    hold on
%end
title('Auditory Left');
axis([0000, rt30, -700, 700])
subplot(6,1,4)
for lfl = 1:size(Vecto_FL3,1)
    plot(Tvec,Vecto_FL3(lfl,:));
    hold on
end
title('Frontal Left');
axis([0, rt30, -700, 700])
subplot(6,1,5)
for lfr = 1:size(Vecto_FR3,1)
    plot(Tvec,Vecto_FR3(lfr,:));
    hold on
end
title('Frontal Right');
axis([0, rt30, -700, 700])
subplot(6,1,6)
for lfr = 1:size(Vecto_EMG3,1)
    plot(Tvec,Vecto_EMG3(lfr,:));
    hold on
end
title('EMG');
axis([0, rt30, -700, 700])

%% Filtering 30_sleep

for k = size(Vecto_AR3,1):-1:1
    for l = size(Vecto_AR3,2):-1:2
    if (Vecto_AR3(k,l)) >= 700;
        %Vecto_AL3(k,:)=[];
        Vecto_FR3(k,:)=[];
        Vecto_FL3(k,:)=[];
        Vecto_AR3(k,:)=[];
        VectoS3(k,:)=[];
        Vecto_EMG3(k,:)=[];
        Stim_order3(k) = [];
    end
    end
end

for k = size(Vecto_AR3,1):-1:1
    for l = size(Vecto_AR3,2):-1:2
    if (Vecto_AR3(k,l)) <= -700;
        %Vecto_AL3(k,:)=[];
        Vecto_FR3(k,:)=[];
        Vecto_FL3(k,:)=[];
        Vecto_AR3(k,:)=[];
        VectoS3(k,:)=[];
        Vecto_EMG3(k,:)=[];
        Stim_order3(k) = [];
    end
    end
end

%size(Vecto_AL3)
size(Vecto_EMG3)
size(VectoS3)

%% Plot 30 sec awake

rt60 = 90000;
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
%for lal = 1:size(Vecto_AL3w,1)
%    plot(Tvec,Vecto_AL3w(lal,:));
%    hold on
%end
%title('Auditory Left');
%axis([0000, rt60, -700, 700])
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


%% Filtering from noise

%% Filtering 30_sleep

for k = size(Vecto_AR3w,1):-1:1
    for l = size(Vecto_AR3w,2):-1:2
    if (Vecto_AR3w(k,l)) >= 800;
       % Vecto_AL3w(k,:)=[];
        Vecto_FR3w(k,:)=[];
        Vecto_FL3w(k,:)=[];
        Vecto_AR3w(k,:)=[];
        VectoS3w(k,:)=[];
        Vecto_EMG3w(k,:)=[];
         Stim_order3w(k) = [];
    end
    end
end


%size(Vecto_AL3w)
size(Vecto_EMG3w)
size(VectoS3w)