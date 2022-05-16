%% Trimming code

ptt = size(VectoS3,1)
pt = size(VectoS3,2)

for s = 1:ptt
 for p = 2: pt
    if (VectoS3(s,p-1) <= 1) && (VectoS3(s,p)>=2);
        
        Vecto_AL3t(s,:) = Vecto_AL3(s,(p-10000):(p+9999));
        Vecto_AR3t(s,:) = Vecto_AR3(s,(p-10000):(p+9999));
        %Vecto_FL3t(s,:) = Vecto_FL3(s,(p-10000):(p+9999));
        Vecto_FR3t(s,:) = Vecto_FR3(s,(p-10000):(p+9999));
        Vecto_EMG3t(s,:) = Vecto_EMG3(s,(p-10000):(p+9999));
        VectoS3t(s,:) = VectoS3(s,(p-10000):(p+9999));
    end
 end
end

%%

Vecto_AL3 = Vecto_AL3t;
Vecto_AR3 = Vecto_AR3t;
%Vecto_FL3 = Vecto_FL3t;
Vecto_FR3 = Vecto_FR3t;
Vecto_EMG3 = Vecto_EMG3t;
VectoS3 = VectoS3t;
%% Plot 30 seconds sleep

rt30 = 25000;
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
for lal = 1:size(Vecto_AL3,1)
    plot(Tvec,Vecto_AL3(lal,:));
    hold on
end
% title('Auditory Left');
% axis([0000, rt30, -700, 700])
% subplot(6,1,4)
% for lfl = 1:size(Vecto_FL3,1)
%     plot(Tvec,Vecto_FL3(lfl,:));
%     hold on
% end
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