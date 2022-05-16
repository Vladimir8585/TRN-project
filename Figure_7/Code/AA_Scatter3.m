


A = ismember(K30f(:,2),K5w(:,3)); %% K30/Vecto
R = find(A == 1);

B = ismember(K5w(:,3),K30f(:,2)); %%K5w (
D = find(B == 1);
Prun = K5w(D,:);
PDrun = K30f(R,:);

C = ismember(Prun(:,1), PDrun(:,1))
K = find(C ==1);
k30 = K5w(C,:);
Vecto_FR3 = Vecto_FR3f(C,:);
VectoS3 = VectoS3f(C,:);
%%

winp = 10;

Del_FR =  D_Delta_FR(winp,:)';
The_FR = D_Theta_FR(winp,:)';
Alp_FR = D_Alph_FR(winp,:)';
Bet_FR = D_Beta_FR(winp,:)';
Gam_FR = D_Gamma_FR(winp,:)';
GamX_FR = D_GammaX_FR(winp,:)';

figure
subplot(2,2,3)
scatter(K30f(:,2),Alp_FR,5,'MarkerEdgeColor',[0.2 0.2 0.7],...
              'MarkerFaceColor',[0.2 0.2 0.7],...
              'LineWidth',0.5)
hold on
s = polyfit(K30f(:,2),Alp_FR,1)
y1 = polyval(s,K30f(:,2));
plot(K30f(:,2),y1,'Color',[0.2 0.2 0.7])
[Z,K] = corrcoef(K30f(:,2),Alp_FR)
%xlim([0 80])
%title('Alpha')

Dot(10,1) = 0;
Dot(11,1) = s(1);
Dot(12,1) = s(2);
Dot(13,1) = Z(1,2);
Dot(14,1) = K(1,2);

%%
winp = 10;

Del_FR =  D_Delta_FR(winp,:)';
The_FR = D_Theta_FR(winp,:)';
Alp_FR = D_Alph_FR(winp,:)';
Bet_FR = D_Beta_FR(winp,:)';
Gam_FR = D_Gamma_FR(winp,:)';
GamX_FR = D_GammaX_FR(winp,:)';

figure
subplot(2,2,3)
scatter(K30f(:,2),Alp_FR,5,'MarkerEdgeColor',[0.8 .3 .6],...
              'MarkerFaceColor',[0.8 .3 .6],...
              'LineWidth',0.5)
hold on
s = polyfit(K30f(:,2),Alp_FR,1)
y1 = polyval(s,K30f(:,2));
plot(K30f(:,1)/1000,y1,'Color',[0.8 .3 .6])
[Z,K] = corrcoef(K30f(:,2),Alp_FR)
%xlim([0 80])
%title('Alpha')

Dot(10,1) = 0;
Dot(11,1) = s(1);
Dot(12,1) = s(2);
Dot(13,1) = Z(1,2);
Dot(14,1) = K(1,2);
%%
figure
subplot(2,2,3)
scatter(K30f(:,2),K30f(:,1)/1000,5,'MarkerEdgeColor',[0.8 .3 .6],...
              'MarkerFaceColor',[0.8 .3 .6],...
              'LineWidth',0.5)
hold on
s = polyfit(K30f(:,2),K30f(:,1)/1000,1)
y1 = polyval(s,K30f(:,2));
plot(K30f(:,2),y1,'Color',[0.8 .3 .6])
[Z,K] = corrcoef(K30f(:,2)/1000,K30f(:,1)/1000)
xlim([0 50])
ylim([0 80])
%title('Alpha')

Dot(10,1) = 0;
Dot(11,1) = s(1);
Dot(12,1) = s(2);
Dot(13,1) = Z(1,2);
Dot(14,1) = K(1,2);

%%

%%
figure
subplot(2,2,3)
scatter(K30f(:,2),K30f(:,1)/1000,5,'MarkerEdgeColor',[0.2 0.2 0.7],...
              'MarkerFaceColor',[0.2 0.2 0.7],...
              'LineWidth',0.5)
hold on
s = polyfit(K30f(:,2),K30f(:,1)/1000,1)
y1 = polyval(s,K30f(:,2));
plot(K30f(:,2),y1,'Color',[0.2 0.2 0.7])
[Z,K] = corrcoef(K30f(:,2)/1000,K30f(:,1)/1000)
xlim([0 50])
ylim([0 80])
%title('Alpha')

Dot(10,1) = 0;
Dot(11,1) = s(1);
Dot(12,1) = s(2);
Dot(13,1) = Z(1,2);
Dot(14,1) = K(1,2);

%%
figure
subplot(2,2,3)
scatter(K30f(:,2),K30f(:,1)/1000,5,'MarkerEdgeColor',[0.2 0.2 0.7],...
              'MarkerFaceColor',[0.2 0.2 0.7],...
              'LineWidth',0.5)
hold on
s = polyfit(K30f(:,2),K30f(:,1)/1000,1)
y1 = polyval(s,K30f(:,2));
plot(K30f(:,2),y1,'Color',[0.2 0.2 0.7])
[Z,K] = corrcoef(K30f(:,2)/1000,K30f(:,1)/1000)
xlim([0 50])
ylim([0 80])
%title('Alpha')

Dot(10,1) = 0;
Dot(11,1) = s(1);
Dot(12,1) = s(2);
Dot(13,1) = Z(1,2);
Dot(14,1) = K(1,2);


%%
[0.8 .3 .6] - magenda
[0.2 0.2 0.7] - blue
[0 1 0] - green
[0.3 0.4 0.7] - blue
%%
%%
winp = 10;

Del_FR =  D_Delta_FR(winp,:)';
The_FR = D_Theta_FR(winp,:)';
Alp_FR = D_Alph_FR(winp,:)';
Bet_FR = D_Beta_FR(winp,:)';
Gam_FR = D_Gamma_FR(winp,:)';
GamX_FR = D_GammaX_FR(winp,:)';

figure
subplot(2,2,3)
scatter(Stim_order3,Bet_FR,5,'MarkerEdgeColor',[0.2 0.2 0.7],...
              'MarkerFaceColor',[0.2 0.2 0.7],...
              'LineWidth',0.5)
hold on
s = polyfit(Stim_order3,Bet_FR,1)
y1 = polyval(s,Stim_order3);
plot(Stim_order3,y1,'Color',[0.2 0.2 0.7])
[Z,K] = corrcoef(Stim_order3,Bet_FR)
xlim([0 50])
ylim([-5 5])
%title('Alpha')

Dot(10,1) = 0;
Dot(11,1) = s(1);
Dot(12,1) = s(2);
Dot(13,1) = Z(1,2);
Dot(14,1) = K(1,2);




%%
winp = 10;

Del_FR =  D_Delta_FR(winp,:)';
The_FR = D_Theta_FR(winp,:)';
Alp_FR = D_Alph_FR(winp,:)';
Bet_FR = D_Beta_FR(winp,:)';
Gam_FR = D_Gamma_FR(winp,:)';
GamX_FR = D_GammaX_FR(winp,:)';

figure
subplot(2,2,3)
scatter(Stim_order3,Del_FR,5,'MarkerEdgeColor',[0.2 0.2 0.7],...
              'MarkerFaceColor',[0.2 0.2 0.7],...
              'LineWidth',0.5)
hold on
s = polyfit(Stim_order3,Del_FR,1)
y1 = polyval(s,Stim_order3);
plot(Stim_order3,y1,'Color',[0.2 0.2 0.7])
[Z,K] = corrcoef(Stim_order3,Del_FR)
xlim([0 50])
ylim([-5 5])
%title('Alpha')

Dot(10,1) = 0;
Dot(11,1) = s(1);
Dot(12,1) = s(2);
Dot(13,1) = Z(1,2);
Dot(14,1) = K(1,2);




%%
winp = 10;

Del_FR =  D_Delta_FR(winp,:)';
The_FR = D_Theta_FR(winp,:)';
Alp_FR = D_Alph_FR(winp,:)';
Bet_FR = D_Beta_FR(winp,:)';
Gam_FR = D_Gamma_FR(winp,:)';
GamX_FR = D_GammaX_FR(winp,:)';

%
figure
subplot(3,2,1)
scatter(K30f(:,1), Del_FR,'MarkerEdgeColor',[0 0.1 0.1],...
              'MarkerFaceColor',[0.3 .7 .7],...
              'LineWidth',1.5)
hold on
s = polyfit(K30f(:,1),Del_FR,1);
y1 = polyval(s,K30f(:,1)/1000);
plot(K30f(:,1),y1)
[Z,K] = corrcoef(K30f(:,1),Del_FR)
%xlim([0 80])
title('Delta')

Dot(1,1) = s(1);
Dot(2,1) = s(2);
Dot(3,1) = Z(1,2);
Dot(4,1) = K(1,2);
%

%
subplot(3,2,2)
scatter(K30f(:,1),The_FR ,'MarkerEdgeColor',[0 0.1 0.1],...
              'MarkerFaceColor',[0.3 .7 .7],...
              'LineWidth',1.5)
hold on
s = polyfit(K30f(:,1),The_FR,1)
y1 = polyval(s,K30f(:,1));
plot(K30f(:,1),y1)
[Z,K] = corrcoef(K30f(:,1),The_FR );
%xlim([0 80])
title('theta')

Dot(5,1) = 0;
Dot(6,1) = s(1);
Dot(7,1) = s(2);
Dot(8,1) = Z(1,2);
Dot(9,1) = K(1,2);

%

subplot(3,2,3)
scatter(K30f(:,1),Alp_FR,'MarkerEdgeColor',[0 0.1 0.1],...
              'MarkerFaceColor',[0.3 .7 .7],...
              'LineWidth',1)
hold on
s = polyfit(K30f(:,1),Alp_FR,1)
y1 = polyval(s,K30f(:,1));
plot(K30f(:,1),y1)
[Z,K] = corrcoef(K30f(:,1),Alp_FR)
%xlim([0 80])
title('Alpha')

Dot(10,1) = 0;
Dot(11,1) = s(1);
Dot(12,1) = s(2);
Dot(13,1) = Z(1,2);
Dot(14,1) = K(1,2);

%

subplot(3,2,4)
scatter(K30f(:,1),Bet_FR,'MarkerEdgeColor',[0 0.1 0.1],...
              'MarkerFaceColor',[0.3 .7 .7],...
              'LineWidth',1.5)
hold on
s = polyfit(K30f(:,1),Bet_FR,1);
y1 = polyval(s,K30f(:,1));
plot(K30f(:,1),y1)
[Z,K] = corrcoef(K30f(:,1),Bet_FR)
%xlim([0 80])
title('Beta')

Dot(15,1) = 0;
Dot(16,1) = s(1);
Dot(17,1) = s(2);
Dot(18,1) = Z(1,2);
Dot(19,1) = K(1,2);

%

subplot(3,2,5)
scatter(K30f(:,1),Gam_FR,'MarkerEdgeColor',[0 0.1 0.1],...
              'MarkerFaceColor',[0.3 .7 .7],...
              'LineWidth',1.5)
hold on
s = polyfit(K30f(:,1),Gam_FR,1);
y1 = polyval(s,K30f(:,1));
plot(K30f(:,1),y1)
[Z,K] = corrcoef(K30f(:,1),Gam_FR);
%xlim([0 80])
title('Gamma')

%xlim([0 60])
Dot(20,1) = 0;
Dot(21,1) = s(1);
Dot(22,1) = s(2);
Dot(23,1) = Z(1,2);
Dot(24,1) = K(1,2);

%

subplot(3,2,6)
scatter(K30f(:,1),GamX_FR,'MarkerEdgeColor',[0 0.1 0.1],...
              'MarkerFaceColor',[0.3 .7 .7],...
              'LineWidth',1.5)
hold on
s = polyfit(K30f(:,1),GamX_FR,1)
y1 = polyval(s,K30f(:,1));
plot(K30f(:,1),y1)
% xlim([0 30])
[Z,K] = corrcoef(K30f(:,1),GamX_FR)
% xlim([0 80])
title('GammaX')

Dot(25,1) = 0;
Dot(26,1) = s(1);
Dot(27,1) = s(2);
Dot(28,1) = Z(1,2);
Dot(29,1) = K(1,2);

Dot
