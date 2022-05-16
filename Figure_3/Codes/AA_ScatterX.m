
K5w1 = K5w;
K5w2 = K5w;
K5w3 = K5w;

K5w = [K5w1; K5w2];
K5w3 = [K5w1; K5w2];

%%

dsh_lx = ones(1,11)*21;
dsh_ly = ones(1,21)*11;

figure

scatter(K5w(:,2)/1000,K5w(:,1)/1000 ,5,'MarkerEdgeColor',[0.2 0.2 0.6],...
              'MarkerFaceColor',[0.2 0.2 0.6],...
              'LineWidth',0.5)
hold on
s = polyfit(K5w(:,2)/1000,K5w(:,1)/1000,1);
y1 = polyval(s,K5w(:,2)/1000);
plot(K5w(:,2)/1000,y1,'Color',[0.2 0.2 0.6])
[Z,K] = corrcoef(K5w(:,2)/1000,K5w(:,1)/1000 );
plot(dsh_lx,1:11,'k:');
plot(1:21,dsh_ly,'k:');

scatter(K5w3(:,2)/1000,K5w3(:,1)/1000 ,5,'MarkerEdgeColor',[0.8 .3 .6],...
              'MarkerFaceColor',[0.8 .3 .6],...
              'LineWidth',0.5)

z = polyfit(K5w3(:,2)/1000,K5w3(:,1)/1000,1);
y2 = polyval(z,K5w3(:,2)/1000);
plot(K5w3(:,2)/1000,y2,'Color',[0.8 .3 .6])
[Zi,Ki] = corrcoef(K5w3(:,2)/1000,K5w3(:,1)/1000 );


xlim([10 40])
ylim([0 80])

%%

Dot(1,1) = 0;
Dot(2,1) = s(1);
Dot(3,1) = s(2);
Dot(4,1) = Z(1,2);
Dot(5,1) = K(1,2);
Dot(6,1) = 0;
Dot(7,1) = z(1);
Dot(8,1) = z(2);
Dot(9,1) = Zi(1,2);
Dot(10,1) = Ki(1,2);
Dot


%%
[0.8 .3 .6] - magenda
[0.2 0.2 0.7] - blue
[0 1 0] - green
[0.3 0.4 0.7] - blue
%%

dsh_lx = ones(1,11)*21;
dsh_ly = ones(1,21)*11;

figure

scatter(K5w(:,2)/1000,K5w(:,1)/1000 ,5,'MarkerEdgeColor',[0.3 0.4 0.7],...
              'MarkerFaceColor',[0.3 0.4 0.7],...
              'LineWidth',0.5)
hold on
s = polyfit(K5w(:,2)/1000,K5w(:,1)/1000,1);
y1 = polyval(s,K5w(:,2)/1000);
plot(K5w(:,2)/1000,y1,'Color',[0.3 0.4 0.7])
[Z,K] = corrcoef(K5w(:,2)/1000,K5w(:,1)/1000 );
plot(dsh_lx,1:11,'k:');
plot(1:21,dsh_ly,'k:');

scatter(K5w3(:,2)/1000,K5w3(:,1)/1000 ,5,'MarkerEdgeColor',[0.1 0.9 0.4],...
              'MarkerFaceColor',[0.1 0.9 0.4],...
              'LineWidth',0.5)
hold on
z = polyfit(K5w3(:,2)/1000,K5w3(:,1)/1000,1);
y2 = polyval(z,K5w3(:,2)/1000);
plot(K5w3(:,2)/1000,y2,'Color',[0.1 0.9 0.4])
[Zi,Ki] = corrcoef(K5w3(:,2)/1000,K5w3(:,1)/1000 );


xlim([10 40])
ylim([0 80])


Dot(1,1) = 0;
Dot(2,1) = s(1);
Dot(3,1) = s(2);
Dot(4,1) = Z(1,2);
Dot(5,1) = K(1,2);
Dot(6,1) = 0;
Dot(7,1) = z(1);
Dot(8,1) = z(2);
Dot(9,1) = Zi(1,2);
Dot(10,1) = Ki(1,2);
Dot



