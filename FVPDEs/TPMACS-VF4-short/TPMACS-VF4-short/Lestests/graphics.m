load('data_Pb33.txt');
load('data_Pb34.txt');
load('data_Pb35.txt');
load('data_Pb36.txt');
load('data_Pb37.txt');
load('data_Pb38.txt');
load('data_Pb39.txt');
load('data_Pb41.txt');
load('data_Pb42.txt');
load('data_Pb43.txt');

linewidth = 3;
figure('Position', get(0, 'Screensize'))%('Position', [10,10,1300,600]);
%h = suptitle('Tumor mass through and after treatment');
%set(h,'FontSize',20,'FontWeight','normal')

ax1=subplot(2,2,1);
hold on;
d33=plot(data_Pb33(:,1),data_Pb33(:,2),'Color','black','LineWidth',linewidth);
d34=plot(data_Pb34(:,1),data_Pb34(:,2),'--','Color','black','LineWidth',linewidth);
d35=plot(data_Pb35(:,1),data_Pb35(:,2),'--','Color',[0.3010 0.7450 0.9330],'LineWidth',linewidth);
d36=plot(data_Pb36(:,1),data_Pb36(:,2),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',linewidth);
xlabel('Time (days)');
ylabel('Tumor proportion');
ylim([1e-6,1e-1]);
title('Surgery only');
legend({'No treatment','Surgery only - threshold 0.03','Surgery only - threshold 0.02','Surgery only - threshold 0.01'},'Location','southeast');

ax2=subplot(2,2,2);
hold on;
d34=plot(data_Pb34(:,1),data_Pb34(:,2),'--','Color','black','LineWidth',linewidth);
d37=plot(data_Pb37(:,1),data_Pb37(:,2),':','Color','red','LineWidth',linewidth);
d38=plot(data_Pb38(:,1),data_Pb38(:,2),':','Color',[0 0.4470 0.7410],'LineWidth',linewidth);
line([14,14],[1e-9,1e-1],'Color',[0.4660 0.6740 0.1880],'LineStyle','--','LineWidth',1);line([56,56],[1e-9,1e-1],'Color',[0.4660 0.6740 0.1880],'LineStyle','--','LineWidth',1);
xlabel('Time (days)');
ylabel('Tumor proportion');
ylim([1e-9,1e-1]);
title('Surgery with Chemotherapy');
legend({'Surgery only','S + chemotherapy : dose 0.2','S + chemotherapy : dose 0.3'},'Location','best')

ax3=subplot(2,2,3);
hold on;
d34=plot(data_Pb34(:,1),data_Pb34(:,2),'--','Color','black','LineWidth',linewidth);
d39=plot(data_Pb39(:,1),data_Pb39(:,2),':','Color','red','LineWidth',linewidth);
d41=plot(data_Pb41(:,1),data_Pb41(:,2),':','Color','green','LineWidth',linewidth);
line([14,14],[1e-6,2e-5],'Color','blue','LineStyle','--','LineWidth',1);line([54,54],[1e-6,2e-5],'Color','blue','LineStyle','--','LineWidth',1);
line([19,19],[1e-6,2e-5],'Color','blue','LineStyle','--','LineWidth',1);line([21,21],[1e-6,2e-5],'Color','blue','LineStyle','--','LineWidth',1);line([26,26],[1e-6,2e-5],'Color','blue','LineStyle','--','LineWidth',1);line([28,28],[1e-6,2e-5],'Color','blue','LineStyle','--','LineWidth',1);line([33,33],[1e-6,2e-5],'Color','blue','LineStyle','--','LineWidth',1);
line([35,35],[1e-6,2e-5],'Color','blue','LineStyle','--','LineWidth',1);line([40,40],[1e-6,2e-5],'Color','blue','LineStyle','--','LineWidth',1);line([42,42],[1e-6,2e-5],'Color','blue','LineStyle','--','LineWidth',1);line([47,47],[1e-6,2e-5],'Color','blue','LineStyle','--','LineWidth',1);line([49,49],[1e-6,2e-5],'Color','blue','LineStyle','--','LineWidth',1);
line([14,19],[2e-5,2e-5],'Color','blue','LineStyle','--','LineWidth',1);line([21,26],[2e-5,2e-5],'Color','blue','LineStyle','--','LineWidth',1);line([28,33],[2e-5,2e-5],'Color','blue','LineStyle','--','LineWidth',1);line([35,40],[2e-5,2e-5],'Color','blue','LineStyle','--','LineWidth',1);line([42,47],[2e-5,2e-5],'Color','blue','LineStyle','--','LineWidth',1);
line([49,54],[2e-5,2e-5],'Color','blue','LineStyle','--','LineWidth',1);
xlabel('Time (days)');
ylabel('Tumor proportion');
ylim([1e-5,1e-1]);
title('Surgery with radiotherapy');
legend({'Surgery only','S + radiotherapy : 30*2Gy/day','S + radiotherapy : 30*1.8Gy/day'},'Location','northeast');


ax4=subplot(2,2,4);
hold on;
d33=plot(data_Pb33(:,1),data_Pb33(:,2),'Color','black','LineWidth',linewidth);
d42=plot(data_Pb42(:,1),data_Pb42(:,2),':','Color','red','LineWidth',linewidth);
d43=plot(data_Pb43(:,1),data_Pb43(:,2),'-.','Color','red','LineWidth',linewidth);
line([14,14],[1e-9,1e-1],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);line([56,56],[1e-9,1e-1],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);
line([54,54],[1e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);
line([19,19],[1e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);line([21,21],[1e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);line([26,26],[1e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);line([28,28],[1e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);line([33,33],[1e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);
line([35,35],[1e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);line([40,40],[1e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);line([42,42],[1e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);line([47,47],[1e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);line([49,49],[1e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);
line([14,19],[2e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);line([21,26],[2e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);line([28,33],[2e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);line([35,40],[2e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);line([42,47],[2e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);
line([49,54],[2e-9,2e-9],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);
xlabel('Time (days)');
ylabel('Tumor proportion');
ylim([1e-9,1e-1]);
legend({'No treatment','Radiotherapy + chemotherapy','Surgery + radiotherapy + chemotherapy'},'Location','best');
title('Chemotherapy and radiotherapy with and without surgery');

ax1.YScale = 'log';ax2.YScale = 'log';ax3.YScale = 'log';ax4.YScale = 'log';

x = [0.62 0.755];    % adjust length and location of arrow 
y = [0.6 0.6];      % adjust hieght and width of arrow
annotation('textarrow',x,y,'Color',[0.4660 0.6740 0.1880],'FontSize',9,'Linewidth',2)
annotation('textbox',[.62 .335 .67 .305],'Color',[0.4660 0.6740 0.1880],'EdgeColor','none','String','Chemotherapy','FontSize',11,'Linewidth',2)

annotation('textarrow',x-[0.44 0.45],y-0.45,'Color','blue','FontSize',9,'Linewidth',2)
annotation('textbox',[.18 .108 .23 .078],'Color','blue','EdgeColor','none','String','Radiotherapy','FontSize',11,'Linewidth',2)

annotation('textarrow',x,y-0.46,'Color',[0.4940 0.1840 0.5560],'FontSize',9,'Linewidth',2)
annotation('textbox',[.618 .105 .668 .075],'Color',[0.4940 0.1840 0.5560],'EdgeColor','none','String','Chemotherapy + Radiotherapy','FontSize',11,'Linewidth',2)

saveas(gcf,'data.png')

%% 
load('data_Pb46.txt');
figure;
hold on;
plot(data_Pb46(:,1),data_Pb46(:,4));
plot(data_Pb46(:,1),data_Pb46(:,8));
plot(data_Pb46(:,1),data_Pb46(:,12));


%% Les données sont dans l'ordre : 
%%temps+dt, masse tumorale totale,dt,Tum(P1),Nut(P1),Endo(P1),VEGF(P1),Tum(P2),...,VEGF(P3)
load('data_Pb50.txt');
load('data_Pb51.txt');
load('data_Pb52.txt');
load('data_Pb53.txt');
load('data_Pb54.txt');
load('data_Pb55.txt');
linewidth = 3;
figure('Position', get(0, 'Screensize'));
hold on;
xlim([0,100]);
plot(data_Pb50(:,1),data_Pb50(:,4),'Color','black','LineWidth',linewidth);
plot(data_Pb50(:,1),data_Pb50(:,8),':','Color','black','LineWidth',linewidth);
plot(data_Pb50(:,1),data_Pb50(:,12),'--','Color','black','LineWidth',linewidth);

plot(data_Pb51(:,1),data_Pb51(:,4),'Color','red','LineWidth',linewidth);
plot(data_Pb51(:,1),data_Pb51(:,8),':','Color','red','LineWidth',linewidth);
plot(data_Pb51(:,1),data_Pb51(:,12),'--','Color','red','LineWidth',linewidth);

plot(data_Pb52(:,1),data_Pb52(:,4),'Color','magenta','LineWidth',linewidth);
plot(data_Pb52(:,1),data_Pb52(:,8),':','Color','magenta','LineWidth',linewidth);
plot(data_Pb52(:,1),data_Pb52(:,12),'--','Color','magenta','LineWidth',linewidth);

plot(data_Pb53(:,1),data_Pb53(:,4),'Color','blue','LineWidth',linewidth);
plot(data_Pb53(:,1),data_Pb53(:,8),':','Color','blue','LineWidth',linewidth);
plot(data_Pb53(:,1),data_Pb53(:,12),'--','Color','blue','LineWidth',linewidth);

plot(data_Pb54(:,1),data_Pb54(:,4),'Color','green','LineWidth',linewidth);
plot(data_Pb54(:,1),data_Pb54(:,8),':','Color','green','LineWidth',linewidth);
plot(data_Pb54(:,1),data_Pb54(:,12),'--','Color','green','LineWidth',linewidth);

plot(data_Pb55(:,1),data_Pb55(:,4),'Color','yellow','LineWidth',linewidth);
plot(data_Pb55(:,1),data_Pb55(:,8),':','Color','yellow','LineWidth',linewidth);
plot(data_Pb55(:,1),data_Pb55(:,12),'--','Color','yellow','LineWidth',linewidth);
%% 
load('data_Pb60.txt');
load('data_Pb61.txt');
load('data_Pb62.txt');
load('data_Pb63.txt');
load('data_Pb64.txt');
load('data_Pb65.txt');
umax = 2.39e8;
figure('Position', get(0, 'Screensize'));
ax1=subplot(2,2,1);
hold on;
xlim([0,100]);
ylim(umax*[1e-6,1]);
plot(data_Pb60(:,1),umax*data_Pb60(:,4),'Color','black','LineWidth',linewidth);
plot(data_Pb60(:,1),umax*data_Pb60(:,8),':','Color','black','LineWidth',linewidth);
plot(data_Pb60(:,1),umax*data_Pb60(:,12),'--','Color','black','LineWidth',linewidth);

plot(data_Pb61(:,1),umax*data_Pb61(:,4),'Color','red','LineWidth',linewidth);
plot(data_Pb61(:,1),umax*data_Pb61(:,8),':','Color','red','LineWidth',linewidth);
plot(data_Pb61(:,1),umax*data_Pb61(:,12),'--','Color','red','LineWidth',linewidth);

plot(data_Pb62(:,1),umax*data_Pb62(:,4),'Color','magenta','LineWidth',linewidth);
plot(data_Pb62(:,1),umax*data_Pb62(:,8),':','Color','magenta','LineWidth',linewidth);
plot(data_Pb62(:,1),umax*data_Pb62(:,12),'--','Color','magenta','LineWidth',linewidth);

plot(data_Pb63(:,1),umax*data_Pb63(:,4),'Color','blue','LineWidth',linewidth);
plot(data_Pb63(:,1),umax*data_Pb63(:,8),':','Color','blue','LineWidth',linewidth);
plot(data_Pb63(:,1),umax*data_Pb63(:,12),'--','Color','blue','LineWidth',linewidth);

plot(data_Pb64(:,1),umax*data_Pb64(:,4),'Color','green','LineWidth',linewidth);
plot(data_Pb64(:,1),umax*data_Pb64(:,8),':','Color','green','LineWidth',linewidth);
plot(data_Pb64(:,1),umax*data_Pb64(:,12),'--','Color','green','LineWidth',linewidth);

plot(data_Pb65(:,1),umax*data_Pb65(:,4),'Color','yellow','LineWidth',linewidth);
plot(data_Pb65(:,1),umax*data_Pb65(:,8),':','Color','yellow','LineWidth',linewidth);
plot(data_Pb65(:,1),umax*data_Pb65(:,12),'--','Color','yellow','LineWidth',linewidth);
xlabel('Time (days)');
ylabel('Cancer cells per cm^2');
%legend({'No treatment - P1','No treatment - P2','No treatment - P3','Surgery - P1','Surgery - P2','Surgery - P3','Chemotherapy - P1','Chemotherapy - P2','Chemotherapy - P3','Radiotherapy-P1','Radiotherapy-P2','Radiotherapy-P3','Chemotherapy+Radiotherapy-P1','Chemotherapy+Radiotherapy-P2','Chemotherapy+Radiotherapy-P3','Surgery+Chemotherapy+Radiotherapy-P1','Surgery+Chemotherapy+Radiotherapy-P2','Surgery+Chemotherapy+Radiotherapy-P3'},'Location','best');
ax1.YScale = 'log';

ax2=subplot(2,2,2);
hold on;
xlim([0,100]);
ylim([2e-2,2e2]);
plot(data_Pb60(:,1),data_Pb60(:,5),'Color','black','LineWidth',linewidth);
plot(data_Pb60(:,1),data_Pb60(:,9),':','Color','black','LineWidth',linewidth);
plot(data_Pb60(:,1),data_Pb60(:,13),'--','Color','black','LineWidth',linewidth);

plot(data_Pb61(:,1),data_Pb61(:,5),'Color','red','LineWidth',linewidth);
plot(data_Pb61(:,1),data_Pb61(:,9),':','Color','red','LineWidth',linewidth);
plot(data_Pb61(:,1),data_Pb61(:,13),'--','Color','red','LineWidth',linewidth);

plot(data_Pb62(:,1),data_Pb62(:,5),'Color','magenta','LineWidth',linewidth);
plot(data_Pb62(:,1),data_Pb62(:,9),':','Color','magenta','LineWidth',linewidth);
plot(data_Pb62(:,1),data_Pb62(:,13),'--','Color','magenta','LineWidth',linewidth);

plot(data_Pb63(:,1),data_Pb63(:,5),'Color','blue','LineWidth',linewidth);
plot(data_Pb63(:,1),data_Pb63(:,9),':','Color','blue','LineWidth',linewidth);
plot(data_Pb63(:,1),data_Pb63(:,13),'--','Color','blue','LineWidth',linewidth);

plot(data_Pb64(:,1),data_Pb64(:,5),'Color','green','LineWidth',linewidth);
plot(data_Pb64(:,1),data_Pb64(:,9),':','Color','green','LineWidth',linewidth);
plot(data_Pb64(:,1),data_Pb64(:,13),'--','Color','green','LineWidth',linewidth);

plot(data_Pb65(:,1),data_Pb65(:,5),'Color','yellow','LineWidth',linewidth);
plot(data_Pb65(:,1),data_Pb65(:,9),':','Color','yellow','LineWidth',linewidth);
plot(data_Pb65(:,1),data_Pb65(:,13),'--','Color','yellow','LineWidth',linewidth);
xlabel('Time (days)');
ylabel('Nutriment (mol.1e-6)/cm^2');
ax2.YScale = 'log';

ax3=subplot(2,2,3);
hold on;
xlim([0,100]);
ylim([1e5,1e7]);
plot(data_Pb60(:,1),umax*data_Pb60(:,6),'Color','black','LineWidth',linewidth);
plot(data_Pb60(:,1),umax*data_Pb60(:,10),':','Color','black','LineWidth',linewidth);
plot(data_Pb60(:,1),umax*data_Pb60(:,14),'--','Color','black','LineWidth',linewidth);

plot(data_Pb61(:,1),umax*data_Pb61(:,6),'Color','red','LineWidth',linewidth);
plot(data_Pb61(:,1),umax*data_Pb61(:,10),':','Color','red','LineWidth',linewidth);
plot(data_Pb61(:,1),umax*data_Pb61(:,14),'--','Color','red','LineWidth',linewidth);

plot(data_Pb62(:,1),umax*data_Pb62(:,6),'Color','magenta','LineWidth',linewidth);
plot(data_Pb62(:,1),umax*data_Pb62(:,10),':','Color','magenta','LineWidth',linewidth);
plot(data_Pb62(:,1),umax*data_Pb62(:,14),'--','Color','magenta','LineWidth',linewidth);

plot(data_Pb63(:,1),umax*data_Pb63(:,6),'Color','blue','LineWidth',linewidth);
plot(data_Pb63(:,1),umax*data_Pb63(:,10),':','Color','blue','LineWidth',linewidth);
plot(data_Pb63(:,1),umax*data_Pb63(:,14),'--','Color','blue','LineWidth',linewidth);

plot(data_Pb64(:,1),umax*data_Pb64(:,6),'Color','green','LineWidth',linewidth);
plot(data_Pb64(:,1),umax*data_Pb64(:,10),':','Color','green','LineWidth',linewidth);
plot(data_Pb64(:,1),umax*data_Pb64(:,14),'--','Color','green','LineWidth',linewidth);

plot(data_Pb65(:,1),umax*data_Pb65(:,6),'Color','yellow','LineWidth',linewidth);
plot(data_Pb65(:,1),umax*data_Pb65(:,10),':','Color','yellow','LineWidth',linewidth);
plot(data_Pb65(:,1),umax*data_Pb65(:,14),'--','Color','yellow','LineWidth',linewidth);
xlabel('Time (days)');
ylabel('Endothelial cells per cm^2');
ax3.YScale = 'log';

ax4=subplot(2,2,4);
hold on;
xlim([0,100]);
ylim([1e-5,1e1]);
plot(data_Pb60(:,1),data_Pb60(:,7),'Color','black','LineWidth',linewidth);
plot(data_Pb60(:,1),data_Pb60(:,11),':','Color','black','LineWidth',linewidth);
plot(data_Pb60(:,1),data_Pb60(:,15),'--','Color','black','LineWidth',linewidth);

plot(data_Pb61(:,1),data_Pb61(:,7),'Color','red','LineWidth',linewidth);
plot(data_Pb61(:,1),data_Pb61(:,11),':','Color','red','LineWidth',linewidth);
plot(data_Pb61(:,1),data_Pb61(:,15),'--','Color','red','LineWidth',linewidth);

plot(data_Pb62(:,1),data_Pb62(:,7),'Color','magenta','LineWidth',linewidth);
plot(data_Pb62(:,1),data_Pb62(:,11),':','Color','magenta','LineWidth',linewidth);
plot(data_Pb62(:,1),data_Pb62(:,15),'--','Color','magenta','LineWidth',linewidth);

plot(data_Pb63(:,1),data_Pb63(:,7),'Color','blue','LineWidth',linewidth);
plot(data_Pb63(:,1),data_Pb63(:,11),':','Color','blue','LineWidth',linewidth);
plot(data_Pb63(:,1),data_Pb63(:,15),'--','Color','blue','LineWidth',linewidth);

plot(data_Pb64(:,1),data_Pb64(:,7),'Color','green','LineWidth',linewidth);
plot(data_Pb64(:,1),data_Pb64(:,11),':','Color','green','LineWidth',linewidth);
plot(data_Pb64(:,1),data_Pb64(:,15),'--','Color','green','LineWidth',linewidth);

plot(data_Pb65(:,1),data_Pb65(:,7),'Color','yellow','LineWidth',linewidth);
plot(data_Pb65(:,1),data_Pb65(:,11),':','Color','yellow','LineWidth',linewidth);
plot(data_Pb65(:,1),data_Pb65(:,15),'--','Color','yellow','LineWidth',linewidth);
xlabel('Time (days)');
ylabel('VEGF concentration (mol.1e-6)/cm^2');
ax4.YScale = 'log';

%% Couleurs

color1 = [0, 0.4470, 0.7410]; %No treatment
color2 = [0.8500, 0.3250, 0.0980];%Surgery
color3 = [0.9290, 0.6940, 0.1250];%Chemotherapy
color4 = [0.4940, 0.1840, 0.5560];%Radiotherapy
color5 = [0.4660, 0.6740, 0.1880];%Chemotherapy + Radiotherapy
color6 = [0.3010, 0.7450, 0.9330];%Surgery + Chemotherapy + Radiotherapy
color7 = [0.25, 0.25, 0.25];%For radiotherapy
color8 = [0.6350, 0.0780, 0.1840];%For chemotherapy
color9 = [0,0.392,0];%For surgery

linewidth = 3;%Epaisseur des traits
linewidthpoint = 2;%Epaisseur che et rad

umax = 2.39e8;
%% Masse tumorale totale
load('data_Pb60.txt');
load('data_Pb61.txt');
load('data_Pb62.txt');
load('data_Pb63.txt');
load('data_Pb64.txt');
load('data_Pb65.txt');
figure('Position', get(groot, 'Screensize'));
ymin = umax*1e-6;
ymax = umax;
ymaxchemo = umax*2e-6;
ax = gca;
set(ax,'FontSize',20)
hold on;
xlim([0,100]);
ylim([ymin,ymax]);
xlabel('Time (days)','FontSize', 30);
ylabel('Number of tumour cells in the brain','FontSize', 30);
plot(data_Pb60(:,1),umax*data_Pb60(:,2),'Color',color1,'LineWidth',linewidth);
plot(data_Pb61(:,1),umax*data_Pb61(:,2),'Color',color2,'LineWidth',linewidth);
plot(data_Pb62(:,1),umax*data_Pb62(:,2),'Color',color3,'LineWidth',linewidth);
plot(data_Pb63(:,1),umax*data_Pb63(:,2),'Color',color4,'LineWidth',linewidth);
plot(data_Pb64(:,1),umax*data_Pb64(:,2),'Color',color5,'LineWidth',linewidth);
plot(data_Pb65(:,1),umax*data_Pb65(:,2),'Color',color6,'LineWidth',linewidth);

line([14,14],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([54,54],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([19,19],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([21,21],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([26,26],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([28,28],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([33,33],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([35,35],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([40,40],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([42,42],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([47,47],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([49,49],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([14,19],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([21,26],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([28,33],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([35,40],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([42,47],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([49,54],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([14,14],[ymin,ymax],'Color',color8,'LineStyle','--','LineWidth',linewidthpoint);line([56,56],[ymin,ymax],'Color',color8,'LineStyle','--','LineWidth',linewidthpoint);
line([54,54],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([1,1],[ymin,2*ymin],'Color',color9,'LineStyle','--','LineWidth',linewidthpoint);
ax.YScale = 'log';
legend({'No treatment','Surgery','TMZ','Radiotherapy','TMZ+Radiotherapy','Surgery+TMZ+Radiotherapy'},'Location','northeast','FontSize', 20);

x = [0.25 0.54];    % adjust length and location of arrow 
y = [0.3 0.3];      % adjust hieght and width of arrow
annotation('textarrow',x,y+0.55,'Color',color8,'FontSize',14,'Linewidth',2)
annotation('textbox',[.32 .478 .37 .448],'Color',color8,'EdgeColor','none','String','Chemotherapy','FontSize',30,'Linewidth',2)

annotation('textarrow',x,y-0.13,'Color',color7,'FontSize',11,'Linewidth',2)
annotation('textbox',[.32 .138 .37 .108],'Color',color7,'EdgeColor','none','String','Radiotherapy','FontSize',30,'Linewidth',2)

annotation('textarrow',[0.19 0.14],[0.17 0.12],'Color',color9,'FontSize',11,'Linewidth',2)
annotation('textbox',[.135 .138 .185 .108],'Color',color9,'EdgeColor','none','String','Surgery','FontSize',30,'Linewidth',2)

saveas(gcf,'Tum_number.png');

%% norm(u_T(t,.),infinity)
load('uT_Pb60_max.txt');uT_Pb60_max(1)=4e-1;
load('uT_Pb61_max.txt');uT_Pb61_max(1)=4e-1;
load('uT_Pb62_max.txt');uT_Pb62_max(1)=4e-1;
load('uT_Pb63_max.txt');uT_Pb63_max(1)=4e-1;
load('uT_Pb64_max.txt');uT_Pb64_max(1)=4e-1;
load('uT_Pb65_max.txt');uT_Pb65_max(1)=4e-1;
figure('Position', get(groot, 'Screensize'));
ymin = 3.5e-1;
ymax = 5.6e-1;
ymaxchemo = 3.6e-1;
set(gca,'FontSize',20)
hold on;
xlim([0,100]);
ylim([ymin,ymax]);%ylim([0.39,0.41]);
xlabel('Time (days)','FontSize', 30);
ylabel('||u_T(t,.)||_{L^{\infty}(\Omega)}','FontSize', 30,'rotation',90);
plot(0:100,uT_Pb60_max,'Color',color1,'LineWidth',linewidth);
plot(0:100,uT_Pb61_max,'Color',color2,'LineWidth',linewidth);
plot(0:100,uT_Pb62_max,'Color',color3,'LineWidth',linewidth);
plot(0:100,uT_Pb63_max,'Color',color4,'LineWidth',linewidth);
plot(0:100,uT_Pb64_max,'Color',color5,'LineWidth',linewidth);
plot(0:100,uT_Pb65_max,'Color',color6,'LineWidth',linewidth);

line([14,14],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([54,54],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([19,19],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([21,21],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([26,26],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([28,28],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([33,33],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([35,35],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([40,40],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([42,42],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([47,47],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([49,49],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([14,19],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([21,26],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([28,33],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([35,40],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([42,47],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([49,54],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([14,14],[ymin,ymax],'Color',color8,'LineStyle','--','LineWidth',linewidthpoint);line([56,56],[ymin,ymax],'Color',color8,'LineStyle','--','LineWidth',linewidthpoint);
line([54,54],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([1,1],[ymin,2*ymin],'Color',color9,'LineStyle','--','LineWidth',linewidthpoint);
%ax.YScale = 'log';
line([50,50],[0.394,0.4005],'Color',color7,'LineStyle','-','LineWidth',linewidthpoint);
line([50,60],[0.4005,0.4005],'Color',color7,'LineStyle','-','LineWidth',linewidthpoint);
line([60,60],[0.394,0.4005],'Color',color7,'LineStyle','-','LineWidth',linewidthpoint);
line([50,60],[0.394,0.394],'Color',color7,'LineStyle','-','LineWidth',linewidthpoint);
legend({'No treatment','Surgery','TMZ','Radiotherapy','TMZ+Radiotherapy','Surgery+TMZ+Radiotherapy'},'Location','northeast','FontSize', 20);

x = [0.25 0.54];    % adjust length and location of arrow 
y = [0.3 0.3];      % adjust hieght and width of arrow
annotation('textarrow',x,y+0.55,'Color',color8,'FontSize',14,'Linewidth',2)
annotation('textbox',[.32 .478 .37 .448],'Color',color8,'EdgeColor','none','String','Chemotherapy','FontSize',30,'Linewidth',2)

annotation('textarrow',x,y-0.13,'Color',color7,'FontSize',11,'Linewidth',2)
annotation('textbox',[.32 .138 .37 .108],'Color',color7,'EdgeColor','none','String','Radiotherapy','FontSize',30,'Linewidth',2)

annotation('textarrow',[0.19 0.14],[0.17 0.12],'Color',color9,'FontSize',11,'Linewidth',2)
annotation('textbox',[.135 .138 .185 .108],'Color',color9,'EdgeColor','none','String','Surgery','FontSize',30,'Linewidth',2)

annotation('arrow',[.52 .48],[.33 .42],'Linewidth',2)


a2 = axes();

a2.Position = [0.3 0.4600 0.20 0.25]; % xlocation, ylocation, xsize, ysize
hold on;axis([50 60 0.396 0.3985]);set(gca,'FontSize',15)
plot(50:60,uT_Pb60_max(50:60),'Color',color1,'LineWidth',linewidth);
plot(50:60,uT_Pb61_max(50:60),'Color',color2,'LineWidth',linewidth);
plot(50:60,uT_Pb62_max(50:60),'Color',color3,'LineWidth',linewidth);
plot(50:60,uT_Pb63_max(50:60),'Color',color4,'LineWidth',linewidth);
plot(50:60,uT_Pb64_max(50:60),'Color',color5,'LineWidth',linewidth);
plot(50:60,uT_Pb65_max(50:60),'Color',color6,'LineWidth',linewidth);
line([56,56],[0.396,0.3985],'Color',color8,'LineStyle','--','LineWidth',linewidthpoint);
%axis tight;
annotation('rectangle',[0.26 0.43 0.27 0.3],'Linewidth',2);


hold off;

saveas(gcf,'u_T_infty.png');



%% Concentration en O2 totale
load('data_Pb60.txt');
load('data_Pb61.txt');
load('data_Pb62.txt');
load('data_Pb63.txt');
load('data_Pb64.txt');
load('data_Pb65.txt');
figure('Position', get(0, 'Screensize'));
linewidth = 3;
ax = gca;
hold on;
xlim([0,100]);
ylim([2e-2,2e2]);
xlabel('Time (days)');
ylabel('Total concentration of O2 in mol.1e-6');
plot(data_Pb60(:,1),data_Pb60(:,16),'Color','black','LineWidth',linewidth);
plot(data_Pb61(:,1),data_Pb61(:,16),'Color','red','LineWidth',linewidth);
plot(data_Pb62(:,1),data_Pb62(:,16),'Color','magenta','LineWidth',linewidth);
plot(data_Pb63(:,1),data_Pb63(:,16),'Color','green','LineWidth',linewidth);
plot(data_Pb64(:,1),data_Pb64(:,16),'Color','blue','LineWidth',linewidth);
plot(data_Pb65(:,1),data_Pb65(:,16),'Color','yellow','LineWidth',linewidth);

line([14,14],[2e-2,2e-1],'Color','blue','LineStyle','--','LineWidth',1);line([54,54],[2e-2,2e-1],'Color','blue','LineStyle','--','LineWidth',1);
line([19,19],[2e-2,2e-1],'Color','blue','LineStyle','--','LineWidth',1);line([21,21],[2e-2,2e-1],'Color','blue','LineStyle','--','LineWidth',1);line([26,26],[2e-2,2e-1],'Color','blue','LineStyle','--','LineWidth',1);line([28,28],[2e-2,2e-1],'Color','blue','LineStyle','--','LineWidth',1);line([33,33],[2e-2,2e-1],'Color','blue','LineStyle','--','LineWidth',1);
line([35,35],[2e-2,2e-1],'Color','blue','LineStyle','--','LineWidth',1);line([40,40],[2e-2,2e-2],'Color','blue','LineStyle','--','LineWidth',1);line([42,42],[2e-2,2e-1],'Color','blue','LineStyle','--','LineWidth',1);line([47,47],[2e-2,2e-1],'Color','blue','LineStyle','--','LineWidth',1);line([49,49],[2e-2,2e-1],'Color','blue','LineStyle','--','LineWidth',1);
line([14,19],[2e-1,2e-1],'Color','blue','LineStyle','--','LineWidth',1);line([21,26],[2e-2,2e-1],'Color','blue','LineStyle','--','LineWidth',1);line([28,33],[2e-1,2e-1],'Color','blue','LineStyle','--','LineWidth',1);line([35,40],[2e-1,2e-1],'Color','blue','LineStyle','--','LineWidth',1);line([42,47],[2e-1,2e-1],'Color','blue','LineStyle','--','LineWidth',1);
line([49,54],[2e-1,2e-1],'Color','blue','LineStyle','--','LineWidth',1);
line([14,14],[2e-2,2e2],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);line([56,56],umax*[2e-2,2e2],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);
line([54,54],[2e-2,2e-1],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1);
ax.YScale = 'log';

x = [0.25 0.54];    % adjust length and location of arrow 
y = [0.3 0.3];      % adjust hieght and width of arrow
annotation('textarrow',x,y+0.55,'Color',[0.4660 0.6740 0.1880],'FontSize',9,'Linewidth',2)
annotation('textbox',[.33 .455 .38 .425],'Color',[0.4660 0.6740 0.1880],'EdgeColor','none','String','Chemotherapy','FontSize',11,'Linewidth',2)

annotation('textarrow',x,y-0.13,'Color','blue','FontSize',9,'Linewidth',2)
annotation('textbox',[.33 .118 .38 .088],'Color','blue','EdgeColor','none','String','Radiotherapy','FontSize',11,'Linewidth',2)

%% Masse tumorale en P1
load('data_Pb60.txt');
load('data_Pb61.txt');
load('data_Pb62.txt');
load('data_Pb63.txt');
load('data_Pb64.txt');
load('data_Pb65.txt');
figure('Position', get(0, 'Screensize'));
ymin = 1;
ymax = 1e6;
ymaxchemo = 2;
ax = gca;
set(ax,'FontSize',20);
hold on;
xlim([0,100]);
ylim([ymin,ymax]);
xlabel('Time (days)','FontSize', 30);
ylabel('Number of tumour cells per cm^2 at P1','FontSize', 30);
plot(data_Pb60(:,1),umax*data_Pb60(:,4),'Color',color1,'LineWidth',linewidth);
%plot(data_Pb61(:,1),umax*data_Pb61(:,4),'Color','red','LineWidth',linewidth);
plot(data_Pb62(:,1),umax*data_Pb62(:,4),'Color',color3,'LineWidth',linewidth);
plot(data_Pb63(:,1),umax*data_Pb63(:,4),'Color',color4,'LineWidth',linewidth);
plot(data_Pb64(:,1),umax*data_Pb64(:,4),'Color',color5,'LineWidth',linewidth);
%plot(data_Pb65(:,1),umax*data_Pb65(:,4),'Color','yellow','LineWidth',linewidth);

line([14,14],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([54,54],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([19,19],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([21,21],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([26,26],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([28,28],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([33,33],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([35,35],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([40,40],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([42,42],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([47,47],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([49,49],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([14,19],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([21,26],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([28,33],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([35,40],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([42,47],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([49,54],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([14,14],[ymin,ymax],'Color',color8,'LineStyle','--','LineWidth',linewidthpoint);line([56,56],[ymin,ymax],'Color',color8,'LineStyle','--','LineWidth',linewidthpoint);
line([54,54],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
ax.YScale = 'log';
legend({'No treatment','TMZ','Radiotherapy','TMZ+Radiotherapy'},'Location','best','FontSize',20);

x = [0.25 0.54];    % adjust length and location of arrow 
y = [0.3 0.3];      % adjust hieght and width of arrow
annotation('textarrow',x,y+0.55,'Color',color8,'FontSize',14,'Linewidth',2)
annotation('textbox',[.32 .478 .37 .448],'Color',color8,'EdgeColor','none','String','Chemotherapy','FontSize',30,'Linewidth',2)

annotation('textarrow',x,y-0.13,'Color',color7,'FontSize',11,'Linewidth',2)
annotation('textbox',[.32 .138 .37 .108],'Color',color7,'EdgeColor','none','String','Radiotherapy','FontSize',30,'Linewidth',2)

saveas(gcf,'Tum_number_P1.png');

%% Masse tumorale en P2
load('data_Pb60.txt');
load('data_Pb61.txt');
load('data_Pb62.txt');
load('data_Pb63.txt');
load('data_Pb64.txt');
load('data_Pb65.txt');
figure('Position', get(0, 'Screensize'));
ymin = 1e3;
ymax = umax;
ymaxchemo = 2e3;
ax = gca;
set(ax,'FontSize',20);
hold on;
xlim([0,100]);
ylim([ymin,ymax]);
xlabel('Time (days)','FontSize', 30);
ylabel('Number of tumour cells per cm^2 at P2','FontSize', 30);
plot(data_Pb60(:,1),umax*data_Pb60(:,8),'Color',color1,'LineWidth',linewidth);
%plot(data_Pb61(:,1),umax*data_Pb61(:,8),'Color','red','LineWidth',linewidth);
plot(data_Pb62(:,1),umax*data_Pb62(:,8),'Color',color3,'LineWidth',linewidth);
plot(data_Pb63(:,1),umax*data_Pb63(:,8),'Color',color4,'LineWidth',linewidth);
plot(data_Pb64(:,1),umax*data_Pb64(:,8),'Color',color5,'LineWidth',linewidth);
%plot(data_Pb65(:,1),umax*data_Pb65(:,8),'Color','yellow','LineWidth',linewidth);

line([14,14],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([54,54],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([19,19],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([21,21],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([26,26],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([28,28],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([33,33],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([35,35],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([40,40],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([42,42],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([47,47],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([49,49],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([14,19],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([21,26],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([28,33],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([35,40],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([42,47],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([49,54],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([14,14],[ymin,ymax],'Color',color8,'LineStyle','--','LineWidth',linewidthpoint);line([56,56],[ymin,ymax],'Color',color8,'LineStyle','--','LineWidth',linewidthpoint);
line([54,54],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
ax.YScale = 'log';
legend({'No treatment','TMZ','Radiotherapy','TMZ+Radiotherapy'},'Location','best','FontSize',20);

x = [0.25 0.54];    % adjust length and location of arrow 
y = [0.3 0.3];      % adjust hieght and width of arrow
annotation('textarrow',x,y+0.55,'Color',color8,'FontSize',14,'Linewidth',2)
annotation('textbox',[.32 .478 .37 .448],'Color',color8,'EdgeColor','none','String','Chemotherapy','FontSize',30,'Linewidth',2)

annotation('textarrow',x,y-0.13,'Color',color7,'FontSize',11,'Linewidth',2)
annotation('textbox',[.32 .138 .37 .108],'Color',color7,'EdgeColor','none','String','Radiotherapy','FontSize',30,'Linewidth',2)
saveas(gcf,'Tum_number_P2.png');

%% Masse tumorale en P3
load('data_Pb60.txt');
load('data_Pb61.txt');
load('data_Pb62.txt');
load('data_Pb63.txt');
load('data_Pb64.txt');
load('data_Pb65.txt');
figure('Position', get(0, 'Screensize'));
ymin = 2e2;
ymax = 8e9;
ymaxchemo = 4e2;
ax = gca;
set(ax,'FontSize',20);
hold on;
xlim([0,100]);
ylim([ymin,ymax]);
xlabel('Time (days)','FontSize', 30);
ylabel('Number of tumour cells per cm^2 at P3','FontSize', 30);
plot(data_Pb60(:,1),umax*data_Pb60(:,12),'Color',color1,'LineWidth',linewidth);
plot(data_Pb61(:,1),umax*data_Pb61(:,12),'Color',color2,'LineWidth',linewidth);
plot(data_Pb62(:,1),umax*data_Pb62(:,12),'Color',color3,'LineWidth',linewidth);
plot(data_Pb63(:,1),umax*data_Pb63(:,12),'Color',color4,'LineWidth',linewidth);
plot(data_Pb64(:,1),umax*data_Pb64(:,12),'Color',color5,'LineWidth',linewidth);
plot(data_Pb65(:,1),umax*data_Pb65(:,12),'Color',color6,'LineWidth',linewidth);

line([14,14],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([54,54],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([19,19],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([21,21],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([26,26],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([28,28],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([33,33],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([35,35],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([40,40],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([42,42],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([47,47],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([49,49],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([14,19],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([21,26],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([28,33],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([35,40],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);line([42,47],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([49,54],[ymaxchemo,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([14,14],[ymin,ymax],'Color',color8,'LineStyle','--','LineWidth',linewidthpoint);line([56,56],[ymin,ymax],'Color',color8,'LineStyle','--','LineWidth',linewidthpoint);
line([54,54],[ymin,ymaxchemo],'Color',color7,'LineStyle','--','LineWidth',linewidthpoint);
line([1,1],[ymin,2*ymin],'Color',color9,'LineStyle','--','LineWidth',linewidthpoint);
ax.YScale = 'log';
legend({'No treatment','Surgery','TMZ','Radiotherapy','TMZ+Radiotherapy','Surgery+TMZ+Radiotherapy'},'Location','northeast');

x = [0.25 0.54];    % adjust length and location of arrow 
y = [0.3 0.3];      % adjust hieght and width of arrow
annotation('textarrow',x,y+0.55,'Color',color8,'FontSize',14,'Linewidth',2)
annotation('textbox',[.32 .478 .37 .448],'Color',color8,'EdgeColor','none','String','Chemotherapy','FontSize',30,'Linewidth',2)

annotation('textarrow',x,y-0.13,'Color',color7,'FontSize',11,'Linewidth',2)
annotation('textbox',[.32 .138 .37 .108],'Color',color7,'EdgeColor','none','String','Radiotherapy','FontSize',30,'Linewidth',2)

annotation('textarrow',[0.19 0.14],[0.17 0.12],'Color',color9,'FontSize',11,'Linewidth',2)
annotation('textbox',[.135 .138 .185 .108],'Color',color9,'EdgeColor','none','String','Surgery','FontSize',30,'Linewidth',2)

saveas(gcf,'Tum_number_P3.png');