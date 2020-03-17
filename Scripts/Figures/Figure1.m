
addpath(genpath(pwd));

saveFigures = true;

%% Figure1D

clearvars -except saveFigures
clc

%Figure1D left
%network 2.1, parameter 95

% load('/Volumes/MELANOMA/Data/Simulations/2nodes/S_outpar1000_2_1_1')
% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Simulations/2nodes/S_outpar1000_2_1_1')
load('./Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
param = 95;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.25;

figure
subplot(1,3,1)
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(976800:977800,S_outpar{param}(1:2,976800:977800),'LineWidth',2);
hold on 
plot(976800:977800,repmat(thres(param),1,length(976800:977800)),':','LineWidth',...
    2, 'Color', 'k')
set(h, {'color'},  {[135 170 222]./255; [22 45 80]./255});
xlim([976800, 977800])
xticks([976800, 977800])
xticklabels({'0', '1000'})
yticks([0, 3500])
yticklabels({'0', '3.5'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product (10^3)')
xlabel('time')
box off

%Figure1D middle
%network 3.2, parameter 95

% load('/Volumes/MELANOMA/Data/Simulations/3nodes/S_outpar1000_3_2_1')
load('./Data/Simulations/3nodes/S_outpar1000_3_2_1')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
param = 95;

subplot(1,3,2)
i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(310700:311700,S_outpar{param}(1:3,310700:311700),'LineWidth',2);
hold on 
plot(310700:311700,repmat(thres(param),1,length(310700:311700)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255});
xlim([310700, 311700])
xticks([310700, 311700])
xticklabels({'0', '1000'})
yticks([])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
xlabel('time')
box off

%Figure1D right
%network 5.4, parameter 95

% load('/Volumes/MELANOMA/Data/Simulations/5nodes/S_outpar1000_5_4_10')
load('./Data/Simulations/5nodes/S_outpar1000_5_4_10')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
param = 5;

subplot(1,3,3)
i = 3;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(310700:311700,S_outpar{param}(1:5,310700:311700),'LineWidth',2);
hold on 
plot(310700:311700,repmat(thres(95),1,length(310700:311700)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255;...
    [44 90 160]./255; [33 68 120]./255});
xlim([310700, 311700])
xticks([310700, 311700])
xticklabels({'0', '1000'})
yticks([])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
xlabel('time')
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure1/Figure1D'])
    print('-dpdf',['./Figures/Figure1/Figure1D'])
end

%% Figure1E

clearvars -except saveFigures
clc

%Figure1E left
%network 2.1, parameter 52

% load('/Volumes/MELANOMA/Data/Simulations/2nodes/S_outpar1000_2_1_1')
% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Simulations/2nodes/S_outpar1000_2_1_1')
load('./Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
param = 52;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.25;

figure
subplot(1,3,1)
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(328000:329000,S_outpar{param}(1:2,328000:329000),'LineWidth',2);
hold on 
plot(328000:329000,repmat(thres(52),1,length(328000:329000)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'},  {[135 170 222]./255; [22 45 80]./255});
xlim([328000 329000])
xticks([328000, 329000])
xticklabels({'0', '1000'})
ylim([0, 500])
yticks([0, 500])
yticklabels({'0', '5'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product (10^2)')
xlabel('time')
box off

%Figure1E middle
%network 3.2, parameter 52

% load('/Volumes/MELANOMA/Data/Simulations/3nodes/S_outpar1000_3_2_1')
% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Simulations/3nodes/S_outpar1000_3_2_1')
load('./Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
param = 52;

i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(855900:856900,S_outpar{param}(1:3,855900:856900),'LineWidth',2);
hold on 
plot(855900:856900,repmat(thres(52),1,length(855900:856900)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255});
xlim([855900, 856900])
xticks([855900, 856900])
xticklabels({'0', '1000'})
ylim([0, 500])
yticks([])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
xlabel('time')
box off

%Figure1E right
%network 5.4, parameter 52

% load('/Volumes/MELANOMA/Data/Simulations/5nodes/S_outpar1000_5_4_6')
% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Simulations/5nodes/S_outpar1000_5_4_6')
load('./Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
param = 2;

i = 3;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(605800:606800,S_outpar{param}(1:5,605800:606800),'LineWidth',2);
hold on 
plot(605800:606800,repmat(thres(52),1,length(605800:606800)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'},{[135 170 222]./255; [22 45 80]./255; [95 141 211]./255;...
    [44 90 160]./255; [33 68 120]./255});
set(gca,'linewidth',2)
box off
xlim([605800, 606800])
xticks([605800, 606800])
xticklabels({'0', '1000'})
ylim([0, 500])
yticks([])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
xlabel('time')
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure1/Figure1E'])
    print('-dpdf',['./Figures/Figure1/Figure1E'])
end

%% Figure1F

clearvars -except saveFigures
clc

%Figure1F left
%network2.1, parameter 48

% load('./Data/Simulations/2nodes/S_outpar1000_2_1_1')
% load('./Data/Data1000')
load('/Volumes/MELANOMA/Data/Simulations/2nodes/S_outpar1000_2_1_1')
load('/Volumes/MELANOMA/Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
param = 48;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.25;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(499000:500000,S_outpar{param}(1:2,499000:500000),'LineWidth',2);
hold on 
plot(499000:500000,repmat(thres(param),1,length(499000:500000)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'},  {[135 170 222]./255; [22 45 80]./255});
xlim([499000, 500000])
xticks([499000, 500000])
xticklabels({'0', '1000'})
ylim([0, 1800])
yticks([0, 1800])
yticklabels({'0', '1.8'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product (10^3)')
xlabel('time')
box off

%Figure1F middle
%network3.2, parameter 48

% load('/Volumes/MELANOMA/Data/Simulations/3nodes/S_outpar1000_3_2_1')
% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Simulations/3nodes/S_outpar1000_3_2_1')
load('./Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
param = 48;

i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(447000:448000,S_outpar{param}(1:3,447000:448000),'LineWidth',2);
hold on 
plot(447000:448000,repmat(thres(param),1,length(447000:448000)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255});
xlim([447000, 448000])
xticks([447000, 448000])
xticklabels({'0', '1000'})
ylim([0, 1800])
yticks([])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
xlabel('time')
box off

%Figure1F right
%network5.4, parameter 48

% load('/Volumes/MELANOMA/Data/Simulations/5nodes/S_outpar1000_5_4_5')
% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Simulations/5nodes/S_outpar1000_5_4_5')
load('./Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
param = 8;

i = 3;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(146800:147800,S_outpar{param}(1:5,146800:147800),'LineWidth',2);
hold on 
plot(146800:147800,repmat(thres(48),1,length(146800:147800)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255;...
    [44 90 160]./255; [33 68 120]./255});
xlim([146800, 147800])
xticks([146800 147800])
xticklabels({'0', '1000'})
ylim([0, 1800])
yticks([])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
xlabel('time')
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure1/Figure1F'])
    print('-dpdf',['./Figures/Figure1/Figure1F'])
end

%% Figure1G

clearvars -except saveFigures
clc

%Figure1G left
%network 2.1, parameter 915

% load('/Volumes/MELANOMA/Data/Simulations/2nodes/S_outpar1000_2_1_10')
% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Simulations/2nodes/S_outpar1000_2_1_10')
load('./Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
param = 15;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.25;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(286500:287500,S_outpar{param}(1:2,286500:287500),'LineWidth',2);
hold on 
plot(286500:287500,repmat(thres(915),1,length(286500:287500)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'},  {[135 170 222]./255; [22 45 80]./255});
xlim([286500, 287500])
xticks([286500, 287500])
xticklabels({'0', '1000'})
ylim([0, 250])
yticks([0, 250])
yticklabels({'0', '2.5'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product (10^2)')
xlabel('time')
box off


%Figure1G middle
%network 3.2, parameter 915

% load('/Volumes/MELANOMA/Data/Simulations/3nodes/S_outpar1000_3_2_10')
% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Simulations/3nodes/S_outpar1000_3_2_10')
load('./Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
param = 15;

i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(411100:412100,S_outpar{param}(1:3,411100:412100),'LineWidth',2);
hold on 
plot(411100:412100,repmat(thres(915),1,length(411100:412100)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'},  {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255});
xlim([411100, 412100])
xticks([411100 412100])
xticklabels({'0', '1000'})
ylim([0, 250])
yticks([])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
xlabel('time')
box off

%Figure1G right
%network 5.4, parameter 915

% load('/Volumes/MELANOMA/Data/Simulations/5nodes/S_outpar1000_5_4_92')
% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Simulations/5nodes/S_outpar1000_5_4_92')
load('./Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
param = 5;

i = 3;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(327300:328300,S_outpar{param}(1:5,327300:328300),'LineWidth',2);
hold on 
plot(327300:328300,repmat(thres(915),1,length(327300:328300)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255;...
    [44 90 160]./255; [33 68 120]./255});
xlim([327300, 328300])
xticks([327300 328300])
xticklabels({'0', '1000'})
ylim([0, 250])
yticks([])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
xlabel('time')
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure1/Figure1G'])
    print('-dpdf',['./Figures/Figure1/Figure1G'])
end
