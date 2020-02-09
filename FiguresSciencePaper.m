% Figure 1 (Plots of global average temperature vs. time for the different
% simulations) Please go to bottom of script for Figure 1. This part
% doesn't work properly


clear;


%load data from the spin-up run
smoothing = 12;
splineT = 8*12;
%{
load startGlobalOAT32;
bA = TempATM(13:end);
%calculate the monthly average from teh spin-up run years 1-31
AnnualAvg = nan(1,12);
for i=1:12
    AnnualAvg(1,i) = nanmean(bA(i:12:end-24));
end
%load the branched runs with constant forcing (generated with function
%TempOA)
load aGlobalOAT_1;
cA = [TempATM];
cA = cA - repmat(AnnualAvg,[1 length(cA)/12]);
ctime = [0:length(TempATM)-1]/12;

%if some data got corrupted, take it to be the average of the simultaneous
%month 1 year before and 1 year after
r = find(isnan(cA));
cA(r) = (cA(r-12)+cA(r+12))/2;

[modes,imf] = ceemdan(cA,.1,1000,100,1);
Modes{1} = modes;

load aGlobalOAT_8;
cA2 = [TempATM(1:end-12)];
cA2 = cA2 - repmat(AnnualAvg,[1 length(cA2)/12]);
ctime2 = [13:length(TempATM)]/12;

%if some data got corrupted, take it to be the average of the simultaneous
%month 1 year before and 1 year after
r = find(isnan(cA2));
cA2(r) = (cA2(r-12)+cA2(r+12))/2;

r = find(isnan(cA2));
cA2(r) = (cA2(r-24)+cA2(r+24))/2;

[modes,imf] = ceemdan(cA2,.1,1000,100,1);
Modes{2} = modes;

%load the branched runs with constant + Nino3.4 forcing (generated with
%function TempOA)
load aGlobalOAT_3;
ocA = [TempATM];
ocA = ocA - repmat(AnnualAvg,[1 length(ocA)/12]);
octime = [0:length(TempATM)-1]/12;

%if some data got corrupted, take it to be the average of the simultaneous
%month 1 year before and 1 year after
r = find(isnan(ocA));
ocA(r) = (ocA(r-12)+ocA(r+12))/2;

[modes,imf] = ceemdan(ocA,.1,1000,100,1);
Modes{3} = modes;

load aGlobalOAT_9;
ocA2 = [TempATM(1:end-12)];
ocA2 = ocA2 - repmat(AnnualAvg,[1 length(ocA2)/12]);
octime2 = [12:length(TempATM)-1]/12;

%if some data got corrupted, take it to be the average of the simultaneous
%month 1 year before and 1 year after
r = find(isnan(ocA2));
ocA2(r) = (ocA2(r-12)+ocA2(r+12))/2;

[modes,imf] = ceemdan(ocA2,.1,1000,100,1);
Modes{4} = modes;

%plot the global average temperatures vs time in a figure
figure;
subplot(2,2,1);
set(gca,'FontSize',16);
plot(ctime,tsmovavg(cA,'s',smoothing),'r');
hold on;
mode = Modes{1};
plot(ctime, mode(end,:),'k');

title('Constant Forcing - Run 1');
xlabel('Years');
ylabel('Temperature Anomaly (C)');
ylim([-.2 1.6]);
xlim([0 120]);
text(3,1.5,'a',...
        'color','k',...
        'FontSize',9);

subplot(2,2,2);
set(gca,'FontSize',16);
plot(octime,tsmovavg(ocA,'s',smoothing),'r');
title('Constant and Nino 3.4 Forcing - Run 1');
hold on;
mode = Modes{3};
plot(octime, mode(end,:),'k');

xlabel('Years');
ylabel('Temperature Anomaly (C)');
ylim([-.2 1.6]);
xlim([0 120]);
text(3,1.5,'c',...
        'color','k',...
        'FontSize',9);
    
subplot(2,2,3);
set(gca,'FontSize',16);
plot(ctime2,tsmovavg(cA2(1:end),'s',smoothing),'r');
hold on;
mode = Modes{2};
plot(ctime2, mode(end,:),'k');

title('Constant Forcing - Run 2');
xlabel('Years');
ylabel('Temperature Anomaly (C)');
ylim([-.2 1.6]);
xlim([0 120]);
text(3,1.5,'b',...
        'color','k',...
        'FontSize',9);
    
subplot(2,2,4);
set(gca,'FontSize',16);
plot(octime2,tsmovavg(ocA2(1:end),'s',smoothing),'r');
hold on;
mode = Modes{4};
plot(octime2, mode(end,:),'k');

title('Constant and Nino 3.4 Forcing - Run 2');
xlabel('Years');
ylabel('Temperature Anomaly (C)');
ylim([-.2 1.6]);
xlim([0 120]);
text(3,1.5,'d',...
        'color','k',...
        'FontSize',9);
%}
%%
%------------------------------------------------------------------------
 
%Figure 2 - plots the leading low frequency PC regressed onto the
%temperature data for each run. EOF analysis was performed using EOFs_CESM
%function and results are simply plotted here.

hspace = .08
height = (1-2*hspace)/2;
%hspace = .01;
hlegend = .1;
width = (1-hlegend-2*hspace)/2;

% horizontal plot
p = [0 1-(hspace+height) width height; hspace+width 1-(hspace+height) width height;...
    0 1-2*(hspace+height) width height; 
    hspace+width 1-2*(hspace+height) width height]; 


f=figure; set(f,'Units','normalized');
    set(f,'Position',[0 0 .8 1]);

versions = [1,3,8,9];
title1 = 'Constant forcing - run 1'
title2 = 'Constant and Nino 3.4 forcing - run 1'
title3 = 'Constant forcing - run 2'
title4 = 'Constant and Nino 3.4 forcing - run 2'

LatB = linspace(-90,90,97);
LonB = linspace(0,360,145);
[LonB,LatB] = meshgrat(LonB,LatB);
lon_data = LonB; 
lat_data = LatB;
load coast;

% width = 45;
% height = 30;
% figure('units','centimeters','Position',[1 1 width height]);
for i=1:4
    subplot('position',p(i,:));
    set(gca,'FontSize',22,'FontName','Myriad');
    load(strcat('EOFs_',num2str(versions(i))),'RegressedTs');
    toplot = RegressedTs;
    toplot(isnan(toplot))=0;
    axesm('miller','Origin',180);
    pcolorm(lat_data,lon_data,toplot);
    load coast;
    hold on;
    plotm(lat,long,'k');
    tightmap;
    if i==1
        title(title1);
    elseif i==2
        title(title2);
    elseif i==3
        title(title3);
    else
        title(title4);
    end
    colormap(b2r_noWhite(-1.5, 1.5));
end
axes('Position', [0.05 0.015 0.89 0.9], 'Visible', 'off'); 
    set(gca,'CLim',[-1.5, 1.5]);
c=colorbar ('FontSize',22,'FontName','Myriad'); ylabel(c,'Sea surface temperature (C)') 

label = {'A','B','C','D'};

axes('Position', [0 0 1 1], 'Visible', 'off'); 
text([.01;.45;.01;.45],[.97;.97;.465;.465],label,...
        'color','k',...
        'FontSize',22,'FontName','Myriad','FontWeight','bold');

% %Add the subplot labels
% uicontrol('style','text','string','A','units','centimeters','position',[5 height-3 .5 .5],'FontSize',12);
% uicontrol('style','text','string','B','units','centimeters','position',[5 height/2-1.25 .5 .5],'FontSize',12);
% uicontrol('style','text','string','C','units','centimeters','position',[width/2+2.5 height-3 .5 .5],'FontSize',12);
% uicontrol('style','text','string','D','units','centimeters','position',[width/2+2.5 height/2-1.25 .5 .5],'FontSize',12);

%%
%------------------------------------------------------------------------
 
%Figure 3 - From gostaChangingCoverage.m in Harvard EPS/ENSO/

%------------------------------------------------------------------------

%Figure 4 - From HContentVideo

%------------------------------------------------------------------------

%Figure S1
%determines the scaling between subduction temperature standard deviations
%and standard deviations in the Nino3.4 region. Then plots the results of
%the monthly ratios of standard deviations

file = 'tos_Omon_CESM1-CAM5_piControl_r1i1p1_000101-031912.nc';
%cmip5.output1.NSF-DOE-NCAR.CESM1-CAM5.piControl.mon.ocean.Omon.r1i1p1.v20140822|tds.ucar.edu

%load pi control data
[V,G] = nc_readPH(file);

%determine the weights using Lat/Lon bounds for the grid
V(7).data = double(V(7).data);
V(8).data = double(V(8).data);
weights = areaquad(V(7).data(1,:,:),V(8).data(1,:,:),V(7).data(3,:,:),V(8).data(2,:,:));
weights = squeeze(weights);

lat_data = double(V(5).data);
lon_data = double(V(6).data);

%bounding region for N subduction region
NsubN = 40;
NsubS = 30;
NsubE = 360-130;
NsubW = 160;

%weights for N subduction region - outside of region weights are 'NaN'
weightN = weights;
weightN(V(5).data>NsubN | V(5).data<NsubS) = NaN;
weightN(V(6).data>NsubE | V(6).data<NsubW) = NaN;

%bounding region for S subduction zone
SsubN = -25;
SsubS = -40;
SsubE = 360-85;
SsubW = 360-160;
%weights for S subduction zone
weightS = weights;
weightS(V(5).data>SsubN | V(5).data<SsubS) = NaN;
weightS(V(6).data>SsubE | V(6).data<SsubW) = NaN;

%bounding region and weights for Nino3.4 region
Nino34N = 5;
Nino34S = -5;
Nino34E = 360-120;
Nino34W = 360-170;
weightNino = weights./weights;
weightNino(V(5).data>Nino34N | V(5).data<Nino34S) = NaN;
weightNino(V(6).data>Nino34E | V(6).data<Nino34W) = NaN;
V(9).data(V(9).data>500)=nan;

%replicate matrix of weights over the whole time series
WeightN = repmat(weightN, [1 1 size(V(9).data,3)]);
WeightS = repmat(weightS, [1 1 size(V(9).data,3)]);
WeightNino = repmat(weightNino, [1 1 size(V(9).data,3)]);

%get the Nino3.4 region temperatures of the time series
NinoT = WeightNino.*V(9).data;

%for each month determine the scaling between subduction temperatures and
%Nino3.4 temperature standard deviations
for i=1:12
    %get the time series of the mean subducting temperature using 1/3 N and
    %2/3 S regions
    SubT = 1/3*nansum(nansum(WeightN(:,:,i:12:end).*V(9).data(:,:,i:12:end)))/nansum(nansum(weightN))...
    +2/3*nansum(nansum(WeightS(:,:,i:12:end).*V(9).data(:,:,i:12:end)))/nansum(nansum(weightS));
    Sub_T(:,i) = squeeze(SubT);
    %average monthly temperature of subducting region
    Clim_SubT(i) = squeeze(nanmean(Sub_T(:,i)));
    %Temperature standard deviation in Nino region associated with 1
    %standard deviation of subducting temperature
    Scale(:,:,i) = std(NinoT(:,:,i:12:end),0,3)/std(Sub_T(:,i),0,1);
end

%only keep the data in the array if it is in the Nino3.4 region
lat_data(:,~any(~isnan(Scale(:,:,1)),1)) = [];
lon_data(:,~any(~isnan(Scale(:,:,1)),1)) = [];
Scale(:,~any(~isnan(Scale(:,:,1)),1),:) = [];
lat_data(~any(~isnan(Scale(:,:,1)),2),:) = [];
lon_data(~any(~isnan(Scale(:,:,1)),2),:) = [];
Scale(~any(~isnan(Scale(:,:,1)),2),:,:) = [];


%Labels fro graphs
months = {'January','February','March','April','May','June','July','August','September','October','November','December'};
%mycolormap = [ ones(1,3); jet(100)];
%for each month, plot the scaling of the ratio of standard deviations in
%the Nino3.4 region to the std of the subducting temperatures

hspace = .05
height = (1-4*hspace)/4;
%hspace = .01;
hlegend = .1;
width = (1-hlegend-2*hspace)/3;

% horizontal plot
p = [0 1-(hspace+height) width height; hspace+width 1-(hspace+height) width height;...
    2*(hspace+width) 1-(hspace+height) width height; 0 1-2*(hspace+height) width height; 
    hspace+width 1-2*(hspace+height) width height; 2*(hspace+width) 1-2*(hspace+height) width height;...
    0 1-3*(hspace+height) width height; hspace+width 1-3*(hspace+height) width height;...
    2*(hspace+width) 1-3*(hspace+height) width height; 0 1-4*(hspace+height) width height; 
    hspace+width 1-4*(hspace+height) width height; 2*(hspace+width) 1-4*(hspace+height) width height]; 


f=figure; set(f,'Units','normalized');
    set(f,'Position',[0 0 1 1]);

for i=1:12
    subplot('Position',p(i,:));
    toplot = squeeze(Scale(:,:,i));
    set(gca,'FontSize',22,'FontName','Myriad');
%    toplot(isnan(toplot))=0;
    axesm('miller','Origin',180,'MapLatLimit',[-15,15],'MapLonLimit',[180,285]);
    pcolorm(lat_data,lon_data,toplot);
    colormap(jet);
    caxis([0 8]);
    load coast;
    hold on;
    plotm(lat,long,'k');
    tightmap;
    gridm;
    title(months{i});
end
%add one colorbar that applies to the whole graph
axes('Position', [0.05 0.025 0.92 0.9], 'Visible', 'off'); 
    set(gca,'CLim',[0, 8]);

c=colorbar('FontSize',22,'FontName','Myriad'); ylabel(c,'Ratio of standard deviations of temperature anomalies')

 %%
%------------------------------------------------------------------------

%Figure S2 - EEMD analysis and plot of the 3 lowest-frequency components
%from the decomposition

titles = {'Constant Forcing - Run 1', 'Constant Forcing - Run 2', 'Constant and Nino 3.4 Step Forcing',...
     'Constant and Nino 3.4 Sine Forcing'};

 %load global temperature data for the constant forcing runs
load aGlobalOAT_1;
constA = NaN(2,length(TempATM));
constA(1,:) = TempATM;
ctime = [0:length(TempATM)-1]/12;

load aGlobalOAT_8;
constA(2,:) = TempATM;

%calculate the monthly mean data over the length of the runs
AnnualAvg = nan(2,12);
for i=1:12
    AnnualAvg(:,i) = nanmean(constA(:,i:12:end),2);
end

AnnualAvg = repmat(AnnualAvg,[1,size(constA,2)/12]);

%remove monthly means from runs
constA = constA-AnnualAvg;

%if and data is missing/corrupted replace it by taking the average of the
%temperature in that month a year before and year after
for i = 1:2
    temp = constA(i,:);
    r = find(isnan(temp))
    temp(r) = (temp(r-12)+temp(r+12))/2;
    if i==2
        r = find(isnan(temp));
        temp(r) = (temp(r-24)+temp(r+24))/2;
    end
    constA(i,:) = temp;
end

hspace = .09;
height = .90;
offset = .06;
p = [offset offset+2*height/3 height/4-hspace/4 height/3-hspace/3;...
    offset+height/4-hspace/16 offset+2*height/3 height/4-hspace/4 height/3-hspace/3;...
    offset+2*height/4-hspace/16 offset+2*height/3 height/4-hspace/4 height/3-hspace/3;...
    offset+3*height/4-hspace/16 offset+2*height/3 height/4-hspace/4 height/3-hspace/3;...
    offset offset+height/3 height/4-hspace/4 height/3-hspace/3;...
    offset+height/4-hspace/16 offset+height/3 height/4-hspace/4 height/3-hspace/3;...
    offset+2*height/4-hspace/16 offset+height/3 height/4-hspace/4 height/3-hspace/3;...
    offset+3*height/4-hspace/16 offset+height/3 height/4-hspace/4 height/3-hspace/3;...
    offset offset height/4-hspace/4 height/3-hspace/3;...
    offset+height/4-hspace/16 offset height/4-hspace/4 height/3-hspace/3;...
    offset+2*height/4-hspace/16 offset height/4-hspace/4 height/3-hspace/3;...
    offset+3*height/4-hspace/16 offset height/4-hspace/4 height/3-hspace/3;];

range = [-.8 .8; -.2 .2; -.15 .15];
%calculate the EEMD using 1000 simulations and plot the 3 lowest frequency 
%components of the EEMD analysis
figure;
for i=1:2;
    [modes,imf] = ceemdan(constA(i,:),.1,1000,100,1);
    for j=1:3    
        subplot('Position',p((j-1)*4+i,:));
        plot((0:size(modes,2)-1)/12, modes(end-j+1,:));
        if j==1
            title(titles{i},'FontSize',16,'FontName','Myriad');
        end
        ylim(range(j,:));
        if i==2
            if j<3
                set(gca,'XTick',[],'YTick',[]);
            else
                set(gca,'YTick',[]);
            end
        end
        if i==1
            if j<3
                set(gca,'XTick',[]);
            end
        end
        set(gca,'FontSize',14,'FontName','Myriad');
    end
    Modes{i} = modes;
%     xlabel('Years');
%     ylabel('Global Surface Temperature Anomaly (C)');
end
 
%load global temperature data for the constant + Nino3.4 forcing runs
%load aGlobalOAT_3;
load GridT30AStep;
TS = TS(1:1440);
oceanA = NaN(2,length(TS));
oceanA(1,:) = TS;
%oceanA(1,:) = TempATM;
otime = [0:length(TempATM)-1]/12;

% load aGlobalOAT_9;
% oceanA(2,13:end) = TempATM(1:end-12);
load GridT50ASine;
TS = TS(1:1440);
oceanA(2,1:length(TS)) = TS;
%calculate the monthly mean data over the length of the runs
AnnualAvg = nan(2,12);
for i=1:12
    AnnualAvg(:,i) = nanmean(oceanA(:,i:12:end),2);
end

%remove monthly means from runs
AnnualAvg = repmat(AnnualAvg,[1,size(oceanA,2)/12]);

oceanA = oceanA-AnnualAvg;

%if and data is missing/corrupted replace it by taking the average of the
%temperature in that month a year before and year after
% for i = 1:2
%     temp = oceanA(i,:);
%     r = find(isnan(temp));
%     if i==2
%         r = r(13:end);
%     end
%     temp(r) = (temp(r-12)+temp(r+12))/2;
%     oceanA(i,:) = temp;
% end

%calculate the EEMD using 1000 simulations and plot the 3 lowest frequency 
%components of the EEMD analysis
for i=1:2;
    if i==2
        [modes,imf] = ceemdan(oceanA(i,1:length(TS)),.1,1000,100,1);
    else
        [modes,imf] = ceemdan(oceanA(i,:),.1,1000,100,1);
    end
    for j=1:3    
        subplot('Position',p((j-1)*4+i+2,:));
        plot((0:size(modes,2)-1)/12, modes(end-j+1,:));
        if j==1
            title(titles{i+2},'FontSize',16,'FontName','Myriad');
        end
        ylim(range(j,:));
        xlim([0 120]);
        if j<3
            set(gca,'XTick',[],'YTick',[]);
        else
            set(gca,'YTick',[]);
        end
        set(gca,'FontSize',14);
    end
    Modes{i+2} = modes;
end
ax = axes('Position', [.0 .0 1 1],'Visible','off');
text(.5, .02, 'Years','FontSize',16,'FontName','Myriad');
text(.02,.33,'Global surface temperature anomaly (C)','FontSize',16,'rotation',90);
axis off;
label = {'A','B','C','D'};
text([.06;.28;.505;.73],[.98;.98;.98;.98],label,...
        'color','k',...
        'FontSize',16,'FontName','Myriad','FontWeight','bold');
%%
%-------------------------------------------------------------------------

%Figure S3 - multiple linear regression of the low-frequency time series
%from EEMD analysis onto the spatial pattern of temperature

r = 15;

hspace = .09;
height = .90;
offset = .03;
p = [offset offset+2*height/3 height/4-hspace/4 height/3-hspace/3;...
    offset+height/4-hspace/16 offset+2*height/3 height/4-hspace/4 height/3-hspace/3;...
    offset+2*height/4-hspace/16 offset+2*height/3 height/4-hspace/4 height/3-hspace/3;...
    offset+3*height/4-hspace/16 offset+2*height/3 height/4-hspace/4 height/3-hspace/3;...
    offset offset+height/3 height/4-hspace/4 height/3-hspace/3;...
    offset+height/4-hspace/16 offset+height/3 height/4-hspace/4 height/3-hspace/3;...
    offset+2*height/4-hspace/16 offset+height/3 height/4-hspace/4 height/3-hspace/3;...
    offset+3*height/4-hspace/16 offset+height/3 height/4-hspace/4 height/3-hspace/3;...
    offset offset height/4-hspace/4 height/3-hspace/3;...
    offset+height/4-hspace/16 offset height/4-hspace/4 height/3-hspace/3;...
    offset+2*height/4-hspace/16 offset height/4-hspace/4 height/3-hspace/3;...
    offset+3*height/4-hspace/16 offset height/4-hspace/4 height/3-hspace/3;];
%load temperature grids from the runs (complied using the function
%TempGrid)
load('GridT_1','TempAGrid');
for i=1:12
    TempAGrid(i:12:end,:,:) = TempAGrid(i:12:end,:,:) - ...
        repmat(nanmean(TempAGrid(i:12:end,:,:)),[size(TempAGrid(i:12:end,:,:),1) 1 1]);
end

%multiple linear regression of 3 lowest frequency signals from EEMD onto
%the timeseries of gridded temperature
RegressedTs = nan([3 (size(squeeze(TempAGrid(1,:,:))))]);
mode = cell2mat(Modes(1));
x1 = squeeze(mode(end,:)); x1 = x1(:);
x2 = squeeze(mode(end-1,:)); x2 = x2(:);
x3 = squeeze(mode(end-2,:)); x3 = x3(:);
for i=1:size(TempAGrid,2)
    for j=1:size(TempAGrid,3)
        data = TempAGrid(:,i,j);
        data = data(:);
        b = regress(data,[x1 x2 x3]);
        RegressedTs(:,i,j) = b;
    end
end

%plot the spatial pattern of results
LatB = linspace(-90,90,97);
LonB = linspace(0,360,145);
f=figure; set(f,'Units','normalized');
set(f,'Position',[0 0 1 1]);
load coast;
for i=1:3
    subplot('Position',p((i-1)*4+1,:))
    axesm('miller','Origin',180);
    pcolorm(LatB,LonB,squeeze(RegressedTs(i,:,:))');
%     r = max(max(abs(squeeze(RegressedTs(i,:,:)))));
    colormap(b2r_noWhite(-r,r));
    hold on;
    plotm(lat,long,'k');
    tightmap;
    if i==1
        title(titles{1},'FontSize',16,'FontName','Myriad');
    end
end

%load temperature grids from the runs (complied using the function
%TempGrid)
load('GridT_8');
for i=1:12
    TempAGrid(i:12:end,:,:) = TempAGrid(i:12:end,:,:) - ...
        repmat(nanmean(TempAGrid(i:12:end,:,:)),[size(TempAGrid(i:12:end,:,:),1) 1 1]);
end

%multiple linear regression of 3 lowest frequency signals from EEMD onto
%the timeseries of gridded temperature
RegressedTs = nan([3 (size(squeeze(TempAGrid(1,:,:))))]);
mode = cell2mat(Modes(2));
x1 = squeeze(mode(end,:)); x1 = x1(:);
x2 = squeeze(mode(end-1,:)); x2 = x2(:);
x3 = squeeze(mode(end-2,:)); x3 = x3(:);
for i=1:size(TempAGrid,2)
    for j=1:size(TempAGrid,3)
        data = TempAGrid(:,i,j);
        data = data(:);
        b = regress(data,[x1 x2 x3]);
        RegressedTs(:,i,j) = b;
    end
end

for i=1:3
    subplot('Position',p((i-1)*4+2,:))
    axesm('miller','Origin',180);
    pcolorm(LatB,LonB,squeeze(RegressedTs(i,:,:))');
%     r = max(max(abs(squeeze(RegressedTs(i,:,:)))));
    colormap(b2r_noWhite(-r,r));
    hold on;
    plotm(lat,long,'k');
    tightmap;
    if i==1
        title(titles{2},'FontSize',16,'FontName','Myriad');
    end
end


%load temperature grids from the runs (complied using the function
%TempGrid)
load('GridT30AStep');
TempAGrid = TempAGrid(1:1440,:,:);
%load('GridT_3');
for i=1:12
    TempAGrid(i:12:end,:,:) = TempAGrid(i:12:end,:,:) - ...
        repmat(nanmean(TempAGrid(i:12:end,:,:)),[size(TempAGrid(i:12:end,:,:),1) 1 1]);
end

%multiple linear regression of 3 lowest frequency signals from EEMD onto
%the timeseries of gridded temperature
RegressedTs = nan([3 (size(squeeze(TempAGrid(1,:,:))))]);
mode = cell2mat(Modes(3));
x1 = squeeze(mode(end,:)); x1 = x1(:);
x2 = squeeze(mode(end-1,:)); x2 = x2(:);
x3 = squeeze(mode(end-2,:)); x3 = x3(:);
for i=1:size(TempAGrid,2)
    for j=1:size(TempAGrid,3)
        data = TempAGrid(:,i,j);
        data = data(:);
        b = regress(data,[x1 x2 x3]);
        RegressedTs(:,i,j) = b;
    end
end
   
for i=1:3
    subplot('Position',p((i-1)*4+3,:))
    axesm('miller','Origin',180);
    pcolorm(LatB,LonB,squeeze(RegressedTs(i,:,:))');
%     r = max(max(abs(squeeze(RegressedTs(i,:,:)))));
    colormap(b2r_noWhite(-r,r));
    hold on;
    plotm(lat,long,'k');
    tightmap;
    if i==1
        title(titles{3},'FontSize',16,'FontName','Myriad');
    end
end

%load temperature grids from the runs (complied using the function
%TempGrid)
load('GridT50ASine');
%load('GridT_9');
TempAGrid = TempAGrid(1:1440,:,:);
for i=1:12
    TempAGrid(i:12:end,:,:) = TempAGrid(i:12:end,:,:) - ...
        repmat(nanmean(TempAGrid(i:12:end,:,:)),[size(TempAGrid(i:12:end,:,:),1) 1 1]);
end

%multiple linear regression of 3 lowest frequency signals from EEMD onto
%the timeseries of gridded temperature
RegressedTs = nan([3 (size(squeeze(TempAGrid(1,:,:))))]);
mode = cell2mat(Modes(4));
x1 = squeeze(mode(end,:)); x1 = x1(:);
x2 = squeeze(mode(end-1,:)); x2 = x2(:);
x3 = squeeze(mode(end-2,:)); x3 = x3(:);
for i=1:size(TempAGrid,2)
    for j=1:size(TempAGrid,3)
        data = TempAGrid(:,i,j);
        data = data(:);
        b = regress(data,[x1 x2 x3]);
        RegressedTs(:,i,j) = b;
    end
end

for i=1:3
    subplot('Position',p((i-1)*4+4,:))
    axesm('miller','Origin',180);
    pcolorm(LatB,LonB,squeeze(RegressedTs(i,:,:))');
%     r = max(max(abs(squeeze(RegressedTs(i,:,:)))));
    colormap(b2r_noWhite(-r,r));
    hold on;
    plotm(lat,long,'k');
    tightmap;
    if i==1
        title(titles{4},'FontSize',16,'FontName','Myriad');
    end
end
axes('Position', [0 0 1 1], 'Visible', 'off');
label = {'A','B','C','D'};
text([.03;.25;.475;.7],[.95;.95;.95;.95],label,...
        'color','k',...
        'FontSize',16,'FontName','Myriad','FontWeight','bold');

axes('Position', [0.05 0.03 0.92 0.88], 'Visible', 'off');
    set(gca,'CLim',[-r, r],'FontSize',16,'FontName','Myriad');

c=colorbar('FontSize',16,'FontName','Myriad'); ylabel(c,'Temperature (C)')

%%
label = {'A','B'};
height = 1; hspace = .09;
width = .4;
f=figure; set(f,'Units','normalized');
set(f,'Position',[0 0 1 1]);

height = .8;
offset = .1;
p = [offset offset width height;...
    offset*1.5+width offset width height];

for i=1:2
    subplot('position',p(i,:));
    %i=i+2;
    if i==1
        cA = constA(1,:);
        k = 1;
    elseif i==2
        cA = oceanA(1,:);
        k = 3;
    elseif i==3
        cA = constA(2,:);
        k = 2;
    else
        cA = oceanA(2,:);
        k = 4;
    end
    
    plot([0:length(cA)-1]/12,tsmovavg(cA,'s',smoothing),'r');
    hold on;
    mode = Modes{k};
    plot([0:length(mode(end,:))-1]/12, mode(end,:),'k');

    title(titles{k}); set(gca,'FontSize',18,'FontName','Myriad');
    %xlabel('Years');
    %ylabel('Temperature Anomaly (C)');
    ylim([-.85 .85]);
    xlim([0 120]);
    if i==1
        set(gca,'YTick',[-.5,0,.5],'XTick',[0,30,60,90,120]);
        ylabel('Temperature anomaly (C)');
        xlabel('Years');

    elseif i==2
        set(gca,'XTick',[0,30,60,90,120],'YTick',[]);
        xlabel('Years');

    elseif i==3
        set(gca,'YTick',[-.5,0,.5],'XTick',[0,30,60,90,120]);
        xlabel('Years');
        ylabel('Temperature anomaly (C)');
    else
        set(gca,'YTick',[],'XTick',[0,30,60,90,120]);
        xlabel('Years');
    end
end

axes('Position', [0 0 1 1], 'Visible', 'off'); 
text([.11;.56],[.86;.86],label,...
        'color','k',...
        'FontSize',16,'FontName','Myriad','FontWeight','bold');
    

%%
label = {'A','B','C','D'};
height = 1; hspace = .09;
width = .55;
f=figure; set(f,'Units','normalized');
set(f,'Position',[0 0 width height]);

height = .94;
offset = .06;
p = [offset*1.5 offset+height/2 height/2-hspace/2 height/2-hspace/2;...
    offset*1.5+height/2-hspace/16 offset+height/2 height/2-hspace/2 height/2-hspace/2;...
    offset*1.5 offset height/2-hspace/2 height/2-hspace/2;...
    offset*1.5+height/2-hspace/16 offset height/2-hspace/2 height/2-hspace/2];

for i=1:4
    subplot('position',p(i,:));
    if i==1
        cA = constA(1,:);
        k = 1;
    elseif i==2
        cA = oceanA(1,:);
        k = 3;
    elseif i==3
        cA = constA(2,:);
        k = 2;
    else
        cA = oceanA(2,:);
        k = 4;
    end
    
    plot([0:length(cA)-1]/12,tsmovavg(cA,'s',smoothing),'r');
    hold on;
    mode = Modes{k};
    plot([0:length(mode(end,:))-1]/12, mode(end,:),'k');

    title(titles{k}); set(gca,'FontSize',18,'FontName','Myriad');
    %xlabel('Years');
    %ylabel('Temperature Anomaly (C)');
    ylim([-.85 .85]);
    xlim([0 120]);
    if i==1
        set(gca,'YTick',[-.5,0,.5],'XTick',[]);
        ylabel('Temperature anomaly (C)');
    elseif i==2
            set(gca,'XTick',[],'YTick',[]);
    elseif i==3
        set(gca,'YTick',[-.5,0,.5],'XTick',[0,30,60,90,120]);
        xlabel('Years');
        ylabel('Temperature anomaly (C)');
    else
        set(gca,'YTick',[],'XTick',[0,30,60,90,120]);
        xlabel('Years');
    end
end

axes('Position', [0 0 1 1], 'Visible', 'off'); 
text([.1;.56;.1;.56],[.94;.94;.47;.47],label,...
        'color','k',...
        'FontSize',14,'FontName','Myriad','FontWeight','bold');
    
%%
%load global temperature data for the constant + Nino3.4 forcing runs
load aGlobalOAT_7;
oceanA = TempATM;
otime = [0:length(TempATM)-1]/12;

%calculate the monthly mean data over the length of the runs
AnnualAvg = nan(1,12);
for i=1:12
    AnnualAvg(:,i) = nanmean(oceanA(:,i:12:end),2);
end

%remove monthly means from runs
AnnualAvg = repmat(AnnualAvg,[1,size(oceanA,2)/12]);

oceanA = oceanA-AnnualAvg;

%if and data is missing/corrupted replace it by taking the average of the
%temperature in that month a year before and year after
temp = oceanA(:);
r = find(isnan(temp));
temp(r) = (temp(r-12)+temp(r+12))/2;
oceanA(:) = temp;

temp = oceanA(:);
r = find(isnan(temp));
temp(r) = (temp(r-1)+temp(r+1))/2;
oceanA(:) = temp;

%calculate the EEMD using 1000 simulations
[modes,imf] = ceemdan(oceanA,.1,1000,100,1);
  
height = 1; hspace = .09;
width = .65;
f=figure; set(f,'Units','normalized');
set(f,'Position',[0 0 width height]);
   
plot([0:length(oceanA)-1]/12,tsmovavg(oceanA,'s',smoothing),'r');
hold on;
plot([0:length(modes(end,:))-1]/12, modes(end,:),'k');

title('Constant and Nino 3.4 forcing - 20yr Lag  g'); set(gca,'FontSize',18,'FontName','Myriad');
xlabel('Years');
ylabel('Temperature Anomaly (C)');
ylim([-.7 .75]);
xlim([0 110]);

set(gca,'YTick',[-.5,0,.5],'XTick',[0,30,60,90]);
