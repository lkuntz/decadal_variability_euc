load coast;

%------------------------------------------------------------------------

%Figure S3 - EEMD analysis and plot of the 3 lowest-frequency components
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

%Figure S4 - multiple linear regression of the low-frequency time series
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

%% Figure 2 and Figure S2

smoothing = 12;
splineT = 8*12;

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
    