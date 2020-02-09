%Figures 1 and S1
%explores changes in the EUC in terms of temperature and velocity using the
%TAO array data

clear;
%Load the temperature data from the TAO array
load('AllT(z)Data');

toPlot = 'UPROFILE';

%make sure all data is double not single
All_temperature = double(All_temperature);
All_depth = double(All_depth);
All_latitude = double(All_latitude);
All_longitude = double(All_longitude);

%latitudes and longitudes along equator which adcp, current profilers, and
%temperature readings are all present
lat = 0;
LON = [220, 250];
%create figure
f=figure; set(f,'units','normalized','position',[0 0 1 1]);
HLocation = nan(length(LON),425);
Transported = nan(length(LON),425);

height = .41;
width = .85;
hspace = .09;
wspace = .05;

p = [wspace 1-hspace-height width height-.15; wspace 1-2*(hspace+height) width height-.15;...
    wspace 1-hspace-.13 width .15; wspace 1-2*(hspace)-height-.13 width .15;];
%for each longitude location
for k = 1:length(LON)

    lon = LON(k);

    La = find(All_latitude==lat);
    Lo = find(All_longitude==lon);

    %extract the temperature record in time and depth at that location
    TLocation = All_temperature(:,:,La,Lo);

    %use a cubic fit to interpolate to depths between upper and lower
    %reading (so lon as there is some data present)
    for j=1:size(TLocation)
        Tdata = squeeze(TLocation(j,:));
        Ddata = All_depth;
        Ddata(isnan(Tdata)) = [];
        Tdata(isnan(Tdata))=[];
        if length(Ddata)>1
            Dmin = find(All_depth==Ddata(1));
            Dmax = find(All_depth==Ddata(end));
            Tinterp = interp1(Ddata,Tdata,All_depth(Dmin:Dmax),'pchip');
            TLocation(j,Dmin:Dmax)=Tinterp;
        end

    end

    %load the current readings
    load('AllCur(z)Data');

    La = find(curLat == lat);
    Lo = find(curLon == lon);
    %keep only data at lat/lon of interest
    cur_loc = Ucur(:,:,La,Lo);

    %load the dopler profiler current readings (these are more reliable)
    load('adcpData');

    La = find(adcp_lat == lat);
    Lo = find(adcp_lon == lon);

    %keep only data at lat/lon of interest
    ADCP_loc = adcp_u(:,:,La,Lo);

    %combine the depth scale of the doppler and current data removing any
    %duplicates
    depthComb = [curDepth;double(adcp_depth)];
    depthComb = sort(depthComb);
    depthComb = unique(depthComb);

    %get the index of the current readings where the adcp readings being
    %and end
    iStartTime = find(abs(curTime-adcp_time(1))==min(abs(curTime-adcp_time(1))));
    iEndTime = find(abs(curTime-adcp_time(end))==min(abs(curTime-adcp_time(end))));

    %create array to hold combined current data
    ULocation = nan(length(curTime),length(depthComb));

    %for each depth of the current readings, add it to the combined data
    c = find(ismember(depthComb,curDepth));
    for i=1:length(c)
        ULocation(:,c(i)) = cur_loc(:,i);
    end

    %so long as the adcp data exists, add it to the combined data
    %(overwriting the current reading if there is one at the same depth and
    %time)
    c = find(ismember(depthComb,adcp_depth));
    for i=1:length(c)
        for n=1:(iEndTime-iStartTime+1)
            if ~isnan(ADCP_loc(n,i))
                ULocation(iStartTime+n-1,c(i)) = ADCP_loc(n,i);
            end
        end
    end

    %use a cubic fit to interpolate to depths between upper and lower
    %reading (so lon as there is some data present)
    for j=1:size(ULocation,1)
        Udata = squeeze(ULocation(j,:));
        Ddata = depthComb;
        Ddata(isnan(Udata)) = [];
        Udata(isnan(Udata))=[];
        if length(Ddata)>1
            Dmin = find(depthComb==Ddata(1));
            Dmax = find(depthComb==Ddata(end));
            Uinterp = interp1(Ddata,Udata,depthComb(Dmin:Dmax),'pchip');
            ULocation(j,Dmin:Dmax)=Uinterp;
        end

    end
    
    %plot either the temperature profiles with velocity contours, velocity
    %profiles with temperature contours, or anomalies with profile contours

    subplot('position',p(k,:));
    contourf(curTime,depthComb,ULocation',[0:10:160],'LineStyle','none');
    colormap('jet');
    shading flat;
    hold on;
    caxis([0 160]);
    set(gca,'YDir','reverse','fontsize',20);
    ylabel('Depth (m)');
    xlim([1980,2015]);
    %if k<length(LON)
        %set(gca,'XTickLabel',[],'XLabel',[]);
    %else
        xlabel('Year');
    %end
    ylim([0 300]);
    
end

%%
%Load the temperature data from the TAO array
%latitudes and longitudes along equator which adcp, current profilers, and
%temperature readings are all present
lat = 0;
LON = [220, 250];

UTOT = [];
UMA = [];
%for each longitude location
for k = 1:length(LON)

    lon = LON(k);
    
    %extract the temperature record in time and depth at that location

    %load the current readings
    load('AllCur(z)Data');

    La = find(curLat == lat);
    Lo = find(curLon == lon);
    %keep only data at lat/lon of interest
    cur_loc = Ucur(:,:,La,Lo);

    %load the dopler profiler current readings (these are more reliable)
    load('adcpData');

    La = find(adcp_lat == lat);
    Lo = find(adcp_lon == lon);

    %keep only data at lat/lon of interest
    ADCP_loc = adcp_u(:,:,La,Lo);

    %combine the depth scale of the doppler and current data removing any
    %duplicates
    depthComb = [curDepth;double(adcp_depth)];
    depthComb = sort(depthComb);
    depthComb = unique(depthComb);

    %get the index of the current readings where the adcp readings being
    %and end
    iStartTime = find(abs(curTime-adcp_time(1))==min(abs(curTime-adcp_time(1))));
    iEndTime = find(abs(curTime-adcp_time(end))==min(abs(curTime-adcp_time(end))));

    %create array to hold combined current data
    ULocation = nan(length(curTime),length(depthComb));

    %for each depth of the current readings, add it to the combined data
    c = find(ismember(depthComb,curDepth));
    for i=1:length(c)
        ULocation(:,c(i)) = cur_loc(:,i);
    end

    %so long as the adcp data exists, add it to the combined data
    %(overwriting the current reading if there is one at the same depth and
    %time)
    c = find(ismember(depthComb,adcp_depth));
    for i=1:length(c)
        for n=1:(iEndTime-iStartTime+1)
            if ~isnan(ADCP_loc(n,i))
                ULocation(iStartTime+n-1,c(i)) = ADCP_loc(n,i);
            end
        end
    end

    Umax = nan(size(ULocation,1),1);
    Umean = nan(size(ULocation,1),1);
    Utotal = nan(size(ULocation,1),1);
    %use a cubic fit to interpolate to depths between upper and lower
    %reading (so lon as there is some data present)
    for j=1:size(ULocation,1)
        Udata = squeeze(ULocation(j,:));
        Ddata = depthComb;
        Ddata(isnan(Udata)) = [];
        Udata(isnan(Udata))=[];
        if length(Ddata)>1
            Dmin = find(depthComb==Ddata(1));
            Dmax = find(depthComb==Ddata(end));
            Uinterp = interp1(Ddata,Udata,depthComb(Dmin:Dmax),'pchip');
            ULocation(j,Dmin:Dmax)=Uinterp;
            Umax(j) = max(Uinterp);
            if Uinterp(1)==Umax(j) || Uinterp(end)==Umax(j)
                Umax(j) = nan;
                Utotal(j) = nan;
                Umean(j) = nan;
            else
                if Uinterp(end)>=80
                    Umean(j) = nan;
                else
                    Umean(j) = nansum(Uinterp(Uinterp>=80));
                end
                if Uinterp(1)<=10 && Uinterp(end)<=10
                    Utotal(j) = nansum(Uinterp);
                else
                    Utotal(j) = nan;
                end
            end
        end

    end
    
    Uclim = [];
    Um_clim = Umean;
    Um_clim(find(curTime==1997):find(curTime<1999,1,'last')) = nan;
    for i=1:12
        Uclim(i) = nanmean(Um_clim(i:12:end));
    end
    Uclim = Uclim-nanmean(Uclim);
    Uclim = repmat(Uclim, [1 40]);
    Uclim = Uclim';
    
    %calculate the temperature and velocity anomalies at each depth by
    %removing the climatology
    Umean_anom = Umean-Uclim(1:length(Umean));
    
    UMA = [UMA Umean_anom(:)];
    UTOT = [UTOT Umean(:)./Utotal(:)];
end

for i=1:length(UMA)-12*5+1
    runmean(i,:) = nanmean(UMA(i:i+5*12-1,:),1);
    timemean(i) = nanmean(curTime(i:i+5*12-1));
    for j=1:2
        if sum(isnan(UMA(i:i+5*12-1,j)))>12*2.5
            runmean(i,j)=NaN;
        end
    end
end

%Remove 97 98 El Nino/La Nina
UMA(find(curTime==1997):find(curTime<1999,1,'last'),:) = nan;
for j=1:length(LON)
    subplot('position',p(2+j,:));
    hold on;
    plot([1980 1999],nanmean(UMA(1:find(curTime==1999)-1,j))*5/100*[1 1],'b','linewidth',5);
    %plot([1990 2000],nanmean(UMA(find(curTime==1990):find(curTime==2000)-1,j))*[1 1]*5/100,'k','linewidth',3);
    plot([1999 2015],nanmean(UMA(find(curTime==1999):end,j))*5/100*[1 1],'r','linewidth',5);
    plot(timemean,runmean(:,j)*5/100,'k','linewidth',2);
    xlim([1980 2015]);
    title(strcat(num2str(LON(j)),' ^oE'),'FontSize',20);
    set(gca,'fontsize',20);
    ylabel({'Flow rate per', 'unit width (m^2/s)'});
    xlim([1980,2015]);
    set(gca,'XTickLabel',[],'XLabel',[]);
end

%add colorbar to plot
axes('Position', [0.05 0.04 0.9 0.901], 'Visible', 'off'); 
set(gca,'CLim',[0, 160]);
colormap('jet');
c=colorbar ('FontSize',20,'FontName','Myraid'); ylabel(c,'Zonal velocity (cm/s)') 


    