%% This script reproduces Figure S8 in Chen et al. (2025). 

% Note that the P-values displayed may vary slightly from the ones in the
% original figure because of the randomness of peak sampling in the null
% model simulations used 

% first load input_Figure_S8.mat in Matlab (version used
% originally was 2024a) to load the necessary input data. First make sure the
% Matlab path is set to the directory in which the .m and
% .mat files were downloaded. 

load input_Figure_S8.mat;

%choose scenario for peak randomization
%sc = 0 : no restrictions on random peak ages chosen
%sc = 1 : distance between randomly chosen peak ages sampled from fitted distribution of distances between observed peak ages

sc=0;

% set x- and y-axis limits for all panels
xlimits = [-130 0;-130 0;-130 0;-130 0;-130 0;-130 0];
ylimits = [-1.5 2;-1.5 2;-1.5 2;-1.5 2;-1.5 2;-1.5 2];

figure(1)

%plot geological events on all panels, minimum line thickness 0.2 Ma, events of longer duration have corresponding thickness
for i=1:1:6
    ax = subplot(3,2,i);

    % OAE1a (Selli Event)
    oae1a=[-120.5 ylimits(i,1) 1.157 ylimits(i,2)-ylimits(i,1)] ;
    rectangle('Position',oae1a,'FaceColor',[1 0.8 0.8],'EdgeColor','none')
    if (i==1)
        t=text(oae1a(1)+ 0.5*oae1a(3), 2+(2-(-1.5))./40, 'OAE1a','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==2)
        t=text(oae1a(1)+ 0.5*oae1a(3), 2+(2-(-0.5))./40, 'OAE1a','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==3)
        t=text(oae1a(1)+ 0.5*oae1a(3), 35+(35-(10))./40, 'OAE1a','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    
    
    % Aptian extinction
    aptian=[-116 ylimits(i,1) 1 ylimits(i,2)-ylimits(i,1)];
    rectangle('Position',aptian,'FaceColor',[1 0.8 0.8],'EdgeColor','none')
    if (i==1)
       t = text(aptian(1)+ 0.5*aptian(3), 2+(2-(-1.5))./40, 'Aptian ext.','Color',[1 0.6 0.6],'FontSize',10);
       t.Rotation = 90;
    end
    if (i==2)
       t = text(aptian(1)+ 0.5*aptian(3), 2+(2-(-0.5))./40, 'Aptian ext.','Color',[1 0.6 0.6],'FontSize',10);
       t.Rotation = 90;
    end
    if (i==3)
       t = text(aptian(1)+ 0.5*aptian(3), 35+(35-(10))./40, 'Aptian ext.','Color',[1 0.6 0.6],'FontSize',10);
       t.Rotation = 90;
    end

    % OAE1b (Jacob, Kilian, Paquier, and Leenhardt events)    
    oae1b=[-113.2-1.55 ylimits(i,1) 4.03 ylimits(i,2)-ylimits(i,1)] ;
    rectangle('Position',oae1b,'FaceColor',[1 0.8 0.8],'EdgeColor','none')
    if (i==1)
        t=text(oae1b(1)+ 0.5*oae1b(3)+2, 2+(2-(-1.5))./40, 'OAE1b','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==2)
        t=text(oae1b(1)+ 0.5*oae1b(3)+2, 2+(2-(-0.5))./40, 'OAE1b','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==3)
        t=text(oae1b(1)+ 0.5*oae1b(3)+2, 35+(35-(10))./40, 'OAE1b','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    
    
    % OAE1c (Amadeus event)
    oae1c=[-106 ylimits(i,1) 0.567 ylimits(i,2)-ylimits(i,1)];
    rectangle('Position',oae1c,'FaceColor',[1 0.8 0.8],'EdgeColor','none')
    if (i==1)
        t=text(oae1c(1)+ 0.5*oae1c(3), 2+(2-(-1.5))./40, 'OAE1c','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==2)
        t=text(oae1c(1)+ 0.5*oae1c(3), 2+(2-(-0.5))./40, 'OAE1c','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==3)
        t=text(oae1c(1)+ 0.5*oae1c(3), 35+(35-(10))./40, 'OAE1c','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    
    
    % OAE1d (Breistroffer event)
    oae1d=[-101 ylimits(i,1) 1 ylimits(i,2)-ylimits(i,1)];
    rectangle('Position',oae1d,'FaceColor',[1 0.8 0.8],'EdgeColor','none')
    if (i==1)
        t=text(oae1d(1)+ 0.5*oae1d(3), 2+(2-(-1.5))./40, 'OAE1d','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==2)
        t=text(oae1d(1)+ 0.5*oae1d(3), 2+(2-(-0.5))./40, 'OAE1d','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==3)
        t=text(oae1d(1)+ 0.5*oae1d(3), 35+(35-(10))./40, 'OAE1d','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end

    % 0AE2 = Cenomanian-Turonian boundary event, 
    oae2=[-93.9 ylimits(i,1) 0.9 ylimits(i,2)-ylimits(i,1)];
    rectangle('Position',oae2,'FaceColor',[1 0.8 0.8],'EdgeColor','none')
    if (i==1)
        t=text(oae2(1)+ 0.5*oae2(3), 2+(2-(-1.5))./40, 'OAE2','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==2)
        t=text(oae2(1)+ 0.5*oae2(3), 2+(2-(-0.5))./40, 'OAE2','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==3)
        t=text(oae2(1)+ 0.5*oae2(3), 35+(35-(10))./40, 'OAE2','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end

    % OAE3    
    oae3=[-87.3 ylimits(i,1) 2.7 ylimits(i,2)-ylimits(i,1)];
    rectangle('Position',oae3,'FaceColor',[1 0.8 0.8],'EdgeColor','none')
    if (i==1)
        t=text(oae3(1)+ 0.5*oae3(3), 2+(2-(-1.5))./40, 'OAE3','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==2)
        t=text(oae3(1)+ 0.5*oae3(3), 2+(2-(-0.5))./40, 'OAE3','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==3)
        t=text(oae3(1)+ 0.5*oae3(3), 35+(35-(10))./40, 'OAE3','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    
    % K-Pg extinction
    kpg=[-66.1 ylimits(i,1) 0.2 ylimits(i,2)-ylimits(i,1)];
    rectangle('Position',kpg,'FaceColor',[1 0.8 0.8],'EdgeColor','none')
    if (i==1)
        t=text(kpg(1)+ 0.5*kpg(3), 2+(2-(-1.5))./40, 'K-Pg','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==2)
        t=text(kpg(1)+ 0.5*kpg(3), 2+(2-(-0.5))./40, 'K-Pg','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==3)
        t=text(kpg(1)+ 0.5*kpg(3), 35+(35-(10))./40, 'K-Pg','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    
    % Paleocene-Eocene Thermal Maximum, PETM   
    petm=[-55.8 ylimits(i,1) 0.2 ylimits(i,2)-ylimits(i,1)];
    rectangle('Position',petm,'FaceColor',[1 0.8 0.8],'EdgeColor','none')
    if (i==1)
        t=text(petm(1)+ 0.5*petm(3), 2+(2-(-1.5))./40, 'PETM','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==2)
        t=text(petm(1)+ 0.5*petm(3), 2+(2-(-0.5))./40, 'PETM','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==3)
        t=text(petm(1)+ 0.5*petm(3), 35+(35-(10))./40, 'PETM','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    
    % Eocene-Oligocene transition (EOT), Eocene-Oligocene extinction event or Grande Coupure 
    eot=[-33.9 ylimits(i,1) 0.5 ylimits(i,2)-ylimits(i,1)];
    rectangle('Position',eot,'FaceColor',[1 0.8 0.8],'EdgeColor','none')
    if (i==1)
        t=text(eot(1)+ 0.5*eot(3), 2+(2-(-1.5))./40, 'EOT','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==2)
        t=text(eot(1)+ 0.5*eot(3), 2+(2-(-0.5))./40, 'EOT','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==3)
        t=text(eot(1)+ 0.5*eot(3), 35+(35-(10))./40, 'EOT','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end

    % Middle Miocene Climatic Transition (MMCT) or Middle Miocene Disruption   
    mmct=[-15 ylimits(i,1) 2 ylimits(i,2)-ylimits(i,1)];
    rectangle('Position',mmct,'FaceColor',[1 0.8 0.8],'EdgeColor','none');
    if (i==1)
        t=text(mmct(1)+ 0.5*mmct(3), 2+(2-(-1.5))./40, 'MMCT','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==2)
        t=text(mmct(1)+ 0.5*mmct(3), 2+(2-(-0.5))./40, 'MMCT','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    if (i==3)
        t=text(mmct(1)+ 0.5*mmct(3), 35+(35-(10))./40, 'MMCT','Color',[1 0.6 0.6],'FontSize',10);
        t.Rotation = 90;
    end
    
    ax.Layer = 'top';

end



%construction of panel A. Gray dots represent the relative residuals of a
%linear model in which the moving average of the WGD establishment rate
%(window size 4 Ma) is modeled as a function of time and the number of
%angiosperm lineages in the 466-species tree through time (Figure 1). The
%blue curve is a LOESS fit with span set to 10% of the total number of data
%points. Peaks in the blue curve that are recovered in the optimal
%exponential-Gaussian mixture model (panel B) are indicated with vertical
%black dashed lines and the associated numbers indicate their ages in Ma.
%The extinction events and OAEs displayed in pink are located significantly
%closer to the WGD peaks than expected under a null model of uniform random
%WGD peak placement across the past 115 Ma, as indicated by the P-value in
%the upper right corner of the panel.

l=1;

subplot(3,2,l);
hold on

% construction of relative residuals profile

%moving window analysis on WGD count data (WGD_Count contains the number of observed WGDs in bins [0,1] Ma, [1,2] Ma,... Window size 4, moving averages instead of moving sums by ./4
%moving window analysis on WGD count data (WGD_Count contains the number of observed WGDs in bins [0,1] Ma, [1,2] Ma,... Window size 4, moving averages instead of moving sums by ./4. 
%Note that first movsum element only incorporates 2 counts, not 4, and second element incorporates 3 counts instead of 4. Last movsum element also incorporates 3 counts only,
%but adjustment not needed as corresponding factor 2 or 4/3 corrections in numerator and denominator of mwWGD_rate_per_lineage cancel each other out. 
%In the same vein, ./4 in numerator and denominator actually also cancel out.
mwWGD_rate = movsum(WGD_Count(1:125),4)./4;
averageLineageCount = movsum(Total_Lineage_Count(1:125),4)./4;
mwWGD_rate_per_lineage = mwWGD_rate./averageLineageCount;
%first element of movsum captures range [-2 2] so average is 0, element 2 of movsum captures range [-1 3] with average 1, element 125 captures range [122 126] with average 124
mwx = [0:1:124];
mwWGD_rate_per_lineageflip = flip(mwWGD_rate_per_lineage(2:125));
mwxflip = -flip(mwx(2:125));
% Total_Lineage_Count is total lineage count at 0,1,2... Ma
Total_Lineage_Count_flip = flip(Total_Lineage_Count(1:125));
% leave out 0 Mya from Total_Lineage_Count_0123_flip, fit linear model to mwWGD_rate_per_lineageflip with independent variables time (mwxflip) and number of lineages (Total_Lineage_Count_flip)
mdl = fitlm([mwxflip',Total_Lineage_Count_flip(2:125)],mwWGD_rate_per_lineageflip);
% calculate relative residuals of model
residuals = mwWGD_rate_per_lineageflip - mdl.Fitted;
relativeresiduals = residuals./mdl.Fitted;

%plot panel A

xlabel( 'Ma', 'Interpreter', 'none' );
ylabel( 'Relative residuals', 'Interpreter', 'tex' );

plot(mwxflip',relativeresiduals,"Marker","o","MarkerSize",3,"MarkerFaceColor",[0.7 0.7 0.7],"MarkerEdgeColor",[0.7 0.7 0.7],"LineStyle","none");
hold on
%plot smoothed relative residuals profile
smoothedrelativeresiduals = smooth(mwxflip',relativeresiduals,0.1,'loess');
plot(mwxflip',smoothedrelativeresiduals,'Color','blue','LineWidth',1.5);
%find peaks in relative residuals profile, these will be plotted as vertical dashed black lines
[pks,locs] = findpeaks(smoothedrelativeresiduals);

x = [-130 0];
y = [0 0];
line(x,y,'Color','black','LineStyle','-')
hold off;

y = [-1.5 2];
for(i=1:1:size(pks,1))
  x = [mwxflip(locs(i)),mwxflip(locs(i))];
  line(x,y,'Color','black','LineStyle','--')
  text(mwxflip(locs(i)), smoothedrelativeresiduals(locs(i))+0.1, num2str(-mwxflip(locs(i)),'%.0f'));
end

xlim([-130 0]);
ylim([-1.5 2]);

%fit Gaussian distribution to distances between peaks in residuals data (to be used as distance distribution to be sampled from for positioning random peaks in scenario 1)
respeakdists = [];
%calculate distances between subsequent peaks
for i=1:1:size(pks,1)-1
    respeakdists(i) = -mwxflip(locs(i)) - (-mwxflip(locs(i+1)));
end
respeakdists = respeakdists';
pd = fitdist(respeakdists,'Normal');

%coupling distance analysis between relative residuals peaks and geological event peaks

%geological events in first 115 Ma
geopeaks =[14 33.65 55.7 66 85.95 93.45 101.5 105.7165 112.735]';
% use couplingDistance function in couplingDistance.m to calculate mean nearest distance between geological events and WGD peaks
cdist = couplingDistance(-mwxflip(locs(1:1:(size(pks,1)))),geopeaks);

%now do the same for 1,000,000 sets of randomized WGD establishment rate per lineage peaks (depends on randomization scenario chosen)
rndcdist=[];
for i=1:1:1000000
    rndpks=[];
    if (sc==0)
        for j=1:1:(size(pks,1))
            rndpks(j)=random('Uniform',0,115);
        end

    elseif (sc==1)
        %sample distances from distribution of residuals peak distances
        rndpkdists=[];
        for j=1:1:(size(respeakdists,1))
            rndpkdists(j)=random(pd);
        end
    
        %sample random starting point in interval 0-30 Ma and position random peaks
        rndpks(1)=random('Uniform',0,30);
        for j=2:1:(size(pks,1))
            rndpks(j)=rndpks(j-1) + rndpkdists(j-1);
        end
    end
    rndcdist(i)=couplingDistance(rndpks,geopeaks);
end


%approximate p-value for peak coupling p=(n+1)/(k+1) with n the number of times the randomized coupling distance is lower or equal to the observed coupling distance 
p1(l) = find(sort(rndcdist) >= cdist, 1)./1000001;

%plot p-value
t=text((xlimits(l,2) - (xlimits(l,2)-xlimits(l,1)).*0.25), (ylimits(l,2) - (ylimits(l,2)-ylimits(l,1)).*0.05), strcat('p=',num2str(p1(l),'%.2e')),'Color','black','FontSize',10);

xticks(-120:20:0);
xticklabels({'120','100','80','60','40','20','0'});
text(xlimits(l,1) - (xlimits(l,2)-xlimits(l,1)).*0.15 ,ylimits(l,2) + (ylimits(l,2)-ylimits(l,1)).*0.10,'A','FontSize',20)










%construction of panel B. Gray dots represent the relative residuals of a
%linear model in which the moving average of the WGD establishment rate
%(window size 4 Ma) is modeled as a function of time and the number of
%angiosperm species across time as obtained from Niklas et al. (1983) Nature
%https://doi.org/10.1038/303614a0. The blue curve is a LOESS fit with span set 
%to 10% of the total number of data points. Peaks in the blue curve that are 
%recovered in the optimal exponential-Gaussian mixture model (panel B) are 
%indicated with vertical black dashed lines and the associated numbers indicate 
%their ages in Ma. The extinction events and OAEs displayed in pink are located 
%significantly closer to the WGD peaks than expected under a null model of uniform
%random WGD peak placement across the past 115 Ma, as indicated by the P-value in
%the upper right corner of the panel.

l=3;

subplot(3,2,l);
hold on

% construction of relative residuals profile

%moving window analysis on WGD count data (WGD_Count contains the number of observed WGDs in bins [0,1] Ma, [1,2] Ma,... Window size 4, moving averages instead of moving sums by ./4. 
%Note that first movsum element only incorporates 2 counts, not 4, and second element incorporates 3 counts instead of 4. Last movsum element also incorporates 3 counts only,
%but adjustment not needed as corresponding factor 2 or 4/3 corrections in numerator and denominator of mwWGD_rate_per_lineage cancel each other out. 
%In the same vein, ./4 in numerator and denominator actually also cancel out.
%Niklas data starts from 3.5 Ma, while WGD count data starts at 0.5 Ma -> start at index 4 instead of 1
mwWGD_rate = movsum(WGD_Count(4:125),4)./4;
averageLineageCount = movsum(Niklas_Species_Count(1:122),4)./4;
mwWGD_rate_per_lineage = mwWGD_rate./averageLineageCount;
%first element of movsum captures range [1 5] so average is 3, element 2 of movsum captures range [2 6] with average 4, element 118 captures range [122 126] with average 124
mwx = [3:1:124];
mwWGD_rate_per_lineageflip = flip(mwWGD_rate_per_lineage(1:122));
mwxflip = -flip(mwx(1:122));
% Warning, moving averages at 3, 4, ... Ma may be 0.5 Ma off from lineage counts at 3.5, 4.5... Ma
Niklas_Species_Count_flip = flip(Niklas_Species_Count(1:122));
% fit linear model to mwWGD_rate_per_lineageflip with independent variables time (mwxflip) and number of lineages (Niklas_Species_Count_flip)
mdl = fitlm([mwxflip',Niklas_Species_Count_flip],mwWGD_rate_per_lineageflip);
% calculate relative residuals of model
residuals = mwWGD_rate_per_lineageflip - mdl.Fitted;
relativeresiduals = residuals./mdl.Fitted;

%plot panel B

xlabel( 'Ma', 'Interpreter', 'none' );
ylabel( 'Relative residuals', 'Interpreter', 'tex' );

plot(mwxflip',relativeresiduals,"Marker","o","MarkerSize",3,"MarkerFaceColor",[0.7 0.7 0.7],"MarkerEdgeColor",[0.7 0.7 0.7],"LineStyle","none");
hold on
%plot smoothed relative residuals profile
smoothedrelativeresiduals = smooth(mwxflip',relativeresiduals,0.1,'loess');
plot(mwxflip',smoothedrelativeresiduals,'Color','blue','LineWidth',1.5);
%find peaks in relative residuals profile, these will be plotted as vertical dashed black lines
[pks,locs] = findpeaks(smoothedrelativeresiduals);

x = [-130 0];
y = [0 0];
line(x,y,'Color','black','LineStyle','-')
hold off;

y = [-1.5 2];
for(i=1:1:size(pks,1))
  x = [mwxflip(locs(i)),mwxflip(locs(i))];
  line(x,y,'Color','black','LineStyle','--')
  text(mwxflip(locs(i)), smoothedrelativeresiduals(locs(i))+0.1, num2str(-mwxflip(locs(i)),'%.0f'));
end

xlim([-130 0]);
ylim([-1.5 2]);

%fit Gaussian distribution to distances between peaks in residuals data (to be used as distance distribution to be sampled from for positioning random peaks in scenario 1)
respeakdists = [];
%calculate distances between subsequent peaks
for i=1:1:size(pks,1)-1
    respeakdists(i) = -mwxflip(locs(i)) - (-mwxflip(locs(i+1)));
end
respeakdists = respeakdists';
pd = fitdist(respeakdists,'Normal');

%coupling distance analysis between relative residuals peaks and geological event peaks

%geological events in first 115 Ma
geopeaks =[14 33.65 55.7 66 85.95 93.45 101.5 105.7165 112.735]';
% use couplingDistance function in couplingDistance.m to calculate mean nearest distance between geological events and WGD peaks
cdist = couplingDistance(-mwxflip(locs(1:1:(size(pks,1)))),geopeaks);

%now do the same for 1,000,000 sets of randomized WGD establishment rate per lineage peaks (depends on randomization scenario chosen)
rndcdist=[];
for i=1:1:1000000
    rndpks=[];
    if (sc==0)
        for j=1:1:(size(pks,1))
            rndpks(j)=random('Uniform',0,115);
        end

    elseif (sc==1)
        %sample distances from distribution of residuals peak distances
        rndpkdists=[];
        for j=1:1:(size(respeakdists,1))
            rndpkdists(j)=random(pd);
        end
    
        %sample random starting point in interval 0-30 Ma and position random peaks
        rndpks(1)=random('Uniform',0,30);
        for j=2:1:(size(pks,1))
            rndpks(j)=rndpks(j-1) + rndpkdists(j-1);
        end
    end
    rndcdist(i)=couplingDistance(rndpks,geopeaks);
end


%approximate p-value for peak coupling p=(n+1)/(k+1) with n the number of times the randomized coupling distance is lower or equal to the observed coupling distance 
p1(l) = find(sort(rndcdist) >= cdist, 1)./1000001;

%plot p-value
t=text((xlimits(l,2) - (xlimits(l,2)-xlimits(l,1)).*0.25), (ylimits(l,2) - (ylimits(l,2)-ylimits(l,1)).*0.05), strcat('p=',num2str(p1(l),'%.2e')),'Color','black','FontSize',10);

xticks(-120:20:0);
xticklabels({'120','100','80','60','40','20','0'});
text(xlimits(l,1) - (xlimits(l,2)-xlimits(l,1)).*0.15 ,ylimits(l,2) + (ylimits(l,2)-ylimits(l,1)).*0.10,'B','FontSize',20)









%construction of panel C. Gray dots represent the relative residuals of a
%linear model in which the moving average of the WGD establishment rate
%(window size 4 Ma) is modeled as a function of time and the number of
%angiosperm families across time as obtained from Niklas et al. (1983) Nature
%https://doi.org/10.1038/303614a0. The blue curve is a LOESS fit with span set 
%to 10% of the total number of data points. Peaks in the blue curve that are 
%recovered in the optimal exponential-Gaussian mixture model (panel B) are 
%indicated with vertical black dashed lines and the associated numbers indicate 
%their ages in Ma. The extinction events and OAEs displayed in pink are located 
%significantly closer to the WGD peaks than expected under a null model of uniform
%random WGD peak placement across the past 115 Ma, as indicated by the P-value in
%the upper right corner of the panel.

l=5;

subplot(3,2,l);
hold on

% construction of relative residuals profile

%moving window analysis on WGD count data (WGD_Count contains the number of observed WGDs in bins [0,1] Ma, [1,2] Ma,... Window size 4, moving averages instead of moving sums by ./4. 
%Note that first movsum element only incorporates 2 counts, not 4, and second element incorporates 3 counts instead of 4. Last movsum element also incorporates 3 counts only,
%but adjustment not needed as corresponding factor 2 or 4/3 corrections in numerator and denominator of mwWGD_rate_per_lineage cancel each other out. 
%In the same vein, ./4 in numerator and denominator actually also cancel out.
%Niklas data starts from 3.5 Ma, while WGD count data starts at 0.5 Ma -> start at index 4 instead of 1
mwWGD_rate = movsum(WGD_Count(4:125),4)./4;
averageLineageCount = movsum(Niklas_Family_Count(1:122),4)./4;
mwWGD_rate_per_lineage = mwWGD_rate./averageLineageCount;
%first element of movsum captures range [1 5] so average is 3, element 2 of movsum captures range [2 6] with average 4, element 118 captures range [122 126] with average 124
mwx = [3:1:124];
mwWGD_rate_per_lineageflip = flip(mwWGD_rate_per_lineage(1:122));
mwxflip = -flip(mwx(1:122));
% Warning, moving averages at 3, 4, ... Ma may be 0.5 Ma off from lineage counts at 3.5, 4.5... Ma
Niklas_Family_Count_flip = flip(Niklas_Family_Count(1:122));
% fit linear model to mwWGD_rate_per_lineageflip with independent variables time (mwxflip) and number of lineages (Niklas_Family_Count_flip)
mdl = fitlm([mwxflip',Niklas_Family_Count_flip],mwWGD_rate_per_lineageflip);
% calculate relative residuals of model
residuals = mwWGD_rate_per_lineageflip - mdl.Fitted;
relativeresiduals = residuals./mdl.Fitted;



%plot panel C

xlabel( 'Ma', 'Interpreter', 'none' );
ylabel( 'Relative residuals', 'Interpreter', 'tex' );

plot(mwxflip',relativeresiduals,"Marker","o","MarkerSize",3,"MarkerFaceColor",[0.7 0.7 0.7],"MarkerEdgeColor",[0.7 0.7 0.7],"LineStyle","none");
hold on
%plot smoothed relative residuals profile
smoothedrelativeresiduals = smooth(mwxflip',relativeresiduals,0.1,'loess');
plot(mwxflip',smoothedrelativeresiduals,'Color','blue','LineWidth',1.5);
%find peaks in relative residuals profile, these will be plotted as vertical dashed black lines
[pks,locs] = findpeaks(smoothedrelativeresiduals);

x = [-130 0];
y = [0 0];
line(x,y,'Color','black','LineStyle','-')
hold off;

y = [-1.5 2];
for(i=1:1:size(pks,1))
  x = [mwxflip(locs(i)),mwxflip(locs(i))];
  line(x,y,'Color','black','LineStyle','--')
  text(mwxflip(locs(i)), smoothedrelativeresiduals(locs(i))+0.1, num2str(-mwxflip(locs(i)),'%.0f'));
end

xlim([-130 0]);
ylim([-1.5 2]);

%fit Gaussian distribution to distances between peaks in residuals data (to be used as distance distribution to be sampled from for positioning random peaks in scenario 1)
respeakdists = [];
%calculate distances between subsequent peaks
for i=1:1:size(pks,1)-1
    respeakdists(i) = -mwxflip(locs(i)) - (-mwxflip(locs(i+1)));
end
respeakdists = respeakdists';
pd = fitdist(respeakdists,'Normal');

%coupling distance analysis between relative residuals peaks and geological event peaks

%geological events in first 115 Ma
geopeaks =[14 33.65 55.7 66 85.95 93.45 101.5 105.7165 112.735]';
% use couplingDistance function in couplingDistance.m to calculate mean nearest distance between geological events and WGD peaks
cdist = couplingDistance(-mwxflip(locs(1:1:(size(pks,1)))),geopeaks);

%now do the same for 1,000,000 sets of randomized WGD establishment rate per lineage peaks (depends on randomization scenario chosen)
rndcdist=[];
for i=1:1:1000000
    rndpks=[];
    if (sc==0)
        for j=1:1:(size(pks,1))
            rndpks(j)=random('Uniform',0,115);
        end

    elseif (sc==1)
        %sample distances from distribution of residuals peak distances
        rndpkdists=[];
        for j=1:1:(size(respeakdists,1))
            rndpkdists(j)=random(pd);
        end
    
        %sample random starting point in interval 0-30 Ma and position random peaks
        rndpks(1)=random('Uniform',0,30);
        for j=2:1:(size(pks,1))
            rndpks(j)=rndpks(j-1) + rndpkdists(j-1);
        end
    end
    rndcdist(i)=couplingDistance(rndpks,geopeaks);
end


%approximate p-value for peak coupling p=(n+1)/(k+1) with n the number of times the randomized coupling distance is lower or equal to the observed coupling distance 
p1(l) = find(sort(rndcdist) >= cdist, 1)./1000001;

%plot p-value
t=text((xlimits(l,2) - (xlimits(l,2)-xlimits(l,1)).*0.25), (ylimits(l,2) - (ylimits(l,2)-ylimits(l,1)).*0.05), strcat('p=',num2str(p1(l),'%.2e')),'Color','black','FontSize',10);

xticks(-120:20:0);
xticklabels({'120','100','80','60','40','20','0'});
text(xlimits(l,1) - (xlimits(l,2)-xlimits(l,1)).*0.15 ,ylimits(l,2) + (ylimits(l,2)-ylimits(l,1)).*0.10,'C','FontSize',20)












%construction of panel D. Gray dots represent the relative residuals of a
%linear model in which the moving average of the WGD establishment rate
%(window size 4 Ma) is modeled as a function of time and the number of
%angiosperm species across time as obtained from Zuntini et al. (2024) Nature
%https://doi.org/10.1038/s41586-024-07324-0. The blue curve is a LOESS fit with span set 
%to 10% of the total number of data points. Peaks in the blue curve that are 
%recovered in the optimal exponential-Gaussian mixture model (panel B) are 
%indicated with vertical black dashed lines and the associated numbers indicate 
%their ages in Ma. The extinction events and OAEs displayed in pink are located 
%significantly closer to the WGD peaks than expected under a null model of uniform
%random WGD peak placement across the past 115 Ma, as indicated by the P-value in
%the upper right corner of the panel.

l=2;

subplot(3,2,l);
hold on

% construction of relative residuals profile

%moving window analysis on WGD count data (WGD_Count contains the number of observed WGDs in bins [0,1] Ma, [1,2] Ma,... Window size 4, moving averages instead of moving sums by ./4. 
%Note that first movsum element only incorporates 2 counts, not 4, and second element incorporates 3 counts instead of 4. Last movsum element also incorporates 3 counts only,
%but adjustment not needed as corresponding factor 2 or 4/3 corrections in numerator and denominator of mwWGD_rate_per_lineage cancel each other out. 
%In the same vein, ./4 in numerator and denominator actually also cancel out.
mwWGD_rate = movsum(WGD_Count(1:125),4)./4;
averageLineageCount = movsum(Zuntini_Species_Count_young(1:125),4)./4;
mwWGD_rate_per_lineage = mwWGD_rate./averageLineageCount;
%first element of movsum captures range [-2 2] so average is 0, element 2 of movsum captures range [-1 3] with average 1, element 125 captures range [122 126] with average 124
mwx = [0:1:124];
mwWGD_rate_per_lineageflip = flip(mwWGD_rate_per_lineage(2:125));
mwxflip = -flip(mwx(2:125));
% Warning, moving averages at 0, 1, ... Ma may be 0.5 Ma off from lineage counts at 0.5, 1.5... Ma
Zuntini_Species_Count_young_flip = flip(Zuntini_Species_Count_young(1:125));
% leave out 0 Mya from Zuntini_Species_Count_young_flip, fit linear model to mwWGD_rate_per_lineageflip with independent variables time (mwxflip) and number of lineages (Zuntini_Species_Count_young_flip)
mdl = fitlm([mwxflip',Zuntini_Species_Count_young_flip(2:125)],mwWGD_rate_per_lineageflip);
% calculate relative residuals of model
residuals = mwWGD_rate_per_lineageflip - mdl.Fitted;
relativeresiduals = residuals./mdl.Fitted;

%plot panel D

xlabel( 'Ma', 'Interpreter', 'none' );
ylabel( 'Relative residuals', 'Interpreter', 'tex' );

plot(mwxflip',relativeresiduals,"Marker","o","MarkerSize",3,"MarkerFaceColor",[0.7 0.7 0.7],"MarkerEdgeColor",[0.7 0.7 0.7],"LineStyle","none");
hold on
%plot smoothed relative residuals profile
smoothedrelativeresiduals = smooth(mwxflip',relativeresiduals,0.1,'loess');
plot(mwxflip',smoothedrelativeresiduals,'Color','blue','LineWidth',1.5);
%find peaks in relative residuals profile, these will be plotted as vertical dashed black lines
[pks,locs] = findpeaks(smoothedrelativeresiduals);

x = [-130 0];
y = [0 0];
line(x,y,'Color','black','LineStyle','-')
hold off;

y = [-1.5 2];
for(i=1:1:size(pks,1))
  x = [mwxflip(locs(i)),mwxflip(locs(i))];
  line(x,y,'Color','black','LineStyle','--')
  text(mwxflip(locs(i)), smoothedrelativeresiduals(locs(i))+0.1, num2str(-mwxflip(locs(i)),'%.0f'));
end

xlim([-130 0]);
ylim([-1.5 2]);

%fit Gaussian distribution to distances between peaks in residuals data (to be used as distance distribution to be sampled from for positioning random peaks in scenario 1)
respeakdists = [];
%calculate distances between subsequent peaks
for i=1:1:size(pks,1)-1
    respeakdists(i) = -mwxflip(locs(i)) - (-mwxflip(locs(i+1)));
end
respeakdists = respeakdists';
pd = fitdist(respeakdists,'Normal');

%coupling distance analysis between relative residuals peaks and geological event peaks

%geological events in first 115 Ma
geopeaks =[14 33.65 55.7 66 85.95 93.45 101.5 105.7165 112.735]';
% use couplingDistance function in couplingDistance.m to calculate mean nearest distance between geological events and WGD peaks
cdist = couplingDistance(-mwxflip(locs(1:1:(size(pks,1)))),geopeaks);

%now do the same for 1,000,000 sets of randomized WGD establishment rate per lineage peaks (depends on randomization scenario chosen)
rndcdist=[];
for i=1:1:1000000
    rndpks=[];
    if (sc==0)
        for j=1:1:(size(pks,1))
            rndpks(j)=random('Uniform',0,115);
        end

    elseif (sc==1)
        %sample distances from distribution of residuals peak distances
        rndpkdists=[];
        for j=1:1:(size(respeakdists,1))
            rndpkdists(j)=random(pd);
        end
    
        %sample random starting point in interval 0-30 Ma and position random peaks
        rndpks(1)=random('Uniform',0,30);
        for j=2:1:(size(pks,1))
            rndpks(j)=rndpks(j-1) + rndpkdists(j-1);
        end
    end
    rndcdist(i)=couplingDistance(rndpks,geopeaks);
end


%approximate p-value for peak coupling p=(n+1)/(k+1) with n the number of times the randomized coupling distance is lower or equal to the observed coupling distance 
p1(l) = find(sort(rndcdist) >= cdist, 1)./1000001;

%plot p-value
t=text((xlimits(l,2) - (xlimits(l,2)-xlimits(l,1)).*0.25), (ylimits(l,2) - (ylimits(l,2)-ylimits(l,1)).*0.05), strcat('p=',num2str(p1(l),'%.2e')),'Color','black','FontSize',10);

xticks(-120:20:0);
xticklabels({'120','100','80','60','40','20','0'});
text(xlimits(l,1) - (xlimits(l,2)-xlimits(l,1)).*0.15 ,ylimits(l,2) + (ylimits(l,2)-ylimits(l,1)).*0.10,'D','FontSize',20)










%construction of panel E. Gray dots represent the relative residuals of a
%linear model in which the moving average of the WGD establishment rate
%(window size 4 Ma) is modeled as a function of time and the number of
%angiosperm genera across time as obtained from Zuntini et al. (2024) Nature
%https://doi.org/10.1038/s41586-024-07324-0. The blue curve is a LOESS fit with span set 
%to 10% of the total number of data points. Peaks in the blue curve that are 
%recovered in the optimal exponential-Gaussian mixture model (panel B) are 
%indicated with vertical black dashed lines and the associated numbers indicate 
%their ages in Ma. The extinction events and OAEs displayed in pink are located 
%significantly closer to the WGD peaks than expected under a null model of uniform
%random WGD peak placement across the past 115 Ma, as indicated by the P-value in
%the upper right corner of the panel.

l=4;

subplot(3,2,l);
hold on

% construction of relative residuals profile

%moving window analysis on WGD count data (WGD_Count contains the number of observed WGDs in bins [0,1] Ma, [1,2] Ma,... Window size 4, moving averages instead of moving sums by ./4. 
%Note that first movsum element only incorporates 2 counts, not 4, and second element incorporates 3 counts instead of 4. Last movsum element also incorporates 3 counts only,
%but adjustment not needed as corresponding factor 2 or 4/3 corrections in numerator and denominator of mwWGD_rate_per_lineage cancel each other out. 
%In the same vein, ./4 in numerator and denominator actually also cancel out.
mwWGD_rate = movsum(WGD_Count(1:125),4)./4;
averageLineageCount = movsum(Zuntini_Genus_Count_young(1:125),4)./4;
mwWGD_rate_per_lineage = mwWGD_rate./averageLineageCount;
%first element of movsum captures range [-2 2] so average is 0, element 2 of movsum captures range [-1 3] with average 1, element 125 captures range [122 126] with average 124
mwx = [0:1:124];
mwWGD_rate_per_lineageflip = flip(mwWGD_rate_per_lineage(2:125));
mwxflip = -flip(mwx(2:125));
% Warning, moving averages at 0, 1, ... Ma may be 0.5 Ma off from lineage counts at 0.5, 1.5... Ma
Zuntini_Genus_Count_young_flip = flip(Zuntini_Genus_Count_young(1:125));
% leave out 0 Mya from Zuntini_Genus_Count_young_flip, fit linear model to mwWGD_rate_per_lineageflip with independent variables time (mwxflip) and number of lineages (Zuntini_Genus_Count_young_flip)
mdl = fitlm([mwxflip',Zuntini_Genus_Count_young_flip(2:125)],mwWGD_rate_per_lineageflip);
% calculate relative residuals of model
residuals = mwWGD_rate_per_lineageflip - mdl.Fitted;
relativeresiduals = residuals./mdl.Fitted;

%plot panel E

xlabel( 'Ma', 'Interpreter', 'none' );
ylabel( 'Relative residuals', 'Interpreter', 'tex' );

plot(mwxflip',relativeresiduals,"Marker","o","MarkerSize",3,"MarkerFaceColor",[0.7 0.7 0.7],"MarkerEdgeColor",[0.7 0.7 0.7],"LineStyle","none");
hold on
%plot smoothed relative residuals profile
smoothedrelativeresiduals = smooth(mwxflip',relativeresiduals,0.1,'loess');
plot(mwxflip',smoothedrelativeresiduals,'Color','blue','LineWidth',1.5);
%find peaks in relative residuals profile, these will be plotted as vertical dashed black lines
[pks,locs] = findpeaks(smoothedrelativeresiduals);

x = [-130 0];
y = [0 0];
line(x,y,'Color','black','LineStyle','-')
hold off;

y = [-1.5 2];
for(i=1:1:size(pks,1))
  x = [mwxflip(locs(i)),mwxflip(locs(i))];
  line(x,y,'Color','black','LineStyle','--')
  text(mwxflip(locs(i)), smoothedrelativeresiduals(locs(i))+0.1, num2str(-mwxflip(locs(i)),'%.0f'));
end

xlim([-130 0]);
ylim([-1.5 2]);

%fit Gaussian distribution to distances between peaks in residuals data (to be used as distance distribution to be sampled from for positioning random peaks in scenario 1)
respeakdists = [];
%calculate distances between subsequent peaks
for i=1:1:size(pks,1)-1
    respeakdists(i) = -mwxflip(locs(i)) - (-mwxflip(locs(i+1)));
end
respeakdists = respeakdists';
pd = fitdist(respeakdists,'Normal');

%coupling distance analysis between relative residuals peaks and geological event peaks

%geological events in first 115 Ma
geopeaks =[14 33.65 55.7 66 85.95 93.45 101.5 105.7165 112.735]';
% use couplingDistance function in couplingDistance.m to calculate mean nearest distance between geological events and WGD peaks
cdist = couplingDistance(-mwxflip(locs(1:1:(size(pks,1)))),geopeaks);

%now do the same for 1,000,000 sets of randomized WGD establishment rate per lineage peaks (depends on randomization scenario chosen)
rndcdist=[];
for i=1:1:1000000
    rndpks=[];
    if (sc==0)
        for j=1:1:(size(pks,1))
            rndpks(j)=random('Uniform',0,115);
        end

    elseif (sc==1)
        %sample distances from distribution of residuals peak distances
        rndpkdists=[];
        for j=1:1:(size(respeakdists,1))
            rndpkdists(j)=random(pd);
        end
    
        %sample random starting point in interval 0-30 Ma and position random peaks
        rndpks(1)=random('Uniform',0,30);
        for j=2:1:(size(pks,1))
            rndpks(j)=rndpks(j-1) + rndpkdists(j-1);
        end
    end
    rndcdist(i)=couplingDistance(rndpks,geopeaks);
end


%approximate p-value for peak coupling p=(n+1)/(k+1) with n the number of times the randomized coupling distance is lower or equal to the observed coupling distance 
p1(l) = find(sort(rndcdist) >= cdist, 1)./1000001;

%plot p-value
t=text((xlimits(l,2) - (xlimits(l,2)-xlimits(l,1)).*0.25), (ylimits(l,2) - (ylimits(l,2)-ylimits(l,1)).*0.05), strcat('p=',num2str(p1(l),'%.2e')),'Color','black','FontSize',10);

xticks(-120:20:0);
xticklabels({'120','100','80','60','40','20','0'});
text(xlimits(l,1) - (xlimits(l,2)-xlimits(l,1)).*0.15 ,ylimits(l,2) + (ylimits(l,2)-ylimits(l,1)).*0.10,'E','FontSize',20)












%construction of panel F. Gray dots represent the relative residuals of a
%linear model in which the moving average of the WGD establishment rate
%(window size 4 Ma) is modeled as a function of time and the number of
%angiosperm families across time as obtained from Zuntini et al. (2024) Nature
%https://doi.org/10.1038/s41586-024-07324-0. The blue curve is a LOESS fit with span set 
%to 10% of the total number of data points. Peaks in the blue curve that are 
%recovered in the optimal exponential-Gaussian mixture model (panel B) are 
%indicated with vertical black dashed lines and the associated numbers indicate 
%their ages in Ma. The extinction events and OAEs displayed in pink are located 
%significantly closer to the WGD peaks than expected under a null model of uniform
%random WGD peak placement across the past 115 Ma, as indicated by the P-value in
%the upper right corner of the panel.

l=6;

subplot(3,2,l);
hold on

% construction of relative residuals profile

%moving window analysis on WGD count data (WGD_Count contains the number of observed WGDs in bins [0,1] Ma, [1,2] Ma,... Window size 4, moving averages instead of moving sums by ./4. 
%Note that first movsum element only incorporates 2 counts, not 4, and second element incorporates 3 counts instead of 4. Last movsum element also incorporates 3 counts only,
%but adjustment not needed as corresponding factor 2 or 4/3 corrections in numerator and denominator of mwWGD_rate_per_lineage cancel each other out. 
%In the same vein, ./4 in numerator and denominator actually also cancel out.

mwWGD_rate = movsum(WGD_Count(1:125),4)./4;
averageLineageCount = movsum(Zuntini_Family_Count_young(1:125),4)./4;
mwWGD_rate_per_lineage = mwWGD_rate./averageLineageCount;
%first element of movsum captures range [-2 2] so average is 0, element 2 of movsum captures range [-1 3] with average 1, element 125 captures range [122 126] with average 124
mwx = [0:1:124];
mwWGD_rate_per_lineageflip = flip(mwWGD_rate_per_lineage(2:125));
mwxflip = -flip(mwx(2:125));
% Warning, moving averages at 0, 1, ... Ma may be 0.5 Ma off from lineage counts at 0.5, 1.5... Ma
Zuntini_Family_Count_young_flip = flip(Zuntini_Family_Count_young(1:125));
% leave out 0 Mya from Zuntini_Family_Count_young_flip, fit linear model to mwWGD_rate_per_lineageflip with independent variables time (mwxflip) and number of lineages (Zuntini_Family_Count_young_flip)
mdl = fitlm([mwxflip',Zuntini_Family_Count_young_flip(2:125)],mwWGD_rate_per_lineageflip);
% calculate relative residuals of model
residuals = mwWGD_rate_per_lineageflip - mdl.Fitted;
relativeresiduals = residuals./mdl.Fitted;

%plot panel F

xlabel( 'Ma', 'Interpreter', 'none' );
ylabel( 'Relative residuals', 'Interpreter', 'tex' );

plot(mwxflip',relativeresiduals,"Marker","o","MarkerSize",3,"MarkerFaceColor",[0.7 0.7 0.7],"MarkerEdgeColor",[0.7 0.7 0.7],"LineStyle","none");
hold on
%plot smoothed relative residuals profile
smoothedrelativeresiduals = smooth(mwxflip',relativeresiduals,0.1,'loess');
plot(mwxflip',smoothedrelativeresiduals,'Color','blue','LineWidth',1.5);
%find peaks in relative residuals profile, these will be plotted as vertical dashed black lines
[pks,locs] = findpeaks(smoothedrelativeresiduals);

x = [-130 0];
y = [0 0];
line(x,y,'Color','black','LineStyle','-')
hold off;

y = [-1.5 2];
for(i=1:1:size(pks,1))
  x = [mwxflip(locs(i)),mwxflip(locs(i))];
  line(x,y,'Color','black','LineStyle','--')
  text(mwxflip(locs(i)), smoothedrelativeresiduals(locs(i))+0.1, num2str(-mwxflip(locs(i)),'%.0f'));
end

xlim([-130 0]);
ylim([-1.5 2]);

%fit Gaussian distribution to distances between peaks in residuals data (to be used as distance distribution to be sampled from for positioning random peaks in scenario 1)
respeakdists = [];
% calculate distances between subsequent peaks
% leave out first peak because outside of [0,115] Ma range
for i=2:1:size(pks,1)-1
    respeakdists(i) = -mwxflip(locs(i)) - (-mwxflip(locs(i+1)));
end
respeakdists = respeakdists';
pd = fitdist(respeakdists,'Normal');

%coupling distance analysis between relative residuals peaks and geological event peaks

%geological events in first 115 Ma
geopeaks =[14 33.65 55.7 66 85.95 93.45 101.5 105.7165 112.735]';
% use couplingDistance function in couplingDistance.m to calculate mean nearest distance between geological events and WGD peaks
% leave out first peak because outside of [0,115] Ma range
cdist = couplingDistance(-mwxflip(locs(2:1:(size(pks,1)))),geopeaks);

% now do the same for 1,000,000 sets of randomized WGD establishment rate per lineage peaks (depends on randomization scenario chosen)
% use size(pks,1)-1 because first peak outside of [0,115] Ma range
rndcdist=[];
for i=1:1:1000000
    rndpks=[];
    if (sc==0)
        for j=1:1:(size(pks,1)-1)
            rndpks(j)=random('Uniform',0,115);
        end

    elseif (sc==1)
        %sample distances from distribution of residuals peak distances
        rndpkdists=[];
        for j=1:1:(size(respeakdists,1))
            rndpkdists(j)=random(pd);
        end
    
        %sample random starting point in interval 0-30 Ma and position random peaks
        rndpks(1)=random('Uniform',0,30);
        for j=2:1:(size(pks,1))
            rndpks(j)=rndpks(j-1) + rndpkdists(j-1);
        end
    end
    rndcdist(i)=couplingDistance(rndpks,geopeaks);
end


%approximate p-value for peak coupling p=(n+1)/(k+1) with n the number of times the randomized coupling distance is lower or equal to the observed coupling distance 
p1(l) = find(sort(rndcdist) >= cdist, 1)./1000001;

%plot p-value
t=text((xlimits(l,2) - (xlimits(l,2)-xlimits(l,1)).*0.25), (ylimits(l,2) - (ylimits(l,2)-ylimits(l,1)).*0.05), strcat('p=',num2str(p1(l),'%.2e')),'Color','black','FontSize',10);

xticks(-120:20:0);
xticklabels({'120','100','80','60','40','20','0'});
text(xlimits(l,1) - (xlimits(l,2)-xlimits(l,1)).*0.15 ,ylimits(l,2) + (ylimits(l,2)-ylimits(l,1)).*0.10,'F','FontSize',20)





















