%% This script reproduces Figure S12 in Chen et al. (2025). 

% Note that the P-values displayed may vary slightly from the ones in the
% original figure because of the randomness of peak sampling in the null
% model simulations used 

% first load input_Figure_3_and_Figure_S12.mat in Matlab (version used
% originally was 2024a) to load the necessary input data. First make sure the
% Matlab path is set to the directory in which the .m and
% .mat files were downloaded

load input_Figure_3_and_Figure_S12.mat;

%choose scenario for peak randomization
%sc = 0 : no restrictions on random peak ages chosen
%sc = 1 : distance between randomly chosen peak ages sampled from fitted distribution of distances between observed peak ages

sc=1;

% set x- and y-axis limits for all panels
xlimits = [-130 0;-130 0;-130 0;-130 0;-130 0;-130 0;-130 0;-130 0; -130 0];
ylimits = [-1.5 2 ; -0.5 2 ; 10 35 ; -0.00001 0.00007; -4 4 ; 5 50 ; -0.00001 0.00007; -6 6; -100 1900];

%p-values for coupling analyses between extrema with peaks in panel A, only used for the panels in which multiple extrema sets are tested (see below)
pp = NaN(9,200);
%minimum p-values for coupling analyses between extrema in the different panels with peaks in panel A
p1 = [];

figure(1)

%plot geological events on all panels, minimum line thickness 0.2 Ma, events of longer duration have corresponding thickness
for i=1:1:9
    ax = subplot(3,3,i);

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

subplot(3,3,l);
hold on

% construction of relative residuals profile

%moving window analysis on WGD count data (WGD_Count contains the number of observed WGDs in bins [0,1] Ma, [1,2] Ma,... Window size 4, moving averages instead of moving sums by ./4
mwWGD_rate = movsum(WGD_Count(1:125),4)./4;
averageLineageCount = movsum(Total_Lineage_Count(1:125),4)./4;
mwWGD_rate_per_lineage = mwWGD_rate./averageLineageCount;
%first element of movsum captures range [-2 2] so average is 0, element 2 of movsum captures range [-1 3] with average 1, element 125 captures range [122 126] with average 124
mwx = [0:1:124];
mwWGD_rate_per_lineageflip = flip(mwWGD_rate_per_lineage(2:125));
mwxflip = -flip(mwx(2:125));
% Total_Lineage_Count_0123 is total lineage count at 0,1,2... Ma
Total_Lineage_Count_0123_flip = flip(Total_Lineage_Count_0123(1:125));
% leave out 0 Mya from Total_Lineage_Count_0123_flip, fit linear model to mwWGD_rate_per_lineageflip with independent variables time (mwxflip) and number of lineages (Total_Lineage_Count_0123_flip)
mdl = fitlm([mwxflip',Total_Lineage_Count_0123_flip(2:125)],mwWGD_rate_per_lineageflip);
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
%find peaks in relative residuals profile, these will be plotted in all panels as vertical dashed black lines
[pks,locs] = findpeaks(smoothedrelativeresiduals);

x = [-130 0];
y = [0 0];
line(x,y,'Color','black','LineStyle','-')
hold off;

% leave out peak at 8 Ma because not significant in exponential-Gaussian mixture model analysis in panel B
y = [-1.5 2];
for(i=1:1:size(pks,1)-1)
  x = [mwxflip(locs(i)),mwxflip(locs(i))];
  line(x,y,'Color','black','LineStyle','--')
  text(mwxflip(locs(i)), smoothedrelativeresiduals(locs(i))+0.1, num2str(-mwxflip(locs(i)),'%.0f'));
end

xlim([-130 0]);
ylim([-1.5 2]);

%fit Gaussian distribution to distances between peaks in residuals data (to be used as distance distribution to be sampled from for positioning random peaks in scenario 1)
respeakdists = [];
%calculate distances between subsequent peaks (leaving out last insignificant peak at 8 Ma)
for i=1:1:size(pks,1)-2
    respeakdists(i) = -mwxflip(locs(i)) - (-mwxflip(locs(i+1)));
end
respeakdists = respeakdists';
pd = fitdist(respeakdists,'Normal');

%coupling distance analysis between relative residuals peaks and geological event peaks

%geological events in first 115 Ma
geopeaks =[14 33.65 55.7 66 85.95 93.45 101.5 105.7165 112.735]';
% use couplingDistance function in couplingDistance.m to calculate mean nearest distance between geological events and significant peaks in panel A
cdist = couplingDistance(-mwxflip(locs(1:1:(size(pks,1)-1))),geopeaks);

%now do the same for 1,000,000 sets of randomized WGD establishment rate per lineage peaks (depends on randomization scenario chosen)
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
        for j=2:1:(size(pks,1)-1)
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












% construction of panel B. Optimal exponential-Gaussian mixture model for
% WGD establishment rate data. Gray dots are kernel density estimates (KDE)
% of the WGD establishment rate (KDE bandwidth = 3.1). The blue line is the
% optimal exponential-Gaussian mixture fit (adjusted R^2 = 0.9986) to the
% KDE data for the past 115 Ma, containing an exponential background and 9
% Gaussian components of which the means are indicated by vertical red
% bars.

l=4;

subplot(3,3,l);
hold on

% calculate kernel density-estimated WGD establishment rate per lineage profile from consensus mean WGD ages in ConsensusMean

WGD_ages=ConsensusMean';

mirror1 = -WGD_ages;
mirror2 = Timemya1(size(Timemya1,1)) + (Timemya1(size(Timemya1,1))-WGD_ages);
WGD_ages = cat(2,mirror2,WGD_ages);
WGD_ages = cat(2,WGD_ages,mirror1);

[kd,xkd,bw] = kde(WGD_ages, Bandwidth=3.1, EvaluationPoints=0.5:1:Timemya1(size(Timemya1,1)));

WGD_rate_per_lineage = kd'./Total_Lineage_Count;
WGD_rate_per_lineage = WGD_rate_per_lineage(1:124,1);
xkd = (xkd(1:124))';

xkdflip = -flip(xkd);
WGD_rate_per_lineageflip = flip(WGD_rate_per_lineage);

%use curveFitter toolbox to fit exponential-Gaussian mixture of the following form to WGD_rate_per_lineage
%a*exp(b*x) -a + a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + a3*exp(-((x-b3)/c3)^2) + a4*exp(-((x-b4)/c4)^2) + a5*exp(-((x-b5)/c5)^2) + a6*exp(-((x-b6)/c6)^2)+ a7*exp(-((x-b7)/c7)^2) + a8*exp(-((x-b8)/c8)^2) + a9*exp(-((x-b9)/c9)^2)
%The optimal exponential-Gaussian mixture model fit was exported as the fit
%code below, and the associated parameters were also saved in fittedmodel_WGDmean_9

[xData, yData] = prepareCurveData( xkdflip, WGD_rate_per_lineageflip );

% Set up fittype and options.
ft = fittype( 'a*exp(b*x)-a+a1*exp(-((x-b1)/c1)^2)+a2*exp(-((x-b2)/c2)^2)+a3*exp(-((x-b3)/c3)^2)+a4*exp(-((x-b4)/c4)^2)+a5*exp(-((x-b5)/c5)^2)+a6*exp(-((x-b6)/c6)^2)+a7*exp(-((x-b7)/c7)^2)+a8*exp(-((x-b8)/c8)^2)+a9*exp(-((x-b9)/c9)^2)', 'independent', 'x', 'dependent', 'y' );
excludedPoints = (xData < -115) | (xData > 0);
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMaxChange = 0.01;
opts.Display = 'Off';
opts.Lower = [-0.1 0 0 0 0 0 0 0 0 0 0 -130 -130 -130 -130 -130 -130 -130 -130 -130 0 0 0 0 0 0 0 1 1];
opts.MaxFunEvals = 1000000000;
opts.MaxIter = 1000000000;
opts.StartPoint = [-0.001 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.01 -105 -100 -85 -75 -65 -55 -35 -25 -15 5 5 5 5 5 5 5 5 5];
opts.TolFun = 1e-18;
opts.TolX = 1e-18;
opts.Upper = [0 Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf 0 0 0 0 0 0 0 0 0 Inf Inf Inf Inf Inf Inf Inf Inf Inf];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
h = plot( fitresult, '-b', xData, yData, '.k', 'predobs', 0.99 );
h(1).Marker="o";
h(1).MarkerSize=3;
h(1).MarkerFaceColor=[0.7 0.7 0.7];
h(1).MarkerEdgeColor=[0.7 0.7 0.7];
h(2).LineWidth=1.5;

% plot peaks of panel A, leave out peak at 8 Ma because not significant in exponential-Gaussian mixture model analysis in panel B
y = [-0.00001 0.00007];
for(i=1:1:size(pks,1)-1)
  x = [mwxflip(locs(i)),mwxflip(locs(i))];
  line(x,y,'Color','black','LineStyle','--')
end


legend off;
% Label axes
xlabel( 'Ma', 'Interpreter', 'none' );
ylabel( 'WGD establishment rate (Ma^{-1})', 'Interpreter', 'tex' );

%plot little red bars at the locations of the Gaussian means in the optimal exponential-Gaussian mixture model
sml = (ylimits(l,2)-ylimits(l,1))./40;

y = [WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b1))+sml WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b1))+4.*sml];
x = [fittedmodel_WGDmean_9.b1 fittedmodel_WGDmean_9.b1];
line(x,y,'Color','red','LineStyle','-','LineWidth',1.5)

y = [WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b2))+sml WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b2))+4.*sml];
x = [fittedmodel_WGDmean_9.b2 fittedmodel_WGDmean_9.b2];
line(x,y,'Color','red','LineStyle','-','LineWidth',1.5)

y = [WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b3))+sml WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b3))+4.*sml];
x = [fittedmodel_WGDmean_9.b3 fittedmodel_WGDmean_9.b3];
line(x,y,'Color','red','LineStyle','-','LineWidth',1.5)

y = [WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b4))+sml WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b4))+4.*sml];
x = [fittedmodel_WGDmean_9.b4 fittedmodel_WGDmean_9.b4];
line(x,y,'Color','red','LineStyle','-','LineWidth',1.5)

y = [WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b5))+sml WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b5))+4.*sml];
x = [fittedmodel_WGDmean_9.b5 fittedmodel_WGDmean_9.b5];
line(x,y,'Color','red','LineStyle','-','LineWidth',1.5)

y = [WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b6))+sml WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b6))+4.*sml];
x = [fittedmodel_WGDmean_9.b6 fittedmodel_WGDmean_9.b6];
line(x,y,'Color','red','LineStyle','-','LineWidth',1.5)

y = [WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b7))+sml WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b7))+4.*sml];
x = [fittedmodel_WGDmean_9.b7 fittedmodel_WGDmean_9.b7];
line(x,y,'Color','red','LineStyle','-','LineWidth',1.5)

y = [WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b8))+sml WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b8))+4.*sml];
x = [fittedmodel_WGDmean_9.b8 fittedmodel_WGDmean_9.b8];
line(x,y,'Color','red','LineStyle','-','LineWidth',1.5)

y = [WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b9))+sml WGD_rate_per_lineage(-round(fittedmodel_WGDmean_9.b9))+4.*sml];
x = [fittedmodel_WGDmean_9.b9 fittedmodel_WGDmean_9.b9];
line(x,y,'Color','red','LineStyle','-','LineWidth',1.5)


xlim([-130 0]);
ylim([-0.00001 0.00007]);
xticks(-120:20:0);
xticklabels({'120','100','80','60','40','20','0'});
ax=gca;
ax.TickDir = 'out';
ax.Layer='top';










% construction of panel C. Comparison of the observed KDE profile for WGD
% establishment rate (black line) with a null model assuming constant
% clade-specific WGD establishment rates (see Methods). The blue line
% represents the average WGD establishment rate profile obtained across
% 10,000 null model simulations, while the dark and light blue bands depict
% the 50% and 95% confidence intervals, respectively. Six of the nine peaks
% identified in panels A and B surpass the 95% confidence threshold.

l=7;

subplot(3,3,l);
hold on

% A WGD establishment rate per lineage profile was calculated for 10,000
% simulated WGD age profiles under a constant rate model (different rates
% for different major clades), simulated ages were imported as
% 'Simulated_WGD_dates_constant_rate'

WGD_rate_per_lineage_sim=[];

for i=1:1:size(Simulated_WGD_dates_constant_rate,1)
    % remove NaNs from parsed Simulated_WGD_dates_constant_rate matrix (not all simulations had the same number of WGDS, rows were extended with NaNs until the maximum number of WGDs obtained in any of the simulations)
    cleanSimWGDages=[];
    cleanSimWGDages=Simulated_WGD_dates_constant_rate(i,:);
    ind4=find(isnan(cleanSimWGDages));
    cleanSimWGDages(ind4)=[];

    % kernel density estimation and construction of simulated WGD establishment rate per lineage profile for simulation replicate i 
    mirror1 = -cleanSimWGDages;
    mirror2 = Timemya1(size(Timemya1,1)) + (Timemya1(size(Timemya1,1))-cleanSimWGDages);
    cleanSimWGDages = cat(2,mirror2,cleanSimWGDages);
    cleanSimWGDages = cat(2,cleanSimWGDages,mirror1);
    [kd_3,xkd_3,bw_3] = kde(cleanSimWGDages, Bandwidth=3.1, EvaluationPoints=0.5:1:Timemya1(size(Timemya1,1)));
    tmp = kd_3'./Total_Lineage_Count;
    WGD_rate_per_lineage_sim(:,i) = tmp(1:124,1);
end
xkd_3 = xkd_3(1,1:124);

% calculate min, max, mean and 50% and 95% confidence intervals on simulated WGD establishment rate per lineage profiles 
WGD_rate_per_lineage_simmin = min(WGD_rate_per_lineage_sim,[],2);
WGD_rate_per_lineage_simmax = max(WGD_rate_per_lineage_sim,[],2);
WGD_rate_per_lineage_simmean = mean(WGD_rate_per_lineage_sim,2,"omitmissing");

P = prctile(WGD_rate_per_lineage_sim, 97.5, 2);
P2 = prctile(WGD_rate_per_lineage_sim, 2.5, 2);

P3 = prctile(WGD_rate_per_lineage_sim, 75, 2);
P4 = prctile(WGD_rate_per_lineage_sim, 25, 2);

x2 = [-xkd_3, fliplr(-xkd_3)];
inBetween1 = [P', fliplr(P2')];
fill(x2, inBetween1, [0.9 0.9 1]);

inBetween2 = [P3', fliplr(P4')];
fill(x2, inBetween2, [0.7 0.7 1]);

plot(-xkd_3,WGD_rate_per_lineage_simmean,'Marker','none','LineStyle','-','Color','blue','LineWidth',1.5);
plot(-xkd_3,WGD_rate_per_lineage,'Marker','none','LineStyle','-','Color','black','LineWidth',1.5);

% plot peaks of panel A, leave out peak at 8 Ma because not significant in exponential-Gaussian mixture model analysis in panel B
y = [-0.00001 0.00007];

for(i=1:1:size(pks,1)-1)
  x = [mwxflip(locs(i)),mwxflip(locs(i))];
  line(x,y,'Color','black','LineStyle','--')
end

xlabel( 'Ma', 'Interpreter', 'none' );
ylabel( 'WGD establishment rate (Ma^{-1})', 'Interpreter', 'tex' );

xlim([-130 0]);
ylim([-0.00001 0.00007]);
xticks(-120:20:0);
xticklabels({'120','100','80','60','40','20','0'});











% construction of panel D. Genus-level angiosperm extinction rate profiles
% for the past 130 Ma. Black lines are extinction rate profiles obtained
% from the angiosperm fossil record in the Paleobiology database using
% three methods of sampling standardization and four macroevolutionary rate
% models in divDyn (v0.8.3). The blue line is the average estimated
% angiosperm extinction rate across the 12 models used, with dark blue and
% light blue bands depicting the interquartile range and total range of the
% values. Five extinction rate peaks are indicated with vertical red bars.
% These angiosperm extinction rate peaks are located significantly closer
% to a WGD peak in panel A (dashed vertical black lines) than expected
% under a null model of uniform random WGD peak placement across the past
% 115 Ma, as indicated by the P-value in the upper right corner of the
% panel.


l=2;

subplot(3,3,l);
hold on

age_ang_extinctions = estimatesqsPCstages.xcoordinate(4:size(estimatesqsPCstages.xcoordinate)-1) ;
ang_extinction_rate(:,1)=estimatecr2f3stages.ycoordinate(4:size(estimatesqsPCstages.xcoordinate)-1);
ang_extinction_rate(:,2)=estimatecrC3tstages.ycoordinate(4:size(estimatesqsPCstages.xcoordinate)-1);
ang_extinction_rate(:,3)=estimatecrGFstages.ycoordinate(4:size(estimatesqsPCstages.xcoordinate)-1);
ang_extinction_rate(:,4)=estimatecrPCstages.ycoordinate(4:size(estimatesqsPCstages.xcoordinate)-1);
ang_extinction_rate(:,5)=estimateraw2f3stages.ycoordinate(4:size(estimatesqsPCstages.xcoordinate)-1);
ang_extinction_rate(:,6)=estimaterawC3tstages.ycoordinate(4:size(estimatesqsPCstages.xcoordinate)-1);
ang_extinction_rate(:,7)=estimaterawGFstages.ycoordinate(4:size(estimatesqsPCstages.xcoordinate)-1);
ang_extinction_rate(:,8)=estimaterawPCstages.ycoordinate(4:size(estimatesqsPCstages.xcoordinate)-1);
ang_extinction_rate(:,9)=estimatesqs2f3stages.ycoordinate(4:size(estimatesqsPCstages.xcoordinate)-1);
ang_extinction_rate(:,10)=estimatesqsC3tstages.ycoordinate(4:size(estimatesqsPCstages.xcoordinate)-1);
ang_extinction_rate(:,11)=estimatesqsGFstages.ycoordinate(4:size(estimatesqsPCstages.xcoordinate)-1);
ang_extinction_rate(:,12)=estimatesqsPCstages.ycoordinate(4:size(estimatesqsPCstages.xcoordinate)-1);
ang_extinction_ratemin = min(ang_extinction_rate,[],2);
ang_extinction_ratemax = max(ang_extinction_rate,[],2);
ang_extinction_ratemean = mean(ang_extinction_rate,2,"omitmissing");

[r,Q]= iqr(ang_extinction_rate,2);

x2 = [-age_ang_extinctions', fliplr(-age_ang_extinctions')];
inBetween1 = [ang_extinction_ratemax', fliplr(ang_extinction_ratemin')];
fill(x2, inBetween1, [0.9 0.9 1]);

inBetween2 = [Q(:,1)', fliplr(Q(:,2)')];
fill(x2, inBetween2, [0.7 0.7 1]);

for(i=1:1:size(ang_extinction_rate,2))
    plot(-age_ang_extinctions,ang_extinction_rate(:,i),'Marker','none','LineStyle','-','Color','black');
end

plot(-age_ang_extinctions,ang_extinction_ratemean,'Marker','none','LineStyle','-','Color','blue','LineWidth',1.5);

% plot peaks of panel A, leave out peak at 8 Ma because not significant in exponential-Gaussian mixture model analysis in panel B
y = [-0.5 2];

for(i=1:1:size(pks,1)-1)
  x = [mwxflip(locs(i)),mwxflip(locs(i))];
  line(x,y,'Color','black','LineStyle','--')
end

xlabel( 'Ma', 'Interpreter', 'none' );
ylabel( 'Angiosperm extinction rate (Ma^{-1})', 'Interpreter', 'tex' );

%find angiosperm extinction rate peaks
[angpeaks, anglocs, angw, angp] = findpeaks(ang_extinction_ratemean);

%plot little red bars at the locations of the angiosperm extinction rate peaks, leave out first peak because at left edge of data
sml = (ylimits(l,2)-ylimits(l,1))./40;
for(i=2:1:size(angpeaks,1))
    x = [-age_ang_extinctions(anglocs(i)),-age_ang_extinctions(anglocs(i))];
    y = [(ang_extinction_ratemax(anglocs(i))+sml), (ang_extinction_ratemax(anglocs(i))+4.*sml)];
    line(x,y,'Color','red','LineStyle','-','LineWidth',1.5)
end

%coupling distance analysis between relative residuals peaks and angiosperm extinction rate peaks

% use couplingDistance function in couplingDistance.m to calculate mean nearest distance between angiosperm extinction rate peaks and significant peaks in panel A
cdist = couplingDistance(-mwxflip(locs(1:1:(size(pks,1)-1))),age_ang_extinctions(anglocs(2:1:size(angpeaks,1))));

%now do the same for 1,000,000 sets of randomized WGD establishment rate per lineage peaks (depends on randomization scenario chosen)
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
        for j=2:1:(size(pks,1)-1)
            rndpks(j)=rndpks(j-1) + rndpkdists(j-1);
        end
    end

    rndcdist(i)=couplingDistance(rndpks,age_ang_extinctions(anglocs(2:1:size(angpeaks,1))));

end

%approximate p-value for peak coupling p=(n+1)/(k+1) with n the number of times the randomized coupling distance is lower or equal to the observed coupling distance 
p1(l) = find(sort(rndcdist) >= cdist, 1)./1000001;

%plot p-value
t=text((xlimits(l,2) - (xlimits(l,2)-xlimits(l,1)).*0.25), (ylimits(l,2) - (ylimits(l,2)-ylimits(l,1)).*0.05), strcat('p=',num2str(p1(l),'%.2e')),'Color','black','FontSize',10);


xlim([-130 0]);
ylim([-0.5 2]);
xticks(-120:20:0);
xticklabels({'120','100','80','60','40','20','0'});









% construction of panel E. Benthic foraminiferal d13C isotopic signature
% over the past 115 Ma, data from Friedrich et al. (2012)
% https://doi.org/10.1130/G32701.1 . Gray dots are d13C measurements, the
% blue line is a smoothing spline with smoothing parameter 0.05.

l=5;

subplot(3,3,l);
hold on
gca.XLim = [-130 0];

d13C_age = Age_frie12;
d13C = d13C_frie12;

plot(-d13C_age,d13C,"Marker","o","MarkerSize",2,"MarkerFaceColor",[0.7 0.7 0.7],"MarkerEdgeColor",[0.7 0.7 0.7],"LineStyle","none");

%sort on age and remove NaN values
[d13C_ageb,ind]=sort(-d13C_age);
d13Cb=d13C(ind);
ind2 = find(isnan(d13C_ageb));
d13C_ageb(ind2)=[];
d13Cb(ind2)=[];
ind3 = find(isnan(d13Cb));
d13C_ageb(ind3)=[];
d13Cb(ind3)=[];

%smoothed d13C profile
fitspline = fit(d13C_ageb,d13Cb,'smoothingspline','SmoothingParam',0.05);
splinex = -115:0.1:0;
smoothedd13C = fitspline(splinex);
plot(splinex ,smoothedd13C,'Color','blue', 'LineWidth',1.5);

% plot peaks of panel A, leave out peak at 8 Ma because not significant in exponential-Gaussian mixture model analysis in panel B
y = [-4 4];
for(i=1:1:size(pks,1)-1)
  x = [mwxflip(locs(i)),mwxflip(locs(i))];
  line(x,y,'Color','black','LineStyle','--')
end

%find maxima in smoothed d13C profile
[peaks13C, locs13C, w13C, p13C] = findpeaks(smoothedd13C);
%find minima in smoothed d13C profile
[peaksmin13C, locsmin13C, wmin13C, pmin13C] = findpeaks(-smoothedd13C);

%combine extrema and sort on prominence p13C, pmin13C
[combined13Cp,II] = sort([p13C' pmin13C']',1) ;
combined13Clocs = [locs13C' locsmin13C']' ;
combined13Clocs = combined13Clocs(II);
combined13Cpeaks = [peaks13C' peaksmin13C']' ;
combined13Cpeaks = combined13Cpeaks(II);
% extremum type = 1 for maximum, -1 for minimum
combined13Cextrtype = [ones(size(p13C')) -ones(size(pmin13C'))]' ;
combined13Cextrtype = combined13Cextrtype(II);

%limit to [0 115] Ma
JJ=[];
for k=1:1:size(combined13Cpeaks,1)
    if(-splinex(combined13Clocs(k))<115)
        JJ = [JJ k];
    end
end

combined13Cpeaks = combined13Cpeaks(JJ);
combined13Clocs = combined13Clocs(JJ);
combined13Cp = combined13Cp(JJ);
combined13Cextrtype = combined13Cextrtype(JJ);


%coupling distance analysis between relative residuals peaks (panel A) and
%d13C extrema, for different prominence thresholds. Sets consisting of the
%top-1 up to the top-18 extrema were tested by comparing the observed
%average nearest distance to a WGD establishement rate peak for any given
%set of extrema to a null distribution of average nearest distances
%obtained by uniformly sampling 9 WGD peaks across the past 115 Ma a
%million times (scenario 0) or to a null distribution of average nearest
%distances obtained by sampling the age of the first peak randomly in the
%interval [0 Ma, 30 Ma] and positioning 8 additional WGD peaks randomly at
%a distance from the previous peak sampled from a Gaussian distribution fit
%to the observed WGD establishment rate inter-peak distances (scenario 1,
%10^6 samples were used to construct the null distribution).


for k=size(combined13Cpeaks,1):-1:(size(combined13Cpeaks,1)-2*(size(pks,1)-1)+1)
    k
    
    % use couplingDistance function in couplingDistance.m to calculate mean nearest distance between selected d13C extrema k:1:size(combined13Cpeaks,1) and significant peaks in panel A
    cdist = couplingDistance(-mwxflip(locs(1:1:(size(pks,1)-1))),-splinex(combined13Clocs(k:1:size(combined13Cpeaks,1))));
    
    %now do the same for 1,000,000 sets of randomized WGD establishment rate per lineage peaks (depends on randomization scenario chosen)
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
            for j=2:1:(size(pks,1)-1)
                rndpks(j)=rndpks(j-1) + rndpkdists(j-1);
            end
        end

        rndcdist(i)=couplingDistance(rndpks,-splinex(combined13Clocs(k:1:size(combined13Cpeaks,1))));
    end
    
    %approximate p-value for peak coupling p=(n+1)/(k+1) with n the number of times the randomized coupling distance is lower or equal to the observed coupling distance 
    ppp = find(sort(rndcdist) >= cdist, 1)./1000001;
    %store in p-value matrix across panels (l) and extremum sets (k)
    if(~isempty(ppp))
        pp(l,k) = ppp;
    else
        pp(l,k) = 1;
    end
end

%find set of extrema with minimum p-value    
ind1=size(combined13Cpeaks,1)-2*(size(pks,1)-1)+1;
if ind1<=0
    ind1=1;
end
ind2=size(combined13Cpeaks,1);
[Ma,I3a] = min(pp(l,ind1:ind2));
if(~isempty(I3a))
    p1(l)=pp(l,ind1+I3a-1);
else
    p1(l)= NaN; 
end

%plot minimum p-value   
t=text((xlimits(l,2) - (xlimits(l,2)-xlimits(l,1)).*0.25), (ylimits(l,2) - (ylimits(l,2)-ylimits(l,1)).*0.05), strcat('p=',num2str(p1(l),'%.2e')),'Color','black','FontSize',10);

%plot little red (maximum) or purple (minimum) bars for the associated extrema
sml = (ylimits(l,2)-ylimits(l,1))./40;

for(i=ind1+I3a-1:1:size(combined13Cpeaks,1))
    if combined13Cextrtype(i) == 1
        x = [splinex(combined13Clocs(i)),splinex(combined13Clocs(i))];
        y = [(smoothedd13C(combined13Clocs(i))+sml),(smoothedd13C(combined13Clocs(i))+4.*sml)];
        line(x,y,'Color','red','LineStyle','-','LineWidth',1.5)
    else
        x = [splinex(combined13Clocs(i)),splinex(combined13Clocs(i))];
        y = [(smoothedd13C(combined13Clocs(i))-sml),(smoothedd13C(combined13Clocs(i))-4.*sml)];
        line(x,y,'Color',[1 0 1],'LineStyle','-','LineWidth',1.5)
    end
end

xlabel( 'Ma', 'Interpreter', 'none' );
ylabel( '\delta^{13}C (‰)', 'Interpreter', 'tex' );


xlim([-130 0]);
ylim([-4 4]);
xticks(-120:20:0);
xticklabels({'120','100','80','60','40','20','0'});









% construction of panel F. Benthic foraminiferal d18O isotopic signature
% over the past 115 Ma, data from Friedrich et al. (2012)
% https://doi.org/10.1130/G32701.1 . Gray dots are d18O measurements, the
% blue line is a smoothing spline with smoothing parameter 0.05.

l=8;

subplot(3,3,l);
hold on

d18O_age = Age_frie12;
d18O = d18O_frie12;

plot(-d18O_age,d18O,"Marker","o","MarkerSize",2,"MarkerFaceColor",[0.7 0.7 0.7],"MarkerEdgeColor",[0.7 0.7 0.7],"LineStyle","none");

%sort on age and remove NaN values
[d18O_ageb,ind]=sort(-d18O_age);
d18Ob=d18O(ind);
ind2 = find(isnan(d18O_ageb));
d18O_ageb(ind2)=[];
d18Ob(ind2)=[];
ind3 = find(isnan(d18Ob));
d18O_ageb(ind3)=[];
d18Ob(ind3)=[];

%smoothed d18O profile
fitspline = fit(d18O_ageb,d18Ob,'smoothingspline','SmoothingParam',0.05);
splinex = -115:0.1:0;
smoothedd18O = fitspline(splinex);

plot(splinex ,smoothedd18O,'Color','blue', 'LineWidth',1.5);

% plot peaks of panel A, leave out peak at 8 Ma because not significant in exponential-Gaussian mixture model analysis in panel B
y = [-6 6];
for(i=1:1:size(pks,1)-1)
  x = [mwxflip(locs(i)),mwxflip(locs(i))];
  line(x,y,'Color','black','LineStyle','--')
end

%find maxima in smoothed d18O profile
[peaks18O, locs18O, w18O, p18O] = findpeaks(smoothedd18O);
%find minima in smoothed d18O profile
[peaksmin18O, locsmin18O, wmin18O, pmin18O] = findpeaks(-smoothedd18O);

%combine extrema and sort on prominence p18O, pmin18O
[combined18Op,II] = sort([p18O' pmin18O']',1) ;
combined18Olocs = [locs18O' locsmin18O']' ;
combined18Olocs = combined18Olocs(II);
combined18Opeaks = [peaks18O' peaksmin18O']' ;
combined18Opeaks = combined18Opeaks(II);
% extremum type = 1 for maximum, -1 for minimum
combined18Oextrtype = [ones(size(p18O')) -ones(size(pmin18O'))]' ;
combined18Oextrtype = combined18Oextrtype(II);

%limit to [0 115] Ma
JJ=[];
for k=1:1:size(combined18Opeaks,1)
    if(-splinex(combined18Olocs(k))<115)
        JJ = [JJ k];
    end
end

combined18Opeaks = combined18Opeaks(JJ);
combined18Olocs = combined18Olocs(JJ);
combined18Op = combined18Op(JJ);
combined18Oextrtype = combined18Oextrtype(JJ);

%coupling distance analysis between relative residuals peaks (panel A) and
%d18O extrema, for different prominence thresholds. Sets consisting of the
%top-1 up to the top-18 extrema were tested by comparing the observed
%average nearest distance to a WGD establishement rate peak for any given
%set of extrema to a null distribution of average nearest distances
%obtained by uniformly sampling 9 WGD peaks across the past 115 Ma a
%million times (scenario 0) or to a null distribution of average nearest
%distances obtained by sampling the age of the first peak randomly in the
%interval [0 Ma, 30 Ma] and positioning 8 additional WGD peaks randomly at
%a distance from the previous peak sampled from a Gaussian distribution fit
%to the observed WGD establishment rate inter-peak distances (scenario 1,
%10^6 samples were used to construct the null distribution).


for k=size(combined18Opeaks,1):-1:(size(combined18Opeaks,1)-2*(size(pks,1)-1)+1)
    k
    % use couplingDistance function in couplingDistance.m to calculate mean nearest distance between selected d18O extrema k:1:size(combined18Opeaks,1) and significant peaks in panel A
    cdist = couplingDistance(-mwxflip(locs(1:1:(size(pks,1)-1))),-splinex(combined18Olocs(k:1:size(combined18Opeaks,1))));

    %now do the same for 1,000,000 sets of randomized WGD establishment rate per lineage peaks (depends on randomization scenario chosen)
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
            for j=2:1:(size(pks,1)-1)
                rndpks(j)=rndpks(j-1) + rndpkdists(j-1);
            end
        end

        rndcdist(i)=couplingDistance(rndpks,-splinex(combined18Olocs(k:1:size(combined18Opeaks,1))));
    end
    
    %approximate p-value for peak coupling p=(n+1)/(k+1) with n the number of times the randomized coupling distance is lower or equal to the observed coupling distance 
    ppp = find(sort(rndcdist) >= cdist, 1)./1000001;  
    %store in p-value matrix across panels (l) and extremum sets (k)
    if(~isempty(ppp))
        pp(l,k) = ppp;
    else
        pp(l,k) = 1;
    end
end

%find set of extrema with minimum p-value 
ind1=size(combined18Opeaks,1)-2*(size(pks,1)-1)+1;
if ind1<=0
    ind1=1;
end
ind2=size(combined18Opeaks,1);
[Ma,I3a] = min(pp(l,ind1:ind2));
if(~isempty(I3a))
    p1(l)=pp(l,ind1+I3a-1);
else
    p1(l)= NaN; 
end

%plot minimum p-value 
t=text((xlimits(l,2) - (xlimits(l,2)-xlimits(l,1)).*0.25), (ylimits(l,2) - (ylimits(l,2)-ylimits(l,1)).*0.05), strcat('p=',num2str(p1(l),'%.2e')),'Color','black','FontSize',10);

%plot little red (maximum) or purple (minimum) bars for the associated extrema
sml = (ylimits(l,2)-ylimits(l,1))./40;

for(i=ind1+I3a-1:1:size(combined18Opeaks,1))
    if combined18Oextrtype(i) == 1
        x = [splinex(combined18Olocs(i)),splinex(combined18Olocs(i))];
        y = [(smoothedd18O(combined18Olocs(i))+sml),(smoothedd18O(combined18Olocs(i))+4.*sml)];
        line(x,y,'Color','red','LineStyle','-','LineWidth',1.5)
    else
        x = [splinex(combined18Olocs(i)),splinex(combined18Olocs(i))];
        y = [(smoothedd18O(combined18Olocs(i))-sml),(smoothedd18O(combined18Olocs(i))-4.*sml)];
        line(x,y,'Color',[1 0 1],'LineStyle','-','LineWidth',1.5)
    end
end


xlabel( 'Ma', 'Interpreter', 'none' );
ylabel( '\delta^{18}O (‰)', 'Interpreter', 'tex' );


xlim([-130 0]);
ylim([-6 6]);
xticks(-120:20:0);
xticklabels({'120','100','80','60','40','20','0'});











% construction of panel G: Low-latitude sea surface temperature in the past
% 130 Ma, data from Song et al. (2019)
% https://doi.org/10.1007/s12583-018-1002-2 . Gray dots represent mean sea
% surface temperature estimates, and the blue curve is a LOESS fit with the
% span set to 2% of the total number of data points.

l=3;

subplot(3,3,l);
hold on

plot(-seasurftemp_AgeMa,seasurftemp_MeanTemperatureC,"Marker","o","MarkerSize",3,"MarkerFaceColor",[0.7 0.7 0.7],"MarkerEdgeColor",[0.7 0.7 0.7],"LineStyle","none");

%smoothed sea surface temperature profile
smoothedseasurftemp = smooth(seasurftemp_AgeMa,seasurftemp_MeanTemperatureC,0.02,'loess');
plot(-seasurftemp_AgeMa,smoothedseasurftemp,'Marker','none','LineStyle','-','Color','blue','LineWidth',1.5);

% plot peaks of panel A, leave out peak at 8 Ma because not significant in exponential-Gaussian mixture model analysis in panel B
y = [10 35];
for(i=1:1:size(pks,1)-1)
  x = [mwxflip(locs(i)),mwxflip(locs(i))];
  line(x,y,'Color','black','LineStyle','--')
end

%find maxima in smoothed sea surface T profile
[seaTpeaks, seaTlocs, seaTw, seaTp] = findpeaks(smoothedseasurftemp);
%find mimima in smoothed sea surface T profile
[seaTpeaksmin, seaTlocsmin, seaTwmin, seaTpmin] = findpeaks(-smoothedseasurftemp);

%combine extrema and sort on prominence seaTp, seaTpmin
[combinedseaTp,II] = sort([seaTp' seaTpmin']',1) ;
combinedseaTlocs = [seaTlocs' seaTlocsmin']' ;
combinedseaTlocs = combinedseaTlocs(II);
combinedseaTpeaks = [seaTpeaks' seaTpeaksmin']' ;
combinedseaTpeaks = combinedseaTpeaks(II);
% extremum type = 1 for maximum, -1 for minimum
combinedseaTextrtype = [ones(size(seaTp')) -ones(size(seaTpmin'))]' ;
combinedseaTextrtype = combinedseaTextrtype(II);

%limit to [0 115] Ma
JJ=[];
for k=1:1:size(combinedseaTpeaks,1)
    if(seasurftemp_AgeMa(combinedseaTlocs(k))<115)
        JJ = [JJ k];
    end
end

combinedseaTpeaks = combinedseaTpeaks(JJ);
combinedseaTlocs = combinedseaTlocs(JJ);
combinedseaTp = combinedseaTp(JJ);
combinedseaTextrtype = combinedseaTextrtype(JJ);

%coupling distance analysis between relative residuals peaks (panel A) and
%sea surface T extrema, for different prominence thresholds. Sets consisting of the
%top-1 up to the top-18 extrema were tested by comparing the observed
%average nearest distance to a WGD establishement rate peak for any given
%set of extrema to a null distribution of average nearest distances
%obtained by uniformly sampling 9 WGD peaks across the past 115 Ma a
%million times (scenario 0) or to a null distribution of average nearest
%distances obtained by sampling the age of the first peak randomly in the
%interval [0 Ma, 30 Ma] and positioning 8 additional WGD peaks randomly at
%a distance from the previous peak sampled from a Gaussian distribution fit
%to the observed WGD establishment rate inter-peak distances (scenario 1,
%10^6 samples were used to construct the null distribution).


for k=size(combinedseaTpeaks,1):-1:(size(combinedseaTpeaks,1)-2*(size(pks,1)-1)+1)
    k
    % use couplingDistance function in couplingDistance.m to calculate mean nearest distance between selected sea surface T extrema k:1:size(combinedseaTpeaks,1) and significant peaks in panel A
    cdist = couplingDistance(-mwxflip(locs(1:1:(size(pks,1)-1))),seasurftemp_AgeMa(combinedseaTlocs(k:1:size(combinedseaTpeaks,1))));

    %now do the same for 1,000,000 sets of randomized WGD establishment rate per lineage peaks (depends on randomization scenario chosen)
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
            for j=2:1:(size(pks,1)-1)
                rndpks(j)=rndpks(j-1) + rndpkdists(j-1);
            end
        end

        rndcdist(i)=couplingDistance(rndpks,seasurftemp_AgeMa(combinedseaTlocs(k:1:size(combinedseaTpeaks,1))));
    end
    
    %approximate p-value for peak coupling p=(n+1)/(k+1) with n the number of times the randomized coupling distance is lower or equal to the observed coupling distance 
    ppp = find(sort(rndcdist) >= cdist, 1)./1000001;
    %store in p-value matrix across panels (l) and extremum sets (k)
    if(~isempty(ppp))
        pp(l,k) = ppp;
    else
        pp(l,k) = 1;
    end
end

%find set of extrema with minimum p-value 
ind1=size(combinedseaTpeaks,1)-2*(size(pks,1)-1)+1;
if ind1<=0
    ind1=1;
end
ind2=size(combinedseaTpeaks,1);
[Ma,I3a] = min(pp(l,ind1:ind2));
if(~isempty(I3a))
    p1(l)=pp(l,ind1+I3a-1);
else
    p1(l)= NaN; 
end

%plot minimum p-value   
t=text((xlimits(l,2) - (xlimits(l,2)-xlimits(l,1)).*0.25), (ylimits(l,2) - (ylimits(l,2)-ylimits(l,1)).*0.05), strcat('p=',num2str(p1(l),'%.2e')),'Color','black','FontSize',10);

%plot little red (maximum) or purple (minimum) bars for the associated extrema
seaTx = -seasurftemp_AgeMa;
sml = (ylimits(l,2)-ylimits(l,1))./40;

for(i=ind1+I3a-1:1:size(combinedseaTpeaks,1))
    if combinedseaTextrtype(i) == 1
        x = [seaTx(combinedseaTlocs(i)),seaTx(combinedseaTlocs(i))];
        y = [(smoothedseasurftemp(combinedseaTlocs(i))+sml),(smoothedseasurftemp(combinedseaTlocs(i))+4.*sml)];
        line(x,y,'Color','red','LineStyle','-','LineWidth',1.5)
    else
        x = [seaTx(combinedseaTlocs(i)),seaTx(combinedseaTlocs(i))];
        y = [(smoothedseasurftemp(combinedseaTlocs(i))-sml),(smoothedseasurftemp(combinedseaTlocs(i))-4.*sml)];
        line(x,y,'Color',[1 0 1],'LineStyle','-','LineWidth',1.5)
    end
end


xlabel( 'Ma', 'Interpreter', 'none' );
ylabel( 'Sea surface T (^{\circ}C)', 'Interpreter', 'tex' );


xlim([-130 0]);
ylim([10 35]);
xticks(-120:20:0);
xticklabels({'120','100','80','60','40','20','0'});










% construction of panel H: Global mean surface temperature in the past 130
% Ma, data from Judd et al. (2024)
% https://www.science.org/doi/abs/10.1126/science.adk3705 . The blue line
% represents the median estimate while the dark and light gray bands
% represent the 68% and 90% confidence intervals, respectively.

l=6;

subplot(3,3,l);
hold on

x3 = [-AverageAge_GMST', fliplr(-AverageAge_GMST')];

%GMST_95 and GMST_05 are upper and lower thresholds of 90% confidence intervals in Judd et al. (2024)
inBetween3 = [GMST_95', fliplr(GMST_05')];
fill(x3, inBetween3, [0.9 0.9 0.9]);

%GMST_84 and GMST_16 are upper and lower thresholds of 68% confidence intervals in Judd et al. (2024)
inBetween4 = [GMST_84', fliplr(GMST_16')];
fill(x3, inBetween4, [0.7 0.7 0.7]);

plot(-AverageAge_GMST,GMST_50,'Marker','none','LineStyle','-','Color','blue','LineWidth',1.5);

% plot peaks of panel A, leave out peak at 8 Ma because not significant in exponential-Gaussian mixture model analysis in panel B
y = [5 50];
for(i=1:1:size(pks,1)-1)
  x = [mwxflip(locs(i)),mwxflip(locs(i))];
  line(x,y,'Color','black','LineStyle','--')
end

xlabel( 'Ma', 'Interpreter', 'none' );
ylabel( 'Global mean surface T (^{\circ}C)', 'Interpreter', 'tex' );

%find maxima in GMST profile
[GMSTpeaks, GMSTlocs, GMSTw, GMSTp] = findpeaks(GMST_50);
%find minima in GMST profile
[GMSTpeaksmin, GMSTlocsmin, GMSTwmin, GMSTpmin] = findpeaks(-GMST_50);

%combine extrema and sort on prominence GMSTp, GMSTpmin
[combinedGMSTp,II] = sort([GMSTp' GMSTpmin']',1) ;
combinedGMSTlocs = [GMSTlocs' GMSTlocsmin']' ;
combinedGMSTlocs = combinedGMSTlocs(II);
combinedGMSTpeaks = [GMSTpeaks' GMSTpeaksmin']' ;
combinedGMSTpeaks = combinedGMSTpeaks(II);
% extremum type = 1 for maximum, -1 for minimum
combinedGMSTextrtype = [ones(size(GMSTp')) -ones(size(GMSTpmin'))]' ;
combinedGMSTextrtype = combinedGMSTextrtype(II);

%limit to [0 115] Ma
JJ=[];
for k=1:1:size(combinedGMSTpeaks,1)
    if(AverageAge_GMST(combinedGMSTlocs(k))<115)
        JJ = [JJ k];
    end
end

combinedGMSTpeaks = combinedGMSTpeaks(JJ);
combinedGMSTlocs = combinedGMSTlocs(JJ);
combinedGMSTp = combinedGMSTp(JJ);
combinedGMSTextrtype = combinedGMSTextrtype(JJ);

%coupling distance analysis between relative residuals peaks (panel A) and
%GMST extrema, for different prominence thresholds. Sets consisting of the
%top-1 up to the top-18 extrema were tested by comparing the observed
%average nearest distance to a WGD establishement rate peak for any given
%set of extrema to a null distribution of average nearest distances
%obtained by uniformly sampling 9 WGD peaks across the past 115 Ma a
%million times (scenario 0) or to a null distribution of average nearest
%distances obtained by sampling the age of the first peak randomly in the
%interval [0 Ma, 30 Ma] and positioning 8 additional WGD peaks randomly at
%a distance from the previous peak sampled from a Gaussian distribution fit
%to the observed WGD establishment rate inter-peak distances (scenario 1,
%10^6 samples were used to construct the null distribution).


for k=size(combinedGMSTpeaks,1):-1:1
    k
    % use couplingDistance function in couplingDistance.m to calculate mean nearest distance between selected GMST extrema k:1:size(combinedGMSTpeaks,1) and significant peaks in panel A
    cdist = couplingDistance(-mwxflip(locs(1:1:(size(pks,1)-1))),AverageAge_GMST(combinedGMSTlocs(k:1:size(combinedGMSTpeaks,1))));

    %now do the same for 1,000,000 sets of randomized WGD establishment rate per lineage peaks (depends on randomization scenario chosen)
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
            for j=2:1:(size(pks,1)-1)
                rndpks(j)=rndpks(j-1) + rndpkdists(j-1);
            end
        end

        rndcdist(i)=couplingDistance(rndpks,AverageAge_GMST(combinedGMSTlocs(k:1:size(combinedGMSTpeaks,1))));
    end
    
    %approximate p-value for peak coupling p=(n+1)/(k+1) with n the number of times the randomized coupling distance is lower or equal to the observed coupling distance 
    ppp = find(sort(rndcdist) >= cdist, 1)./1000001;
    %store in p-value matrix across panels (l) and extremum sets (k)
    if(~isempty(ppp))
        pp(l,k) = ppp;
    else
        pp(l,k) = 1;
    end
end

%find set of extrema with minimum p-value 
ind1=size(combinedGMSTpeaks,1)-2*(size(pks,1)-1)+1;
if ind1<=0
    ind1=1;
end
ind2=size(combinedGMSTpeaks,1);
[Ma,I3a] = min(pp(l,ind1:ind2));
if(~isempty(I3a))
    p1(l)=pp(l,ind1+I3a-1);
else
    p1(l)= NaN; 
end

%plot minimum p-value   
t=text((xlimits(l,2) - (xlimits(l,2)-xlimits(l,1)).*0.25), (ylimits(l,2) - (ylimits(l,2)-ylimits(l,1)).*0.05), strcat('p=',num2str(p1(l),'%.2e')),'Color','black','FontSize',10);

%plot little red (maximum) or purple (minimum) bars for the associated extrema
GMSTx = -AverageAge_GMST;
sml = (ylimits(l,2)-ylimits(l,1))./40;

for(i=ind1+I3a-1:1:size(combinedGMSTpeaks,1))
    if combinedGMSTextrtype(i) == 1
        x = [GMSTx(combinedGMSTlocs(i)),GMSTx(combinedGMSTlocs(i))];
        y = [(GMST_95(combinedGMSTlocs(i))+sml),(GMST_95(combinedGMSTlocs(i))+4.*sml)];
        line(x,y,'Color','red','LineStyle','-','LineWidth',1.5)
    else
        x = [GMSTx(combinedGMSTlocs(i)),GMSTx(combinedGMSTlocs(i))];
        y = [(GMST_05(combinedGMSTlocs(i))-sml),(GMST_05(combinedGMSTlocs(i))-4.*sml)];
        line(x,y,'Color',[1 0 1],'LineStyle','-','LineWidth',1.5)
    end
end


xlim([-130 0]);
ylim([5 50]);
xticks(-120:20:0);
xticklabels({'120','100','80','60','40','20','0'});










% construction of panel I: Atmospheric CO2 concentration in the past 130
% Ma, data from Foster et al. (2017) https://doi.org/10.1038/ncomms14845 .
% The dark blue solid curve represents the probability maximum while the
% dark and light gray bands represent the 68% and 95% confidence intervals,
% respectively.

l=9;

subplot(3,3,l);
hold on

x3 = [-CO2_AgeMa', fliplr(-CO2_AgeMa')];
%up95 and lw95 are upper and lower thresholds of 95% confidence intervals in Foster et al. (2017)
inBetween3 = [CO2_up95', fliplr(CO2_lw95')];
fill(x3, inBetween3, [0.9 0.9 0.9]);

%up68 and lw68 are upper and lower thresholds of 68% confidence intervals in Foster et al. (2017)
inBetween4 = [CO2_up68', fliplr(CO2_lw68')];
fill(x3, inBetween4, [0.7 0.7 0.7]);

plot(-CO2_AgeMa,pCO2ProbabilityMaximum,'Marker','none','LineStyle','-','Color','blue','LineWidth',1.5);

% plot peaks of panel A, leave out peak at 8 Ma because not significant in exponential-Gaussian mixture model analysis in panel B
y = [-100 1900];
for(i=1:1:size(pks,1)-1)
  x = [mwxflip(locs(i)),mwxflip(locs(i))];
  line(x,y,'Color','black','LineStyle','--')
end

xlabel( 'Ma', 'Interpreter', 'none' );
ylabel( 'CO_2 (ppm)', 'Interpreter', 'tex' );

%find maxima in CO2 profile
[CO2peaks, CO2locs, CO2w, CO2p] = findpeaks(pCO2ProbabilityMaximum);
%find minima in CO2 profile
[CO2peaksmin, CO2locsmin, CO2wmin, CO2pmin] = findpeaks(-pCO2ProbabilityMaximum);

%combine extrema and sort on prominence CO2p, CO2pmin
[combinedCO2p,II] = sort([CO2p' CO2pmin']',1) ;
combinedCO2locs = [CO2locs' CO2locsmin']' ;
combinedCO2locs = combinedCO2locs(II);
combinedCO2peaks = [CO2peaks' CO2peaksmin']' ;
combinedCO2peaks = combinedCO2peaks(II);
% extremum type = 1 for maximum, -1 for minimum
combinedCO2extrtype = [ones(size(CO2p')) -ones(size(CO2pmin'))]' ;
combinedCO2extrtype = combinedCO2extrtype(II);

%limit to [0 115] Ma
JJ=[];
for k=1:1:size(combinedCO2peaks,1)
    if(CO2_AgeMa(combinedCO2locs(k))<115)
        JJ = [JJ k];
    end
end

combinedCO2peaks = combinedCO2peaks(JJ);
combinedCO2locs = combinedCO2locs(JJ);
combinedCO2p = combinedCO2p(JJ);
combinedCO2extrtype = combinedCO2extrtype(JJ);

%coupling distance analysis between relative residuals peaks (panel A) and
%CO2 extrema, for different prominence thresholds. Sets consisting of the
%top-1 up to the top-18 extrema were tested by comparing the observed
%average nearest distance to a WGD establishement rate peak for any given
%set of extrema to a null distribution of average nearest distances
%obtained by uniformly sampling 9 WGD peaks across the past 115 Ma a
%million times (scenario 0) or to a null distribution of average nearest
%distances obtained by sampling the age of the first peak randomly in the
%interval [0 Ma, 30 Ma] and positioning 8 additional WGD peaks randomly at
%a distance from the previous peak sampled from a Gaussian distribution fit
%to the observed WGD establishment rate inter-peak distances (scenario 1,
%10^6 samples were used to construct the null distribution).


for k=size(combinedCO2peaks,1):-1:(size(combinedCO2peaks,1)-2*(size(pks,1)-1)+1)
    k
    % use couplingDistance function in couplingDistance.m to calculate mean nearest distance between selected CO2 extrema k:1:size(combinedCO2peaks,1) and significant peaks in panel A
    cdist = couplingDistance(-mwxflip(locs(1:1:(size(pks,1)-1))),CO2_AgeMa(combinedCO2locs(k:1:size(combinedCO2peaks,1))));

    %now do the same for 1,000,000 sets of randomized WGD establishment rate per lineage peaks (depends on randomization scenario chosen)
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
            for j=2:1:(size(pks,1)-1)
                rndpks(j)=rndpks(j-1) + rndpkdists(j-1);
            end
        end

        rndcdist(i)=couplingDistance(rndpks,CO2_AgeMa(combinedCO2locs(k:1:size(combinedCO2peaks,1))));
    end
    
    %approximate p-value for peak coupling p=(n+1)/(k+1) with n the number of times the randomized coupling distance is lower or equal to the observed coupling distance 
    ppp = find(sort(rndcdist) >= cdist, 1)./1000001;
    %store in p-value matrix across panels (l) and extremum sets (k)
    if(~isempty(ppp))
        pp(l,k) = ppp;
    else
        pp(l,k) = 1;
    end
end

%find set of extrema with minimum p-value 
ind1=size(combinedCO2peaks,1)-2*(size(pks,1)-1)+1;
if ind1<=0
    ind1=1;
end
ind2=size(combinedCO2peaks,1);
[Ma,I3a] = min(pp(l,ind1:ind2));
if(~isempty(I3a))
    p1(l)=pp(l,ind1+I3a-1);
else
    p1(l)= NaN; 
end

%plot minimum p-value   
t=text((xlimits(l,2) - (xlimits(l,2)-xlimits(l,1)).*0.25), (ylimits(l,2) - (ylimits(l,2)-ylimits(l,1)).*0.05), strcat('p=',num2str(p1(l),'%.2e')),'Color','black','FontSize',10);

%plot little red (maximum) or purple (minimum) bars for the associated extrema
CO2x = -CO2_AgeMa;
sml = (ylimits(l,2)-ylimits(l,1))./40;

for(i=ind1+I3a-1:1:size(combinedCO2peaks,1))
    if combinedCO2extrtype(i) == 1
        x = [CO2x(combinedCO2locs(i)),CO2x(combinedCO2locs(i))];
        y = [(CO2_up95(combinedCO2locs(i))+sml),(CO2_up95(combinedCO2locs(i))+4.*sml)];
        line(x,y,'Color','red','LineStyle','-','LineWidth',1.5)
    else
        x = [CO2x(combinedCO2locs(i)),CO2x(combinedCO2locs(i))];
        y = [(CO2_lw95(combinedCO2locs(i))-sml),(CO2_lw95(combinedCO2locs(i))-4.*sml)];
        line(x,y,'Color',[1 0 1],'LineStyle','-','LineWidth',1.5)
    end
end

xlim([-130 0]);
ylim([-100 1900]);
xticks(-120:20:0);
xticklabels({'120','100','80','60','40','20','0'});




















