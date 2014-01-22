function [Runlength_mean,Runlength_std,Velocity_mean,Velocity_std]=Calc_Velocity_RunLength(varargin)
% Calc_Velocity_RunLength calculates the average velocity and run length for a
% dataset of single molecule runs. Inputs can be either velocity (nm/s) or run
% length (um), or both. User can also change input optional parameters for
% setting X-axis limit for velocity and run length.
%
%SYNOPSIS [Runlength_mean,Runlength_std,Velocity_mean,Velocity_std] = ...
%    Calc_Velocity_RunLength(velocity,[VelMax])
%    Calc_Velocity_RunLength(runlength,[MaxRun])
%    Calc_Velocity_RunLength(velocity,runlength,[MaxRun],[VelMax])
%       
%INPUT  velocity - velocity measurements (nm/s)
%       runlength - run length measurements (nm/s) 
%       
%       Optional: 
%           MaxRun - maximum run length distance (um) used for plotting [Default = 30 um]
%           VelMax - maximum velocity (nm/s) used for plotting [DEFAULT = 200 nm/s]
%           
%DEFAULTS   MaxRun = 30 um
%           VelMax = 200 nm/s
%
%OUTPUT Runlength_mean = average run length determined from exponential fit
%       Runlength_std = standard deviation of run length from expon. fit
%       Velocity_mean = average velocity
%       Velocity_std = standard deviation of velocity       
%
%       Three plots:
%            Velocity histogram
%            Exponential distribution of run lengths (fraction of runs)
%            Exponential distribution of run lengths (number of runs)
%           
%Written and used by the Reck-Peterson lab (http://reck-peterson.med.harvard.edu)
%using programs originally written by the Danuser lab (http://lccb.hms.harvard.edu). 

%% Read inputs and return errors
close all;                              %Close all open plots
optargs = {30 200};                     %Optional inputs. 30 um and 200 nm/s

%Flags for run length and velocity calculations
RunLengthFlag = 0;                      %0 - inactive; 1 - active
VelocityFlag = 0;

%Check that there are args given; error if > 4 or < 1: 
numvarargs = length(varargin);
if numvarargs > 4
    error('myApp:argChk','Too many inputs; requires at most 4 inputs')
end
if numvarargs < 1
    error('myApp:argChk','No inputs provided; requires at least 1 input')
end

%%Input options for a single input variable

%Parse arguments for all possible combinations.
%Set defaults: MaxRun = 30 um; MaxVel = 200 nm/s if not provided

if numvarargs == 1
    %If there is only one input, this means that default parameters will be used:
    MaxRun = cell2mat(optargs(:,1)); 
    VelMax = cell2mat(optargs(:,2));
    %Determine if input parameter is run length or velocity
    if mean(cell2mat(varargin(:,1))) > 40      %If mean is > 40, then input is velocity
        VelocityFlag = 1;
        disp('Velocity provided, calculating histogram')
        velocity = cell2mat(varargin(:,1));
    else
        RunLengthFlag = 1;
        disp('Run length provided, calculating histogram and exponential fit')
        runlength = cell2mat(varargin(:,1));
    end
end

%%Input options for two input options
if numvarargs == 2
    %Read input options
    test1 = cell2mat(varargin(:,1));
    test2 = cell2mat(varargin(:,2));
    
    %If both are velocity and run length
    if size(test1,1) > 1 && size(test2,1) > 1      
        velocity = test1;
        runlength = test2;
        MaxRun = cell2mat(optargs(:,1)); 
        VelMax = cell2mat(optargs(:,2));
        RunLengthFlag = 1;
        VelocityFlag = 1;
        disp('Velocity & run length data provided')
        if size(velocity,1) ~= size(runlength,1)
            error('myApp:argChk','-->Velocity and run length are not the same size')
        end
    end
    
    %If one is velocity OR run length and the other is optional input
    if size(test1,1) > 1 && size(test2,1) == 1 
        %Determine if test1 is velocity or run length
        if mean(test1) > 40      %If mean is > 40, then input is velocity
            VelocityFlag = 1;
            disp('Velocity provided, calculating histogram')
            velocity = cell2mat(varargin(:,1));
        else
            RunLengthFlag = 1;
            disp('Run length provided, calculating histogram and exponential fit')
            runlength = cell2mat(varargin(:,1));
        end
        %Determine if test2 is optional input for run length or velocity:
        if test2 > 100
            VelMax = test2;
            MaxRun = cell2mat(optargs(:,1));   
        else
            MaxRun = test2;
            VelMax = cell2mat(optargs(:,2));
        end
        
    end
    
end

if numvarargs == 3
    %Read input options
    test1 = cell2mat(varargin(:,1));
    test2 = cell2mat(varargin(:,2));
    test3 = cell2mat(varargin(:,3));
    
    %Determine if first two are velocity and run length
    if size(test1,1) > 1 && size(test2,1) > 1      
        velocity = test1;
        runlength = test2;
        if size(velocity,1) ~= size(runlength,1)
            error('myApp:argChk','-->Velocity and run length are not the same size')
        end
        disp('Velocity & run length data provided')
        %Determine if test3 is run length or velocity input option
        RunLengthFlag = 1;
        VelocityFlag = 1;
        if test3 > 100
            VelMax = test3;
            MaxRun = cell2mat(optargs(:,1));   
        else
            MaxRun = test3;
            VelMax = cell2mat(optargs(:,2));
        end
    end
    %Return error if they gave velocity or run length and TWO input options
    if size(test1,1) > 1 && size(test2,1) == 1 && size(test3,1) == 1
        error('myApp:argChk','Cannot use two default parameters if only calculating a single plot')
    end
end

if numvarargs == 4
    %Read input options
    velocity = cell2mat(varargin(:,1));
    runlength = cell2mat(varargin(:,2));
    MaxRun = cell2mat(varargin(:,3));
    VelMax = cell2mat(varargin(:,4));
    disp('Velocity & run length data provided')
    if size(velocity,1) ~= size(runlength,1)
            error('myApp:argChk','-->Velocity and run length are not the same size')
    end
    RunLengthFlag = 1;
    VelocityFlag = 1;
end

%% Run plotting and fitting
%Inputs
RunRange = 0:1:MaxRun;                  %Range for run length plot
if VelocityFlag == 1
    N=size(velocity,1);                     %Determine number of entries
    Velocity_mean=mean(velocity);           %Calculate mean velocity
    Velocity_std=std(velocity);             %Calculate std. dev. of velocity
    bin_size=3.5*Velocity_std/N^(1/3);      %Determine bin size
end

%Plot histogram of velocity
if VelocityFlag == 1
    plotFitVel(velocity,1,bin_size/2:bin_size:VelMax,0,Velocity_mean,Velocity_std); 
end

%Determine exponential decay of run lengths
if RunLengthFlag == 1
    [RL_mean,RL_std]=fitSingleExp(runlength,RunRange,4);
end

%Outputs
if RunLengthFlag == 1
    Runlength_mean=RL_mean;                 
    Runlength_std=RL_std;                   
end

%Plot histogram of run lengths along with fitted exponential decay
if RunLengthFlag == 1
    fancyHistRL(runlength,0,RunRange,Runlength_mean,Runlength_std); 
end

%Print average run length and velocity to console
if RunLengthFlag == 1
    disp(['Average run length = ',num2str(Runlength_mean),' +/- ',num2str(Runlength_std)])  
end
if VelocityFlag == 1
    disp(['Average velocity = ',num2str(Velocity_mean),' +/- ',num2str(Velocity_std)])
end

%% Subfunction for fitting exponential distribution

function plotFitVel(data,numComponents,numBins,normalized,vel_mean, vel_std)

options = statset('Display','final');
if isempty(numComponents)
    BIC = zeros(1,6);
    for k = 1:6
        a = gmdistribution.fit(data(:),k,'Options',options);
        BIC(k)= a.BIC;
        clear a;
    end
    [minBIC,numComponents] = min(BIC);
end


fit = gmdistribution.fit(data(:),numComponents,'Options',options);
%figure; hist(data(:),numBins);

[histVal,Xvals]=hist(data(:),numBins);

if normalized==1
    histVal=histVal/(sum(histVal)*(Xvals(2)-Xvals(1)));
    count=1;
else
    count=(sum(histVal)*(Xvals(2)-Xvals(1)));
end
figure1 = figure('Color',[1 1 1]); 
axes1 = axes('Parent',figure1,'FontName','Times New Roman');
bar(Xvals,histVal,'FaceColor',[0.729411780834198 0.831372559070587 0.95686274766922]);
normBinWidth=(Xvals(2)-Xvals(1))/10;
hold on;
for i=1:fit.NComponents
    mmgY(i,:)=fit.PComponents(i)*normpdf(min(Xvals):normBinWidth:max(Xvals),fit.mu(i),sqrt(fit.Sigma(1,1,i)));
    plot(min(Xvals):normBinWidth:max(Xvals),fit.PComponents(i)*count*normpdf(min(Xvals):normBinWidth:max(Xvals),fit.mu(i),sqrt(fit.Sigma(1,1,i))),'LineWidth',2,'Color',[0 0 0]);
end
mmgYcomb=sum(mmgY);
plot(min(Xvals):normBinWidth:max(Xvals),count*mmgYcomb','LineWidth',2,'Color',[0 1 0]);

set(gcf,'color','w');
axis([min(Xvals)-(Xvals(2)-Xvals(1)) max(Xvals)*1.1 0 max(histVal)*1.1]);
hold off;

if ~isempty(numComponents)
    BIC=fit.BIC;
end

% Create title
titleTxt = sprintf('Average velocity = %5.2f +/- %5.2f nm/s',vel_mean,vel_std);
title(titleTxt,'FontWeight','bold','FontSize',16,...
    'FontName','Times New Roman');

% Create xlabel
xlabel('Velocity (nm/s)','FontWeight','demi','FontSize',14,...
    'FontName','Times New Roman');

% Create ylabel
ylabel('Number','FontWeight','demi','FontSize',14,...
    'FontName','Times New Roman');

end
%% Plotting single exponential
   
function [mu,mu_std] = fitSingleExp(data,range,plotDistr)
% mu - fitted parameter
% mu_std = error of the fitted parameter
% data = list of the data to be fit with a single exponential function
%range = either number of bins or min:increment:max as described in hist
%function
% plotDistr dummy variable if 1-plots 1-CDF
%                             if 2-plots pdf
%                             if 3 plots both 1-cdf and pdf
%
opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

[ECDF,t_ECDF] = ecdf(data);
%get the parameter error
[mu,resnorm,~,~,~,~,J] = lsqnonlin(@costFct, mean(data), [], [], opts, t_ECDF, ECDF);

J = full(J);
mu_std = sqrt(resnorm/(numel(t_ECDF)-2)*inv(J'*J));

if plotDistr==4
    
    return
    
end


%plotting
if plotDistr==1 || plotDistr==3
    T = 1-exp(-1./mu.*t_ECDF(1));
    ECDF1 = 1 - exp(-1./mu.*t_ECDF);
    ECDF1 = (ECDF1 - T) / (1-T);
    
    figure; hold on; title('1-CDF');
    plot(t_ECDF,1-ECDF1,'r');
    plot(t_ECDF,1-ECDF,'b*');
    hold off;
end

if plotDistr==2 || plotDistr==3
    figure; hold on; title('PDF');
    [Y,X]=hist(data,range);
    dx=(X(2)-X(1))/10;
    t=0:dx:max(X);
    PDF1 = 1/mu * exp(-1/mu * t);
    PDF1 = PDF1*sum(Y)*dx*10;
    bar(X,Y); plot(t,PDF1,'r','lineWidth',2);
end


end


function v = costFct(mu, t_ECDF, ECDF)
T = 1-exp(-1./mu.*t_ECDF(1));
CDF = 1 - exp(-1./mu.*t_ECDF);
CDF = (CDF - T) / (1-T);
v = CDF - ECDF;
end

%% Plot histogram and fit

function fancyHistRL(data, lowerCutOff,binRange,RLmean,RLstd)

[mustats(1),mustats(2)] = fitSingleExp(data,binRange,4);
yval=hist(data,binRange);
% area under curve above cutoff minRunLength:
a = exp(-lowerCutOff/mustats(1));
% normalize to this area:
yval2 = yval*a/sum(yval*(binRange(2)-binRange(1)));
% Create figure
figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.149805447470817 0.130337078651685 0.755194552529183 0.794662921348315],...
    'FontName','Times New Roman');
bar(binRange,yval2,'FaceColor',[0.925490200519562 0.839215695858002 0.839215695858002]); hold on;
(binRange(2)-binRange(1))/100;
binRange2=binRange(1):(binRange(2)-binRange(1))/100:binRange(end);

plot(binRange2,exp(-binRange2/mustats(1))/mustats(1),'LineWidth',2,'Color',[0 0 0]);
axis([binRange(1) max(binRange)*1.1 0 max(yval2)*1.1]);

titleTxt = sprintf('Average run length = %5.2f +/- %5.2f \\mum',RLmean,RLstd);
% Create title
title(titleTxt,'FontWeight','bold',...
    'FontSize',16,...
    'FontName','Times New Roman');

% Create xlabel
xlabel('Run length (\mum)','FontWeight','demi','FontSize',14,...
    'FontName','Times New Roman');

% Create ylabel
ylabel('Number of runs','FontWeight','demi','FontSize',14,...
    'FontName','Times New Roman');
%['Average run length = ',num2str(Runlength_mean),' +/- ',num2str(Runlength_std)])
end
end
