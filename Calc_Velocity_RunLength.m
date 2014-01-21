function [Runlength_mean,Runlength_std,Velocity_mean,Velocity_std]=Calc_Velocity_RunLength(velocity,runlength,varargin)
% Calc_Velocity_RunLength calculates the average velocity and run length for a
% dataset of single molecule runs. Inputs are velocity (nm/s) and run
% length (um).
%
%SYNOPSIS [Runlength_mean,Runlength_std,Velocity_mean,Velocity_std] = ...
%    Calc_Velocity_RunLength(velocity,runlength,MaxRun,VelMax)
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
%Written and used by the Reck-Peterson lab (reck-peterson.med.harvard.edu)
%using programs originally written by the Danuser lab: 
%
%           (Danuser command)   -->     (Sub-function command)
%           plotMMGFit          -->     plotFit
%           fitSingleExpCDF     -->     fitSingleExp
%           pltExpFitHist       -->     pltExpFit

%% Master function that will run fitting and plotting
close all;                              %Close all open plots

%Check number of optional inputs. If > 2 send error. 
numvarargs = length(varargin);
if numvarargs > 2
    error('myApp:argChk','Too many inputs; requires at most 2 optional inputs')
end

%Check that velocitues and run lengths are the same size
if size(velocity,1) ~= size(runlength,1)
    error('myApp:argChk','Velocity and run length are not the same size')
end

%Set defaults: MaxRun = 30 um; MaxVel = 200 nm/s
optargs = {30 200};                     %Set defaults in optargs                 
optargs(1:numvarargs) = varargin;       %For each emptpy value in varagin put in default
MaxRun = cell2mat(optargs(:,1));        %Rename optarg(:,1) as MaxRun
VelMax = cell2mat(optargs(:,2));        %Rename optarg(:,2) as VelMax

%Inputs
RunRange = 0:1:MaxRun;                  %Range for run length plot
N=size(velocity,1);                     %Determine number of entries
Velocity_mean=mean(velocity);           %Calculate mean velocity
Velocity_std=std(velocity);             %Calculate std. dev. of velocity
bin_size=3.5*Velocity_std/N^(1/3);      %Determine bin size

%Plot histogram of velocity
plotFit(velocity,1,bin_size/2:bin_size:VelMax,0);  

%Determine exponential decay of run lengths
[RL_mean,RL_std]=fitSingleExp(runlength,RunRange,4);

%Plot histogram of run lengths along with fitted exponential decay
pltExpFit(runlength,0,RunRange);     

%Outputs
Runlength_mean=RL_mean;                 
Runlength_std=RL_std;                   

%Print average run length and velocity to console
disp(['Average run length = ',num2str(Runlength_mean),' +/- ',num2str(Runlength_std)])  
disp(['Average velocity = ',num2str(Velocity_mean),' +/- ',num2str(Velocity_std)])      
end

%% Subfunction for fitting exponential distribution

function plotFit(data,numComponents,numBins,normalized)

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
figure; set(0,'DefaultAxesFontSize',16); bar(Xvals,histVal);
normBinWidth=(Xvals(2)-Xvals(1))/10;
hold on;
for i=1:fit.NComponents
    mmgY(i,:)=fit.PComponents(i)*normpdf(min(Xvals):normBinWidth:max(Xvals),fit.mu(i),sqrt(fit.Sigma(1,1,i)));
    plot(min(Xvals):normBinWidth:max(Xvals),fit.PComponents(i)*count*normpdf(min(Xvals):normBinWidth:max(Xvals),fit.mu(i),sqrt(fit.Sigma(1,1,i))),'r','LineWidth',2);
end
mmgYcomb=sum(mmgY);
plot(min(Xvals):normBinWidth:max(Xvals),count*mmgYcomb','g','LineWidth',2);

set(gcf,'color','w');
axis([min(Xvals)-(Xvals(2)-Xvals(1)) max(Xvals)*1.1 0 max(histVal)*1.1]);
hold off;

if ~isempty(numComponents)
    BIC=fit.BIC;
end

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

function pltExpFit(data,lowerCutOff,binRange)

[mustats(1),mustats(2)] = fitSingleExp(data,binRange,4);
yval=hist(data,binRange);
% area under curve above cutoff minRunLength:
a = exp(-lowerCutOff/mustats(1));
% normalize to this area:
yval2 = yval*a/sum(yval*(binRange(2)-binRange(1)));
figure; set(0,'DefaultAxesFontSize',16);
bar(binRange,yval2); hold on;
(binRange(2)-binRange(1))/100;
binRange2=binRange(1):(binRange(2)-binRange(1))/100:binRange(end);

plot(binRange2,exp(-binRange2/mustats(1))/mustats(1),'r','linewidth',2);
axis([binRange(1) max(binRange)*1.1 0 max(yval2)*1.1]);
set(gcf,'color','w');

hold off;
figure; set(0,'DefaultAxesFontSize',16);
bar(binRange,yval); hold on;
plot(binRange2,exp(-binRange2/mustats(1))/mustats(1) / a*sum(yval*(binRange(2)-binRange(1))),'r','linewidth',2);
axis([0 max(binRange)*1.1 0 max(yval)*1.1]);
set(gcf,'color','w');
hold off;
end
