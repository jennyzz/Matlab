function [Runlength_mean,Runlength_std,Velocity_mean,Velocity_std]=Velocity_RunLength(velocity,runlength)
% DISCONTINUED. This script has been renamed as Calc_Velocity_RunLength.
% Velocity_RunLength calculates the average velocity and run length for a
% dataset. 
%
%SYNOPSIS [Runlength_mean,Runlength_std,Velocity_mean,Velocity_std] = ...
%    Velocity_RunLength(velocity,runlength)
%
%INPUT  velocity - velocity in nm/s
%       runlength - run length in um 

%% Script that will run fitting and plotting
close all;

error('myApp:argChk','This script has been discontinued. Use Calc_Velocity_Runlength instead.')

RunRange = 0:1:30;                         %Range for histogram
VelMax = 200;

N=size(velocity,1);                     %Determine number of entries

Velocity_mean=mean(velocity);           %Calculate mean velocity
Velocity_std=std(velocity);             %Calculate std. dev. of velocity
bin_size=3.5*Velocity_std/N^(1/3);      %Determine bin size
plotFit(velocity,1,bin_size/2:bin_size:VelMax,0);   %Plot histogram of velocity 
[RL_mean,RL_std]=fitSingleExp(runlength,RunRange,4);      %Determine exponential decay of run lengths
pltExpFitHist(runlength,0,RunRange);     %Plot histogram of run lengths
Runlength_mean=RL_mean;                 %Save output
Runlength_std=RL_std;                   %Save output

disp(['Average run length = ',num2str(Runlength_mean),' +/- ',num2str(Runlength_std)])  %Print to console run length +/- std. dev
disp(['Average velocity = ',num2str(Velocity_mean),' +/- ',num2str(Velocity_std)])      %Print to console velocity +/- std. dev
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

function pltExpFitHist(data,lowerCutOff,binRange)

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
