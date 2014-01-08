function [mu, mu_std] = fitSingleExpCDF(data,range,plotDistr)
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


function v = costFct(mu, t_ECDF, ECDF)
T = 1-exp(-1./mu.*t_ECDF(1));
CDF = 1 - exp(-1./mu.*t_ECDF);
CDF = (CDF - T) / (1-T);
v = CDF - ECDF;
