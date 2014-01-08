function [Runlength_mean,Runlength_std,Velocity_mean,Velocity_std]=Velocity_RunLength(velocity,runlength)
% Velocity_RunLength calculates the average velocity and run length for a
% dataset. 
%
%SYNOPSIS [Runlength_mean,Runlength_std,Velocity_mean,Velocity_std] = ...
%    Velocity_RunLength(velocity,runlength)
%
%INPUT  velocity - velocity in nm/s
%       runlength - run length in um

close all;

N=size(velocity,1);                     %Determine number of entries

Velocity_mean=mean(velocity);           %Calculate mean velocity
Velocity_std=std(velocity);             %Calculate std. dev. of velocity
bin_size=3.5*Velocity_std/N^(1/3);      %Determine bin size
plotMMGfit(velocity,1,bin_size/2:bin_size:250,0);   %Plot histogram of velocity 
[RL_mean,RL_std] = fitSingleExpCDF(runlength,0:1:50,2);      %Determine exponential decay of run lengths
plotExpFitHist(runlength,0,0:1:50);     %Plot histogram of run lengths
Runlength_mean=RL_mean;                 %Save output
Runlength_std=RL_std;                   %Save output

disp(['Average run length = ',num2str(Runlength_mean),' +/- ',num2str(Runlength_std)])  %Print to console run length +/- std. dev
disp(['Average velocity = ',num2str(Velocity_mean),' +/- ',num2str(Velocity_std)])      %Print to console velocity +/- std. dev
end