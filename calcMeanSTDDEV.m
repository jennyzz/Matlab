function out = calcMeanSTDDEV(input,framerate)
%calcMeanSTDDEV takes a matrix of individual mean squared displacement 
%calculations [m x n], where each row (m) is a specific timepoint and 
%each column (n) is a separately tracked particle.  
%
%SYNOPSIS outputMatrix = ...
%    calcMeanSTDDEV(input,framerate)
%
%INPUT  input - matrix of mean squared displacement calculations, where
%each column is a separate particle.
%       framerate - length of time between exposures (ie framerate)

out = zeros(size(input,1),3);
time = framerate;

for i=1:size(input,1)
    
    avg = mean(input(i,:));
    stddev = std(input(i,:));
    
    out(i,1) = time;
    out(i,2) = avg;
    out(i,3) = stddev;
    
    time = time + framerate;
end

