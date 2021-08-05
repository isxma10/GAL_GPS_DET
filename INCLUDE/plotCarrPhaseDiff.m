clear all;close all;
CH1=load('trackingResultsPRN21_CH1_RHC.mat');
CH2=load('trackingResultsPRN21_CH2_RHC.mat');

diff = zeros(1,length(CH1.trackResults.carrPhase));
for i=1:length(CH1.trackResults.carrPhase)
    diff(i) = CH1.trackResults.carrPhase(i) - CH2.trackResults.carrPhase(i);
    if diff(i) < -2%-pi
       diff(i) = diff(i) + 2*pi; 
    end
    if diff(i) > 2%pi
       diff(i) = diff(i) - 2*pi; 
    end
   
end
logLengthSec = 9;
startPoint = 1;
endPoint = 89;
time = linspace(0.1,logLengthSec, endPoint)
diffNoMean = diff(startPoint:endPoint) - mean(diff(30:endPoint));

figure(5)
plot(time,diffNoMean)

mean(CH1.trackResults.CNo(3000:end))

mean(CH2.trackResults.CNo(3000:end))
rad2deg(mean(diff(30:endPoint)))
rad2deg(std(diffNoMean(30:endPoint)))