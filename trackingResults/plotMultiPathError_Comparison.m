%% This function Generates and Compares the Multipath Errror Envelopes for 
%% GPS, E1B & DET 
% Created by Matthew Alcock and Dr Paul Blunt

%% Clean the environment ==================================================

close all;clear all;

%% Load Tracking Results ==================================================

% Load GPS Tracking Results
multipathGPS = load('trackResultsMultipath_GPS_1_Amp_0p5_0_Inc1_25s_EL_0p16.mat');
noMultipathGPS = load('trackResultsNoMultipath_GPS_1_25s_EL_0p16.mat');

% Load GAL Tracking Results
multipathGAL = load('trackResultsMultipath_GAL_1_Amp_0p5_0_Inc1_25s_EL_0p16.mat');
noMultipathGAL = load('trackResultsNoMultipath_GAL_1_25s_EL_0p16.mat');

% Load DET Tracking Results
multipathDET = load('trackResultsMultipath_DET_1_Amp_0p5_0_Inc1_25s_EL_0p16.mat');
noMultipathDET = load('trackResultsNoMultipath_DET_1_25s_EL_0p16.mat');

%% Initialize tracking variables ==========================================

secondsTracked = 25;

samplingFreq = 53e6;
codeFreq = 1.023e6;
sampleShiftsPerInc = 5;
incrementTime = 1;
chipsInSeconds = incrementTime*((samplingFreq/codeFreq)/sampleShiftsPerInc);

chipsTracked = secondsTracked/chipsInSeconds;

chipMetres = 299792458/1.023e6;
samples = secondsTracked*10;
MpathShift = floor(samples/40);

%% Initialize result structure ============================================

Results.meanErrorGPS = zeros(1, MpathShift);
Results.stdErrorGPS = zeros(1 ,MpathShift);

Results.meanErrorGAL = zeros(1, MpathShift);
Results.stdErrorGAL = zeros(1 ,MpathShift);

Results.meanErrorDET = zeros(1, MpathShift);
Results.stdErrorDET = zeros(1 ,MpathShift);

%% Calculate GPS Code Phase Difference ====================================

% Calculate GPS Multipath Error
codePhaseDiffGPS = chipMetres.*(noMultipathGPS.trackResults.codePhase - multipathGPS.trackResults.codePhase);
codePhaseDiffGPS(1) = 0;

% Average the Error every 4 seconds

for n = 1:1:MpathShift
    meanCodeErrorGPS = mean(codePhaseDiffGPS(((n-1)*40)+25:((n-1)*40)+35));
    stdCodeErrorGPS = std(codePhaseDiffGPS(((n-1)*40)+25:((n-1)*40)+35))/sqrt(n);
    Results.meanCodeErrorGPS(n) = meanCodeErrorGPS;
    Results.stdCodeErrorGPS(n) = stdCodeErrorGPS;
    err1 = Results.stdCodeErrorGPS;
end

%% Calculate Galileo Code Phase Difference ================================

% Calculate GAL Multipath Error
codePhaseDiffGAL = chipMetres.*(noMultipathGAL.trackResults.codePhase - multipathGAL.trackResults.codePhase);
codePhaseDiffGAL(1) = 0;

% Average the Error every 4 seconds

for n = 1:1:MpathShift
    meanCodeErrorGAL = mean(codePhaseDiffGAL(((n-1)*40)+25:((n-1)*40)+35));
    stdCodeErrorGAL = std(codePhaseDiffGAL(((n-1)*40)+25:((n-1)*40)+35))/sqrt(n);
    Results.meanCodeErrorGAL(n) = meanCodeErrorGAL;
    Results.stdCodeErrorGAL(n) = stdCodeErrorGAL;
    err2 = Results.stdCodeErrorGAL;
end

%% Calculate DET Code Phase Difference ====================================

% Calculate DET Multipath Error
codePhaseDiffDET = chipMetres.*(noMultipathDET.trackResults.codePhase - multipathDET.trackResults.codePhase);
codePhaseDiffDET(1) = 0;

% Average the Error every 4 seconds

for n = 1:1:MpathShift
    meanCodeErrorDET = mean(codePhaseDiffDET(((n-1)*40)+25:((n-1)*40)+35));
    stdCodeErrorDET = std(codePhaseDiffDET(((n-1)*40)+25:((n-1)*40)+35))/sqrt(n);
    Results.meanCodeErrorDET(n) = meanCodeErrorDET;
    Results.stdCodeErrorDET(n) = stdCodeErrorDET;
    err3 = Results.stdCodeErrorDET;
end

%% Plot the Multipath Error Envelopes =====================================

t = linspace(0,chipsTracked,MpathShift);

figure(79)
subplot(2,2,1)
plot(codePhaseDiffGPS,'r')
hold on
plot(codePhaseDiffGAL,'g')
hold on
plot(codePhaseDiffDET,'b')
title('Raw Multipath Envelopes')
xlabel('Time ms)')
ylabel('Range Error in Metres')

subplot(2,1,2)
plot(t,Results.meanCodeErrorGPS,'r');
hold on
errorbar(t,Results.meanCodeErrorGPS,err1,'.r');
hold on
plot(t,Results.meanCodeErrorGAL,'g');
hold on
errorbar(t,Results.meanCodeErrorGAL,err2,'.g');
hold on
plot(t,Results.meanCodeErrorDET,'b');
hold on
errorbar(t,Results.meanCodeErrorDET,err3,'.b');
title('Filtered Multipath Error and Loop Differences')
xlabel('Multipath in Chips')
ylabel('Range Error in Metres')
