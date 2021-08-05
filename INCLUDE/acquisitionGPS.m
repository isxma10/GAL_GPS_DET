function acqResults = acquisitionGPS(longSignal, settings)
%Function performs cold start acquisition on the collected "data". It
%searches for GPS signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - 11 ms of raw signal from the front-end 
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the 
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number. 
 
%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis and Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on Peter Rinder and Nicolaj Bertelsen
% Adapted by Matthew Alcock and Dr Paul Blunt
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%CVS record:
%$Id: acquisition.m,v 1.1.2.12 2006/08/14 12:08:03 dpl Exp $

%% Initialization =========================================================

% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));

% Create two 1msec vectors of data to correlate with and one with zero DC
signal1 = longSignal(1 : 10*samplesPerCode);
signal2 = longSignal(10*samplesPerCode+1 : 20*samplesPerCode);

signal0DC = longSignal - mean(longSignal); 

% Find sampling period
ts = 1 / settings.samplingFreq;

% Find phase points of the local carrier wave 
phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;

% Number of the frequency bins for the given acquisition band (50Hz steps)
numberOfFrqBins = round(settings.acqSearchBand * 20) + 1;

% Generate all C/A codes and sample them according to the sampling freq.
caCodesTable = makeCaTable(settings);


%--- Initialize arrays to speed up the code -------------------------------
% Search results of all frequency bins and code shifts (for one satellite)
results     = zeros(numberOfFrqBins, samplesPerCode);

% Carrier frequencies of the frequency bins
frqBins     = zeros(1, numberOfFrqBins);


%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 32);
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 32);
% noise floor values for SNR estimation
acqResults.noiseValue = 0;
%---- initialise counter for noise floor values to zero
NoiseNum = 0;
% intialise noise monitors
TempNoiseValue = zeros(1,numberOfFrqBins);
noiseValue = zeros(1, 32);

fprintf('(');

% Perform search for all listed PRN numbers ...
for PRN = settings.acqSatelliteList

%% Correlate signals ======================================================   
    %--- Perform DFT of C/A code ------------------------------------------
    caCodeFreqDom = conj(fft(caCodesTable(PRN, :)));

    %--- Make the correlation for whole frequency band (for all freq. bins)
    for frqBinIndex = 1:numberOfFrqBins
        
        acqRes1 = zeros(1,samplesPerCode);
        acqRes2 = zeros(1,samplesPerCode);
        
        for timeIndex = 1:10

            %--- Generate carrier wave frequency grid (0.05kHz step) -----------
            frqBins(frqBinIndex) = settings.IF - ...
                                   (settings.acqSearchBand/2) * 1000 + ...
                                   0.05e3 * (frqBinIndex - 1);

            %--- Generate local sine and cosine -------------------------------
            sinCarr = sin(frqBins(frqBinIndex) * phasePoints);
            cosCarr = cos(frqBins(frqBinIndex) * phasePoints);

            %--- "Remove carrier" from the signal -----------------------------
            I1      = sinCarr .* signal1(((timeIndex-1)*samplesPerCode)+1:timeIndex*samplesPerCode);
            Q1      = cosCarr .* signal1(((timeIndex-1)*samplesPerCode)+1:timeIndex*samplesPerCode);
            I2      = sinCarr .* signal2(((timeIndex-1)*samplesPerCode)+1:timeIndex*samplesPerCode);
            Q2      = cosCarr .* signal2(((timeIndex-1)*samplesPerCode)+1:timeIndex*samplesPerCode);
            %--- Convert the baseband signal to frequency domain --------------
            IQfreqDom1 = fft(I1 + j*Q1);
            IQfreqDom2 = fft(I2 + j*Q2);

            %--- Multiplication in the frequency domain (correlation in time
            %domain)
            convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
            convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;

            %--- Perform inverse DFT and store correlation results ------------
            acqRes1 = acqRes1 + (abs(ifft(convCodeIQ1)) .^ 2);
            acqRes2 = acqRes2 + (abs(ifft(convCodeIQ2)) .^ 2);
        end
        
        %--- Check which msec had the greater power and save that, will
        %"blend" 1st and 2nd msec but will correct data bit issues
        if (max(acqRes1) > max(acqRes2))
            results(frqBinIndex, :) = acqRes1;
        else
            results(frqBinIndex, :) = acqRes2;
        end
        
       TempNoiseValue(frqBinIndex) = mean(acqRes1);
                
    end % frqBinIndex = 1:numberOfFrqBins
    
    noiseValue(PRN) = mean(TempNoiseValue);
    
    %--- Plot FFTs of the signal acquistions if plotFFTs is high ---------
    if settings.plotFFTs == 1
        
        yrange = linspace(- 0.05e3*(numberOfFrqBins/2),0.05e3*(numberOfFrqBins/2),numberOfFrqBins);
        xrange = linspace(-settings.codeLength/2,settings.codeLength/2,samplesPerCode);
        
        figure(PRN)
        surf(xrange,yrange,results);
        shading INTERP;
         
        title ('Acquisition results');
        xlabel('Code delay (chips)');
        ylabel('Frequency offset (Hz)');
        axis  ([-settings.codeLength/2 settings.codeLength/2 ...
            -0.5e3*settings.acqSearchBand 0.5e3*settings.acqSearchBand ...
            min(min(results, [], 2)) max(max(results, [], 2))]);
    end %if

    
%% Look for correlation peaks in the results ==============================
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak
    
    %--- Find the correlation peak and the carrier frequency --------------
    [peakSize frequencyBinIndex] = max(max(results, [], 2));

    %--- Find code phase of the same correlation peak ---------------------
    [peakSize codePhase] = max(max(results));
    

    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
    excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    excludeRangeIndex2 = codePhase + samplesPerCodeChip;

    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < 2
        codePhaseRange = excludeRangeIndex2 : ...
                         (samplesPerCode + excludeRangeIndex1);
                         
    elseif excludeRangeIndex2 >= samplesPerCode
        codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
                         excludeRangeIndex1;
    else
        codePhaseRange = [1:excludeRangeIndex1, ...
                          excludeRangeIndex2 : samplesPerCode];
    end

    %--- Find the second highest correlation peak in the same freq. bin ---
    secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));

    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(PRN) = peakSize/secondPeakSize;
    
    % If the result is above threshold, then there is a signal ...
    if (peakSize/secondPeakSize) > settings.acqThreshold

%% Satellites Detected ====================================================
         
         %--- Indicate PRN number of the detected signal -------------------
         fprintf('%02d ', PRN);

        %--- Save properties of the detected satellite signal -------------
       
        acqResults.carrFreq(PRN)  = frqBins(frequencyBinIndex);
        acqResults.codePhase(PRN) = codePhase;
        %--- flag to remove the PRN from the noise floor calculation
        presentFlag = 1; 
    
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
        %--- flag to keep the PRN in the noise floor calculation
        presentFlag = 0; 
            
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
    
    % --- remove datasets with signals present from the noise floor
    % estimation
    if presentFlag == 1
        noiseValue(PRN) = 0;
    else
        % --- count the number of valid noise measurements
        NoiseNum = NoiseNum + 1;
    end %if
    
end    % for PRN = satelliteList

% ---- Calculate a measured noise floor value for the correlations with no
    % signal present
% acqResults.noiseValue = sum(noiseValue(:)) / (sqrt(2)*NoiseNum);
% find PRN with lowest noise value
% acqResults.noisePRN = find(noiseValue == min(noiseValue(find(noiseValue ~=0))));

    


%=== Acquisition is over ==================================================
fprintf(')\n');
