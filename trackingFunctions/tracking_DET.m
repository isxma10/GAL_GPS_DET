function [trackResults, channel]= tracking_DET(fid, channel, settings, acqResults)
% Performs code and carrier tracking for all channels.
%
%[trackResults, channel]= GIOVEtracking(fid, channel,
%settings, acqResults)
%
%   Inputs:
%       fid             - file identifier of the signal record.
%       channel    - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs:
%       trackResults - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute spreading
%                       code's starting positions, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on code by DMAkos Oct-1999
% Adapted by Matthew Alcock and Dr Paul Blunt
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% CVS record:
% $Id: tracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

%% Initialize result structure ============================================

% only track across code intervals
codePeriods = floor(settings.msToProcess/4);     % For GIOVE one L1B code is 4ms

% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock

% The absolute sample in the record of the L1B code start:
trackResults.absoluteSample = zeros(1, codePeriods);

% Freq of the code:
trackResults.codeFreq       = inf(1, settings.msToProcess);
% Freq of the subcarrier:
trackResults.subCarrFreq       = inf(1, settings.msToProcess);

% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, settings.msToProcess);

% Outputs from the correlators (In-phase):
trackResults.I_I_P            = zeros(1, codePeriods);
trackResults.I_I_E            = zeros(1, codePeriods);
trackResults.I_I_L            = zeros(1, codePeriods);
trackResults.I_E_P            = zeros(1, codePeriods);
trackResults.I_L_P            = zeros(1, codePeriods);
trackResults.I_Q_P            = zeros(1, codePeriods);

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_I_E            = zeros(1, codePeriods);
trackResults.Q_I_P            = zeros(1, codePeriods);
trackResults.Q_I_L            = zeros(1, codePeriods);
trackResults.Q_E_P            = zeros(1, codePeriods);
trackResults.Q_L_P            = zeros(1, codePeriods);
trackResults.Q_Q_P            = zeros(1, codePeriods);

% Loop discriminators
trackResults.dllDiscr       = inf(1, codePeriods);
trackResults.dllDiscrFilt   = inf(1, codePeriods);
trackResults.sllDiscr       = inf(1, codePeriods);
trackResults.sllDiscrFilt   = inf(1, codePeriods);
trackResults.pllDiscr       = inf(1, codePeriods);
trackResults.pllDiscrFilt   = inf(1, codePeriods);

trackResults.codePhase = zeros(1, length(settings.measurementPoints));
trackResults.deltaCarrPhase = zeros(1, length(settings.measurementPoints));
trackResults.carrCycleCount = zeros(1, length(settings.measurementPoints));
trackResults.msOfMeasurement = zeros(1, length(settings.measurementPoints));
trackResults.carrPhase = zeros(1, length(settings.measurementPoints));

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================

%--- DLL variables --------------------------------------------------------
% Define early-late offset (in chips)
% Code tracking loop parameters
settings.dllDampingRatio         = 0.7;
settings.dllNoiseBandwidth       = 2;       %[Hz]
settings.dllCorrelatorSpacing    = 0.5;     %[chips]
% Subcarrier tracking loop parameters
settings.sllDampingRatio         = 0.7;
settings.sllNoiseBandwidth       = 2;       %[Hz]
settings.sllCorrelatorSpacing    = 1;     %[subchips]
% Carrier tracking loop parameters
settings.pllDampingRatio         = 0.7;
settings.pllNoiseBandwidth       = 25;      %[Hz]
dllEarlyLateSpc = settings.dllCorrelatorSpacing;


% Summation interval
PDIcode = ((settings.codeLength/ settings.codeLengthCA)/1000);

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
                                    settings.dllDampingRatio, ...
                                    1.0);
%--- SLL variables --------------------------------------------------------
sllEarlyLateSpc = settings.sllCorrelatorSpacing;
% Summation interval
PDIsubCarr = ((settings.codeLength/ settings.codeLengthCA)/1000);

% Calculate filter coefficient values
[tau1SubCarr, tau2SubCarr] = calcLoopCoef(settings.sllNoiseBandwidth, ...
                                    settings.sllDampingRatio, ...
                                    1.0);
%--- PLL variables --------------------------------------------------------
% Summation interval
PDIcarr = ((settings.codeLength/ settings.codeLengthCA)/1000);

% Calculate filter coefficient values
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
                                    settings.pllDampingRatio, ...
                                    0.25);
hwb = waitbar(0,'Tracking...');

%% Start processing channels ==============================================
  
    % Only process if PRN is non zero (acquisition was successful)
for channelNr = 1:settings.numberOfChannels
    
    % Only process if PRN is non zero (acquisition was successful)
    if (channel(channelNr).PRN ~= 0)
        % Save additional information - each channel's tracked PRN
        trackResults(channelNr).PRN     = channel(channelNr).PRN;
        
        % Move the starting point of processing. Can be used to start the
        % signal processing at any point in the data record (e.g. for long
        % records). In addition skip through that data file to start at the
        % appropriate sample (corresponding to code phase). Assumes sample
        % type is schar (or 1 byte per sample) 
        fseek(fid, ...
              settings.skipNumberOfBytes + channel.codePhase-1, ...
              'bof');

        % DE tracking
        % Get a vector with the E1B code sampled 1x/chip with no BOC
        E1BCode = generateE1Bcode(channel(channelNr).PRN,0)';
        % Then make it possible to do early and late versions
        E1BCode = [E1BCode(4092) E1BCode E1BCode(1)];
        
              
        %--- Perform various initializations ------------------------------

        % define initial code frequency basis of NCO
        codeFreq      = settings.codeFreqBasis;
        % define initial subcarrier frequency basis of NCO
        bocRatio = settings.subFreqBasis/settings.codeFreqBasis;
        subFreq      = bocRatio*settings.codeFreqBasis;
        % define residual code phase (in chips)
        remCodePhase  = 0.0;
        % define residual subcarrier phase (in chips)
        remSubCarrPhase  = 0.0;
        % define carrier frequency which is used over whole tracking period
        carrFreq      = channel.acquiredFreq;
        carrFreqBasis = channel.acquiredFreq;
        % define residual carrier phase
        remCarrPhase  = 0.0;

        %code tracking loop parameters
        oldCodeNco   = 0.0;
        oldCodeError = 0.0;
        
        %subCarr tracking loop parameters
        oldSubCarrNco   = 0.0;
        oldSubCarrError = 0.0;

        %carrier/Costas loop parameters
        oldCarrNco   = 0.0;
        oldCarrError = 0.0;
        
        % running total of samples for the measurements
        totalSamplesRead = 0;
        measurementNumber = 1;
        carrierCycleCount = 0;

        %=== Process the number of specified code periods =================
        for loopCnt =  1:codePeriods
            
%% GUI update -------------------------------------------------------------
            % The GUI is updated every 50ms. This way Matlab GUI is still
            % responsive enough. At the same time Matlab is not occupied
            % all the time with GUI task.
            if (rem(loopCnt, 50) == 0)
                try
                   waitbar(loopCnt/codePeriods, ...
                            hwb, ...
                            ['Tracking: Ch ', int2str(channelNr), ...
                            ' of ', int2str(settings.numberOfChannels), ...
                            '; PRN#', int2str(channel(channelNr).PRN), ...
                            '; Completed ',int2str(loopCnt*4), ...
                            ' of ', int2str(codePeriods*4), ' msec']);                      
                catch
                    % The progress bar was closed. It is used as a signal
                    % to stop, "cancel" processing. Exit.
                    disp('Progress bar closed, exiting...');
                    return
                end
            end

%% Read next block of data ------------------------------------------------            
            % Find the size of a "block" or code period in whole samples
            
            % Update the phasestep based on code freq (variable) and
            % sampling frequency (fixed)
            codePhaseStep = codeFreq / settings.samplingFreq;
            
            blksize = ceil(((settings.codeLength)-remCodePhase) / codePhaseStep);
            
            % Read in the appropriate number of samples to process this
            % interation 
            [rawSignal, samplesRead] = fread(fid, ...
                                             blksize, settings.dataType);
            rawSignal = rawSignal';  %transpose vector
            
                        
            % If did not read in enough samples, then could be out of 
            % data - better exit 
            if (samplesRead ~= blksize)
                disp('Not able to read the specified number of samples  for tracking, exiting!')
                fclose(fid);
                return
            end

%% Set up all the code phase tracking information -------------------------
            % Define index into early code vector
            tcode       = ((remCodePhase-dllEarlyLateSpc)) : ...
                          codePhaseStep : ...
                          (((blksize-1)*codePhaseStep+remCodePhase-dllEarlyLateSpc));
            tcode2      = ceil(tcode) + 1;
            earlyCode   = E1BCode(tcode2);
            
            % Define index into late code vector
            tcode       = (remCodePhase+dllEarlyLateSpc) : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase+dllEarlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            lateCode    = E1BCode(tcode2);
            
            % Define index into prompt code vector
            tcode       = remCodePhase : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase);
            tcode2      = ceil(tcode) + 1;
            promptCode  = E1BCode(tcode2);
           
            %DE tracking
            remCodePhase = (tcode(blksize) + codePhaseStep) - 4092.0;
          

%% Set up all the subcarrier phase tracking information -------------------------

            subTime    = (0:blksize) ./ settings.samplingFreq;
            
            % Get the argument to sin/cos functions
            trigSubArg = ((subFreq * 2.0 * pi) .* subTime) + remSubCarrPhase;
            remSubCarrPhase = rem(trigSubArg(blksize+1), (2 * pi));
            
            trigSubArgEarly = ((subFreq * 2.0 * pi) .* subTime) + remSubCarrPhase - sllEarlyLateSpc;
            
            trigSubArgLate = ((subFreq * 2.0 * pi) .* subTime) + remSubCarrPhase + sllEarlyLateSpc;
                                   
            % Finally compute the signal to mix the collected data to baseband
            earlySubCarr = sign(sin(trigSubArgEarly(1:blksize)));
            lateSubCarr = sign(sin(trigSubArgLate(1:blksize)));
            promptSubCarr = sign(sin(trigSubArg(1:blksize)));
            
%% Generate the carrier frequency to mix the signal to baseband -----------
            time    = (0:blksize) ./ settings.samplingFreq;
            
            % Get the argument to sin/cos functions
            trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
            remCarrPhase = rem(trigarg(blksize+1), (2 * pi));
            
            % Finally compute the signal to mix the collected data to baseband
            carrCos = cos(trigarg(1:blksize));
            carrSin = sin(trigarg(1:blksize));

%% Generate the 12 standard accumulated values ---------------------------
            % First mix to baseband
            qBasebandSignal = carrCos .* rawSignal;
            iBasebandSignal = carrSin .* rawSignal;
                       
            % DE tracking
            % Now get early, late, and prompt values for each
            I_I_E = sum(promptSubCarr.* earlyCode  .* iBasebandSignal);
            Q_I_E = sum(promptSubCarr.* earlyCode  .* qBasebandSignal);
            I_I_P = sum(promptSubCarr.* promptCode .* iBasebandSignal);
            I_E_P = sum(earlySubCarr.* promptCode .* iBasebandSignal);
            I_L_P = sum(lateSubCarr.* promptCode .* iBasebandSignal);
            Q_I_P = sum(promptSubCarr.* promptCode .* qBasebandSignal);
            Q_E_P = sum(earlySubCarr.* promptCode .* qBasebandSignal);
            Q_L_P = sum(lateSubCarr.* promptCode .* qBasebandSignal);
            I_I_L = sum(promptSubCarr.* lateCode   .* iBasebandSignal);
            Q_I_L = sum(promptSubCarr.* lateCode   .* qBasebandSignal);
            
%% Find PLL error and update carrier NCO ----------------------------------

            % Implement carrier loop discriminator (phase detector)
            carrError = atan(Q_I_P / I_I_P) / (2.0 * pi);
            
            % Implement carrier loop filter and generate NCO command
            carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
                (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
            oldCarrNco   = carrNco;
            oldCarrError = carrError;

            % Modify carrier freq based on NCO command
            carrFreq = carrFreqBasis + carrNco;

            trackResults.carrFreq(loopCnt) = carrFreq;

%% Find DLL error and update code NCO -------------------------------------
            codeError = (sqrt(I_I_E * I_I_E + Q_I_E * Q_I_E) - sqrt(I_I_L * I_I_L + Q_I_L * Q_I_L)) / ...
                (sqrt(I_I_E * I_I_E + Q_I_E * Q_I_E) + sqrt(I_I_L * I_I_L + Q_I_L * Q_I_L));
            
            % Implement code loop filter and generate NCO command
            codeNco = oldCodeNco + (tau2code/tau1code) * ...
                (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;
            
            % Modify code freq based on NCO command
            codeFreq = settings.codeFreqBasis - codeNco;
            
            trackResults.codeFreq(loopCnt) = codeFreq;
                        
            
%% Find SLL error and update SLL NCO -------------------------------------
            subError = (sqrt(I_E_P * I_E_P + Q_E_P * Q_E_P) - sqrt(I_L_P * I_L_P + Q_L_P * Q_L_P)) / ...
                (sqrt(I_E_P * I_E_P + Q_E_P * Q_E_P) + sqrt(I_L_P * I_L_P + Q_L_P * Q_L_P));
            % Implement code loop filter and generate NCO command
            subNco = oldSubCarrNco + (tau2SubCarr/tau1SubCarr) * ...
                (subError - oldSubCarrError) + subError * (PDIsubCarr/tau1SubCarr);
            oldSubCarrNco   = subNco;
            oldSubCarrError = subError;
            
            % Modify code freq based on NCO command
            subFreq = settings.codeFreqBasis -(subNco);
%             subFreq = settings.codeFreqBasis;
            trackResults(channelNr).subFreq(loopCnt) = subFreq;
            
%% read mesurment data on the same sample ---------------------------------

            % accumulate the samples read
            totalSamplesRead = totalSamplesRead + samplesRead;
            
            % check if a measurement point is in the current block 
            if (settings.measurementPoints(measurementNumber) < totalSamplesRead)
                
                % increment the measurement point
                samplePoint = blksize-(totalSamplesRead - settings.measurementPoints(measurementNumber))+1;
                
                % record code and carrier phase
                trackResults(channelNr).codePhase(measurementNumber) = tcode(samplePoint);
                if measurementNumber == 1
                    trackResults(channelNr).deltaCarrPhase(measurementNumber) = rem(trigarg(samplePoint), (2 * pi));
                else
                    trackResults(channelNr).deltaCarrPhase(measurementNumber) = rem(trigarg(samplePoint), (2 * pi)) - trackResults(channelNr).deltaCarrPhase(measurementNumber - 1);
                end
                trackResults(channelNr).carrPhase(measurementNumber) = rem(trigarg(samplePoint), (2 * pi));
                % record millisecond of the measurement
                trackResults(channelNr).msOfMeasurement(measurementNumber) =loopCnt;
                
                % count the carrier cycles to the sample point
                for i = 2:samplePoint
                    if (rem(trigarg(i-1), 2*pi) < pi)&&(rem(trigarg(i), 2*pi)>= pi)
                         carrierCycleCount = carrierCycleCount + 1; 
                    end
                end
                
                % record carrier and code cycle count
                trackResults(channelNr).carrCycleCount(measurementNumber) = carrierCycleCount;
                
                % reset cycyle count 
                carrierCycleCount = 0;
                
                % count to the end of the block
                for i = samplePoint+1:length(trigarg)
                    if (rem(trigarg(i-1), 2*pi) < pi)&&(rem(trigarg(i), 2*pi)>= pi)
                        carrierCycleCount = carrierCycleCount + 1; 
                    end
                end
                measurementNumber = measurementNumber + 1;
            else
                % count the carrier cycles
                for i = 2:length(trigarg)
                    if (rem(trigarg(i-1), 2*pi) < pi)&&(rem(trigarg(i), 2*pi)>= pi)
                        carrierCycleCount = carrierCycleCount + 1; 
                    end
                end
            end
            
%% Record various measures to show in postprocessing ----------------------
            % Record sample number (based on 8bit samples)
            trackResults(channelNr).absoluteSample(loopCnt) = ftell(fid);

            trackResults(channelNr).dllDiscr(loopCnt)       = codeError;
            trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
            trackResults(channelNr).sllDiscr(loopCnt)       = subError;
            trackResults(channelNr).sllDiscrFilt(loopCnt)   = subNco;
            trackResults(channelNr).pllDiscr(loopCnt)       = carrError;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;

            trackResults(channelNr).I_I_E(loopCnt) = I_I_E;
            trackResults(channelNr).I_I_P(loopCnt) = I_I_P;
            trackResults(channelNr).I_E_P(loopCnt) = I_E_P;
            trackResults(channelNr).I_L_P(loopCnt) = I_L_P;
            trackResults(channelNr).I_I_L(loopCnt) = I_I_L;
            
            trackResults(channelNr).Q_I_E(loopCnt) = Q_I_E;
            trackResults(channelNr).Q_I_P(loopCnt) = Q_I_P;
            trackResults(channelNr).Q_E_P(loopCnt) = Q_E_P;
            trackResults(channelNr).Q_L_P(loopCnt) = Q_L_P;
            trackResults(channelNr).Q_I_L(loopCnt) = Q_I_L;
            
        end % for loopCnt       
        
    end % if a PRN is assigned
end % for channelNr 


% Close the waitbar
close(hwb)
