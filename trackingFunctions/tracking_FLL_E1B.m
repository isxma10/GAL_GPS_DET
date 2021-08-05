function [trackResults, channel]= tracking_FLL_E1B(fid, channel, settings)
% Performs code and carrier tracking for all channels.
%
%[trackResults, channel] = tracking_FLL_2ndOrderUnknownData(fid, channel, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
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

%CVS record:
%$Id: tracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

%% Initialize result structure ============================================

% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock

% The absolute sample in the record of the C/A code start:
trackResults.absoluteSample = zeros(1, settings.msToProcess);

% Freq of the C/A code:
trackResults.codeFreq       = inf(1, settings.msToProcess);

% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, settings.msToProcess);

% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, settings.msToProcess);
trackResults.I_E            = zeros(1, settings.msToProcess);
trackResults.I_L            = zeros(1, settings.msToProcess);

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_E            = zeros(1, settings.msToProcess);
trackResults.Q_P            = zeros(1, settings.msToProcess);
trackResults.Q_L            = zeros(1, settings.msToProcess);

% Loop discriminators
trackResults.dllDiscr       = inf(1, settings.msToProcess);
trackResults.dllDiscrFilt   = inf(1, settings.msToProcess);
trackResults.pllDiscr       = inf(1, settings.msToProcess);
trackResults.pllDiscrFilt   = inf(1, settings.msToProcess);
trackResults.carrierDoppler = inf(1, settings.msToProcess);
trackResults.codeDoppler   = inf(1, settings.msToProcess);
trackResults.I_P_D          = zeros(1, settings.msToProcess);
trackResults.Q_P_D          = zeros(1, settings.msToProcess);
trackResults.CdLi           = zeros(1, settings.msToProcess);
trackResults.CrLi           = zeros(1, settings.msToProcess);
trackResults.CNo            = zeros(1, settings.msToProcess);
trackResults.NavBits        = zeros(1, settings.msToProcess);
trackResults.bitSyncCnt     = zeros(1, 20);
trackResults.bitSync        = zeros(1, settings.msToProcess);
trackResults.remCodePhase   = zeros(1, settings.msToProcess);

trackResults.codePhase = zeros(1, length(settings.measurementPoints));
trackResults.deltaCarrPhase = zeros(1, length(settings.measurementPoints));
trackResults.carrCycleCount = zeros(1, length(settings.measurementPoints));
trackResults.msOfMeasurement = zeros(1, length(settings.measurementPoints));
trackResults.carrPhase = zeros(1, length(settings.measurementPoints));

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================

codePeriods = settings.msToProcess/4;     % For GPS one C/A code is one ms

hwb = waitbar(0,'Tracking...');

%% Start processing channels ==============================================
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
              settings.skipNumberOfBytes + channel(channelNr).codePhase-1, ...
              'bof');


        % Get a vector with the C/A code sampled 1x/chip
        E1BCode = generateE1Bcode(channel(channelNr).PRN,1);
        % Then make it possible to do early and late versions
        E1BCode = [E1BCode(8183) E1BCode(8184) E1BCode E1BCode(1) E1BCode(2)];

        %--- Perform various initializations ------------------------------

        % define initial code frequency basis of NCO
        codeFreq      = settings.codeFreqBasis;
        % define residual code phase (in chips)
        remCodePhase  = 0.0;
        % define carrier frequency which is used over whole tracking period
        carrFreq      = channel(channelNr).acquiredFreq;
        carrFreqBasis = channel(channelNr).acquiredFreq;
        % define residual carrier phase
        remCarrPhase  = 0.0;

        %carrier/Costas loop parameters
        I_P_D = 0;
        Q_P_D = 0;
        I_E_S = 0;
        Q_E_S = 0;
        I_P_S = 0;
        Q_P_S = 0;
        I_L_S = 0;
        Q_L_S = 0;
        carrNco = 0;
        codeLockInd = 0;
        carrLockInd = 0;
        NavBit = 1;
        accmCount = 1;
        
        % Define early-late offset (in chips)
        earlyLateSpc = 0.5/3;

        % Summation interval
        coherentAccmNum = 1;
        PDIcode = ((settings.codeLength/ settings.codeLengthCA)/1000)*coherentAccmNum;
        PDIcarr = ((settings.codeLength/ settings.codeLengthCA)/1000)*coherentAccmNum;
        %bit sync counters
        bitSync = 0;
        bitSyncCnt = zeros(1, 20);
 
        % inital Doppler shift estimate
        deltaFcaMinus_k = 0;
        % inital rate of change of Doppler shift estimate
        deltaFcaDotMinus_k = 0;
        % Carrier frequency locked loop bandwidth
        Bl_ca = 5;
        % FLL gains
        Kcf1 = 3.4 *(Bl_ca*PDIcarr);
        Kcf2 = 2.04 *((Bl_ca*PDIcarr)^2);
        % code frequency
        Fco = settings.codeFreqBasis;
        % intial code Doppler shift estimate
        deltaFco_k = settings.spectrumInversion*carrFreq/1540;
        TstMinus_k = (((Fco + deltaFco_k)/Fco)*PDIcode);
        Bl_co = 1;
        Kco = 4 * Bl_co * PDIcode;
        
        % Initial values for Cno estimation
        Pw = 0;
        % set to loss of lock threshold
        cdLi = 10;
        CNo = 10*log10(cdLi);
        
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
                            '; Completed ',int2str(loopCnt), ...
                            ' of ', int2str(codePeriods), ' msec']);                       
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
                        
            blksize = ceil((settings.codeLengthE1B-remCodePhase) / codePhaseStep);
            
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
            tcode       = (remCodePhase-earlyLateSpc) : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
            tcode2      = ceil(tcode*2) + 2;
            earlyCode   = E1BCode(tcode2);
            
            % Define index into late code vector
            tcode       = (remCodePhase+earlyLateSpc) : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
            tcode2      = ceil(tcode*2) + 2;
            lateCode    = E1BCode(tcode2);
            
            % Define index into prompt code vector
            tcode       = remCodePhase : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase);
            tcode2      = ceil(tcode*2) + 2;
            promptCode  = E1BCode(tcode2);
            
            remCodePhase = (tcode(blksize) + codePhaseStep) - 4092.0;

%% Generate the carrier frequency to mix the signal to baseband -----------
            time    = (0:blksize) ./ settings.samplingFreq;
            
            % Get the argument to sin/cos functions
            trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
            remCarrPhase = rem(trigarg(blksize+1), (2 * pi));
            
            % Finally compute the signal to mix the collected data to bandband
            carrCos = cos(trigarg(1:blksize));
            carrSin = sin(trigarg(1:blksize));

%% Generate the six standard accumulated values ---------------------------
            % First mix to baseband
            qBasebandSignal = carrCos .* rawSignal;
            iBasebandSignal = carrSin .* rawSignal;

            % Now get early, late, and prompt values for each
            
            % accumulate 1ms
            I_E = sum(earlyCode  .* iBasebandSignal);
            Q_E = sum(earlyCode  .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);
            Q_P = sum(promptCode .* qBasebandSignal);
            I_L = sum(lateCode   .* iBasebandSignal);
            Q_L = sum(lateCode   .* qBasebandSignal);
            
            
            % accumulate coherently further
            I_E_S = I_E_S + I_E;
            Q_E_S = Q_E_S + Q_E;
            I_P_S = I_P_S + I_P;
            Q_P_S = Q_P_S + Q_P;
            I_L_S = I_L_S + I_L;
            Q_L_S = Q_L_S + Q_L;
            
                    
            
%% update when accumulations are ready ----------------------------------
            % Only update at the end of accumulation
            if accmCount == coherentAccmNum   
%% Carrier to noise estimation 
                % Square up the correlations
                I2_plus_Q2 = (I_P_S * I_P_S) + (Q_P_S * Q_P_S);
                % Compute the narrow band power 
                Pn = I2_plus_Q2;
                
                % Accumulate the wide band power
                Pw = Pw + (I_P * I_P) + (Q_P * Q_P);
                
%                 if (coherentAccmNum == 20)
%                    % Use a moving average
%                    cdLi = cdLi + (((Pn/Pw) - cdLi)/100);
%                    cnoEst = (coherentAccmNum*(cdLi - 1))/(PDIcode*(coherentAccmNum - cdLi));
%                    CNo = 10*log10(cnoEst);
%                 end
%% Find PLL error and update carrier NCO         
                    
                % Carrier frequency discriminator
                useFilteredDenominator = 1;
                
                cross = (Q_P_S * I_P_D) - (I_P_S * Q_P_D);
                
                dot = (I_P_S * I_P_D) + (Q_P_S * Q_P_D);
                
                 if (dot == 0)
                    F = 0;
                else
                    F = atan(cross / dot);
                end	
                
                if(dot < 0)

                    NavBit = -NavBit;
                    cross = -cross;
                    dot = -dot;
                end

                carrLockInd = carrLockInd + (dot - carrLockInd)/(3000/coherentAccmNum);
                codeLockInd = codeLockInd + (I2_plus_Q2 - codeLockInd)/(3000/coherentAccmNum);

                if (useFilteredDenominator == 1)
                    if (carrLockInd == 0)
                        F = 0;
                    else
                        F = atan(cross / carrLockInd);
                    end
                end
                 
                % Carrier frequency discriminator normalisation function
                Ncod = (1/(2*pi*PDIcarr));

                % Carrier frequency error
                deltaFca_k = Ncod * F;

                % Doppler shift
                deltaFcaPlus_k = deltaFcaMinus_k + (Kcf1*deltaFca_k);

                % Rate of change of Doppler shift
                deltaFcaDotPlus_k = deltaFcaDotMinus_k + ((Kcf2/PDIcarr)*deltaFca_k); 

                % Doppler shift predicted forward to next iteration
                deltaFcaMinus_kplus1 = deltaFcaPlus_k + (deltaFcaDotPlus_k*PDIcarr);

                % Rate of change of Doppler shift predicted forward to next iteration
                deltaFcaDotMinus_kplus1 = deltaFcaDotPlus_k;

                % NCO update
                carrNco = deltaFcaMinus_kplus1;

                % Store Doppler shift for next iteration
                deltaFcaMinus_k = deltaFcaMinus_kplus1;

                % Store rate of change of Doppler for next iteration
                deltaFcaDotMinus_k = deltaFcaDotMinus_kplus1;
                
                I_P_D = I_P_S;
                Q_P_D = Q_P_S;
             
%% Find DLL error and update code NCO -------------------------------------
                
                %Code discriminator normalisation function
                Nele = sqrt(codeLockInd);
                %Code phase error
                if Nele == 0
                   X_k = 0; 
                else
                   X_k = (sqrt(I_E_S * I_E_S + Q_E_S * Q_E_S) - sqrt(I_L_S * I_L_S + Q_L_S * Q_L_S)) / Nele;
                end
                
                % code phase estimate
                TstPlus_k = TstMinus_k - ((Kco * X_k)/Fco);
                
                % Calculate code Doppler shift from the carrier
                deltaFco_k = settings.spectrumInversion*((carrFreq - settings.IF)/1540);
                 
                % code phase estimate predicted forward
                TstMinus_kPlus1 = TstPlus_k + (((Fco + deltaFco_k)/Fco)*PDIcode); 
                
                % set the NCO frequency
                codeNcoGroves = (((TstMinus_kPlus1 - TstMinus_k)/PDIcode)*Fco);
                
                % store the code Doppler
                codeDoppler = TstMinus_kPlus1;
                
                % store code phase for next iteration
                TstMinus_k = TstMinus_kPlus1;
                                
                % restart the accumulations
                Pw = 0;
                I_E_S = 0;
                Q_E_S = 0;
                I_P_S = 0;
                Q_P_S = 0;
                I_L_S = 0;
                Q_L_S = 0;
                
                accmCount = 1;
            else
                % increment the accumulation counter 
                accmCount = accmCount + 1;
                
                % Accumulate the wide band power
                Pw = Pw + (I_P * I_P) + (Q_P * Q_P);
            end 
                         
                %%Update the NCOs      
                carrFreq = carrFreqBasis + carrNco;
                
                codeFreq = codeNcoGroves;
           
%% read measurement data on the same sample ---------------------------------

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
                
                % reset cycle count 
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
            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq;
            trackResults(channelNr).absoluteSample(loopCnt) = ftell(fid);

            trackResults(channelNr).dllDiscr(loopCnt)       = X_k;
            trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNcoGroves - settings.codeFreqBasis;
            trackResults(channelNr).pllDiscr(loopCnt)       = deltaFca_k;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;
            trackResults(channelNr).carrierDoppler(loopCnt) = deltaFcaPlus_k;
            trackResults(channelNr).codeDoppler(loopCnt) = codeDoppler;
            
            
            trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
            trackResults(channelNr).I_L(loopCnt) = I_L;
            trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            trackResults(channelNr).Q_L(loopCnt) = Q_L;
            trackResults(channelNr).I_P_D(loopCnt) = I_P_D;
            trackResults(channelNr).Q_P_D(loopCnt) = Q_P_D;
            trackResults(channelNr).CdLi(loopCnt) = codeLockInd;
            trackResults(channelNr).CrLi(loopCnt) = carrLockInd;
            trackResults(channelNr).CNo(loopCnt) = CNo;
            trackResults(channelNr).NavBits(loopCnt) = NavBit;
            trackResults(channelNr).bitSync(loopCnt) = bitSync;
            trackResults(channelNr).bitSyncCnt = bitSyncCnt;
            trackResults(channelNr).remCodePhase(loopCnt) = remCodePhase;
        
        end % for loopCnt
        trackResults(channelNr).status  = 'T';
        % C/No based lock detector
%         if CNo > 6.0
%             trackResults(channelNr).status  = 'T';
%         else
%             trackResults(channelNr).status  = '-';
%         end
%         
    end % if a PRN is assigned
end % for channelNr 

% Close the waitbar
close(hwb)
