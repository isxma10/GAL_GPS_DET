% Script postProcessing.m processes the raw signal from the specified data
% file (in settings) operating on blocks of 37 seconds of data.
%
% It runs acquisition code identifying the satellites in the file,

disp ('Starting processing...');

[fid, message] = fopen(settings.fileName, 'rb');

    %If success, then process the data
    if (fid > 0)

        % Move the starting point of processing. Can be used to start the
        % signal processing at any point in the data record (e.g. good for long
        % records or for signal processing in blocks).
        fseek(fid, settings.skipNumberOfBytes, 'bof');
        % Find number of samples per spreading code
        samplesPerCode = round(settings.samplingFreq / ...
                           (settings.codeFreqBasis / settings.codeLength));

        % Read data for acquisition. 11ms of signal are needed for the fine
        % frequency estimation
        data = fread(fid, 21*samplesPerCode, settings.dataType)';

        %--- Do the acquisition -------------------------------------------
        disp ('   Acquiring satellites...');
        acqResults = acquisitionGPS(data, settings);
        plotAcquisition(acqResults);
    end 

%% Initialize channels and prepare for the run ============================

    % Start further processing only if a GNSS signal was acquired (the
    % field FREQUENCY will be set to 0 for all not acquired signals)
    if (any(acqResults.carrFreq))
        channel = preRun(acqResults, settings);
        showChannelStatus(channel, settings);
    else
        % No satellites to track, exit
        disp('No GNSS signals detected, signal processing finished.');
        trackResults = [];
        return;
    end