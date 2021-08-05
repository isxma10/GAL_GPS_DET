function L5codesTable = makeL5Table(settings)
%Function generates L5 codes for all 37 satellites based on the settings
%provided in the structure "settings". The codes are digitized at the
%sampling frequency specified in the settings structure.

%L5codesTable = makeL5Table(settings)
%
%   Inputs:
%       settings        - receiver settings
%   Outputs:
%       L5codesTable    - an array of arrays (matrix) containing L5 codes
%                       for all satellite PRN-s

% load L5I_codes.mat L5I_codes;
load L5I_codesFromSDR L5I_codes;
load L5Q_codesFromSDR L5Q_codes;
%--- Find number of samples per spreading code ----------------------------
samplesPerCode = round(settings.samplingFreq / ...
                           (settings.codeFreqBasis / settings.codeLength));

%--- Prepare the output matrix to speed up function -----------------------
L5codesTable = zeros(37, samplesPerCode);
 
%--- Find time constants --------------------------------------------------
ts = 1/settings.samplingFreq;   % Sampling period in sec
tc = 1/settings.codeFreqBasis;  % L5 chip period in sec
 
%=== For all satellite PRN-s ...
for PRN = 1:37
    %--- Generate L5 code for given PRN -----------------------------------
%     L5code = L5I_codes(PRN,:);
    L5code = L5Q_codes(PRN,:);
    %=== Digitizing =======================================================
    
    %--- Make index array to read L5 code values -------------------------
    % The length of the index array depends on the sampling frequency -
    % number of samples per millisecond (because one L5 code period is one
    % millisecond).
    codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);
    
    %--- Correct the last index (due to number rounding issues) -----------
    codeValueIndex(end) = 10230;
    
    %--- Make the digitized version of the L5 code -----------------------
    % The "upsampled" code is made by selecting values form the L5 code
    % chip array (L5code) for the time instances of each sample.
    L5codesTable(PRN, :) = L5code(codeValueIndex);
    
end % for PRN = 1:37
