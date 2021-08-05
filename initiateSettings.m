%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis and Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Adapted by Matthew Alcock and Dr Paul Blunt
% CVS record:
% $Id: init.m,v 1.14.2.21 2006/08/22 13:46:00 dpl Exp $
%--------------------------------------------------------------------------
%% Script initializes settings and environment of the software receiver
%% Clean up the environment first =========================================
clear; close all; clc;

format ('compact');
format ('long', 'g');

%--- Include folders with functions ---------------------------------------
addpath include             % The software receiver functions
addpath geoFunctions        % Position calculation related functions
addpath trackingFunctions   % Galileo and GPS Tracking Functions
%% Print startup ==========================================================
fprintf(['\n',...
    'Welcome to:  softGNSS\n\n', ...
    'An open source GNSS SDR software project initiated by:\n\n', ...
    '              Danish GPS Center/Aalborg University\n\n', ...
    'The code was improved by GNSS Laboratory/University of Colorado.\n\n',...
    'Enhancements to acquisition and navigation and completely new tracking \n',...
    'algorithms have been added by EPFL.\n\n',...
    'The software receiver softGNSS comes with ABSOLUTELY NO WARRANTY;\n',...
    'for details please read license details in the file license.txt. This\n',...
    'is free software, and  you  are  welcome  to  redistribute  it under\n',...
    'the terms described in the license.\n\n']);
fprintf('                   -------------------------------\n\n');

%% Initialize constants, settings =========================================
[settings] = initialSettings();

%% Generate plot of raw data and ask if ready to start processing =========
try
    fprintf('Probing data (%s)...\n', settings.fileName)
    probeData(settings);
catch
    % There was an error, print it and exit
    errStruct = lasterror;
    disp(errStruct.message);
    disp('  (change settings in "initSettingsNSL_26MHz.m" to reconfigure)')    
    return;
end
    
disp('  Raw IF data plotted ')
disp('  (change settings in "initSettingsNSL_26MHz.m" to reconfigure)')
disp(' ');
disp('  Processing is now split into three stages;  Acquisition, Tracking and Navigation')
disp('  Use runAcquisition and runTracking_...(many different types) and runNav to perform each stage ')

