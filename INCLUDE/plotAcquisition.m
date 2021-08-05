function plotAcquisition(acqResults)
%Functions plots bar plot of acquisition results (acquisition metrics). No
%bars are shown for the satellites not included in the acquisition list (in
%structure SETTINGS). 
%
%plotAcquisition(acqResults)
%
%   Inputs:
%       acqResults    - Acquisition results from function acquisition.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis
%--------------------------------------------------------------------------
%% Plot all results =======================================================
figure(101);

hAxes = newplot();

bar(hAxes, acqResults.peakMetric);

title (hAxes, 'Acquisition results');
xlabel(hAxes, 'PRN number (no bar - SV is not in the acquisition list)');
ylabel(hAxes, 'Acquisition Metric');

oldAxis = axis(hAxes);
axis  (hAxes, [0, 33, 0, oldAxis(4)]);
set   (hAxes, 'XMinorTick', 'on');
set   (hAxes, 'YGrid', 'on');

%% Mark acquired signals ==================================================

acquiredSignals = acqResults.peakMetric .* (acqResults.carrFreq > 0);

hold(hAxes, 'on');
bar (hAxes, acquiredSignals, 'FaceColor', [0 0.8 0]);
hold(hAxes, 'off');

legend(hAxes, 'Not acquired signals', 'Acquired signals');
