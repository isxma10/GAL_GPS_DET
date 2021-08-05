multipath = load('trackResultsMulitpath_Amp_0p5_Phase_0_60s.mat');
noMultipath = load('trackingResultsNoMultpath.mat');

chipMetres = 299792458/1.023e6;

codePhaseDiff = chipMetres.*(noMultipath.trackResults.codePhase - multipath.trackResults.codePhase);

figure(69)
plot(codePhaseDiff)

figure(70)
plot(noMultipath.trackResults.CNo)
hold on
plot(multipath.trackResults.CNo)
