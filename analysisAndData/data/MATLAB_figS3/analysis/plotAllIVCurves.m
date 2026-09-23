% plot IV curve for every cell on one plot

function [steadyVmAll, deltaVmAll] = plotAllIVCurves(allCellVm, stepCurrents, Fs, fileNames)

baselineDuration = 0.5;   % seconds
stepDuration     = 0.5;   % seconds
analysisDuration = 0.1;   % final 100 ms

baselineSamples = round(baselineDuration * Fs);
stepSamples     = round(stepDuration * Fs);
analysisSamples = round(analysisDuration * Fs);

numSteps = numel(stepCurrents);
numCells = size(allCellVm,3);

steadyVmAll = nan(numCells, numSteps);
deltaVmAll  = nan(numCells, numSteps);

baselineStart = baselineSamples - analysisSamples + 1;
baselineEnd   = baselineSamples;

steadyStart = baselineSamples + stepSamples - analysisSamples + 1;
steadyEnd   = baselineSamples + stepSamples;

figure;
hold on;

for cellNum = 1:numCells

    for stepNum = 1:numSteps

        trace = allCellVm(:,stepNum,cellNum);

        baselineVm = mean(trace(baselineStart:baselineEnd),'omitnan');

        steadyVmAll(cellNum,stepNum) = mean(trace(steadyStart:steadyEnd), 'omitnan');

        deltaVmAll(cellNum,stepNum) = steadyVmAll(cellNum,stepNum) - baselineVm;

    end

    plot(stepCurrents, steadyVmAll(cellNum,:),'o-','LineWidth', 1.2, 'DisplayName', char(fileNames(cellNum)));

end

xlabel('Injected Current (nA)');
ylabel('Membrane Potential (mV)');
title('I-V Curves for All Cells');

legend('Location','best', 'Interpreter','none');
yticks(-100:25:50)
grid off;
box off;
hold off;

exportgraphics(gcf,'best5IVcurves.pdf','ContentType','vector','BackgroundColor','white');

end