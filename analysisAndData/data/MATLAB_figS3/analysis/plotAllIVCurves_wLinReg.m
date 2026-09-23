% plot IV curve for every cell on one plot + linear regression

function [steadyVmAll, deltaVmAll, slopes, rSquared] = plotAllIVCurves_wLinReg(allCellVm, stepCurrents, Fs, fileNames)

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

slopes = nan(numCells,1);
rSquared = nan(numCells,1);

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

    %% Plot raw I-V relationship for each cell
    %plot(stepCurrents, steadyVmAll(cellNum,:),'o-','LineWidth',1.2,'DisplayName',char(fileNames(cellNum)));
    
    %% Linear regression for EACH cell
    p = polyfit(stepCurrents, steadyVmAll(cellNum,:),1);

    slope = p(1);
    intercept = p(2);

    slopes(cellNum) = slope;

    % Predicted values
    fittedVm = polyval(p, stepCurrents);

    % Calculate R^2
    SSres = sum((steadyVmAll(cellNum,:) - fittedVm).^2);
    SStot = sum((steadyVmAll(cellNum,:) - mean(steadyVmAll(cellNum,:))).^2);

    rSquared(cellNum) = 1 - (SSres/SStot);

    % % Overlay regression line for EACH cell
    % plot(stepCurrents, fittedVm,'--','LineWidth',2,'HandleVisibility','off');
    fprintf('%s: slope = %.3f mV/nA, R^2 = %.3f\n', char(fileNames(cellNum)), slope, rSquared(cellNum));
end
%% Mean I-V curve + regression

% Mean I-V curve
meanVm = mean(steadyVmAll,1,'omitnan');
semVm = std(steadyVmAll,[],1,'omitnan') ./ sqrt(numCells);

errorbar(stepCurrents, meanVm, semVm,'ko-','LineWidth',2,'MarkerFaceColor','k','DisplayName','Mean ± SEM');

% Regression through mean
pMean = polyfit(stepCurrents, meanVm,1);
meanFitVm = polyval(pMean, stepCurrents);

% R2 for mean curve
SSres = sum((meanVm - meanFitVm).^2);
SStot = sum((meanVm - mean(meanVm)).^2);

R2_mean = 1 - SSres/SStot;

fprintf('Mean I-V fit R^2 = %.3f\n', R2_mean);

% plot(stepCurrents, meanFitVm, k--','LineWidth',3,'DisplayName','Linear fit');

xlabel('Injected Current (nA)');
ylabel('Membrane Potential (mV)');
title('I-V Curves for All Cells');

legend('off');

ylim([-100 50])
yticks(-100:25:50)
grid off;
box off;

hold off;

exportgraphics(gcf,'best5IVcurves_withRegression.pdf','ContentType','vector','BackgroundColor','white');

end