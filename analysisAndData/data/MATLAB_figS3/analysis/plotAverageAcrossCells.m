%% plot average current steps across all cells

function plotAverageAcrossCells(groupAvgVm, allCellVm, t, stepCurrents, numTrialsPerCell)

% Calculate sample sizes
numCells = size(allCellVm,3);
numSteps = length(stepCurrents);
totalTrials = sum(numTrialsPerCell);

% Create one color for each current step
colors = lines(numSteps);

figure
hold on

for step = 1:numSteps

    plot(t, groupAvgVm(:,step),'Color', colors(step,:),'LineWidth', 1.5);

end

xlabel('Time (s)')
ylabel('Vm (mV)')

title(sprintf('Average Current Steps Across Cells (N = %d cells, n = %d trials)', numCells, totalTrials))

legend(compose('%g pA', stepCurrents * 1000), 'Location', 'eastoutside')

box off
hold off

end