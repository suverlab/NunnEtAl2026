%% OMN 11/18/25
    % calculates the average Laa, Raa, Lvm, Rvm across trials for each experiment
%% instructions:
% 1) data = load('yourMAT.mat') % make sure you youve labeled ROIs with GUI
% 2) groupedTrials = group_trials_by_date_exp(data)
% 3) RUN THIS FXN with: comboAvg = compute_combo_averages(data,groupedTrials, buffer) 
        % buffer = 20 (amount of samples to plot before and after ROI)
%%

function comboAvg = compute_combo_averages(S, groupedTrials, buffer)
    if nargin < 3
        display('input 3 arguements: data, groupedTrials (run function), and buffer size (20 is default)')
    end
    
    comboAvg = struct([]);
    
    for g = 1:numel(groupedTrials)
        idx = groupedTrials(g).trialIndices;
        comboAvg(g).key = groupedTrials(g).key;
        comboAvg(g).Laa = {};
        comboAvg(g).Raa = {};
        comboAvg(g).Lvm = {};
        comboAvg(g).Rvm = {};

        num_L_ROIs = 0;
        num_R_ROIs = 0;
        
        for t = idx
            Laa = S.trials(t).Laa - S.trials(t).Laa(1);
            Raa = S.trials(t).Raa - S.trials(t).Raa(1);
            Vm  = S.trials(t).Vm  - S.trials(t).Vm(1);
            N = length(Laa);
            
            L_ROIs = S.ROIs_x(t).Laa;
            R_ROIs = S.ROIs_x(t).Raa;
            
            % LEFT ROIs
            for kROI = 1:size(L_ROIs,1)
                pos = L_ROIs(kROI,:);
                if isempty(pos) || numel(pos) ~= 2, continue; end
                x1 = max(1, round(pos(1)) - buffer);
                x2 = min(N, round(pos(2)) + buffer);
                if x2 < x1, continue; end
                comboAvg(g).Laa{end+1} = (Laa(x1:x2) - Laa(x1));
                comboAvg(g).Lvm{end+1} = (Vm(x1:x2)  - Vm(x1));
                %% absolute value
                %comboAvg(g).Laa{end+1} = abs(Laa(x1:x2) - Laa(x1));
                %comboAvg(g).Lvm{end+1} = abs(Vm(x1:x2)  - Vm(x1));
                num_L_ROIs = num_L_ROIs + 1;
            end
            
            % RIGHT ROIs
            for kROI = 1:size(R_ROIs,1)
                pos = R_ROIs(kROI,:);
                if isempty(pos) || numel(pos) ~= 2, continue; end
                x1 = max(1, round(pos(1)) - buffer);
                x2 = min(N, round(pos(2)) + buffer);
                if x2 < x1, continue; end
                comboAvg(g).Raa{end+1} = (Raa(x1:x2) - Raa(x1));
                comboAvg(g).Rvm{end+1} = (Vm(x1:x2)  - Vm(x1));
                %% absolute value
                %comboAvg(g).Raa{end+1} = abs(Raa(x1:x2) - Raa(x1));
                %comboAvg(g).Rvm{end+1} = abs(Vm(x1:x2)  - Vm(x1));
                num_R_ROIs = num_R_ROIs + 1;
            end
        end
        % Store ROI counts for plotting/summary
        comboAvg(g).num_L_ROIs = num_L_ROIs;
        comboAvg(g).num_R_ROIs = num_R_ROIs;
        
        % Compute mean for this combo
        comboAvg(g).Laa_mean = mean_pad(comboAvg(g).Laa);
        comboAvg(g).Raa_mean = mean_pad(comboAvg(g).Raa);
        comboAvg(g).Lvm_mean = mean_pad(comboAvg(g).Lvm);
        comboAvg(g).Rvm_mean = mean_pad(comboAvg(g).Rvm);
    end
end


