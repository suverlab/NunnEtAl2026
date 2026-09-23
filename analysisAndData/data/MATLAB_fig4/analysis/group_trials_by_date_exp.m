function groupedTrials = group_trials_by_date_exp(S)
    % Convert trial identifiers
    dates = string({S.trials.date});
    expNums = [S.trials.expNum];
    
    % Combine date + expNum to form a unique key
    keys = strcat(dates, "_", string(expNums));
    uniqueKeys = unique(keys);
    
    % Initialize struct
    groupedTrials = struct([]);
    
    for k = 1:numel(uniqueKeys)
        groupedTrials(k).key = uniqueKeys(k);
        groupedTrials(k).trialIndices = find(keys == uniqueKeys(k));
    end
end
