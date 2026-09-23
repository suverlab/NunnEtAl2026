
%% NT composition of all direct inputs to APN2 
% produces a pie chart of excitatory vs. inhibitory
% uses the connections table to extract these inputs

%%
% extract APN2 input connections
apnInputConns  = connections(ismember(connections.post_root_id,apn2IDs),:);

% extract columns from connections table into their own array (optimizes indexing and increases speed)
apnInputIDs = string(apnInputConns.pre_root_id);
apnInputNTs = apnInputConns.nt_type;

apnInputNtMap = containers.Map(apnInputIDs,apnInputNTs);

% categorize synapse types
apnInputSynType(size(apnInputIDs,1)) = string;
for n = 1:size(apnInputIDs,1)  
   currNt = apnInputNtMap(apnInputIDs(n));
    if strcmp(currNt,'ACH')
        apnInputSynType(n) = 'excitatory';
    elseif strcmp(currNt,'GABA') | strcmp(currNt,'GLUT')
        apnInputSynType(n) = 'inhibitory';
    else 
        apnInputSynType(n) = 'other';
    end
end

% count number of inputs in each group
numExcInputs = length(find(strcmp(apnInputSynType,"excitatory")));
numInhInputs = length(find(strcmp(apnInputSynType,"inhibitory")));

% plot
excVsInhPiechart = figure;
piechart([numExcInputs,numInhInputs],{'Excitatory','Inhibitory'})
title("Direct inputs to APN2")

cd(savedFiguresPath)
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputs_piechart",'fig')
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputs_piechart",'png')
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputs_piechart",'pdf')