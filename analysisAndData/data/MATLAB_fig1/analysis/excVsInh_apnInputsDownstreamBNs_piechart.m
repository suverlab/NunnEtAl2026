%% NT composition of direct inputs to APN2 that are downstream BNs
% use the path connections table to extract these inputs

if loadVars == 1
    cd(savedVariablesPath)
    if ~isempty(bnApnPathConns_fileName)
        load(bnApnPathConns_fileName)
    elseif isempty(bnApnPathConns_fileName) & ~exist("bnApnPathConns",'var')
        error("necessary variable does not exist. please run previous section or upload proper variable.")
    end
end

% extract columns from path connections table into their own array (optimizes indexing and increases speed)
bnApnPreIDs  = string(bnApnPathConns.preRootID);
bnApnPostIDs = string(bnApnPathConns.postRootID);
bnApnNTs = string(bnApnPathConns.NTtype);

% extract direct inputs to APN2 and their NTs
bApnInputIDs = bnApnPreIDs(ismember(bnApnPostIDs,apn2IDs));
bApnInputNTs = bnApnNTs(ismember(bnApnPostIDs,apn2IDs));

bApnInputNtMap = containers.Map(bApnInputIDs,bApnInputNTs);

% categorize synapse types
bApnInputSynType(size(bApnInputIDs,1)) = string;
for n = 1:size(bApnInputIDs,1)  
   currNt = bApnInputNtMap(bApnInputIDs(n));
    if strcmp(currNt,'ACH')
        bApnInputSynType(n) = 'excitatory';
    elseif strcmp(currNt,'GABA') | strcmp(currNt,'GLUT')
        bApnInputSynType(n) = 'inhibitory';
    else 
        bApnInputSynType(n) = 'other';
    end
end

% count number of inputs in each group
numExcInputs = length(find(strcmp(bApnInputSynType,"excitatory")));
numInhInputs = length(find(strcmp(bApnInputSynType,"inhibitory")));

% plot
excVsInhPiechart = figure;
piechart([numExcInputs,numInhInputs],{'Excitatory','Inhibitory'})
title("Direct inputs to APN2 that are downstream BNs")

% save figure
cd(savedFiguresPath)
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputsDowstreamBNs_piechart",'fig')
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputsDowstreamBNs_piechart",'png')
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputsDowstreamBNs_piechart",'pdf')