%% NT composition of direct inputs to APN2 that are downstream JONs
% use the path connections table to extract these inputs

if loadVars == 1
    cd(savedVariablesPath)
    if ~isempty(jonApnPathConns_fileName)
        load(jonApnPathConns_fileName)
    elseif isempty(jonApnPathConns_fileName) & ~exist("jonApnPathConns",'var')
        error("necessary variable does not exist. please run previous section or upload proper variable.")
    end
end

% extract columns from path connections table into their own array (optimizes indexing and increases speed)
jonApnPreIDs  = string(jonApnPathConns.preRootID);
jonApnPostIDs = string(jonApnPathConns.postRootID);
jonApnNTs = jonApnPathConns.NTtype;

% extract direct inputs to APN2 and their NTs
jApnInputIDs = jonApnPreIDs(ismember(jonApnPostIDs,apn2IDs));
jApnInputNTs = jonApnNTs(ismember(jonApnPostIDs,apn2IDs));

jApnInputNtMap = containers.Map(jApnInputIDs,jApnInputNTs);

% categorize synapse types
jApnInputSynType(size(jApnInputIDs,1)) = string;
for n = 1:size(jApnInputIDs,1)  
   currNt = jApnInputNtMap(jApnInputIDs(n));
    if strcmp(currNt,'ACH')
        jApnInputSynType(n) = 'excitatory';
    elseif strcmp(currNt,'GABA') | strcmp(currNt,'GLUT')
        jApnInputSynType(n) = 'inhibitory';
    else 
        jApnInputSynType(n) = 'other';
    end
end

% count number of inputs in each group
numExcInputs = length(find(strcmp(jApnInputSynType,"excitatory")));
numInhInputs = length(find(strcmp(jApnInputSynType,"inhibitory")));

% plot
excVsInhPiechart = figure;
piechart([numExcInputs,numInhInputs],{'Excitatory','Inhibitory'})
title("Direct inputs to APN2 that are downstream JONs")

cd(savedFiguresPath)
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputsDowstreamJons_piechart",'fig')
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputsDowstreamJons_piechart",'png')
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputsDowstreamJons_piechart",'pdf')