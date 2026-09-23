%% Heatmap of adjacency matrix showing APN2 -> AMN connectivity

% extract columns from classifications and connections and table into their own array (optimizes indexing and increases speed) -> repeated code so I can run only this section if I want
allRootIDs = string(classifications.root_id);
allSides = string(classifications.side);
allPreIDs = connections.pre_root_id;
allPostIDs = connections.post_root_id;
allNTs = string(connections.nt_type);

% -- Upload DFS paths data and path connections -- % 
cd(savedVariablesPath)
apnToAmnPathsCSV = readcell("260223_APN2_AMN_paths.csv");
load("260324_categorizedConnections_apnAmnPaths");

% -- Split groups of paths into 1 path/row -- %

apnToAmnPaths = {};

numPathGroups = numel(apnToAmnPathsCSV);
step = max(1, round(numPathGroups/10)); % guarantees mod function below doesn't break if there are less than 10 rows

for p = 1:numPathGroups
    groupedPaths = apnToAmnPathsCSV{p};
    if strcmp(groupedPaths,'None')
        continue
    end

    if mod(p, step) == 0
        fprintf('%.0f %% of grouped paths have been separated\n',(p/numPathGroups)*100)
    end

    singlePaths = regexp(groupedPaths, '\[([^\[\]]+)\]', 'tokens');

    for k = 1:numel(singlePaths)
        apnToAmnPaths{end+1,1} = strsplit(singlePaths{k}{1}, ', ');
    end
end 

% -- preallocate vars -- %
apnAmnInputID(size(apnToAmnPaths,1),1) = string; % IDs of AMN inputs that are along the APN -> AMN path
dsAmnID(size(apnToAmnPaths,1),1) = string; % IDs of AMNs that are downstream APN
apnAmnNumHops(size(apnToAmnPaths,1),1) = nan;
apnAmnNumSyn(size(apnToAmnPaths,1),1) = nan;

% -- lookup table -- %

edgeKey = strcat(apnAmnPathConns.preRootID,"->",apnAmnPathConns.postRootID);

% sum duplicate edges
[G,key] = findgroups(edgeKey);
syn = splitapply(@sum,apnAmnPathConns.numSynapses,G);

% build lookup table
synLookup = containers.Map(cellstr(key),syn);

for n = 1:length(apnToAmnPaths)

    path = apnToAmnPaths{n};

    if numel(path)==1
        disp(path)
    end

    apnAmnInputID(n,1) = path{1};
    dsAmnID(n,1) = path{end};
    apnAmnNumSyn(n,1) = pathSynapses(path,synLookup);
    apnAmnNumHops(n,1) = numel(path)-1;
end

apnAmnResults = table(apnAmnInputID,dsAmnID,apnAmnNumSyn,apnAmnNumHops);
apnAmnResults = rmmissing(apnAmnResults);

% -- organize the data so they can be plotted in the desired order -- %

% sort APN2 columns by side
apnLeft  = apn2IDs(ismember(apn2IDs,...
    allRootIDs(strcmp(allSides,"left"))));

apnRight = apn2IDs(ismember(apn2IDs,...
    allRootIDs(strcmp(allSides,"right"))));

apnOrder = [apnLeft; apnRight];

% sort AMN by pair and side
ALL1 = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"left"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"anterior levator") & strcmp(amnGroups.AMN,"1"))));

ADL2 = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"left"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"anterior depressor") & strcmp(amnGroups.AMN,"2"))));

PDL3 = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"left"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"posterior depressor") & strcmp(amnGroups.AMN,"3"))));

PLL4a = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"left"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"posterior levator") & strcmp(amnGroups.AMN,"4a"))));

PLL4b = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"left"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"posterior levator") & strcmp(amnGroups.AMN,"4b"))));

ALR1 = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"right"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"anterior levator") & strcmp(amnGroups.AMN,"1"))));

ADR2 = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"right"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"anterior depressor") & strcmp(amnGroups.AMN,"2"))));

PDR3 = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"right"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"posterior depressor") & strcmp(amnGroups.AMN,"3"))));

PLR4a = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"right"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"posterior levator") & strcmp(amnGroups.AMN,"4a"))));

PLR4b = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"right"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"posterior levator") & strcmp(amnGroups.AMN,"4b"))));

amnOrder = [unique(ALL1);
                unique(ADL2);
                unique(PDL3);
                unique(PLL4a);
                unique(PLL4b);
                unique(ALR1);
                unique(ADR2);
                unique(PDR3);
                unique(PLR4a);
                unique(PLR4b)];

% count the size of each group
nALL1 = numel(unique(ALL1));
nADL2 = numel(unique(ADL2));
nPDL3 = numel(unique(PDL3));
nPLL4a = numel(unique(PLL4a));
nPLL4b = numel(unique(PLL4b));

nALR1 = numel(unique(ALR1));
nADR2 = numel(unique(ADR2));
nPDR3 = numel(unique(PDR3));
nPLR4a = numel(unique(PLR4a));
nPLR4b = numel(unique(PLR4b));

groupSizes = [nALL1 nADL2 nPDL3 nPLL4a nPLL4b ...
              nALR1 nADR2 nPDR3 nPLR4a nPLR4b];

nApnLeft = numel(apnLeft);
nApnRight = numel(apnRight);

% -- create matrices -- %

[G,apnID,amnID] = findgroups(apnAmnResults.apnAmnInputID,apnAmnResults.dsAmnID);

minHop = splitapply(@min,apnAmnResults.apnAmnNumHops,G);

apnAmnResults.group = G;
apnAmnResults.minHop = minHop(G);

shortestResults = apnAmnResults(apnAmnResults.apnAmnNumHops == apnAmnResults.minHop,:);

[G2,apnID,amnID] = findgroups(shortestResults.apnAmnInputID,...
                                shortestResults.dsAmnID);

hop = splitapply(@min,shortestResults.apnAmnNumHops,G2);
syn = splitapply(@sum,shortestResults.apnAmnNumSyn,G2);

synMat1 = nan(length(amnOrder),length(apnOrder));
synMat2 = nan(length(amnOrder),length(apnOrder));

for i = 1:length(apnID)

    r = find(amnOrder == amnID(i));
    c = find(apnOrder == apnID(i));

    if isempty(r) || isempty(c)
        continue
    end

    if hop(i) == 2
        synMat1(r,c) = syn(i);

    elseif hop(i) == 3
        synMat2(r,c) = syn(i);

    end
end

synMat1(isnan(synMat1)) = 0;
synMat2(isnan(synMat2)) = 0;

% -- plot 2-hop heatmap -- %
apnToAmnHeatmap_2hop = figure;
imagesc(synMat1)
colormap(parula)
colorbar
clim([0 50])
title("APN2 -> AMN (2 hops)",'FontSize',14)

hold on

yBreaks = cumsum(groupSizes);
for k = 1:length(yBreaks)-1
    yline(yBreaks(k)+0.5,'k','LineWidth',1.5)
end
xline(nApnLeft+0.5,'k','LineWidth',2)

% add group labels
groupLabels = {...
    'AMN1,Left,AL',...
    'AMN2,Left,AD',...
    'AMN3,Left,PD',...
    'AMN4a,Left,PL',...
    'AMN4b,Left,PL',...
    'AMN1,Right,AL',...
    'AMN2,Right,AD',...
    'AMN3,Right,PD',...
    'AMN4a,Right,PL',...
    'AMN4b,Right,PL'};

groupCenters = cumsum(groupSizes) - groupSizes/2; % calc center of each group

yticks(groupCenters)
yticklabels(groupLabels)

% apnLabels = {"A,left","B,left","C,left","D,left","E,left","F,left","unmatched,left","A,right","B,right","C,right","D,right","E,right","F,right","unmatched,right"};

set(gca,'XTick',1:1:14,'XTickLabel',apnOrder,'YAxisLocation','left');

text(nApnLeft/2, 0.25, 'Left APN2', ...
    'HorizontalAlignment','center','FontSize',12)

text(nApnLeft+nApnRight/2, 0.25, 'Right APN2', ...
    'HorizontalAlignment','center','FontSize',12)

% -- plot 3-hops heatmap -- %
apnToAmnHeatmap_3hop = figure;
imagesc(synMat2)
colormap(parula)
colorbar
clim([0 500])
title("APN2 -> AMN (3 hops)",'FontSize',14)

hold on

yBreaks = cumsum(groupSizes);
for k = 1:length(yBreaks)-1
    yline(yBreaks(k)+0.5,'k','LineWidth',1.5)
end
xline(nApnLeft+0.5,'k','LineWidth',2)

% add group labels
groupLabels = {...
    'AMN1,Left,AL',...
    'AMN2,Left,AD',...
    'AMN3,Left,PD',...
    'AMN4a,Left,PL',...
    'AMN4b,Left,PL',...
    'AMN1,Right,AL',...
    'AMN2,Right,AD',...
    'AMN3,Right,PD',...
    'AMN4a,Right,PL',...
    'AMN4b,Right,PL'};

groupCenters = cumsum(groupSizes) - groupSizes/2; % calc center of each group

yticks(groupCenters)
yticklabels(groupLabels)

set(gca,'XTick',1:1:14,'XTickLabel',apnOrder,'YAxisLocation','left');

text(nApnLeft/2, 0.25, 'Left APN2', ...
    'HorizontalAlignment','center','FontSize',12)

text(nApnLeft+nApnRight/2, 0.25, 'Right APN2', ...
    'HorizontalAlignment','center','FontSize',12)

% save figures
cd(savedFiguresPath)

saveas(apnToAmnHeatmap_2hop,todaysDate+"_apnToAmnConnMatrix_2hop",'fig')
saveas(apnToAmnHeatmap_2hop,todaysDate+"_apnToAmnConnMatrix_2hop",'png')
saveas(apnToAmnHeatmap_2hop,todaysDate+"_apnToAmnConnMatrix_2hop",'pdf')

saveas(apnToAmnHeatmap_3hop,todaysDate+"_apnToAmnConnMatrix_3hop",'fig')
saveas(apnToAmnHeatmap_3hop,todaysDate+"_apnToAmnConnMatrix_3hop",'png')
saveas(apnToAmnHeatmap_3hop,todaysDate+"_apnToAmnConnMatrix_3hop",'pdf')
