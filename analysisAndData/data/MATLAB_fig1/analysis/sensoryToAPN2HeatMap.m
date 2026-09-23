%% Heatmap of adjacency matrix showing APN2 input connectivity

% extract columns from classifications and connections and table into their own array (optimizes indexing and increases speed) -> repeated code so I can run only this section if I want
allRootIDs = string(classifications.root_id);
allSides = string(classifications.side);
allPreIDs = connections.pre_root_id;
allPostIDs = connections.post_root_id;
allNTs = string(connections.nt_type);

% -- Upload DFS paths data -- % 
cd(savedVariablesPath)
jToApnPathsCSV = readcell("260320_JON_APN2_paths.csv");
bToApnPathsCSV = readcell("260223_BN_APN2_paths.csv");
hToApnPathsCSV = readcell("260730_hBN_APN2_paths.csv");
dToApnPathsCSV = readcell("260730_DNx01_APN2_paths.csv");

% -- Split groups of paths into 1 path/row -- %

% JON
jToApnPaths = {};

numPathGroups = numel(jToApnPathsCSV);
step = max(1, round(numPathGroups/10)); % guarantees mod function below doesn't break if there are less than 10 rows

for p = 1:numPathGroups
    groupedPaths = jToApnPathsCSV{p};
    if strcmp(groupedPaths,'None')
        continue
    end

    if mod(p, step) == 0
        fprintf('%.0f %% of grouped paths have been separated\n',(p/numPathGroups)*100)
    end

    singlePaths = regexp(groupedPaths, '\[([^\[\]]+)\]', 'tokens');

    for k = 1:numel(singlePaths)
        jToApnPaths{end+1,1} = strsplit(singlePaths{k}{1}, ', ');
    end
end 

% ant BN
bToApnPaths = {};

numPathGroups = numel(bToApnPathsCSV);
step = max(1, round(numPathGroups/10)); % guarantees mod function below doesn't break if there are less than 10 rows

for p = 1:numPathGroups
    groupedPaths = bToApnPathsCSV{p};
    if strcmp(groupedPaths,'None')
        continue
    end

    if mod(p, step) == 0
        fprintf('%.0f %% of grouped paths have been separated\n',(p/numPathGroups)*100)
    end

    singlePaths = regexp(groupedPaths, '\[([^\[\]]+)\]', 'tokens');

    for k = 1:numel(singlePaths)
        bToApnPaths{end+1,1} = strsplit(singlePaths{k}{1}, ', ');
    end
end 

% head BN
hToApnPaths = {};

numPathGroups = numel(hToApnPathsCSV);
step = max(1, round(numPathGroups/10)); % guarantees mod function below doesn't break if there are less than 10 rows

for p = 1:numPathGroups
    groupedPaths = hToApnPathsCSV{p};
    if strcmp(groupedPaths,'None')
        continue
    end

    if mod(p, step) == 0
        fprintf('%.0f %% of grouped paths have been separated\n',(p/numPathGroups)*100)
    end

    singlePaths = regexp(groupedPaths, '\[([^\[\]]+)\]', 'tokens');

    for k = 1:numel(singlePaths)
        hToApnPaths{end+1,1} = strsplit(singlePaths{k}{1}, ', ');
    end
end 

% DNx01
dToApnPaths = {};

numPathGroups = numel(dToApnPathsCSV);
step = max(1, round(numPathGroups/10)); % guarantees mod function below doesn't break if there are less than 10 rows

for p = 1:numPathGroups
    groupedPaths = dToApnPathsCSV{p};
    if strcmp(groupedPaths,'None')
        continue
    end

    if mod(p, step) == 0
        fprintf('%.0f %% of grouped paths have been separated\n',(p/numPathGroups)*100)
    end

    singlePaths = regexp(groupedPaths, '\[([^\[\]]+)\]', 'tokens');

    for k = 1:numel(singlePaths)
        dToApnPaths{end+1,1} = strsplit(singlePaths{k}{1}, ', ');
    end
end 

% -- Compile path connections table for head bristles and DNx01 -- %
hbnApnPathConns = splitPaths(hbnApnPathsCSV_fileName,savedVariablesPath,connections); 
dnxApnPathConns = splitPaths(dnxApnPathsCSV_fileName,savedVariablesPath,connections); 

% -- write function for calculate total number of synapses along paths -- %
function totalSyn = pathSynapses(path,synLookup)

totalSyn = 0;

for k = 1:numel(path)-1
    key = char(strcat(string(path{k}),"->",string(path{k+1})));
    if isKey(synLookup,key)
        totalSyn = totalSyn + synLookup(key);
    end
end
end

% -- preallocate vars -- %
jApnInputID(size(jToApnPaths,1),1) = string;
jDsApnID(size(jToApnPaths,1),1) = string;
jToApnNumHops(size(jToApnPaths,1),1) = nan;
jToApnNumSyn(size(jToApnPaths,1),1) = nan;

bApnInputID(size(bToApnPaths,1),1) = string; 
bDsApnID(size(bToApnPaths,1),1) = string;
bToApnNumHops(size(bToApnPaths,1),1) = nan;
bToApnNumSyn(size(bToApnPaths,1),1) = nan;

hApnInputID(size(hToApnPaths,1),1) = string; 
hDsApnID(size(hToApnPaths,1),1) = string;
hToApnNumHops(size(hToApnPaths,1),1) = nan;
hToApnNumSyn(size(hToApnPaths,1),1) = nan;

dApnInputID(size(dToApnPaths,1),1) = string; 
dDsApnID(size(dToApnPaths,1),1) = string;
dToApnNumHops(size(dToApnPaths,1),1) = nan;
dToApnNumSyn(size(dToApnPaths,1),1) = nan;

% -- ant BN lookup table -- %

edgeKey = strcat(bnApnPathConns.preRootID,"->",bnApnPathConns.postRootID);

% sum duplicate edges
[G,key] = findgroups(edgeKey);
syn = splitapply(@sum,bnApnPathConns.numSynapses,G);

% build lookup table
synLookup = containers.Map(cellstr(key),syn);

for n = 1:length(bToApnPaths)

    path = bToApnPaths{n};

    bApnInputID(n,1) = path{1};
    bDsApnID(n,1) = path{end};
    bToApnNumSyn(n,1) = pathSynapses(path,synLookup);
    bToApnNumHops(n,1) = numel(path)-1;
end

% -- JON lookup table -- %

jEdgeKey = strcat(jonApnPathConns.preRootID,"->",jonApnPathConns.postRootID);

% sum duplicate edges
[G,key] = findgroups(jEdgeKey);
jSyn = splitapply(@sum,jonApnPathConns.numSynapses,G);

% build lookup table
jSynLookup = containers.Map(cellstr(key),jSyn);

for n = 1:length(jToApnPaths)

    path = jToApnPaths{n};

    if ismember(path{1},bApnInputID)
        continue
    else
        jApnInputID(n,1) = path{1};
        jDsApnID(n,1) = path{end};
        jToApnNumSyn(n,1) = pathSynapses(path,jSynLookup);
        jToApnNumHops(n,1) = numel(path)-1;
    end
end

% -- head BN lookup table -- %

hEdgeKey = strcat(hbnApnPathConns.preRootID,"->",hbnApnPathConns.postRootID);

% sum duplicate edges
[G,key] = findgroups(hEdgeKey);
hSyn = splitapply(@sum,hbnApnPathConns.numSynapses,G);

% build lookup table
hSynLookup = containers.Map(cellstr(key),hSyn);

for n = 1:length(hToApnPaths)

    path = hToApnPaths{n};

    if ismember(path{1},bApnInputID)
        continue
    else
        hApnInputID(n,1) = path{1};
        hDsApnID(n,1) = path{end};
        hToApnNumSyn(n,1) = pathSynapses(path,hSynLookup);
        hToApnNumHops(n,1) = numel(path)-1;
    end
end

% -- DNx01 lookup table -- %

dEdgeKey = strcat(dnxApnPathConns.preRootID,"->",dnxApnPathConns.postRootID);

% sum duplicate edges
[G,key] = findgroups(dEdgeKey);
dSyn = splitapply(@sum,dnxApnPathConns.numSynapses,G);

% build lookup table
dSynLookup = containers.Map(cellstr(key),dSyn);

for n = 1:length(dToApnPaths)

    path = dToApnPaths{n};

    if ismember(path{1},bApnInputID)
        continue
    else
        dApnInputID(n,1) = path{1};
        dDsApnID(n,1) = path{end};
        dToApnNumSyn(n,1) = pathSynapses(path,dSynLookup);
        dToApnNumHops(n,1) = numel(path)-1;
    end
end

% -- concatanate all input vars -- %
sensInputID = [jApnInputID;bApnInputID;hApnInputID;dApnInputID];
dsApnID = [jDsApnID;bDsApnID;hDsApnID;dDsApnID];
numHops = [jToApnNumHops;bToApnNumHops;hToApnNumHops;dToApnNumHops];
numSyn = [jToApnNumSyn;bToApnNumSyn;hToApnNumSyn;dToApnNumSyn];

results = table(sensInputID,dsApnID,numSyn,numHops);
results = rmmissing(results);

% -- organize the data so they can be plotted in the desired order -- %

% sort APN2 columns
apnLeft  = apn2IDs(ismember(apn2IDs,...
    allRootIDs(strcmp(allSides,"left"))));

apnRight = apn2IDs(ismember(apn2IDs,...
    allRootIDs(strcmp(allSides,"right"))));

apnOrder = [apnLeft; apnRight];

% sort sensory inputs
jLeft = jApnInputID(ismember(jApnInputID,...
    allRootIDs(strcmp(allSides,"left"))));

bLeft = bApnInputID(ismember(bApnInputID,...
    allRootIDs(strcmp(allSides,"left"))));

hLeft = hApnInputID(ismember(hApnInputID,...
    allRootIDs(strcmp(allSides,"left"))));

dLeft = dApnInputID(ismember(dApnInputID,...
    allRootIDs(strcmp(allSides,"left"))));

jRight = jApnInputID(ismember(jApnInputID,...
    allRootIDs(strcmp(allSides,"right"))));

bRight = bApnInputID(ismember(bApnInputID,...
    allRootIDs(strcmp(allSides,"right"))));

hRight = hApnInputID(ismember(hApnInputID,...
    allRootIDs(strcmp(allSides,"right"))));

dRight = dApnInputID(ismember(dApnInputID,...
    allRootIDs(strcmp(allSides,"right"))));

inputOrder = [unique(jLeft);
                unique(bLeft);
                unique(hLeft);
                unique(dLeft);
                unique(jRight);
                unique(bRight);
                unique(hRight);
                unique(dRight)];

% count the size of each group
njLeft = numel(unique(jLeft));
nbLeft = numel(unique(bLeft));
nhLeft = numel(unique(hLeft));
ndLeft = numel(unique(dLeft));

njRight = numel(unique(jRight));
nbRight = numel(unique(bRight));
nhRight = numel(unique(hRight));
ndRight = numel(unique(dRight));

groupSizes = [njLeft nbLeft nhLeft ndLeft ...
              njRight nbRight nhRight ndRight];

nApnLeft = numel(apnLeft);
nApnRight = numel(apnRight);

% -- create matrices -- %

[G,sInputID,apnID] = findgroups(results.sensInputID,results.dsApnID);

minHop = splitapply(@min,results.numHops,G);

results.group = G;
results.minHop = minHop(G);

shortestResults = results(results.numHops == results.minHop,:);

[G2,sInputID,apnID] = findgroups(shortestResults.sensInputID,...
                                shortestResults.dsApnID);

hop = splitapply(@min,shortestResults.numHops,G2);
syn = splitapply(@sum,shortestResults.numSyn,G2);

synMat1 = nan(length(inputOrder),length(apnOrder));
synMat2 = nan(length(inputOrder),length(apnOrder));

for i = 1:length(sInputID)

    r = find(inputOrder == sInputID(i));
    c = find(apnOrder == apnID(i));

    if isempty(r) || isempty(c)
        continue
    end

    if hop(i) == 1
        synMat1(r,c) = syn(i);

    elseif hop(i) == 2
        synMat2(r,c) = syn(i);

    end
end

synMat1(isnan(synMat1)) = 0;
synMat2(isnan(synMat2)) = 0;

% -- plot 1-hop heatmap -- %
sensToApnHeatmap_1hop = figure;
imagesc(synMat1)
colormap(parula)
colorbar
clim([0 25])
title("Sensory -> APN2 (1 hop)",'FontSize',14)

hold on

yBreaks = cumsum(groupSizes);
for k = 1:length(yBreaks)-1
    yline(yBreaks(k)+0.5,'k','LineWidth',1.5)
end
xline(nApnLeft+0.5,'k','LineWidth',2)

% add group labels
groupLabels = {...
    'JON Left',...
    'Ant BN Left',...
    'Head BN Left',...
    'DNx01 Left',...
    'JON Right',...
    'Ant BN Right',...
    'Head BN Right',...
    'DNx01 Right'};

groupCenters = cumsum(groupSizes) - groupSizes/2; % calc center of each group

yticks(groupCenters)
yticklabels(groupLabels)

set(gca,'XTick',1:1:14,'XTickLabel',apnOrder,'YAxisLocation','left');

text(nApnLeft/2, -30, 'Left APN2', ...
    'HorizontalAlignment','center','FontSize',12)

text(nApnLeft+nApnRight/2, -30, 'Right APN2', ...
    'HorizontalAlignment','center','FontSize',12)

% -- plot 2-hops heatmap -- %
sensToApnHeatmap_2hop = figure;
imagesc(synMat2)
colormap(parula)
colorbar
clim([0 25])
title("Sensory -> APN2 (2 hops)",'FontSize',14)

hold on

yBreaks = cumsum(groupSizes);
for k = 1:length(yBreaks)-1
    yline(yBreaks(k)+0.5,'k','LineWidth',1.5)
end
xline(nApnLeft+0.5,'k','LineWidth',2)

% add group labels
groupLabels = {...
    'JON Left',...
    'Ant BN Left',...
    'Head BN Left',...
    'DNx01 Left',...
    'JON Right',...
    'Ant BN Right',...
    'Head BN Right',...
    'DNx01 Right'};

groupCenters = cumsum(groupSizes) - groupSizes/2; % calc center of each group

yticks(groupCenters)
yticklabels(groupLabels)

set(gca,'XTick',1:1:14,'XTickLabel',apnOrder,'YAxisLocation','left');

text(nApnLeft/2, -30, 'Left APN2', ...
    'HorizontalAlignment','center','FontSize',12)

text(nApnLeft+nApnRight/2, -30, 'Right APN2', ...
    'HorizontalAlignment','center','FontSize',12)

% save figures
cd(savedFiguresPath)

saveas(sensToApnHeatmap_1hop,todaysDate+"_sensInputToApnConnMatrix_1hop",'fig')
saveas(sensToApnHeatmap_1hop,todaysDate+"_sensInputToApnConnMatrix_1hop",'png')
saveas(sensToApnHeatmap_1hop,todaysDate+"_sensInputToApnConnMatrix_1hop",'pdf')

saveas(sensToApnHeatmap_2hop,todaysDate+"_sensInputToApnConnMatrix_2hop",'fig')
saveas(sensToApnHeatmap_2hop,todaysDate+"_sensInputToApnConnMatrix_2hop",'png')
saveas(sensToApnHeatmap_2hop,todaysDate+"_sensInputToApnConnMatrix_2hop",'pdf')