%% Network graph of APN2 outputs to premotor neurons (1-2 hops) where each APN2 neuron is a different node

% -- Extract/compile info and create adjacency matrix -- %

% extract columns from classifications and connections and table into their own array (optimizes indexing and increases speed) -> repeated code so I can run only this section if I want
allSides = string(classifications.side);
allPreIDs = connections.pre_root_id;
allPostIDs = connections.post_root_id;
allNTs = string(connections.nt_type);

% extract premotor neurons (direct inputs to AMN neurons) and their NTs -> repeated code so I can run only this section if I want
pmnIDs = allPreIDs(ismember(allPostIDs,amnIDs));

% index out APN2 to premotor connections (max 2 hops)
apnPmnConns = connections(ismember(allPreIDs,apn2IDs) & ismember(allPostIDs,pmnIDs),:);
apnPmnInIDs = unique(allPreIDs(ismember(allPreIDs,allPostIDs(ismember(allPreIDs,apn2IDs))) & ismember(allPostIDs,pmnIDs)));
% apnPmnConns_2hop = connections(ismember(allPreIDs,apn2IDs) & ismember(allPostIDs,apnPmnInIDs) | ismember(allPreIDs,apnPmnInIDs) & ismember(allPostIDs,pmnIDs),:);
% apnPmnConns = [apnPmnConns_1hop;apnPmnConns_2hop];

% find AMN IDs that are only 3 hops downstream APN2
% pmnDsApnIDs = unique(apnPmnConns.post_root_id(ismember(apnPmnConns.post_root_id,pmnIDs)));
pmnDsApnIDs = unique(apnPmnConns.post_root_id);
pmnAmnDsApnConns = connections(ismember(allPreIDs,pmnDsApnIDs) & ismember(allPostIDs,amnIDs),:);
amnDsApnIDs = unique(pmnAmnDsApnConns.post_root_id);
amnDsApnGroups = amnGroups(ismember(amnGroups.amnID,amnDsApnIDs),:);

% sort APN2 to PMN connections based on which AMN pair they project to
pmnIDs1 = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.AMN,"1")))));
pmnIDs2 = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.AMN,"2")))));
pmnIDs3 = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.AMN,"3")))));
pmnIDs4a = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.AMN,"4a")))));
pmnIDs4b = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.AMN,"4b")))));

apnPmnConns.targetAmnPair(ismember(apnPmnConns.post_root_id,pmnIDs1),1) = "1";
apnPmnConns.targetAmnPair(ismember(apnPmnConns.post_root_id,pmnIDs2),1) = "2";
apnPmnConns.targetAmnPair(ismember(apnPmnConns.post_root_id,pmnIDs3),1) = "3";
apnPmnConns.targetAmnPair(ismember(apnPmnConns.post_root_id,pmnIDs4a),1) = "4a";
apnPmnConns.targetAmnPair(ismember(apnPmnConns.post_root_id,pmnIDs4b),1) = "4b";

% sort APN2 to PMN connections based on which muscle the downstream AMNs project to
PDpmnIDs = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.muscleName,"posterior depressor")))));
ADpmnIDs = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.muscleName,"anterior depressor")))));
PLpmnIDs = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.muscleName,"posterior levator")))));
ALpmnIDs = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.muscleName,"anterior levator")))));

apnPmnConns.targetMuscle(ismember(apnPmnConns.post_root_id,PDpmnIDs),1) = "posterior depressor";
apnPmnConns.targetMuscle(ismember(apnPmnConns.post_root_id,ADpmnIDs),1) = "anterior depressor";
apnPmnConns.targetMuscle(ismember(apnPmnConns.post_root_id,PLpmnIDs),1) = "posterior levator";
apnPmnConns.targetMuscle(ismember(apnPmnConns.post_root_id,ALpmnIDs),1) = "anterior levator";

% sum synapses for identical pre/post pairs
[Gpairs, uniqPreID, uniqPostID] = findgroups(apnPmnConns.pre_root_id, apnPmnConns.post_root_id);
totalSyn = splitapply(@sum, apnPmnConns.syn_count, Gpairs);

% node list and indices
nodeIDs = unique([uniqPreID; uniqPostID],'stable');

[~,preIdx]  = ismember(uniqPreID,nodeIDs);
[~,postIdx] = ismember(uniqPostID,nodeIDs);

% adjacency matrix
apnPmnAdjMatrix = sparse(preIdx, postIdx, totalSyn, ...
           numel(nodeIDs), numel(nodeIDs));

% directed graph
apnPmnDirGraph = digraph(apnPmnAdjMatrix,nodeIDs);



% -- Format network graph and plot -- %

% label output (premotor) nodes based on which AMN pair and muscle group they target downstream + label side of brain of each node
nodeNames = apnPmnDirGraph.Nodes.Name;
for n = 1:size(nodeNames,1)
    currentAmnPair = unique(apnPmnConns.targetAmnPair(apnPmnConns.post_root_id==nodeNames(n)));
    if isempty(currentAmnPair)
        apnPmnDirGraph.Nodes.targetAmnPair(n,1) = "";
    else
        apnPmnDirGraph.Nodes.targetAmnPair(n,1) = currentAmnPair;
    end

    currentMuscle = unique(apnPmnConns.targetMuscle(apnPmnConns.post_root_id==nodeNames(n)));
    if isempty(currentMuscle)
        apnPmnDirGraph.Nodes.targetMuscle(n,1) = "";
    else
        apnPmnDirGraph.Nodes.targetMuscle(n,1) = currentMuscle;
    end

    apnPmnDirGraph.Nodes.side(n,1) = unique(classifications.side(ismember(classifications.root_id,nodeNames(n))));
end

% color output nodes based on NT type
nodeColor = zeros(numnodes(apnPmnDirGraph),3); % default node color = black

ntColors = [0 0 0 % black = excitatory (ACH)
            0.7 0.7 0.7]; % light gray = inhibitory (GLUT/GABA)

outputNodes = unique(string(apnPmnConns.post_root_id));

for o = 1:length(outputNodes)
    output = outputNodes(o);

    nt = unique(string(apnPmnConns.nt_type(apnPmnConns.post_root_id == output)));  % get this input's nt
    idx = strcmp(nodeIDs, output); % find the corresponding node in the graph

    % assign color
    if any(nt == "ACH")
        nodeColor(idx,:) = ntColors(1,:); 
    elseif any(nt == "GLUT") || any(nt == "GABA")
        nodeColor(idx,:) = ntColors(2,:);
    end
end

% make nodes of each APN2 mirror pair a different color and make the associated input edges match
nodeColors = [           % create array of colors thats the size of num mirror pairs
    0.85 0.20 0.55   % pink
    0.10 0.75 0.75   % cyan
    0.55 0.20 0.75   % purple
    0.95 0.40 0.15   % orange
    0.95 0.70 0.10   % yellow
    0.10 0.70 0.20   % green
    0.90 0.10 0.10   % red
    0.10 0.45 0.85   % blue 
];

edgeColor = zeros(numedges(apnPmnDirGraph),3); % default = black

for k = 1:size(apnMirrorPairs,1) 
    pairs = string(apnMirrorPairs(k,:));

    nidx = ismember(nodeIDs,pairs); % find nodes of current pair
    nodeColor(nidx,:) = repmat(nodeColors(k,:),sum(nidx),1); % assign color

    eidx = ismember(string(apnPmnDirGraph.Edges.EndNodes(:,1)),pairs); % find edges whose origin is one of these posts
    edgeColor(eidx,:) = repmat(nodeColors(k,:),sum(eidx),1); % assign color
end

% make only the APN2 nodes large
nodeSize = 5*ones(numnodes(apnPmnDirGraph),1); % set default node size
nodeSize(ismember(nodeNames, apn2IDs)) = 15; % make APN2 nodes larger

% plot
apnPmnNetworkGraph = figure;
plot(apnPmnDirGraph,'Layout','force', ...
    'MarkerSize',nodeSize,'NodeColor',nodeColor,'EdgeColor',edgeColor, ...
    'ArrowSize',0, ...
    'NodeLabel',apnPmnDirGraph.Nodes.targetMuscle+" "+apnPmnDirGraph.Nodes.targetAmnPair+" "+apnPmnDirGraph.Nodes.side);
    

hold on

legendNames = ["Pair A","Pair B","Pair C","Pair D","Pair E","Pair F","Unmatched L","Unmatched R","Excitatory Input","Inhibitory Input"];
colors = [nodeColors;ntColors];

h = gobjects(size(colors,1),1);
for k = 1:size(colors,1)
    h(k) = scatter(nan,nan,100,...
        'filled',...
        'MarkerFaceColor',colors(k,:),...
        'MarkerEdgeColor',colors(k,:));
end

legend(h,legendNames,'Location','eastoutside')
title("Individual APN2 input to premotor neurons (1 hop), categorized by downstream muscle target")

% save figure
cd(savedFiguresPath)
saveas(apnPmnNetworkGraph,todaysDate+"_indivApnToPmnNetworkGraph",'fig')
saveas(apnPmnNetworkGraph,todaysDate+"_indivApnToPmnNetworkGraph",'png')
saveas(apnPmnNetworkGraph,todaysDate+"_indivApnToPmnNetworkGraph",'pdf')