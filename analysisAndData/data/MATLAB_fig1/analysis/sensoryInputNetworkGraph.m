%% Network graph of sensory inputs to APN2 where each APN2 neuron is a different node

if loadVars == 1
    cd(savedVariablesPath)
    load(categorized_jonApnPathConns_fileName)
    load(categorized_bnApnPathConns_fileName)
end

% -- Extract/compile info and create adjacency matrix -- %

% extract columns from classifications and connections and tables into their own array (optimizes indexing and increases speed) -> repeated code so I can run only this section if I want
allRootIDs = classifications.root_id;
allSides = string(classifications.side);
allSuperClass = classifications.super_class;

jonApnPreIDs  = string(jonApnPathConns.preRootID);
jonApnPostIDs = string(jonApnPathConns.postRootID);
jonApnNTs = jonApnPathConns.NTtype;

bnApnPreIDs  = string(bnApnPathConns.preRootID);
bnApnPostIDs = string(bnApnPathConns.postRootID);
bnApnNTs = string(bnApnPathConns.NTtype);

% extract direct inputs to APN2 and their NTs
apnAmnInputID = jonApnPreIDs(ismember(jonApnPostIDs,apn2IDs));
bApnInputID = bnApnPreIDs(ismember(bnApnPostIDs,apn2IDs));

sensInputIDs = unique([apnAmnInputID;bApnInputID]);

% index out sensory inputs to APN2 connections -> some repeated code so I can run only this section if I want
apnInputConns  = connections(ismember(connections.post_root_id,apn2IDs),:); % extract APN2 input connections
sensInputConns = apnInputConns(ismember(apnInputConns.pre_root_id,sensInputIDs),:);


% sum synapses for identical pre/post pairs
[Gpairs, uniqPreID, uniqPostID] = findgroups(sensInputConns.pre_root_id, sensInputConns.post_root_id);
totalSyn = splitapply(@sum, sensInputConns.syn_count, Gpairs);

% node list and indices
nodeIDs = unique([uniqPreID; uniqPostID],'stable');

[~,preIdx]  = ismember(uniqPreID,nodeIDs);
[~,postIdx] = ismember(uniqPostID,nodeIDs);

% adjacency matrix
sensApnInputAdjMatrix = sparse(preIdx, postIdx, totalSyn, ...
           numel(nodeIDs), numel(nodeIDs));

% directed graph
sensInputDirGraph = digraph(sensApnInputAdjMatrix,nodeIDs);


% -- Format network graph and plot -- %

nodeColor = zeros(numnodes(sensInputDirGraph),3); % default node color = black

% color input nodes based on NT type
ntColors = [0 0 0 % black = excitatory (ACH)
            % 196/255 164/255 132/255]; % light brown
            % 150/255 121/255 105/255]; % mocha
            0.6 0.6 0.6]; % gray = inhibitory (GLUT/GABA)

inputNodes = unique(string(sensInputConns.pre_root_id));

for i = 1:length(inputNodes)
    output = inputNodes(i);

    nt = unique(string(sensInputConns.nt_type(sensInputConns.pre_root_id == output)));  % get this input's nt
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

edgeColor = zeros(numedges(sensInputDirGraph),3); % default = black

for k = 1:size(apnMirrorPairs,1) 
    pairs = string(apnMirrorPairs(k,:));

    nidx = ismember(nodeIDs,pairs); % find nodes of current pair
    nodeColor(nidx,:) = repmat(nodeColors(k,:),sum(nidx),1); % assign color

    eidx = ismember(string(sensInputDirGraph.Edges.EndNodes(:,2)),pairs); % find edges whose destination is one of these posts
    edgeColor(eidx,:) = repmat(nodeColors(k,:),sum(eidx),1); % assign color
end

% make only the APN2 nodes large
nodeNames = sensInputDirGraph.Nodes.Name;
apnNodes = unique(string(sensInputConns.post_root_id)); % unique APN2 neurons
nodeSize = 2*ones(numnodes(sensInputDirGraph),1); % set default node size
nodeSize(ismember(nodeNames, apnNodes)) = 8; % make APN2 nodes larger

% plot
sensInputNetworkGraph = figure;
plot(sensInputDirGraph,'Layout','force','MarkerSize',nodeSize,'NodeColor',nodeColor,'EdgeColor',edgeColor,'EdgeAlpha',0.4,'ArrowSize',0);

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
title("Direct sensory inputs to individual APN2 neurons")

% save figure
cd(savedFiguresPath)
saveas(sensInputNetworkGraph,todaysDate+"_sensInputToIndivApnNetworkGraph",'fig')
saveas(sensInputNetworkGraph,todaysDate+"_sensInputToIndivApnNetworkGraph",'png')
saveas(sensInputNetworkGraph,todaysDate+"_sensInputToIndivApnNetworkGraph",'pdf')

%%
nodes = table2array(sensInputDirGraph.Nodes);
test = pdist(nodes);
stest = squareform(nodes);
