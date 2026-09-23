%% Network graph of inhibitory inputs to APN2 that are independent of both JONs and BNs

% extract columns from classifications and connections and table into their own array (optimizes indexing and increases speed) -> repeated code so I can run only this section if I want
allRootIDs = classifications.root_id;
allSuperClass = classifications.super_class;
allPreIDs = connections.pre_root_id;
allPostIDs = connections.post_root_id;
allNTs = string(connections.nt_type);

% index out JON-independent and BN-independent inhibitory APN2 inputs and their cell classes -> some repeated code so I can run only this section if I want
inhInputIDs = allPreIDs(ismember(allPostIDs,apn2IDs) & (strcmp(allNTs,"GLUT") | strcmp(allNTs,"GABA")));

jDepInputIDs = inhInputIDs(ismember(inhInputIDs,jApnInputIDs));
jDepInputClass = string(allSuperClass(ismember(allRootIDs,jDepInputIDs)));
jIndInputIDs = inhInputIDs(~ismember(inhInputIDs,jDepInputIDs));
jIndInputClass = string(allSuperClass(ismember(allRootIDs,jIndInputIDs)));

bDepInputIDs = inhInputIDs(ismember(inhInputIDs,jApnInputIDs));
bDepInputClass = string(allSuperClass(ismember(allRootIDs,bDepInputIDs)));
bIndInputIDs = inhInputIDs(~ismember(inhInputIDs,bDepInputIDs));
bIndInputClass = string(allSuperClass(ismember(allRootIDs,bIndInputIDs)));

allIndInputIDs = [jIndInputIDs;bIndInputIDs];

% extract the specific connections for all independent inhibitory inputs
indInputConns = connections(ismember(allPreIDs,allIndInputIDs) & (ismember(allPostIDs,amnIDs) | ismember(allPostIDs,apn2IDs)),:);
indInputConns = renamevars(indInputConns,["pre_root_id","post_root_id","nt_type","syn_count"],["preRootID","postRootID","NTtype","numSynapses"]);

% create directed graph
indInputNodeNames = ["inh.indInput_L","inh.indInput_R","APN2_L","APN2_R"];

[indInputConns,indInputDirGraph] = digraph4n(indInputConns,indInputNodeNames,allIndInputIDs,apn2IDs,classifications);

% save directed graph and categorized path connections variable
cd(savedVariablesPath)
save(todaysDate+"_digraph_indInputConns.mat","indInputDirGraph")
save(todaysDate+"_categorizedConnections_indInputConns.mat","indInputConns")

% plot weighted network graph
indInputNetworkGraph = plotDigraph(indInputDirGraph,indInputConns,indInputNodeNames,0.25);

% save figure
cd(savedFiguresPath)
saveas(indInputNetworkGraph,todaysDate+"_indInputNetworkGraph",'fig')
saveas(indInputNetworkGraph,todaysDate+"_indInputNetworkGraph",'png')
saveas(indInputNetworkGraph,todaysDate+"_indInputNetworkGraph",'pdf')