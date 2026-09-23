%% APN2 to AMN - split paths

apnAmnPathConns = splitPaths(apnAmnPathsCSV_fileName,savedVariablesPath,connections);

% save the paths connections table as a variable so you dont have to keep running this (it takes several minutes)
save(todaysDate+"_connections_apnAmnPaths.mat","apnAmnPathConns")

%% APN2 to AMN - directed graph
if loadVars == 1
    cd(savedVariablesPath)
    if ~isempty(apnAmnPathConns_fileName)
        load(apnAmnPathConns_fileName)
    elseif isempty(apnAmnPathConns_fileName) & ~exist("apnAmnPathConns",'var')
        error("necessary variable does not exist. please run previous section or upload proper variable.")
    end
end

apnAmnNodeNames = ["APN2_L","APN2_R","APN2_C", ...
                    "Exc.Interneuron_L","Exc.Interneuron_R","Exc.Interneuron_C", ...
                    "Inh.Interneuron_L","Inh.Interneuron_R","Inh.Interneuron_C", ...
                    "AMN_L","AMN_R","AMN_C"];

[apnAmnPathConns,apnAmnDirGraph] = digraph12n(apnAmnPathConns,apnAmnNodeNames,apn2IDs,amnIDs,connections,classifications);

% save directed graph and categorized path connections variable
save(todaysDate+"_digraph_apnAmnPaths.mat","apnAmnDirGraph")
save(todaysDate+"_categorizedConnections_apnAmnPaths.mat","apnAmnPathConns")

%% APN2 to AMN - plot weighted network graph
if loadVars == 1
    cd(savedVariablesPath)
    load(categorized_apnAmnPathConns_fileName)
    load(apnAmnDirGraph_fileName)
end

apnAmnNodeNames = ["APN2_L","APN2_R","APN2_C", ...
                    "Exc.Interneuron_L","Exc.Interneuron_R","Exc.Interneuron_C", ...
                    "Inh.Interneuron_L","Inh.Interneuron_R","Inh.Interneuron_C", ...
                    "AMN_L","AMN_R","AMN_C"];

apnAmnNetworkGraph = plotDigraph(apnAmnDirGraph,apnAmnPathConns,apnAmnNodeNames);

% save figure
cd(savedFiguresPath)
saveas(apnAmnNetworkGraph,todaysDate+"_apnAmnNetworkGraph",'fig')
saveas(apnAmnNetworkGraph,todaysDate+"_apnAmnNetworkGraph",'png')
saveas(apnAmnNetworkGraph,todaysDate+"_apnAmnNetworkGraph",'pdf')