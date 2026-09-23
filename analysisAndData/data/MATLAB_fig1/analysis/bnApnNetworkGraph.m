%% BN to APN2 - split paths

bnApnPathConns = splitPaths(bnApnPathsCSV_fileName,savedVariablesPath,connections); 

% save the paths connections table as a variable so you dont have to keep running this (it takes several minutes)
save(todaysDate+"_connections_bnApnPaths.mat","bnApnPathConns")

%% BN to APN2 - directed graph
if loadVars == 1
    cd(savedVariablesPath)
    if ~isempty(bnApnPathConns_fileName)
        load(bnApnPathConns_fileName)
    elseif isempty(bnApnPathConns_fileName) & ~exist("bnApnPathConns",'var')
        error("necessary variable does not exist. please run previous section or upload proper variable.")
    end
end

bnApnNodeNames = ["BN_L","BN_R","BN_C", ...
                    "Exc.Interneuron_L","Exc.Interneuron_R","Exc.Interneuron_C", ...
                    "Inh.Interneuron_L","Inh.Interneuron_R","Inh.Interneuron_C", ...
                    "APN2_L","APN2_R","APN2_C"];

[bnApnPathConns,bnApnDirGraph] = digraph12n(bnApnPathConns,bnApnNodeNames,bnIDs,apn2IDs,connections,classifications);

% save directed graph and categorized path connections variable
save(todaysDate+"_digraph_bnApnPaths.mat","bnApnDirGraph")
save(todaysDate+"_categorizedConnections_bnApnPaths.mat","bnApnPathConns")

%% BN to APN2 - plot weighted network graph
if loadVars == 1
    cd(savedVariablesPath)
    load(categorized_bnApnPathConns_fileName)
    load(bnApnDirGraph_fileName)
end

bnApnNodeNames = ["BN_L","BN_R","BN_C", ...
                    "Exc.Interneuron_L","Exc.Interneuron_R","Exc.Interneuron_C", ...
                    "Inh.Interneuron_L","Inh.Interneuron_R","Inh.Interneuron_C", ...
                    "APN2_L","APN2_R","APN2_C"];

bnApnNetworkGraph = plotDigraph(bnApnDirGraph,bnApnPathConns,bnApnNodeNames,0.005);

% save figure
cd(savedFiguresPath)
saveas(bnApnNetworkGraph,todaysDate+"_bnApnNetworkGraph",'fig')
saveas(bnApnNetworkGraph,todaysDate+"_bnApnNetworkGraph",'png')
saveas(bnApnNetworkGraph,todaysDate+"_bnApnNetworkGraph",'pdf')