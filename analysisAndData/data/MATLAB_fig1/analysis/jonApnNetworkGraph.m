%% JON to APN2 - split paths

jonApnPathConns = splitPaths(jonApnPathsCSV_fileName,savedVariablesPath,connections); 

% save the paths connections table as a variable so you dont have to keep running this (it takes several minutes)
save(todaysDate+"_connections_jonApnPaths.mat","jonApnPathConns")

%% JON to APN2 - directed graph
if loadVars == 1
    cd(savedVariablesPath)
    if ~isempty(jonApnPathConns_fileName)
        load(jonApnPathConns_fileName)
    elseif isempty(jonApnPathConns_fileName) & ~exist("jonApnPathConns",'var')
        error("necessary variable does not exist. please run previous section or upload proper variable.")
    end
end

jonApnNodeNames = ["JON_L","JON_R","JON_C", ...
                    "Exc.Interneuron_L","Exc.Interneuron_R","Exc.Interneuron_C", ...
                    "Inh.Interneuron_L","Inh.Interneuron_R","Inh.Interneuron_C", ...
                    "APN2_L","APN2_R","APN2_C"];

[jonApnPathConns,jonApnDirGraph] = digraph12n(jonApnPathConns,jonApnNodeNames,jonIDs,apn2IDs,connections,classifications);

% save directed graph and categorized path connections variable
save(todaysDate+"_digraph_jonApnPaths.mat","jonApnDirGraph")
save(todaysDate+"_categorizedConnections_jonApnPaths.mat","jonApnPathConns")

%% JON to APN2 - plot weighted network graph
if loadVars == 1
    cd(savedVariablesPath)
    load(categorized_jonApnPathConns_fileName)
    load(jonApnDirGraph_fileName)
end

jonApnNodeNames = ["JON_L","JON_R","JON_C", ...
                    "Exc.Interneuron_L","Exc.Interneuron_R","Exc.Interneuron_C", ...
                    "Inh.Interneuron_L","Inh.Interneuron_R","Inh.Interneuron_C", ...
                    "APN2_L","APN2_R","APN2_C"];

jonApnNetworkGraph = plotDigraph(jonApnDirGraph,jonApnPathConns,jonApnNodeNames,0.0005);

% save figure
cd(savedFiguresPath)
saveas(jonApnNetworkGraph,todaysDate+"_jonApnNetworkGraph",'fig')
saveas(jonApnNetworkGraph,todaysDate+"_jonApnNetworkGraph",'png')
saveas(jonApnNetworkGraph,todaysDate+"_jonApnNetworkGraph",'pdf')