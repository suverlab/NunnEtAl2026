%%
% plots the directed graph with the desired parameters to yield a weighted network graph 
% that contains neuron counts, synapse counts, and shows synaptic weights all on the graph.
%%
% 'arguments' are the inputs to the function. 
% If the argument doesn't equal anything, you must give the function an input 
% If the argument is set equal to something, that argument is optional. 
% If the user puts nothing in the funciton for that input, it will default to the value assigned to it in the arguments section. 
% If the user writes something different in the function inputs, it will use whatever value the user gave

%%
function networkGraph = plotDigraph(directedGraph,categorized_pathConnections,nodeNames,edgeWeightScale,nodeSize,arrowSize)
%   creates weighted network graph

    arguments 
        directedGraph
        categorized_pathConnections
        nodeNames
        edgeWeightScale = 0.005
        nodeSize = 20
        arrowSize = 20
    end

    numNeuronsPerNode = table;
    numNeuronsPerNode.nodeNames = nodeNames';
    numNeuronsPerNode.numNeurons = nan(size(nodeNames,2),1);
    for n = 1:size(nodeNames,2)
        currNodeName = nodeNames(n);
        tempGroupTable = categorized_pathConnections(categorized_pathConnections.preClassification==currNodeName | categorized_pathConnections.postClassification==currNodeName,:);
        [~,uniqueCellIDs,~] = groupcounts([tempGroupTable.preRootID(tempGroupTable.preClassification==currNodeName);tempGroupTable.postRootID(tempGroupTable.postClassification==currNodeName)]);
        numNeuronsPerNode.numNeurons(n) = size(uniqueCellIDs,1);
    end
    
    for n = 1:height(directedGraph.Nodes)
         directedGraph.Nodes.numNeurons(n) = numNeuronsPerNode.numNeurons(numNeuronsPerNode.nodeNames==directedGraph.Nodes.Name(n));
    end
    
    networkGraph = figure;
    plot(directedGraph,'Layout','layered', ...
        'LineWidth',directedGraph.Edges.Weight*edgeWeightScale, ...
        'ArrowSize',arrowSize, ...
        'MarkerSize',nodeSize, ...
        'NodeLabel',string(directedGraph.Nodes.Name)+" ("+string(directedGraph.Nodes.numNeurons)+" neurons)", ...
        'NodeFontWeight','bold', ...
        'EdgeLabel',string(directedGraph.Edges.Weight))
end