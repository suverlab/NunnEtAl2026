%% Interneurons that project to both premotor neurons and APN2

% extract columns from classifications and connections and table into their own array (optimizes indexing and increases speed) -> repeated code so I can run only this section if I want
allRootIDs = classifications.root_id;
allSuperClass = classifications.super_class;
allSides = string(classifications.side);
allPreIDs = connections.pre_root_id;
allPostIDs = connections.post_root_id;
allNTs = string(connections.nt_type);

% extract direct inputs to premotor neurons and their NTs
amnInputIDs = allPreIDs(ismember(allPostIDs,amnIDs));
premInputIDs = allPreIDs(ismember(allPostIDs,amnInputIDs));
premInputNTs = allNTs(ismember(allPreIDs,premInputIDs));

% extract direct inputs to APN2 and find overlap IDs
apnInputIDs = allPreIDs(ismember(allPostIDs,apn2IDs));
inputPremApnInputIDs = unique(premInputIDs(ismember(premInputIDs,apnInputIDs)));

% extract the specific connections for both premotor and APN2
inputsPremApnConns = connections(ismember(allPreIDs,inputPremApnInputIDs) & (ismember(allPostIDs,amnInputIDs) | ismember(allPostIDs,apn2IDs)),:);
inputsPremApnConns = renamevars(inputsPremApnConns,["pre_root_id","post_root_id","nt_type","syn_count"],["preRootID","postRootID","NTtype","numSynapses"]);

% rewrite digraph12n function to have 10 nodes and 2 endCellTypes
function [categorized_pathConnections,directedGraph] = digraph10n(pathConnections,nodeNames,startCellTypeIDs,endCellTypeIDs1,endCellTypeIDs2,codexConnections,codexClassifications)

% adapted from digraph2 function (just has 10 nodes instead of 12)

    % extract columns from codex classifications and connections table into their own array (optimizes indexing and increases speed)
    allRootIDs = codexClassifications.root_id;
    sides = string(codexClassifications.side);
    allPreRootIDs = codexConnections.pre_root_id;
    allNTs = string(codexConnections.nt_type);

    % extract columns from path connections table into their own array (optimizes indexing and increases speed)
    preIDs  = string(pathConnections.preRootID);
    postIDs = string(pathConnections.postRootID);

    % create a container map for the side of the brain that each cell is on
    sideMap = containers.Map(allRootIDs, sides);

    % create new variable in table for cell types
    pathConnections.preClassification = strings(height(pathConnections),1); 
    pathConnections.postClassification = strings(height(pathConnections),1); 

    % find cell types of each connection (based on NT type and side of brain)
    numConnections = height(pathConnections);
    step = max(1, round(numConnections/10)); % guarantees mod function below doesn't break if there are less than 10 connections

    for n = 1:numConnections
        if mod(n, step) == 0
            fprintf('Directed graph is %.0f %% completed\n',(n/numConnections)*100)
        end

        preID  = preIDs(n);
        postID = postIDs(n);
    
        % determine which side of the brain the pre and post nodes are on
        preSideOfBrain  = sideMap(preID);
        postSideOfBrain = sideMap(postID);

        % determine what NT type the current connection is and if its exc/inh
        currConnNT = unique(allNTs(allPreRootIDs==preID)); % NT of current synapse will be based on pre ID (presynaptic neuron)
        nextConnNT = unique(allNTs(allPreRootIDs==postID)); % post ID needs to be classified based on what NT it will SEND, not RECEIVE
        
        if strcmp(currConnNT,'ACH')
            currConnType = 'excitatory';
        elseif strcmp(currConnNT,'GABA') || strcmp(currConnNT,'GLUT')
            currConnType = 'inhibitory';
        else 
            currConnType = 'other';
        end
                
        if isempty(nextConnNT)
            nextConnType = 'no output';
        elseif strcmp(nextConnNT,'ACH')
            nextConnType = 'excitatory';
        elseif strcmp(nextConnNT,'GABA') || strcmp(nextConnNT,'GLUT')
            nextConnType = 'inhibitory';
        else 
            nextConnType = 'other';
        end

        % determine if current pre node is part of the starting cell type and if current post node is part of ending cell type
        startCellTypeCheck = ismember(preID,startCellTypeIDs); 
        endCellTypeCheck1 = ismember(postID,endCellTypeIDs1);
        endCellTypeCheck2 = ismember(postID,endCellTypeIDs2);

        % classify pre ID and store info
        switch preSideOfBrain
            case 'left',   preSideIdx = 1;
            case 'right',  preSideIdx = 2;
        end

        if strcmp(currConnType,'other')
            pathConnections.preClassification(n) = "other";
        else
            if startCellTypeCheck & strcmp(currConnType,'excitatory') 
                offset = 0;
            elseif startCellTypeCheck & strcmp(currConnType,'inhibitory') 
                offset = 2;         
            end
            pathConnections.preClassification(n) = nodeNames(offset + preSideIdx);
        end

        % classify post ID and store info
        switch postSideOfBrain
            case 'left',   postSideIdx = 1;
            case 'right',  postSideIdx = 2;
        end

        if strcmp(nextConnType,'other')
            pathConnections.postClassification(n) = "other";
        else
            if endCellTypeCheck1 & strcmp(currConnType,'excitatory') 
                offset = 4;
            elseif endCellTypeCheck1 & strcmp(currConnType,'inhibitory') 
                offset = 6;
            elseif endCellTypeCheck2
                offset = 8;
            end        
            pathConnections.postClassification(n) = nodeNames(offset + postSideIdx);
        end
  
    end
    fprintf('Directed graph is %.0f %% completed\n',100) % state that the loop is done 

    % delete connections that are other NT types (ex/ SER or OCT) b/c we are only interested in exc vs inh, not neuromodulators
    pathConnections(strcmp(pathConnections.preClassification,'other'),:) = [];
    pathConnections(strcmp(pathConnections.postClassification,'other'),:) = [];  

    % create categorized output variable
    categorized_pathConnections = pathConnections;

    % make directed graph
    tempDirGraph = digraph(categorized_pathConnections.preClassification(:),categorized_pathConnections.postClassification(:),categorized_pathConnections.numSynapses);
    
    % make directed graph
    endNodes = string(tempDirGraph.Edges.EndNodes);

    [G,~,idx] = unique(endNodes,'rows');
    weights = accumarray(idx, tempDirGraph.Edges.Weight);

    directedGraph = digraph(G(:,1), G(:,2), weights);
end

% create directed graph
inputsPremApnNodeNames = ["exc.interneuron_L","exc.interneuron_R","inh.interneuron_L","inh.interneuron_R","exc.premotor_L","exc.premotor_R","inh.premotor_L","inh.premotor_R","APN2_L","APN2_R"];

[inputsPremApnConns,inputsPremApnDirGraph] = digraph10n(inputsPremApnConns,inputsPremApnNodeNames,inputPremApnInputIDs,amnInputIDs,apn2IDs,connections,classifications);

% save directed graph and categorized path connections variable
cd(savedVariablesPath)
save(todaysDate+"_digraph_inputsPremApnPaths.mat","inputsPremApnDirGraph")
save(todaysDate+"_categorizedConnections_premApnAmnConns.mat","inputsPremApnConns")

% plot weighted network graph
inputsPremApnNetworkGraph = plotDigraph(inputsPremApnDirGraph,inputsPremApnConns,inputsPremApnNodeNames,0.005);

% save figure
cd(savedFiguresPath)
saveas(inputsPremApnNetworkGraph,todaysDate+"_inputsPremApnNetworkGraph",'fig')
saveas(inputsPremApnNetworkGraph,todaysDate+"_inputsPremApnNetworkGraph",'png')
saveas(inputsPremApnNetworkGraph,todaysDate+"_inputsPremApnNetworkGraph",'pdf')