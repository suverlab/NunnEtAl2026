%%
% adapted from digraph2 function (6 nodes instead of 12)
% 2 endCellTypes
% specifically used to analyze premotor inputs to both APN2 and AMN
% as well as cells that input onto both premotor neurons and APN2

%%
function [categorized_pathConnections,directedGraph] = digraph8n(pathConnections,nodeNames,startCellTypeIDs,endCellTypeIDs1,endCellTypeIDs2,codexConnections,codexClassifications)

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

        if strcmp(currConnNT,'ACH')
            currConnType = 'excitatory';
        elseif strcmp(currConnNT,'GABA') || strcmp(currConnNT,'GLUT')
            currConnType = 'inhibitory';
        else 
            currConnType = 'other';
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

        if endCellTypeCheck1
            offset = 4;
        elseif endCellTypeCheck2
            offset = 6;
        end        
        pathConnections.postClassification(n) = nodeNames(offset + postSideIdx);
  
    end
    fprintf('Directed graph is %.0f %% completed\n',100) % state that the loop is done 

    % delete connections that are other NT types (ex/ SER or OCT) b/c we are only interested in exc vs inh, not neuromodulators
    pathConnections(strcmp(pathConnections.preClassification,'other'),:) = [];

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