%%
% adapted from digraph12n (4 nodes instead of 12)
% used for analyzing direct inhibitory inputs to APN2 that are independent of both JONs and BNs
% only classfies by left and right sides
% doesn't classify by NT type because we are only looking at inhibitory connections
% since there are no interneurons (just looking at direct connections) there are only start and end cell types (4 nodes)

%%
function [categorized_pathConnections,directedGraph] = digraph4n(pathConnections,nodeNames,codexClassifications)

    % extract columns from codex classifications and connections table into their own array (optimizes indexing and increases speed)
    allRootIDs = codexClassifications.root_id;
    sides = string(codexClassifications.side);

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

        % classify pre and post ID and store info
        switch preSideOfBrain
            case 'left',   preSideIdx = 1;
            case 'right',  preSideIdx = 2;
        end
     
        switch postSideOfBrain
            case 'left',   postSideIdx = 3;
            case 'right',  postSideIdx = 4;
        end
       
        pathConnections.preClassification(n) = nodeNames(preSideIdx);
        pathConnections.postClassification(n) = nodeNames(postSideIdx);


    end
    fprintf('Directed graph is %.0f %% completed\n',100) % state that the loop is done  

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