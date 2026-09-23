%%
% There is a built in digraph() function in Matlab. 
% However, there is some pre-processing of the pathConnections table that needs to be done, 
% as well as some post-processing, to get correct weighted network graph we are shooting for. 
% The processing involves categorizes the pre- and post-IDs of each connection by which node they will belong to. 
% Then I run the digraph() on these classifications, not on their ID numbers 
% (otherwise each cell would have a different node and it would be pure chaos).
%%
% This function creates a directed graph with 12 nodes, 
% which works for any network between 2 cell types where you want to sort 
% by left and right sides and 
% categorize the interneurons into inhibitory vs. excitatory and left vs. right. 
% There are also nodes that are central instead of left or right. 
% Rarely, there are neurons classified as central in Codex, which is why I included those nodes
% If there are no central neurons on the current path, those nodes won't show up on the graph.
%%
function [categorized_pathConnections,directedGraph] = digraph12n(pathConnections,nodeNames,startCellTypeIDs,endCellTypeIDs,codexConnections,codexClassifications)
%   1. creates a digraph with 12 nodes
%   2. adds 2 columns to pathConnectionsTable that classify pre and post IDs for each connection

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
        endCellTypeCheck = ismember(postID,endCellTypeIDs);
         
        % classify pre ID and store info
        switch preSideOfBrain
            case 'left',   preSideIdx = 1;
            case 'right',  preSideIdx = 2;
            case 'center', preSideIdx = 3;
        end

        if strcmp(currConnType,'other')
            pathConnections.preClassification(n) = "other";
        else
            if startCellTypeCheck 
                offset = 0;
            elseif strcmp(currConnType,'excitatory') 
                offset = 3;
            elseif strcmp(currConnType,'inhibitory')
                offset = 6;           
            end
            pathConnections.preClassification(n) = nodeNames(offset + preSideIdx);
        end

        % classify post ID and store info
        switch postSideOfBrain
            case 'left',   postSideIdx = 1;
            case 'right',  postSideIdx = 2;
            case 'center', postSideIdx = 3;
        end

        if strcmp(nextConnType,'other')
            pathConnections.postClassification(n) = "other";
        else
            if endCellTypeCheck
                offset = 9;
            elseif strcmp(nextConnType,'excitatory') 
                offset = 3;
            elseif strcmp(nextConnType,'inhibitory')
                offset = 6;
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

