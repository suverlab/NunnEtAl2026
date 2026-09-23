%% APN2 Connectomics
% run this code after running the DFS algorithm in Python
% By Kaylee Odum
% Created 260110
% Last Update 260721

clearvars, close all

%% user inputs

todaysDate = "260916";

% names of ID files
apn2IDs_fileName = "APN2_IDs.txt";
jonIDs_fileName = "260320_JON_IDs.txt";
bnIDs_fileName = "BN_IDs_antNv.txt";
hbnIDs_fileName = "BN_IDs_head.txt";
dnxIDs_fileName = "DNx01_IDs.txt";
amnIDs_fileName = "AMN_IDs.txt";

% names of files obtained from running the DFS on python
apnAmnPathsCSV_fileName = "260223_APN2_AMN_paths.csv";
jonApnPathsCSV_fileName = "260320_JON_APN2_paths.csv";
bnApnPathsCSV_fileName = "260223_BN_APN2_paths.csv";
hbnApnPathsCSV_fileName = "260730_hBN_APN2_paths.csv";
dnxApnPathsCSV_fileName = "260730_DNx01_APN2_paths.csv";

% Toggle to load the following variables into each section (=1), or just run the script fresh (=0)
loadVars = 1;

% variables obtained from running code in this script that separates CSV file paths 
% (if not run yet, variable should = [])
apnAmnPathConns_fileName = "260324_connections_apnAmnPaths.mat";
jonApnPathConns_fileName = "260324_connections_jonApnPaths.mat";
bnApnPathConns_fileName = "260324_connections_bnApnPaths.mat";

% variables obtained from running code in this script that creates digraph and categorizes each connection in the path connections variables
% (if not run yet, variable should = [])
categorized_apnAmnPathConns_fileName = "260324_categorizedConnections_apnAmnPaths.mat";
categorized_jonApnPathConns_fileName = "260324_categorizedConnections_jonApnPaths.mat";
categorized_bnApnPathConns_fileName = "260324_categorizedConnections_bnApnPaths.mat";

apnAmnDirGraph_fileName = "260324_digraph_apnAmnPaths.mat";
jonApnDirGraph_fileName = "260324_digraph_jonApnPaths.mat";
bnApnDirGraph_fileName = "260324_digraph_bnApnPaths.mat";

% windows OS paths
connectomicsFolderPath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Suver Lab\APN2 Connectomics";
dataUploadPath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Suver Lab\APN2 Connectomics\data";
savedVariablesPath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Suver Lab\APN2 Connectomics\savedVariables";
savedFiguresPath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Suver Lab\APN2 Connectomics\savedFigures";

% mac OS paths
% connectomicsFolderPath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Suver Lab/APN2 Connectomics";
% dataUploadPath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Suver Lab/APN2 Connectomics/data";
% savedVariablesPath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Suver Lab/APN2 Connectomics/savedVariables";
% savedFiguresPath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Suver Lab/APN2 Connectomics/savedFigures";

%% load the data to be analyzed (obtained from various sources)

% set path
cd(dataUploadPath) 

% upload cell type data of all cells in connectome
opts = detectImportOptions("cellTypes_princeton.csv"); % prevents matlab from shortening cell ID numbers
    opts = setvartype(opts,'root_id','string');
cellTypes = readtable("cellTypes_princeton.csv",opts);

% upload classifications data of all cells in connectome
opts = detectImportOptions("classification_princeton.csv"); 
    opts = setvartype(opts,'root_id','string');
classifications = readtable("classification_princeton.csv",opts);

% upload entire connectome
opts = detectImportOptions("connections_princeton.csv");
    opts = setvartype(opts,'pre_root_id','string');
    opts = setvartype(opts,'post_root_id','string');
connections = readtable("connections_princeton.csv",opts);

% upload APN2 cell IDs 
apn2IDs = readlines(apn2IDs_fileName,"LineEnding",',');

% upload JON cell IDs 
jonIDs = readlines(jonIDs_fileName,"LineEnding",',');

% upload bristle cell IDs 
bnIDs = readlines(bnIDs_fileName,"LineEnding",',');
hbnIDs = readlines(hbnIDs_fileName,"LineEnding",',');

% upload DNx01 cell IDs
dnxIDs = readlines(dnxIDs_fileName,"LineEnding",',');

% upload AMN cell IDs (cell names obtained from Marie and IDs downloaded from codex)
amnIDs = readlines(amnIDs_fileName,"LineEnding",',');

% APN2 mirror pairs (column 1 = left; column 2 = right)
apnMirrorPairs = ["720575940614190242","720575940635960292";
                  "720575940629552170","720575940627027087";
                  "720575940606919602","720575940625952002";
                  "720575940629394232","720575940613641794";
                  "720575940650920569","720575940618772719";
                  "720575940603930046","720575940620833974";
                  "720575940630214564","no match";
                  "no match","720575940625597317"];

% Categorize AMNs by which muscle they project to
amnGroups = table;
amnGroups.amnID = amnIDs;
amnGroups.AMN(:,1) = ["3","3","2","2","4a","4a","1","1","4b","4b"]';
amnGroups.muscleName(:,1) = ["posterior depressor","posterior depressor","anterior depressor","anterior depressor","posterior levator","posterior levator","anterior levator","anterior levator", "posterior levator", "posterior levator"];

%% Functions

% 3 main functions
function pathConnections = splitPaths(csvFileName,csvPathName,codexConnections)
%   1. uploads DFS results from Python
%   2. splits the paths in the paths CSV file (obtained from running the Python DFS code) so that each cell ID in the path is its own element
%   3. splits paths up into pre/post root id connections

    % load the paths data 
    cd(csvPathName)
    pathsCSV = readcell(csvFileName);

    % split groups of paths into 1 path/row
    paths = {};

    numPathGroups = numel(pathsCSV);
    step = max(1, round(numPathGroups/10)); % guarantees mod function below doesn't break if there are less than 10 rows

    for p = 1:numPathGroups
        groupedPaths = pathsCSV{p};
        if strcmp(groupedPaths,'None')
            continue
        end

        if mod(p, step) == 0
            fprintf('%.0f %% of grouped paths have been separated\n',(p/numPathGroups)*100)
        end

        singlePaths = regexp(groupedPaths, '\[([^\[\]]+)\]', 'tokens');
    
        for k = 1:numel(singlePaths)
            paths{end+1,1} = strsplit(singlePaths{k}{1}, ', ');
        end
    end 

    % create lookup maps for extracting info from connections dataset
    keys = strcat(string(codexConnections.pre_root_id), "_", ...
              string(codexConnections.post_root_id));

    synMap = containers.Map(keys, codexConnections.syn_count);
    ntMap  = containers.Map(keys, codexConnections.nt_type);

    % split each path into a pre and post id
    preIDs = {};
    postIDs = {};
    numSyn = [];
    NT = {};

    idx = 1;
    
    numPaths = size(paths,1);
    step = max(1, round(numPaths/10)); % guarantees mod function below doesn't break if there are less than 10 paths

    for p = 1:numPaths       
        if mod(p, step) == 0
            fprintf('%.0f %% of paths have been split into connections\n',(p/numPaths)*100)
        end

        currPath = paths{p};
        pathLength = length(currPath) - 1;

        for n = 1:pathLength
            currConn = string(currPath{n}) + "_" + string(currPath{n+1});
          
            preIDs(idx,1)  = currPath(n);
            postIDs(idx,1) = currPath(n+1);
            numSyn(idx,1) = sum(synMap(currConn));
            NT{idx,1} = ntMap(currConn);         

            idx = idx+1;
        end
    end

    pathConnections = table(preIDs,postIDs,numSyn,NT, ...
    'VariableNames',{'preRootID','postRootID','numSynapses','NTtype'});

    % make sure there are no duplicates
    pathConnections = unique(pathConnections);
end

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

% Adapted functions
function [categorized_pathConnections,directedGraph] = digraph4n(pathConnections,nodeNames,codexClassifications)
% adapted from digraph12n function (just has 4 nodes instead of 12)
% doesn't classify NT type b/c we are only looking at inhibitory connections
% specifically used for analyzing direct inhibitory inputs to APN2 that are independent of both JONs and BNs

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

function [categorized_pathConnections,directedGraph] = digraph8n(pathConnections,nodeNames,startCellTypeIDs,endCellTypeIDs1,endCellTypeIDs2,codexConnections,codexClassifications)
% adapted from digraph2 function (just has 8 nodes instead of 12)
% rewrite digraph12n function to have 8 nodes and 2 endCellTypes
% specifically used to analyze premotor inputs to both APN2 and AMN

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

function [categorized_pathConnections,directedGraph] = digraph10n(pathConnections,nodeNames,startCellTypeIDs,endCellTypeIDs1,endCellTypeIDs2,codexConnections,codexClassifications)
% adapted from digraph2 function (just has 10 nodes instead of 12)
% rewrite digraph12n function to have 10 nodes and 2 endCellTypes
% specifically used to analyze neurons that input onto both APN2 and premotor neurons

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

%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% THE FOLLOWING ARE ANALYSES OTHER THAN THE MAIN PATHWAYS
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%
%% NT composition of all direct inputs to APN2 

% extract APN2 input connections
apnInputConns  = connections(ismember(connections.post_root_id,apn2IDs),:);

% extract columns from connections table into their own array (optimizes indexing and increases speed)
apnInputIDs = unique(string(apnInputConns.pre_root_id));
apnInputNTs = apnInputConns.nt_type;

% apnInputNtMap = containers.Map(apnInputIDs,apnInputNTs);

% categorize synapse types
apnInputType = strings(size(apnInputIDs,1));
apnInputNumSyn = nan(size(apnInputIDs,1));
for n = 1:size(apnInputIDs,1)  
   currNt = unique(apnInputConns.nt_type(strcmp(apnInputConns.pre_root_id,apnInputIDs(n))));
    if strcmp(currNt,'ACH')
        apnInputType(n) = 'excitatory';
    elseif strcmp(currNt,'GABA') | strcmp(currNt,'GLUT')
        apnInputType(n) = 'inhibitory';
    else 
        apnInputType(n) = 'other';
    end

    apnInputNumSyn(n) = sum(apnInputConns.syn_count(strcmp(apnInputConns.pre_root_id,apnInputIDs(n))));
end

% count number of inputs in each group
numExcInputs = length(find(strcmp(apnInputType,"excitatory")));
numInhInputs = length(find(strcmp(apnInputType,"inhibitory")));

numExcSyn = sum(apnInputNumSyn(strcmp(apnInputType,"excitatory")));
numInhSyn = sum(apnInputNumSyn(strcmp(apnInputType,"inhibitory")));

% plot
excVsInhPiechart = figure;
piechart([numExcInputs,numInhInputs],{'Excitatory','Inhibitory'})
title("Direct inputs to APN2 (num neurons)")

excVsInhSynPiechart = figure;
piechart([numExcSyn,numInhSyn],{'Excitatory','Inhibitory'})
title("Direct inputs to APN2 (num synapses)")

cd(savedFiguresPath)
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputs_piechart",'fig')
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputs_piechart",'png')
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputs_piechart",'pdf')

saveas(excVsInhSynPiechart,todaysDate+"_excVsInh_apnInputsSynCount_piechart",'fig')
saveas(excVsInhSynPiechart,todaysDate+"_excVsInh_apnInputsSynCount_piechart",'png')
saveas(excVsInhSynPiechart,todaysDate+"_excVsInh_apnInputsSynCount_piechart",'pdf')

%% NT composition of direct inputs to APN2 that are downstream JONs
% use the path connections table to extract these inputs

if loadVars == 1
    cd(savedVariablesPath)
    if ~isempty(jonApnPathConns_fileName)
        load(jonApnPathConns_fileName)
    elseif isempty(jonApnPathConns_fileName) & ~exist("jonApnPathConns",'var')
        error("necessary variable does not exist. please run previous section or upload proper variable.")
    end
end

% extract columns from path connections table into their own array (optimizes indexing and increases speed)
jonApnPreIDs  = string(jonApnPathConns.preRootID);
jonApnPostIDs = string(jonApnPathConns.postRootID);
jonApnNTs = jonApnPathConns.NTtype;

% extract direct inputs to APN2 and their NTs
jApnInputIDs  = jonApnPreIDs(ismember(jonApnPostIDs,apn2IDs));
jApnInputNTs = jonApnNTs(ismember(jonApnPostIDs,apn2IDs));

jApnInputNtMap = containers.Map(jApnInputIDs ,jApnInputNTs);

% categorize synapse types
jApnInputSynType(size(jApnInputIDs ,1)) = string;
for n = 1:size(jApnInputIDs ,1)  
   currNt = jApnInputNtMap(jApnInputIDs (n));
    if strcmp(currNt,'ACH')
        jApnInputSynType(n) = 'excitatory';
    elseif strcmp(currNt,'GABA') | strcmp(currNt,'GLUT')
        jApnInputSynType(n) = 'inhibitory';
    else 
        jApnInputSynType(n) = 'other';
    end
end

% count number of inputs in each group
numExcInputs = length(find(strcmp(jApnInputSynType,"excitatory")));
numInhInputs = length(find(strcmp(jApnInputSynType,"inhibitory")));

% plot
excVsInhPiechart = figure;
piechart([numExcInputs,numInhInputs],{'Excitatory','Inhibitory'})
title("Direct inputs to APN2 that are downstream JONs")

cd(savedFiguresPath)
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputsDowstreamJons_piechart",'fig')
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputsDowstreamJons_piechart",'png')
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputsDowstreamJons_piechart",'pdf')

%% NT composition of direct inputs to APN2 that are downstream BNs
% use the path connections table to extract these inputs

if loadVars == 1
    cd(savedVariablesPath)
    if ~isempty(bnApnPathConns_fileName)
        load(bnApnPathConns_fileName)
    elseif isempty(bnApnPathConns_fileName) & ~exist("bnApnPathConns",'var')
        error("necessary variable does not exist. please run previous section or upload proper variable.")
    end
end

% extract columns from path connections table into their own array (optimizes indexing and increases speed)
bnApnPreIDs  = string(bnApnPathConns.preRootID);
bnApnPostIDs = string(bnApnPathConns.postRootID);
bnApnNTs = string(bnApnPathConns.NTtype);

% extract direct inputs to APN2 and their NTs
bApnInputID = bnApnPreIDs(ismember(bnApnPostIDs,apn2IDs));
bToApnNumHops = bnApnNTs(ismember(bnApnPostIDs,apn2IDs));

bApnInputNtMap = containers.Map(bApnInputID,bToApnNumHops);

% categorize synapse types
bApnInputSynType(size(bApnInputID,1)) = string;
for n = 1:size(bApnInputID,1)  
   currNt = bApnInputNtMap(bApnInputID(n));
    if strcmp(currNt,'ACH')
        bApnInputSynType(n) = 'excitatory';
    elseif strcmp(currNt,'GABA') | strcmp(currNt,'GLUT')
        bApnInputSynType(n) = 'inhibitory';
    else 
        bApnInputSynType(n) = 'other';
    end
end

% count number of inputs in each group
numExcInputs = length(find(strcmp(bApnInputSynType,"excitatory")));
numInhInputs = length(find(strcmp(bApnInputSynType,"inhibitory")));

% plot
excVsInhPiechart = figure;
piechart([numExcInputs,numInhInputs],{'Excitatory','Inhibitory'})
title("Direct inputs to APN2 that are downstream BNs")

% save figure
cd(savedFiguresPath)
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputsDowstreamBNs_piechart",'fig')
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputsDowstreamBNs_piechart",'png')
saveas(excVsInhPiechart,todaysDate+"_excVsInh_apnInputsDowstreamBNs_piechart",'pdf')

%% Bar graph of all inhibitory inputs to APN2, sorted by super class and by JON-dependent vs independent

% extract columns from classifications and connections and table into their own array (optimizes indexing and increases speed)
allRootIDs = classifications.root_id;
allSuperClass = classifications.super_class;
allPreIDs = connections.pre_root_id;
allPostIDs = connections.post_root_id;
allNTs = string(connections.nt_type);

% index out JON-dependent and JON-independent inhibitory APN2 inputs and their cell classes
inhInputIDs = allPreIDs(ismember(allPostIDs,apn2IDs) & (strcmp(allNTs,"GLUT") | strcmp(allNTs,"GABA")));

jDepInputIDs = inhInputIDs(ismember(inhInputIDs,jApnInputIDs));
jDepInputClass = string(allSuperClass(ismember(allRootIDs,jDepInputIDs)));

jIndInputIDs = inhInputIDs(~ismember(inhInputIDs,jDepInputIDs));
jIndInputClass = string(allSuperClass(ismember(allRootIDs,jIndInputIDs)));

% create the bar graph
class = ["ascending","central","descending","sensory"];
dependency = ["JO-dependent","JO-independent"];
jInputsClassAndDep = nan(length(class),length(dependency));

for i = 1:length(class)
    jInputsClassAndDep(i,1) = length(jDepInputClass(strcmp(jDepInputClass,class(i)))); % number of JO-dependent inputs in the current class (column 1 = dep)
    jInputsClassAndDep(i,2) = length(jIndInputClass(strcmp(jIndInputClass,class(i)))); % number of JO-independent inputs in the current class (column 2 = ind)
end

% plot 
bargraph = figure; 
sgtitle("JO-dependent and JO-independent inhibitory APN2 inputs in by class")
subplot(1,2,1) %count
bar(jInputsClassAndDep,'stacked')
xticks(1:1:length(class))
xticklabels(class)
xlabel("Cell Type")
ylabel("num inputs")
legend(dependency)
title("Number of APN2 inputs")

subplot(1,2,2) %percentage
bar((jInputsClassAndDep/sum(jInputsClassAndDep,"all"))*100,'stacked')
xticks(1:1:length(class))
xticklabels(class)
xlabel("Cell Type")
ylabel("% inputs")
legend(dependency)
title("Percent of APN2 inputs")

% save figure
cd(savedFiguresPath)
saveas(bargraph,todaysDate+"_apnInhibInputs_JonDepVsInd_bargraph",'fig')
saveas(bargraph,todaysDate+"_apnInhibInputs_JonDepVsInd_bargraph",'png')
saveas(bargraph,todaysDate+"_apnInhibInputs_JonDepVsInd_bargraph",'pdf')

%% Network graph of inhibitory inputs to APN2 that are independent of both JONs and BNs

% extract columns from classifications and connections and table into their own array (optimizes indexing and increases speed) -> repeated code so I can run only this section if I want
allPreIDs = connections.pre_root_id;
allPostIDs = connections.post_root_id;
allNTs = string(connections.nt_type);

% index out JON-independent and BN-independent inhibitory APN2 inputs -> some repeated code so I can run only this section if I want
inhInputIDs = allPreIDs(ismember(allPostIDs,apn2IDs) & (strcmp(allNTs,"GLUT") | strcmp(allNTs,"GABA")));

jDepInputIDs = inhInputIDs(ismember(inhInputIDs,jApnInputIDs));
jIndInputIDs = inhInputIDs(~ismember(inhInputIDs,jDepInputIDs));

bDepInputIDs = inhInputIDs(ismember(inhInputIDs,bApnInputID));
bIndInputIDs = inhInputIDs(~ismember(inhInputIDs,bDepInputIDs));

allIndInputIDs = [jIndInputIDs;bIndInputIDs]; 
allIndInputIDs = unique(allIndInputIDs);

% extract cell types of the independent inputs
indInputCellTypes = cellTypes.primary_type(ismember(cellTypes.root_id,allIndInputIDs)); 

% extract the specific connections for all independent inhibitory inputs
indInputConns = connections(ismember(allPreIDs,allIndInputIDs) & (ismember(allPostIDs,amnIDs) | ismember(allPostIDs,apn2IDs)),:);
indInputConns = renamevars(indInputConns,["pre_root_id","post_root_id","nt_type","syn_count"],["preRootID","postRootID","NTtype","numSynapses"]);

% create directed graph
indInputNodeNames = ["inh.indInput_L","inh.indInput_R","APN2_L","APN2_R"];

[indInputConns,indInputDirGraph] = digraph4n(indInputConns,indInputNodeNames,classifications);

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

%% Premotor neurons (AMN inputs) that also input onto APN2

% extract columns from classifications and connections and table into their own array (optimizes indexing and increases speed) -> repeated code so I can run only this section if I want
allPreIDs = connections.pre_root_id;
allPostIDs = connections.post_root_id;
allNTs = string(connections.nt_type);

% extract direct inputs to AMN and their NTs
amnInputIDs = allPreIDs(ismember(allPostIDs,amnIDs));
amnInputNTs = allNTs(ismember(allPreIDs,amnInputIDs));

% extract direct inputs to APN2 and find overlap IDs
apnInputIDs = allPreIDs(ismember(allPostIDs,apn2IDs));
premApnAmnInputIDs = unique(amnInputIDs(ismember(amnInputIDs,apnInputIDs)));

% extract the specific connections for both AMN and APN2
premApnAmnConns = connections(ismember(allPreIDs,premApnAmnInputIDs) & (ismember(allPostIDs,amnIDs) | ismember(allPostIDs,apn2IDs)),:);
premApnAmnConns = renamevars(premApnAmnConns,["pre_root_id","post_root_id","nt_type","syn_count"],["preRootID","postRootID","NTtype","numSynapses"]);

% create directed graph
premApnAmnNodeNames = ["exc.premotor_L","exc.premotor_R","inh.premotor_L","inh.premotor_R","APN2_L","APN2_R","AMN_L","AMN_R"];

[premApnAmnConns,premApnAmnDirGraph] = digraph8n(premApnAmnConns,premApnAmnNodeNames,premApnAmnInputIDs,apn2IDs,amnIDs,connections,classifications);

% save directed graph and categorized path connections variable
cd(savedVariablesPath)
save(todaysDate+"_digraph_premApnAmnPaths.mat","premApnAmnDirGraph")
save(todaysDate+"_categorizedConnections_premApnAmnConns.mat","premApnAmnConns")

% plot weighted network graph
premApnAmnNetworkGraph = plotDigraph(premApnAmnDirGraph,premApnAmnConns,premApnAmnNodeNames,0.25);

% save figure
cd(savedFiguresPath)
saveas(premApnAmnNetworkGraph,todaysDate+"_premApnAmnNetworkGraph",'fig')
saveas(premApnAmnNetworkGraph,todaysDate+"_premApnAmnNetworkGraph",'png')
saveas(premApnAmnNetworkGraph,todaysDate+"_premApnAmnNetworkGraph",'pdf')

%% Interneurons that project to both premotor neurons and APN2

% extract columns from classifications and connections and table into their own array (optimizes indexing and increases speed) -> repeated code so I can run only this section if I want
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
inputPremApnInputIDs(ismember(inputPremApnInputIDs,apn2IDs) | ismember(inputPremApnInputIDs,jonIDs)) = [];

% extract the specific connections for both premotor and APN2
inputsPremApnConns = connections(ismember(allPreIDs,inputPremApnInputIDs) & (ismember(allPostIDs,amnInputIDs) | ismember(allPostIDs,apn2IDs)),:);
inputsPremApnConns = renamevars(inputsPremApnConns,["pre_root_id","post_root_id","nt_type","syn_count"],["preRootID","postRootID","NTtype","numSynapses"]);

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

%% Network graph of APN2 outputs to premotor neurons (1-2 hops) where each APN2 neuron is a different node

% -- Extract/compile info and create adjacency matrix -- %

% extract columns from classifications and connections and table into their own array (optimizes indexing and increases speed) -> repeated code so I can run only this section if I want
allSides = string(classifications.side);
allPreIDs = connections.pre_root_id;
allPostIDs = connections.post_root_id;
allNTs = string(connections.nt_type);

% extract premotor neurons (direct inputs to AMN neurons) and their NTs -> repeated code so I can run only this section if I want
pmnIDs = allPreIDs(ismember(allPostIDs,amnIDs));

% index out APN2 to premotor connections (max 2 hops)
apnPmnConns = connections(ismember(allPreIDs,apn2IDs) & ismember(allPostIDs,pmnIDs),:);
apnPmnInIDs = unique(allPreIDs(ismember(allPreIDs,allPostIDs(ismember(allPreIDs,apn2IDs))) & ismember(allPostIDs,pmnIDs)));
% apnPmnConns_2hop = connections(ismember(allPreIDs,apn2IDs) & ismember(allPostIDs,apnPmnInIDs) | ismember(allPreIDs,apnPmnInIDs) & ismember(allPostIDs,pmnIDs),:);
% apnPmnConns = [apnPmnConns_1hop;apnPmnConns_2hop];

% find AMN IDs that are only 3 hops downstream APN2
% pmnDsApnIDs = unique(apnPmnConns.post_root_id(ismember(apnPmnConns.post_root_id,pmnIDs)));
pmnDsApnIDs = unique(apnPmnConns.post_root_id);
pmnAmnDsApnConns = connections(ismember(allPreIDs,pmnDsApnIDs) & ismember(allPostIDs,amnIDs),:);
amnDsApnIDs = unique(pmnAmnDsApnConns.post_root_id);
amnDsApnGroups = amnGroups(ismember(amnGroups.amnID,amnDsApnIDs),:);

% sort APN2 to PMN connections based on which AMN pair they project to
pmnIDs1 = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.AMN,"1")))));
pmnIDs2 = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.AMN,"2")))));
pmnIDs3 = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.AMN,"3")))));
pmnIDs4a = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.AMN,"4a")))));
pmnIDs4b = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.AMN,"4b")))));

apnPmnConns.targetAmnPair(ismember(apnPmnConns.post_root_id,pmnIDs1),1) = "1";
apnPmnConns.targetAmnPair(ismember(apnPmnConns.post_root_id,pmnIDs2),1) = "2";
apnPmnConns.targetAmnPair(ismember(apnPmnConns.post_root_id,pmnIDs3),1) = "3";
apnPmnConns.targetAmnPair(ismember(apnPmnConns.post_root_id,pmnIDs4a),1) = "4a";
apnPmnConns.targetAmnPair(ismember(apnPmnConns.post_root_id,pmnIDs4b),1) = "4b";

% sort APN2 to PMN connections based on which muscle the downstream AMNs project to
PDpmnIDs = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.muscleName,"posterior depressor")))));
ADpmnIDs = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.muscleName,"anterior depressor")))));
PLpmnIDs = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.muscleName,"posterior levator")))));
ALpmnIDs = unique(pmnAmnDsApnConns.pre_root_id(ismember(pmnAmnDsApnConns.post_root_id,amnDsApnGroups.amnID(strcmp(amnDsApnGroups.muscleName,"anterior levator")))));

apnPmnConns.targetMuscle(ismember(apnPmnConns.post_root_id,PDpmnIDs),1) = "posterior depressor";
apnPmnConns.targetMuscle(ismember(apnPmnConns.post_root_id,ADpmnIDs),1) = "anterior depressor";
apnPmnConns.targetMuscle(ismember(apnPmnConns.post_root_id,PLpmnIDs),1) = "posterior levator";
apnPmnConns.targetMuscle(ismember(apnPmnConns.post_root_id,ALpmnIDs),1) = "anterior levator";

% sum synapses for identical pre/post pairs
[Gpairs, uniqPreID, uniqPostID] = findgroups(apnPmnConns.pre_root_id, apnPmnConns.post_root_id);
totalSyn = splitapply(@sum, apnPmnConns.syn_count, Gpairs);

% node list and indices
nodeIDs = unique([uniqPreID; uniqPostID],'stable');

[~,preIdx]  = ismember(uniqPreID,nodeIDs);
[~,postIdx] = ismember(uniqPostID,nodeIDs);

% adjacency matrix
apnPmnAdjMatrix = sparse(preIdx, postIdx, totalSyn, ...
           numel(nodeIDs), numel(nodeIDs));

% directed graph
apnPmnDirGraph = digraph(apnPmnAdjMatrix,nodeIDs);



% -- Format network graph and plot -- %

% label output (premotor) nodes based on which AMN pair and muscle group they target downstream + label side of brain of each node
nodeNames = apnPmnDirGraph.Nodes.Name;
for n = 1:size(nodeNames,1)
    currentAmnPair = unique(apnPmnConns.targetAmnPair(apnPmnConns.post_root_id==nodeNames(n)));
    if isempty(currentAmnPair)
        apnPmnDirGraph.Nodes.targetAmnPair(n,1) = "";
    else
        apnPmnDirGraph.Nodes.targetAmnPair(n,1) = currentAmnPair;
    end

    currentMuscle = unique(apnPmnConns.targetMuscle(apnPmnConns.post_root_id==nodeNames(n)));
    if isempty(currentMuscle)
        apnPmnDirGraph.Nodes.targetMuscle(n,1) = "";
    else
        apnPmnDirGraph.Nodes.targetMuscle(n,1) = currentMuscle;
    end

    apnPmnDirGraph.Nodes.side(n,1) = unique(classifications.side(ismember(classifications.root_id,nodeNames(n))));
end

% color output nodes based on NT type
nodeColor = zeros(numnodes(apnPmnDirGraph),3); % default node color = black

ntColors = [0 0 0 % black = excitatory (ACH)
            0.7 0.7 0.7]; % light gray = inhibitory (GLUT/GABA)

outputNodes = unique(string(apnPmnConns.post_root_id));

for o = 1:length(outputNodes)
    output = outputNodes(o);

    nt = unique(string(apnPmnConns.nt_type(apnPmnConns.post_root_id == output)));  % get this input's nt
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

edgeColor = zeros(numedges(apnPmnDirGraph),3); % default = black

for k = 1:size(apnMirrorPairs,1) 
    pairs = string(apnMirrorPairs(k,:));

    nidx = ismember(nodeIDs,pairs); % find nodes of current pair
    nodeColor(nidx,:) = repmat(nodeColors(k,:),sum(nidx),1); % assign color

    eidx = ismember(string(apnPmnDirGraph.Edges.EndNodes(:,1)),pairs); % find edges whose origin is one of these posts
    edgeColor(eidx,:) = repmat(nodeColors(k,:),sum(eidx),1); % assign color
end

% make only the APN2 nodes large
nodeSize = 5*ones(numnodes(apnPmnDirGraph),1); % set default node size
nodeSize(ismember(nodeNames, apn2IDs)) = 15; % make APN2 nodes larger

% plot
apnPmnNetworkGraph = figure;
plot(apnPmnDirGraph,'Layout','force', ...
    'MarkerSize',nodeSize,'NodeColor',nodeColor,'EdgeColor',edgeColor, ...
    'ArrowSize',0, ...
    'NodeLabel',apnPmnDirGraph.Nodes.targetMuscle+" "+apnPmnDirGraph.Nodes.targetAmnPair+" "+apnPmnDirGraph.Nodes.side);
    

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
title("Individual APN2 input to premotor neurons (1 hop), categorized by downstream muscle target")

% save figure
cd(savedFiguresPath)
saveas(apnPmnNetworkGraph,todaysDate+"_indivApnToPmnNetworkGraph",'fig')
saveas(apnPmnNetworkGraph,todaysDate+"_indivApnToPmnNetworkGraph",'png')
saveas(apnPmnNetworkGraph,todaysDate+"_indivApnToPmnNetworkGraph",'pdf')

%% Heatmap of adjacency matrix showing APN2 input connectivity

% extract columns from classifications and connections and table into their own array (optimizes indexing and increases speed) -> repeated code so I can run only this section if I want
allRootIDs = string(classifications.root_id);
allSides = string(classifications.side);
allPreIDs = connections.pre_root_id;
allPostIDs = connections.post_root_id;
allNTs = string(connections.nt_type);

% -- Upload DFS paths data -- % 
cd(savedVariablesPath)
jToApnPathsCSV = readcell("260320_JON_APN2_paths.csv");
bToApnPathsCSV = readcell("260223_BN_APN2_paths.csv");
hToApnPathsCSV = readcell("260730_hBN_APN2_paths.csv");
dToApnPathsCSV = readcell("260730_DNx01_APN2_paths.csv");

% -- Split groups of paths into 1 path/row -- %

% JON
jToApnPaths = {};

numPathGroups = numel(jToApnPathsCSV);
step = max(1, round(numPathGroups/10)); % guarantees mod function below doesn't break if there are less than 10 rows

for p = 1:numPathGroups
    groupedPaths = jToApnPathsCSV{p};
    if strcmp(groupedPaths,'None')
        continue
    end

    if mod(p, step) == 0
        fprintf('%.0f %% of grouped paths have been separated\n',(p/numPathGroups)*100)
    end

    singlePaths = regexp(groupedPaths, '\[([^\[\]]+)\]', 'tokens');

    for k = 1:numel(singlePaths)
        jToApnPaths{end+1,1} = strsplit(singlePaths{k}{1}, ', ');
    end
end 

% ant BN
bToApnPaths = {};

numPathGroups = numel(bToApnPathsCSV);
step = max(1, round(numPathGroups/10)); % guarantees mod function below doesn't break if there are less than 10 rows

for p = 1:numPathGroups
    groupedPaths = bToApnPathsCSV{p};
    if strcmp(groupedPaths,'None')
        continue
    end

    if mod(p, step) == 0
        fprintf('%.0f %% of grouped paths have been separated\n',(p/numPathGroups)*100)
    end

    singlePaths = regexp(groupedPaths, '\[([^\[\]]+)\]', 'tokens');

    for k = 1:numel(singlePaths)
        bToApnPaths{end+1,1} = strsplit(singlePaths{k}{1}, ', ');
    end
end 

% head BN
hToApnPaths = {};

numPathGroups = numel(hToApnPathsCSV);
step = max(1, round(numPathGroups/10)); % guarantees mod function below doesn't break if there are less than 10 rows

for p = 1:numPathGroups
    groupedPaths = hToApnPathsCSV{p};
    if strcmp(groupedPaths,'None')
        continue
    end

    if mod(p, step) == 0
        fprintf('%.0f %% of grouped paths have been separated\n',(p/numPathGroups)*100)
    end

    singlePaths = regexp(groupedPaths, '\[([^\[\]]+)\]', 'tokens');

    for k = 1:numel(singlePaths)
        hToApnPaths{end+1,1} = strsplit(singlePaths{k}{1}, ', ');
    end
end 

% DNx01
dToApnPaths = {};

numPathGroups = numel(dToApnPathsCSV);
step = max(1, round(numPathGroups/10)); % guarantees mod function below doesn't break if there are less than 10 rows

for p = 1:numPathGroups
    groupedPaths = dToApnPathsCSV{p};
    if strcmp(groupedPaths,'None')
        continue
    end

    if mod(p, step) == 0
        fprintf('%.0f %% of grouped paths have been separated\n',(p/numPathGroups)*100)
    end

    singlePaths = regexp(groupedPaths, '\[([^\[\]]+)\]', 'tokens');

    for k = 1:numel(singlePaths)
        dToApnPaths{end+1,1} = strsplit(singlePaths{k}{1}, ', ');
    end
end 

% -- Compile path connections table for head bristles and DNx01 -- %
hbnApnPathConns = splitPaths(hbnApnPathsCSV_fileName,savedVariablesPath,connections); 
dnxApnPathConns = splitPaths(dnxApnPathsCSV_fileName,savedVariablesPath,connections); 

% -- write function for calculate total number of synapses along paths -- %
function totalSyn = pathSynapses(path,synLookup)

totalSyn = 0;

for k = 1:numel(path)-1
    key = char(strcat(string(path{k}),"->",string(path{k+1})));
    if isKey(synLookup,key)
        totalSyn = totalSyn + synLookup(key);
    end
end
end

% -- preallocate vars -- %
jApnInputID(size(jToApnPaths,1),1) = string;
jDsApnID(size(jToApnPaths,1),1) = string;
jToApnNumHops(size(jToApnPaths,1),1) = nan;
jToApnNumSyn(size(jToApnPaths,1),1) = nan;

bApnInputID(size(bToApnPaths,1),1) = string; 
bDsApnID(size(bToApnPaths,1),1) = string;
bToApnNumHops(size(bToApnPaths,1),1) = nan;
bToApnNumSyn(size(bToApnPaths,1),1) = nan;

hApnInputID(size(hToApnPaths,1),1) = string; 
hDsApnID(size(hToApnPaths,1),1) = string;
hToApnNumHops(size(hToApnPaths,1),1) = nan;
hToApnNumSyn(size(hToApnPaths,1),1) = nan;

dApnInputID(size(dToApnPaths,1),1) = string; 
dDsApnID(size(dToApnPaths,1),1) = string;
dToApnNumHops(size(dToApnPaths,1),1) = nan;
dToApnNumSyn(size(dToApnPaths,1),1) = nan;

% -- ant BN lookup table -- %

edgeKey = strcat(bnApnPathConns.preRootID,"->",bnApnPathConns.postRootID);

% sum duplicate edges
[G,key] = findgroups(edgeKey);
syn = splitapply(@sum,bnApnPathConns.numSynapses,G);

% build lookup table
synLookup = containers.Map(cellstr(key),syn);

for n = 1:length(bToApnPaths)

    path = bToApnPaths{n};

    bApnInputID(n,1) = path{1};
    bDsApnID(n,1) = path{end};
    bToApnNumSyn(n,1) = pathSynapses(path,synLookup);
    bToApnNumHops(n,1) = numel(path)-1;
end

% -- JON lookup table -- %

jEdgeKey = strcat(jonApnPathConns.preRootID,"->",jonApnPathConns.postRootID);

% sum duplicate edges
[G,key] = findgroups(jEdgeKey);
jSyn = splitapply(@sum,jonApnPathConns.numSynapses,G);

% build lookup table
jSynLookup = containers.Map(cellstr(key),jSyn);

for n = 1:length(jToApnPaths)

    path = jToApnPaths{n};

    if ismember(path{1},bApnInputID)
        continue
    else
        jApnInputID(n,1) = path{1};
        jDsApnID(n,1) = path{end};
        jToApnNumSyn(n,1) = pathSynapses(path,jSynLookup);
        jToApnNumHops(n,1) = numel(path)-1;
    end
end

% -- head BN lookup table -- %

hEdgeKey = strcat(hbnApnPathConns.preRootID,"->",hbnApnPathConns.postRootID);

% sum duplicate edges
[G,key] = findgroups(hEdgeKey);
hSyn = splitapply(@sum,hbnApnPathConns.numSynapses,G);

% build lookup table
hSynLookup = containers.Map(cellstr(key),hSyn);

for n = 1:length(hToApnPaths)

    path = hToApnPaths{n};

    if ismember(path{1},bApnInputID)
        continue
    else
        hApnInputID(n,1) = path{1};
        hDsApnID(n,1) = path{end};
        hToApnNumSyn(n,1) = pathSynapses(path,hSynLookup);
        hToApnNumHops(n,1) = numel(path)-1;
    end
end

% -- DNx01 lookup table -- %

dEdgeKey = strcat(dnxApnPathConns.preRootID,"->",dnxApnPathConns.postRootID);

% sum duplicate edges
[G,key] = findgroups(dEdgeKey);
dSyn = splitapply(@sum,dnxApnPathConns.numSynapses,G);

% build lookup table
dSynLookup = containers.Map(cellstr(key),dSyn);

for n = 1:length(dToApnPaths)

    path = dToApnPaths{n};

    if ismember(path{1},bApnInputID)
        continue
    else
        dApnInputID(n,1) = path{1};
        dDsApnID(n,1) = path{end};
        dToApnNumSyn(n,1) = pathSynapses(path,dSynLookup);
        dToApnNumHops(n,1) = numel(path)-1;
    end
end

% -- concatanate all input vars -- %
sensInputID = [jApnInputID;bApnInputID;hApnInputID;dApnInputID];
dsApnID = [jDsApnID;bDsApnID;hDsApnID;dDsApnID];
numHops = [jToApnNumHops;bToApnNumHops;hToApnNumHops;dToApnNumHops];
numSyn = [jToApnNumSyn;bToApnNumSyn;hToApnNumSyn;dToApnNumSyn];

results = table(sensInputID,dsApnID,numSyn,numHops);
results = rmmissing(results);

% -- organize the data so they can be plotted in the desired order -- %

% sort APN2 columns
apnLeft  = apn2IDs(ismember(apn2IDs,...
    allRootIDs(strcmp(allSides,"left"))));

apnRight = apn2IDs(ismember(apn2IDs,...
    allRootIDs(strcmp(allSides,"right"))));

apnOrder = [apnLeft; apnRight];

% sort sensory inputs
jLeft = jApnInputID(ismember(jApnInputID,...
    allRootIDs(strcmp(allSides,"left"))));

bLeft = bApnInputID(ismember(bApnInputID,...
    allRootIDs(strcmp(allSides,"left"))));

hLeft = hApnInputID(ismember(hApnInputID,...
    allRootIDs(strcmp(allSides,"left"))));

dLeft = dApnInputID(ismember(dApnInputID,...
    allRootIDs(strcmp(allSides,"left"))));

jRight = jApnInputID(ismember(jApnInputID,...
    allRootIDs(strcmp(allSides,"right"))));

bRight = bApnInputID(ismember(bApnInputID,...
    allRootIDs(strcmp(allSides,"right"))));

hRight = hApnInputID(ismember(hApnInputID,...
    allRootIDs(strcmp(allSides,"right"))));

dRight = dApnInputID(ismember(dApnInputID,...
    allRootIDs(strcmp(allSides,"right"))));

inputOrder = [unique(jLeft);
                unique(bLeft);
                unique(hLeft);
                unique(dLeft);
                unique(jRight);
                unique(bRight);
                unique(hRight);
                unique(dRight)];

% count the size of each group
njLeft = numel(unique(jLeft));
nbLeft = numel(unique(bLeft));
nhLeft = numel(unique(hLeft));
ndLeft = numel(unique(dLeft));

njRight = numel(unique(jRight));
nbRight = numel(unique(bRight));
nhRight = numel(unique(hRight));
ndRight = numel(unique(dRight));

groupSizes = [njLeft nbLeft nhLeft ndLeft ...
              njRight nbRight nhRight ndRight];

nApnLeft = numel(apnLeft);
nApnRight = numel(apnRight);

% -- create matrices -- %

[G,sInputID,apnID] = findgroups(results.sensInputID,results.dsApnID);

minHop = splitapply(@min,results.numHops,G);

results.group = G;
results.minHop = minHop(G);

shortestResults = results(results.numHops == results.minHop,:);

[G2,sInputID,apnID] = findgroups(shortestResults.sensInputID,...
                                shortestResults.dsApnID);

hop = splitapply(@min,shortestResults.numHops,G2);
syn = splitapply(@sum,shortestResults.numSyn,G2);

synMat1 = nan(length(inputOrder),length(apnOrder));
synMat2 = nan(length(inputOrder),length(apnOrder));

for i = 1:length(sInputID)

    r = find(inputOrder == sInputID(i));
    c = find(apnOrder == apnID(i));

    if isempty(r) || isempty(c)
        continue
    end

    if hop(i) == 1
        synMat1(r,c) = syn(i);

    elseif hop(i) == 2
        synMat2(r,c) = syn(i);

    end
end

synMat1(isnan(synMat1)) = 0;
synMat2(isnan(synMat2)) = 0;

% -- plot 1-hop heatmap -- %
sensToApnHeatmap_1hop = figure;
imagesc(synMat1)
colormap(parula)
colorbar
clim([0 25])
title("Sensory -> APN2 (1 hop)",'FontSize',14)

hold on

yBreaks = cumsum(groupSizes);
for k = 1:length(yBreaks)-1
    yline(yBreaks(k)+0.5,'k','LineWidth',1.5)
end
xline(nApnLeft+0.5,'k','LineWidth',2)

% add group labels
groupLabels = {...
    'JON Left',...
    'Ant BN Left',...
    'Head BN Left',...
    'DNx01 Left',...
    'JON Right',...
    'Ant BN Right',...
    'Head BN Right',...
    'DNx01 Right'};

groupCenters = cumsum(groupSizes) - groupSizes/2; % calc center of each group

yticks(groupCenters)
yticklabels(groupLabels)

set(gca,'XTick',1:1:14,'XTickLabel',apnOrder,'YAxisLocation','left');

text(nApnLeft/2, -30, 'Left APN2', ...
    'HorizontalAlignment','center','FontSize',12)

text(nApnLeft+nApnRight/2, -30, 'Right APN2', ...
    'HorizontalAlignment','center','FontSize',12)

% -- plot 2-hops heatmap -- %
sensToApnHeatmap_2hop = figure;
imagesc(synMat2)
colormap(parula)
colorbar
clim([0 25])
title("Sensory -> APN2 (2 hops)",'FontSize',14)

hold on

yBreaks = cumsum(groupSizes);
for k = 1:length(yBreaks)-1
    yline(yBreaks(k)+0.5,'k','LineWidth',1.5)
end
xline(nApnLeft+0.5,'k','LineWidth',2)

% add group labels
groupLabels = {...
    'JON Left',...
    'Ant BN Left',...
    'Head BN Left',...
    'DNx01 Left',...
    'JON Right',...
    'Ant BN Right',...
    'Head BN Right',...
    'DNx01 Right'};

groupCenters = cumsum(groupSizes) - groupSizes/2; % calc center of each group

yticks(groupCenters)
yticklabels(groupLabels)

set(gca,'XTick',1:1:14,'XTickLabel',apnOrder,'YAxisLocation','left');

text(nApnLeft/2, -30, 'Left APN2', ...
    'HorizontalAlignment','center','FontSize',12)

text(nApnLeft+nApnRight/2, -30, 'Right APN2', ...
    'HorizontalAlignment','center','FontSize',12)

% save figures
cd(savedFiguresPath)

saveas(sensToApnHeatmap_1hop,todaysDate+"_sensInputToApnConnMatrix_1hop",'fig')
saveas(sensToApnHeatmap_1hop,todaysDate+"_sensInputToApnConnMatrix_1hop",'png')
saveas(sensToApnHeatmap_1hop,todaysDate+"_sensInputToApnConnMatrix_1hop",'pdf')

saveas(sensToApnHeatmap_2hop,todaysDate+"_sensInputToApnConnMatrix_2hop",'fig')
saveas(sensToApnHeatmap_2hop,todaysDate+"_sensInputToApnConnMatrix_2hop",'png')
saveas(sensToApnHeatmap_2hop,todaysDate+"_sensInputToApnConnMatrix_2hop",'pdf')

%% APN2 Input Cosine Similarity with Weights (Using Synapse Counts) -> code adapted from Marie's code

apn2IDs_sorted = [apnMirrorPairs(:,1);apnMirrorPairs(:,2)];
apn2IDs_sorted(strcmp(apn2IDs_sorted,"no match")) = [];

% Create weighted matrix
weighted_matrix = zeros(numel(apn2IDs_sorted), numel(sInputID));

% Populate the weighted matrix using synapse counts
for ii = 1:numel(apn2IDs_sorted)
    currentAPN2 = apn2IDs_sorted(ii);
    APN2conn = apnInputConns(apnInputConns.post_root_id == currentAPN2, :);

    for jj = 1:height(APN2conn)
        neuron_idx = find(strcmp(sInputID, APN2conn.pre_root_id{jj}));
        weighted_matrix(ii, neuron_idx) = weighted_matrix(ii, neuron_idx) + APN2conn.syn_count(jj);
    end
end

% Compute cosine similarity
cosine_similarity_weighted = zeros(numel(apn2IDs_sorted));
for ii = 1:numel(apn2IDs_sorted)
    for jj = 1:numel(apn2IDs_sorted)
        if norm(weighted_matrix(ii, :)) > 0 && norm(weighted_matrix(jj, :)) > 0
            cosine_similarity_weighted(ii, jj) = dot(weighted_matrix(ii, :), weighted_matrix(jj, :)) / ...
                                               (norm(weighted_matrix(ii, :)) * norm(weighted_matrix(jj, :)));
        else
            cosine_similarity_weighted(ii, jj) = 0; % Handle zero norm cases
        end
    end
end

% Visualize cosine similarity matrix
figCosine = figure('Color', 'white', 'Position', [500, 50, 700, 650]);
imagesc(cosine_similarity_weighted);
colormap(flipud(gray));  % Apply grayscale colormap
colorbar;
xticks(1:numel(apn2IDs_sorted));
yticks(1:numel(apn2IDs_sorted));
xticklabels(apn2IDs_sorted);
yticklabels(apn2IDs_sorted);
title('Weighted Cosine Similarity of APN2s and Sensory Inputs');
set(gca, 'XTickLabelRotation', 90);

% Add cosine similarity values to the cells
for ii = 1:numel(apn2IDs_sorted)
    for j = 1:numel(apn2IDs_sorted)
        value = cosine_similarity_weighted(ii, j);
        text_color = 'w';  % default black
        if value < 0.5     % if background is dark, switch to white
            text_color = 'k';
        end
        text(j, ii, sprintf('%.2f', value), ...
            'HorizontalAlignment', 'center', 'Color', text_color, 'FontSize', 10);
    end
end

% save figures
cd(savedFiguresPath)

saveas(figCosine,todaysDate+"_apnCosineSimilarityMatrix",'fig')
saveas(figCosine,todaysDate+"_apnCosineSimilarityMatrix",'png')
saveas(figCosine,todaysDate+"_apnCosineSimilarityMatrix",'pdf')

% Perform Pairwise Comparisons using Mann-Whitney U Test
pairwise_p_values = zeros(numel(apn2IDs_sorted));
for ii = 1:numel(apn2IDs_sorted)
    for jj = ii+1:numel(apn2IDs_sorted)
        p_mannwhitney = ranksum(cosine_similarity_weighted(ii, :), cosine_similarity_weighted(j, :));
        pairwise_p_values(ii, jj) = p_mannwhitney;
        pairwise_p_values(jj, ii) = p_mannwhitney; % Symmetric matrix
    end
end

% Display pairwise p-values
disp('Pairwise Mann-Whitney U Test p-values:');
disp(array2table(pairwise_p_values, 'VariableNames',apn2IDs_sorted, 'RowNames', apn2IDs_sorted));

%% Heatmap of adjacency matrix showing APN2 -> AMN connectivity

% extract columns from classifications and connections and table into their own array (optimizes indexing and increases speed) -> repeated code so I can run only this section if I want
allRootIDs = string(classifications.root_id);
allSides = string(classifications.side);
allPreIDs = connections.pre_root_id;
allPostIDs = connections.post_root_id;
allNTs = string(connections.nt_type);

% -- Upload DFS paths data and path connections -- % 
cd(savedVariablesPath)
apnToAmnPathsCSV = readcell("260223_APN2_AMN_paths.csv");
load("260324_categorizedConnections_apnAmnPaths");

% -- Split groups of paths into 1 path/row -- %

apnToAmnPaths = {};

numPathGroups = numel(apnToAmnPathsCSV);
step = max(1, round(numPathGroups/10)); % guarantees mod function below doesn't break if there are less than 10 rows

for p = 1:numPathGroups
    groupedPaths = apnToAmnPathsCSV{p};
    if strcmp(groupedPaths,'None')
        continue
    end

    if mod(p, step) == 0
        fprintf('%.0f %% of grouped paths have been separated\n',(p/numPathGroups)*100)
    end

    singlePaths = regexp(groupedPaths, '\[([^\[\]]+)\]', 'tokens');

    for k = 1:numel(singlePaths)
        apnToAmnPaths{end+1,1} = strsplit(singlePaths{k}{1}, ', ');
    end
end 

% -- preallocate vars -- %
apnAmnInputID(size(apnToAmnPaths,1),1) = string; % IDs of AMN inputs that are along the APN -> AMN path
dsAmnID(size(apnToAmnPaths,1),1) = string; % IDs of AMNs that are downstream APN
apnAmnNumHops(size(apnToAmnPaths,1),1) = nan;
apnAmnNumSyn(size(apnToAmnPaths,1),1) = nan;

% -- lookup table -- %

edgeKey = strcat(apnAmnPathConns.preRootID,"->",apnAmnPathConns.postRootID);

% sum duplicate edges
[G,key] = findgroups(edgeKey);
syn = splitapply(@sum,apnAmnPathConns.numSynapses,G);

% build lookup table
synLookup = containers.Map(cellstr(key),syn);

for n = 1:length(apnToAmnPaths)

    path = apnToAmnPaths{n};

    if numel(path)==1
        disp(path)
    end

    apnAmnInputID(n,1) = path{1};
    dsAmnID(n,1) = path{end};
    apnAmnNumSyn(n,1) = pathSynapses(path,synLookup);
    apnAmnNumHops(n,1) = numel(path)-1;
end

apnAmnResults = table(apnAmnInputID,dsAmnID,apnAmnNumSyn,apnAmnNumHops);
apnAmnResults = rmmissing(apnAmnResults);

% -- organize the data so they can be plotted in the desired order -- %

% sort APN2 columns by side
apnLeft  = apn2IDs(ismember(apn2IDs,...
    allRootIDs(strcmp(allSides,"left"))));

apnRight = apn2IDs(ismember(apn2IDs,...
    allRootIDs(strcmp(allSides,"right"))));

apnOrder = [apnLeft; apnRight];

% sort AMN by pair and side
ALL1 = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"left"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"anterior levator") & strcmp(amnGroups.AMN,"1"))));

ADL2 = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"left"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"anterior depressor") & strcmp(amnGroups.AMN,"2"))));

PDL3 = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"left"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"posterior depressor") & strcmp(amnGroups.AMN,"3"))));

PLL4a = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"left"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"posterior levator") & strcmp(amnGroups.AMN,"4a"))));

PLL4b = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"left"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"posterior levator") & strcmp(amnGroups.AMN,"4b"))));

ALR1 = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"right"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"anterior levator") & strcmp(amnGroups.AMN,"1"))));

ADR2 = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"right"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"anterior depressor") & strcmp(amnGroups.AMN,"2"))));

PDR3 = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"right"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"posterior depressor") & strcmp(amnGroups.AMN,"3"))));

PLR4a = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"right"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"posterior levator") & strcmp(amnGroups.AMN,"4a"))));

PLR4b = dsAmnID(ismember(dsAmnID,allRootIDs(strcmp(allSides,"right"))) ...
    & ismember(dsAmnID,amnGroups.amnID(strcmp(amnGroups.muscleName,"posterior levator") & strcmp(amnGroups.AMN,"4b"))));

amnOrder = [unique(ALL1);
                unique(ADL2);
                unique(PDL3);
                unique(PLL4a);
                unique(PLL4b);
                unique(ALR1);
                unique(ADR2);
                unique(PDR3);
                unique(PLR4a);
                unique(PLR4b)];

% count the size of each group
nALL1 = numel(unique(ALL1));
nADL2 = numel(unique(ADL2));
nPDL3 = numel(unique(PDL3));
nPLL4a = numel(unique(PLL4a));
nPLL4b = numel(unique(PLL4b));

nALR1 = numel(unique(ALR1));
nADR2 = numel(unique(ADR2));
nPDR3 = numel(unique(PDR3));
nPLR4a = numel(unique(PLR4a));
nPLR4b = numel(unique(PLR4b));

groupSizes = [nALL1 nADL2 nPDL3 nPLL4a nPLL4b ...
              nALR1 nADR2 nPDR3 nPLR4a nPLR4b];

nApnLeft = numel(apnLeft);
nApnRight = numel(apnRight);

% -- create matrices -- %

[G,apnID,amnID] = findgroups(apnAmnResults.apnAmnInputID,apnAmnResults.dsAmnID);

minHop = splitapply(@min,apnAmnResults.apnAmnNumHops,G);

apnAmnResults.group = G;
apnAmnResults.minHop = minHop(G);

shortestResults = apnAmnResults(apnAmnResults.apnAmnNumHops == apnAmnResults.minHop,:);

[G2,apnID,amnID] = findgroups(shortestResults.apnAmnInputID,...
                                shortestResults.dsAmnID);

hop = splitapply(@min,shortestResults.apnAmnNumHops,G2);
syn = splitapply(@sum,shortestResults.apnAmnNumSyn,G2);

synMat1 = nan(length(amnOrder),length(apnOrder));
synMat2 = nan(length(amnOrder),length(apnOrder));

for i = 1:length(apnID)

    r = find(amnOrder == amnID(i));
    c = find(apnOrder == apnID(i));

    if isempty(r) || isempty(c)
        continue
    end

    if hop(i) == 2
        synMat1(r,c) = syn(i);

    elseif hop(i) == 3
        synMat2(r,c) = syn(i);

    end
end

synMat1(isnan(synMat1)) = 0;
synMat2(isnan(synMat2)) = 0;

% -- plot 2-hop heatmap -- %
apnToAmnHeatmap_2hop = figure;
imagesc(synMat1)
colormap(parula)
colorbar
clim([0 50])
title("APN2 -> AMN (2 hops)",'FontSize',14)

hold on

yBreaks = cumsum(groupSizes);
for k = 1:length(yBreaks)-1
    yline(yBreaks(k)+0.5,'k','LineWidth',1.5)
end
xline(nApnLeft+0.5,'k','LineWidth',2)

% add group labels
groupLabels = {...
    'AMN1,Left,AL',...
    'AMN2,Left,AD',...
    'AMN3,Left,PD',...
    'AMN4a,Left,PL',...
    'AMN4b,Left,PL',...
    'AMN1,Right,AL',...
    'AMN2,Right,AD',...
    'AMN3,Right,PD',...
    'AMN4a,Right,PL',...
    'AMN4b,Right,PL'};

groupCenters = cumsum(groupSizes) - groupSizes/2; % calc center of each group

yticks(groupCenters)
yticklabels(groupLabels)

% apnLabels = {"A,left","B,left","C,left","D,left","E,left","F,left","unmatched,left","A,right","B,right","C,right","D,right","E,right","F,right","unmatched,right"};

set(gca,'XTick',1:1:14,'XTickLabel',apnOrder,'YAxisLocation','left');

text(nApnLeft/2, 0.25, 'Left APN2', ...
    'HorizontalAlignment','center','FontSize',12)

text(nApnLeft+nApnRight/2, 0.25, 'Right APN2', ...
    'HorizontalAlignment','center','FontSize',12)

% -- plot 3-hops heatmap -- %
apnToAmnHeatmap_3hop = figure;
imagesc(synMat2)
colormap(parula)
colorbar
clim([0 500])
title("APN2 -> AMN (3 hops)",'FontSize',14)

hold on

yBreaks = cumsum(groupSizes);
for k = 1:length(yBreaks)-1
    yline(yBreaks(k)+0.5,'k','LineWidth',1.5)
end
xline(nApnLeft+0.5,'k','LineWidth',2)

% add group labels
groupLabels = {...
    'AMN1,Left,AL',...
    'AMN2,Left,AD',...
    'AMN3,Left,PD',...
    'AMN4a,Left,PL',...
    'AMN4b,Left,PL',...
    'AMN1,Right,AL',...
    'AMN2,Right,AD',...
    'AMN3,Right,PD',...
    'AMN4a,Right,PL',...
    'AMN4b,Right,PL'};

groupCenters = cumsum(groupSizes) - groupSizes/2; % calc center of each group

yticks(groupCenters)
yticklabels(groupLabels)

set(gca,'XTick',1:1:14,'XTickLabel',apnOrder,'YAxisLocation','left');

text(nApnLeft/2, 0.25, 'Left APN2', ...
    'HorizontalAlignment','center','FontSize',12)

text(nApnLeft+nApnRight/2, 0.25, 'Right APN2', ...
    'HorizontalAlignment','center','FontSize',12)

% save figures
cd(savedFiguresPath)

saveas(apnToAmnHeatmap_2hop,todaysDate+"_apnToAmnConnMatrix_2hop",'fig')
saveas(apnToAmnHeatmap_2hop,todaysDate+"_apnToAmnConnMatrix_2hop",'png')
saveas(apnToAmnHeatmap_2hop,todaysDate+"_apnToAmnConnMatrix_2hop",'pdf')

saveas(apnToAmnHeatmap_3hop,todaysDate+"_apnToAmnConnMatrix_3hop",'fig')
saveas(apnToAmnHeatmap_3hop,todaysDate+"_apnToAmnConnMatrix_3hop",'png')
saveas(apnToAmnHeatmap_3hop,todaysDate+"_apnToAmnConnMatrix_3hop",'pdf')
