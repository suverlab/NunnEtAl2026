%%
% The CSV table obtained from running the DFS algorithm on Python is structured weird. 
% Each row contains all the paths for 1 start node and 1 end node combination. 
% So, each row has a bunch of bracketed lists within a giant bracket ('nested' brackets).
%%
% There are 2 main parts to this function:
	%Part 1: separates the nested brackets so that each row contains only 1 path
	%Part 2: separates each path into a list of connections, containing pre- and post-IDs (structured just like the Codex connections dataset)
%%
% When this is finished, all of the pre- and post-IDs for all paths are saved 
% in a table called "pathConnections" (which can be named however you wish when you call the function forward). 
% The table contains the pre- and post-ID of each connection as well as 
% their associated meta data that is found in the Codex connections dataset (e.g. NT type and synapse count).

%%
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