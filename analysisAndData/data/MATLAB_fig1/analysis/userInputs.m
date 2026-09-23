%% user inputs

% a string of the date, used for naming figures and variables when they are saved
todaysDate = "260916";

% strings of the names of the cell ID files
apn2IDs_fileName = "APN2_IDs.txt";
jonIDs_fileName = "260320_JON_IDs.txt";
bnIDs_fileName = "BN_IDs_antNv.txt";
hbnIDs_fileName = "BN_IDs_head.txt";
dnxIDs_fileName = "DNx01_IDs.txt";
amnIDs_fileName = "AMN_IDs.txt";

% strings of the file names of the CSV files obtained from running the python code that runs the depth first search
apnAmnPathsCSV_fileName = "260223_APN2_AMN_paths.csv";
jonApnPathsCSV_fileName = "260320_JON_APN2_paths.csv";
bnApnPathsCSV_fileName = "260223_BN_APN2_paths.csv";
hbnApnPathsCSV_fileName = "260730_hBN_APN2_paths.csv";
dnxApnPathsCSV_fileName = "260730_DNx01_APN2_paths.csv";

% Toggle to load the following variables into each section (=1), or just run the script fresh (=0)
loadVars = 1;

% strings of the file names of the variables saved after running the code (later in the script) 
% that separates the CSV file paths. 
% if you haven't obtained these variables yet, set the variable to empty brackets with no quotation marks ([]);

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

% strings of the file names of the variables saved after running the code (later in the script)
% creates the directed graph between to two cell types
apnAmnDirGraph_fileName = "260324_digraph_apnAmnPaths.mat";
jonApnDirGraph_fileName = "260324_digraph_jonApnPaths.mat";
bnApnDirGraph_fileName = "260324_digraph_bnApnPaths.mat";

% mac OS paths
% connectomicsFolderPath = "/Users/kayleeyeager/Library/CloudStorage/OneDrive-Vanderbilt/Suver Lab/APN2 Connectomics";
% dataUploadPath = "/Users/kayleeyeager/Library/CloudStorage/OneDrive-Vanderbilt/Suver Lab/APN2 Connectomics/data";
% savedVariablesPath = "/Users/kayleeyeager/Library/CloudStorage/OneDrive-Vanderbilt/Suver Lab/APN2 Connectomics/savedVariables";
% savedFiguresPath = "/Users/kayleeyeager/Library/CloudStorage/OneDrive-Vanderbilt/Suver Lab/APN2 Connectomics/savedFigures";


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