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