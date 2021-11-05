% Must have COBRA Toolbox initialized to run
% Can switch LPSolver  to another solver (e.g., 'gurobi') if ibm_cplex is
% not installed

%TO DO: capping food limits does not work properly when food items are
%weighted. Fix this next (see Test 5)
%TO DO: Investigate 'Any' functionality

clear options
load('salarecon_bigg_curated.mat');
ref=pwd;
ref=strsplit(ref,'Salmon_Analysis');
addpath(ref{1});
LPSolver = 'ibm_cplex';
if ~exist('solverInstalled','var')
    [solverOK, solverInstalled] = changeCobraSolver(LPSolver, 'LP',0,1);
end

ExRxnNum=find(contains(model.rxns, 'EX_'));
ExRxns= strcat('Diet_',model.rxns(ExRxnNum));
ExRxns = regexprep(ExRxns,'_e','\[e\]');
ExRxns{12}='Diet_EX_etoh[e]';
model.rxns(ExRxnNum)=ExRxns;
model.lb(ExRxnNum)=0;
model=changeRxnBounds(model,{'Biomass'},1,'b');
% [modelO,pointsModel,sl2,foodMenu]=findViableDiet(model,'Met',{},'A');
model=changeRxnBounds(model,{'Biomass'},1000,'u');

for i=1:length(ExRxns)
    strings=strsplit(ExRxns{i},'Diet_EX_');
    [MW, Ematrix] = computeMW(model, strings{end});
    if Ematrix(2)>0
        model.lb(ExRxnNum(i))=0;
    else
        model.lb(ExRxnNum(i))=-1000;
    end
end
model=changeRxnBounds(model,'Diet_EX_chol[e]',-1000,'l');
model=changeRxnBounds(model,'Diet_EX_nh4[e]',0,'l');
model=changeRxnBounds(model,'Diet_EX_so4[e]',0,'l');
model=changeRxnBounds(model,'Diet_EX_cysi__L[e]',0,'l');

feed=readtable('feed_AA_composition.csv');
feeds=feed.Properties.VariableNames;

for mol=1:length(feed{:,1})
    molecule=strcat(feed{mol,1},'[c]');
    if strcmp(feed{mol,1},'weight')
        continue
    end
    [MW, ~] = computeMW(model, molecule);
    for i=2:length(feeds)
        feed{mol,i}=feed{mol,i}/MW;
    end
end

writetable(feed,'feed_composition_mmolpergram.csv');
feedRxnso=feed.Properties.VariableNames(2:end);
feedRxns=strcat('Diet_EX_',feedRxnso);
feedRxns=strcat(feedRxns,'[c]');
feedMets=strcat(feed{1:end-1,1},'[c]');
feedMatrix=-1*[feed.(feedRxnso{1}),feed.(feedRxnso{2}),feed.(feedRxnso{3})];
feedMatrix(end,:)=[];
model = addMultipleReactions(model, feedRxns, feedMets, feedMatrix,...
    'lb', [0,0,0], 'ub', [0,0,0]);

test=0;
%% Test 1
test=test+1;
disp(['Test ' num2str(test) '->'])
options.slnType='Quick';
options.roiWeights=[1];
% options.targetedFoodItems={'Diet_EX_FM[c]',1.54;'Diet_EX_SBM[c]',0.37;'Diet_EX_BSF[c]',3} ;
options.targetedFoodItems={'All'};
% options.targetedFoodItems={'Diet_EX_but[e]',1;'Diet_EX_chol[e]',1;'Diet_EX_so4[e]',1};
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(model,{'Biomass'},{'max'},options);

%% Test 2
test=test+1;
disp(['Test ' num2str(test) '->'])

options.slnType='Quick';
options.roiWeights=[10];
% options.targetedFoodItems={'Diet_EX_FM[c]',1.54;'Diet_EX_SBM[c]',0.37;'Diet_EX_BSF[c]',3} ;
options.targetedFoodItems={'All'};
% options.targetedFoodItems={'Diet_EX_but[e]',1;'Diet_EX_chol[e]',1;'Diet_EX_so4[e]',1};
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(model,{'Biomass'},{'max'},options);

%% Test 3
test=test+1;
disp(['Test ' num2str(test) '->'])

options.slnType='Quick';
options.roiWeights=[10];
options.targetedFoodItems={'Diet_EX_FM[c]',1;'Diet_EX_SBM[c]',1;'Diet_EX_BSF[c]',1} ;
% options.foodAddedLimit=40;
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(model,{'Biomass'},{'max'},options);

%% Test 4
test=test+1;
disp(['Test ' num2str(test) '->'])

options.slnType='Quick';
options.roiWeights=[10];
options.targetedFoodItems={'Diet_EX_FM[c]',1;'Diet_EX_SBM[c]',1;'Diet_EX_BSF[c]',1} ;
options.foodAddedLimit=30;
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(model,{'Biomass'},{'max'},options);

%% Test 5

test=test+1;
disp(['Test ' num2str(test) '->'])

clear options
options.slnType='Quick';
options.roiWeights=[10];
options.foodAddedLimit=30;
options.targetedFoodItems={'Diet_EX_FM[c]',1.54;'Diet_EX_SBM[c]',0.37;'Diet_EX_BSF[c]',3} ;
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(model,{'Biomass'},{'max'},options);

%TO DO: capping food limits does not work properly when food items are
%weighted. Fix this in future

%% Test 6
test=test+1;
disp(['Test ' num2str(test) '->'])

clear options
options.slnType='Quick';
options.roiWeights=[10];
% options.foodAddedLimit=30;
fdItems=find(contains(model.rxns,'Diet_EX_'));
fdItems=fdItems(find(~contains(model.rxns(fdItems),'[c]')));

options.targetedFoodItems=[model.rxns(fdItems),num2cell(1*ones(length(fdItems),1))];
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(model,{'Biomass'},{'max'},options);

%TO DO: capping food limits does not work properly when food items are
%weighted. Have to fix this

%% Test 7
test=test+1;
disp(['Test ' num2str(test) '->'])

clear options
options.slnType='Quick';
options.roiWeights=[10];
% options.foodAddedLimit=30;
fdItems=find(contains(model.rxns,'Diet_EX_'));
fdItems=fdItems(find(~contains(model.rxns(fdItems),'[c]')));
options.targetedFoodItems=[model.rxns(fdItems),num2cell(30*ones(length(fdItems),1))];
model.lb(find(contains(model.rxns,'Diet_EX_FM[c]')))=-40;
model.ub(find(contains(model.rxns,'Diet_EX_FM[c]')))=-39;
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(model,{'Biomass'},{'max'},options);

%% Test 8
test=test+1;
disp(['Test ' num2str(test) '->'])

clear options
options.slnType='Quick';
options.roiWeights=[10];
% options.foodAddedLimit=30;
fdItems=find(contains(model.rxns,'Diet_EX_'));
fdItems=fdItems(find(~contains(model.rxns(fdItems),'[c]')));
options.targetedFoodItems=[model.rxns(fdItems),num2cell(20*ones(length(fdItems),1))];
model.lb(find(contains(model.rxns,'Diet_EX_FM[c]')))=-40;
model.ub(find(contains(model.rxns,'Diet_EX_FM[c]')))=-39;
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(model,{'Biomass'},{'max'},options);

%% Test 9
test=test+1;
disp(['Test ' num2str(test) '->'])
clear options
options.slnType='Quick';
options.roiWeights=[10];
% options.foodAddedLimit=30;
fdItems=find(contains(model.rxns,'Diet_EX_'));
fdItems=fdItems(find(~contains(model.rxns(fdItems),'[c]')));
options.targetedFoodItems=[model.rxns(fdItems),num2cell(13.44*ones(length(fdItems),1))];
model.lb(find(contains(model.rxns,'Diet_EX_FM[c]')))=-40;
model.ub(find(contains(model.rxns,'Diet_EX_FM[c]')))=-39;
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(model,{'Biomass'},{'max'},options);
model.lb(find(contains(model.rxns,'Diet_EX_FM[c]')))=0;
model.ub(find(contains(model.rxns,'Diet_EX_FM[c]')))=0;
%% Test 10
test=test+1;
disp(['Test ' num2str(test) '->'])

clear options
options.slnType='Quick';
options.roiWeights=[5];
% options.foodAddedLimit=30;
fdItems=find(contains(model.rxns,'Diet_EX_'));
fdItems=fdItems(find(~contains(model.rxns(fdItems),'[c]')));
options.targetedFoodItems=[model.rxns(fdItems),num2cell(13.44*ones(length(fdItems),1))];
model.lb(find(contains(model.rxns,'Diet_EX_BSF[c]')))=-40;
model.ub(find(contains(model.rxns,'Diet_EX_BSF[c]')))=-39;
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(model,{'Biomass'},{'max'},options);

%% Test 11
test=test+1;
disp(['Test ' num2str(test) '->'])
clear options
options.slnType='Quick';
options.roiWeights=[10];
options.targetedFoodItems={'Diet_EX_FM[c]',1.54;'Diet_EX_SBM[c]',0.37;'Diet_EX_BSF[c]',3} ;
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(model,{'Biomass'},{'max'},options);

