function [newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges] = nutritionAlgorithm(model,rois,roisMinMax,options)
%TO DO
%Work on food items removal point system

% Identifies the minimal changes to a diet necessary to get a desired
% change in one or more reactions of interest. One may enter a metabolite
% of the pointsModel instead of a reaction and the algorithm will optimize the
% diet with a sink or demand reaction for the corresponding metabolite of
% interest.
%
% USAGE:
%
%    [newDietModel,pointsModel,slnMin,slnMax,pointsModelSln,itemsRemoved,itemsAdded] = optimizeDiet4ROIs(pointsModel,obj,objMinMax,rois,roisMinMax,options)
%
%     Example: [newDietModel,pointsModel,roiFlux,pointsModelSln,itemsRemoved,itemsAdded] = optimizeDiet4ROIs(WBmodel,'Whole_body_objective_rxn','max',{},{})
%
% INPUTS:
%    pointsModel:          COBRA pointsModel structure with the fields:
%                      * .S
%                      * .b
%                      * .ub
%                      * .ub
%                      * .mets  (required if pointsModel.SIntRxnBool absent)
%                      * .rxns  (required if pointsModel.SIntRxnBool absent)
%
%   obj:           organism's objective function
%
%   objMinMax:     minimize ('min') or maximize ('max') objective function
%
%   rois:          cell array of all reactions of interest
%
%   roisMinMax:    cell array of 'min'/'max' entries for rois
%
% OPTIONAL INPUTS:
%   options:  Structure containing the optional specifications:
%
%       * .foodOrMets:  dictates if the algorithm adds individual
%       metabolites to the diet or food items. Default is food items.
%       "Food Cat" adjust algorithm to identify categories of food rather 
%       than specific items."AllMets" allows any dietary metabolite into 
%       the solution and "FoodMets" only allows metabolites that are in 
%       the fdTable spreadsheet into the solution.
%       Possible inputs are: "Food Items", "Food Cat", "AllMets", "FoodMets".
%
%       * .roiWeights:   a vector of weights for each reaction of interest
%       default is equal to 1
%
%       * .targetedFoodItems: A nx2 cell vector that specifies any food items 
%       to target and the corresponding weight. 
%
%       * .initObjSln: provide an initial solution for the objective
%       function. Output from optimizeWBmodel.
%
%       * .caloricRange: 1x2 vector defining boundries for diet calories
%
%       * .slnType: Specify if solution should be 'Detailed' or 'Quick'.
%                   Default setting is 'Detailed'
%
%       * .roiBound: 'Unbounded' or 'Bounded'. Default is 'Bounded'.
%
%       * .foodAddedLimit: Specify a limit for the units of food that can
%                          be added to the diet
%
%       * .foodRemovedLimit: Specify a limit for the units of food that can
%                          be removed from the diet
%
%       * .optObj: optimize current objective function as well as rois...
%       Aceptable inputs are 'True' and 'False'. Default is 'False';
%
% OUTPUT:
%    solution:       Structure containing the following fields:
%
% relaxedModel       pointsModel structure that admits a flux balance solution
%
% .. Authors: - Bronson R. Weston   2021


disp('_____________________________________________________')


% Determine if any rois are metabolites
metRois=[];
for i=1:length(rois)
    if any(strcmp(model.mets,rois{i}))
        metRois=[metRois,i];
        %         if strcmp(slnType,'Detailed')
        %             slnType='Quick';
        %             disp('slnType changed to Quick because one or more rois defined as a metabolite')
        %         end
        if strcmp(roisMinMax{i},'max')
            model=addDemandReaction(model,rois{i}); %adds demand reaction as 'DM_metabolite'
            rois{i}=['DM_',rois{i}];
            model=changeRxnBounds(model,rois{i},100000,'u');
        else
            model=addSinkReactions(model,rois(i),-100000,0);
            rois{i}=['sink_',rois{i}];
        end
    end
end

%initialize optional variables
roiWeights=10*ones(1,length(rois));
targetedFoodItems={};
slnType='Detailed';
roiBound='Bounded';
foodAddedLimit=1000000;
foodRemovedLimit=1000000;
optObj='False';

if exist('options','var')
    fn = fieldnames(options);
    for k=1:numel(fn)
        %         if( isnumeric(options.(fn{k})) )
        %             % do stuff
        %         end
        if strcmp(fn{k},'roiWeights')
            roiWeights=options.roiWeights;
            if length(roiWeights)~=length(rois)
                error('length of roiWeights vector must be the same as items rois')
            end
        elseif strcmp(fn{k},'optObj')
            optObj=options.optObj;
        elseif strcmp(fn{k},'roiBound')
            roiBound=options.roiBound;
            if ~strcmp(roiBound,'Unbounded') && ~strcmp(roiBound,'Bounded')
                error('Invalid roiBound input. Must be "Unbounded" or "Bounded"')
            end
        elseif strcmp(fn{k},'foodAddedLimit')
            foodAddedLimit=options.foodAddedLimit;
        elseif strcmp(fn{k},'foodRemovedLimit')
            foodRemovedLimit=options.foodRemovedLimit;
        elseif strcmp(fn{k},'slnType')
            slnType=options.slnType;
            if ~strcmp(slnType,'Detailed') && ~strcmp(slnType,'Quick')
                error('Invalid slnType input. Must be "Detailed" or "Quick"')
            end
        elseif strcmp(fn{k},'targetedFoodItems')
            targetedFoodItems=options.targetedFoodItems;
        else
            error(['Invalid "options" field entered: ', fn{k}])
        end
    end
end

if any(roiWeights<=0)
    error('"roiWeights" can not contain entries less than or equal to zero')
end

%adjust ub and lb if roiBound specifies 'Unbound'
if strcmp(roiBound, 'Unbounded')
    for i=1:length(rois)
        f=find(strcmp(model.rxns,rois{i}));
        if strcmp(roisMinMax{i},'max')
            if model.ub(f)~=0
                model.ub(f)=100000;
            end
        else
            if model.lb(f)~=0
                model.lb(f)=-100000;
            end
        end
    end
end


newDietModel=model; %Copy original instance of model for new diet pointsModel
pointsModel=model; %Copy original instance of model for points pointsModel


%Calculate newDietModel objective function and restrict obj in main pointsModel
objIndex=find(model.c==1);

roiIndexO=zeros(1,length(rois));
for i=1:length(roiIndexO)
    roiIndexO(i)=find(strcmp(newDietModel.rxns,rois{i}));
    disp(['Reaction of Interest ', num2str(i),' = ', newDietModel.rxns{roiIndexO(i)}])
end

% get flux of objective function

if model.ub(objIndex)~=model.lb(objIndex) && strcmp(optObj,'True')
    model_Obj = optimizeWBModel(newDietModel);
    f1=model_Obj.f;
    initRoiFlux=model_Obj.v(roiIndexO);
else
    f1=model.ub(objIndex);
    initRoiFlux=NaN(1,length(rois));
end




% If sln type is detailed, check if roi is already min or maxed out and
% if not, define min max range for roi

if strcmp(slnType,'Detailed')
    OroiFluxMin=[];
    OroiFluxMax=[];
    for i=1:length(rois)
        if initRoiFlux(i)==newDietModel.lb(roiIndexO)
            OroiFluxMin(i)=newDietModel.lb(roiIndexO(i));
        else
            pointsModel = changeObjective(pointsModel,rois{i});
            pointsModel.osenseStr = 'min';
            sln = optimizeCbModel(pointsModel);
            if isnan(sln.f)
                warning('Not a viable initial solution, recommend adding nutrients to initial diet');
                OroiFluxMin(i)=NaN;
                continue
            end
            OroiFluxMin(i)=sln.v(roiIndexO(i));
        end
        if initRoiFlux(i)==newDietModel.ub(roiIndexO)
            OroiFluxMax(i)=newDietModel.ub(roiIndexO(i));
        else
            pointsModel = changeObjective(pointsModel,rois{i});
            pointsModel.osenseStr = 'max';
            sln = optimizeCbModel(pointsModel);
            if isnan(sln.f)
                warning('Not a viable initial solution, recommend adding nutrients to initial diet');
                OroiFluxMax(i)=NaN;
                continue
            end
            OroiFluxMax(i)=sln.v(roiIndexO(i));
        end
    end
end

pointsModel=addMetabolite(pointsModel, 'unitOfFoodAdded[dP]');
pointsModel=addMetabolite(pointsModel, 'unitOfFoodRemoved[dP]');
pointsModel=addMetabolite(pointsModel, 'unitOfFoodChange[dP]');
pointsModel=addMetabolite(pointsModel, 'roiPoint[roiP]');
pointsModel=addMetabolite(pointsModel, 'point[P]');

%If necessary, add all diet exchange reactions to targetedFoodItems
if isempty(targetedFoodItems)
    %Add all Diet_EX reactions, and set weight to 1
    dietRxns=find(contains(pointsModel.rxns,'Diet_EX'));
    targetedFoodItems=[pointsModel.rxns(dietRxns),num2cell(ones(length(dietRxns),1))];
elseif any(strcmp(targetedFoodItems(:,1),'All')) && length(targetedFoodItems(:,1))>1
    dietRxns=find(contains(pointsModel.rxns,'Diet_EX'));
    targetedFoodItemsTemp=[pointsModel.rxns(dietRxns), ... 
        num2cell(cell2mat(targetedFoodItems(strcmp(targetedFoodItems(:,1),'All'),2))*ones(length(dietRxns),1))];
    [~,ai,bi]=intersect(targetedFoodItemsTemp(:,1),targetedFoodItems(:,1));
    targetedFoodItemsTemp(ai,2)=targetedFoodItems(bi,2);
    targetedFoodItems=targetedFoodItemsTemp;
elseif any(strcmp(targetedFoodItems(:,1),'All'))
    dietRxns=find(contains(pointsModel.rxns,'Diet_EX'));
    targetedFoodItems=[pointsModel.rxns(dietRxns),num2cell(ones(length(dietRxns),1))];    
end

%Add "Food_Added" reactions to points model

Mets= pointsModel.mets;
[~,ai,bi]=intersect(pointsModel.rxns,targetedFoodItems(:,1));
foodRxns=pointsModel.rxns(ai);
foodRxns=regexprep(foodRxns,'Diet_EX_','Food_Added_EX_');
sMatrix=pointsModel.S(:,ai);
f=find(strcmp(pointsModel.mets,'unitOfFoodAdded[dP]'));
A=-1*cell2mat(targetedFoodItems(bi,2)).';
sMatrix(f,:)=A;
pointsModel = addMultipleReactions(pointsModel, foodRxns, Mets, sMatrix, 'lb', -100000*ones(1,length(foodRxns)), 'ub', zeros(1,length(foodRxns)));
pointsModel = addMultipleReactions(pointsModel, {'Point_EX_unitOfFoodRemoved2Change[dp]','Point_EX_unitOfFoodAdded2Change[dp]','Point_EX_unitOfFoodChange[dP]_[P]','Point_EX_Point[P]'}, {'unitOfFoodRemoved[dP]','unitOfFoodAdded[dP]','unitOfFoodChange[dP]','point[P]'}, [-1 0 0 0;0 -1 0 0;1 1 -1 0;0 0 1 -1], 'lb', [-1000000,-1000000, -1000000,-1000000], 'ub', [foodRemovedLimit,foodAddedLimit,1000000,1000000]);
foodAddedLimit
%Add "Food Removed" reactions to points model
foodDietIndex=find(contains(pointsModel.rxns,'Diet_EX_'));
foodDietIndex=foodDietIndex(pointsModel.lb(foodDietIndex)<0);
foodItems=pointsModel.rxns(foodDietIndex);
foodItems=regexprep(foodItems,'Diet_EX_','Food_Removed_EX_');

sMatrix=-1*pointsModel.S(:,foodDietIndex);
f=find(contains(pointsModel.mets,'unitOfFoodRemoved[dP]'));
A=-1*ones(1,length(foodDietIndex));
sMatrix(f,:)=A;
pointsModel = addMultipleReactions(pointsModel, foodItems, pointsModel.mets, sMatrix, 'lb', -100000*ones(1,length(foodDietIndex)), 'ub', zeros(1,length(foodDietIndex)));


%Get roi Indexes
for i=1:length(rois)
    roiIndexP(i)=find(strcmp(pointsModel.rxns,rois{i}));
end
roiUB=pointsModel.ub(roiIndexP);
roiLB=pointsModel.lb(roiIndexP);

%replace roi function
stoich=pointsModel.S(:,roiIndexP);
if length(roiIndexP)>1
    metInd=find(any(stoich.'~=0));
    metsRoi=pointsModel.mets(any(stoich.'~=0)).';
    metsStoich=full(stoich(any(stoich.'~=0),:));
else
    metsRoi=pointsModel.mets(find(stoich~=0)).';
    metsStoich=full(stoich(find(stoich~=0)));
end

weightVector=zeros(1,length(roiIndexP));
weightVector(contains(roisMinMax,'max'))=-1;
weightVector(contains(roisMinMax,'min'))=1;
metsStoich=[metsStoich;weightVector.*roiWeights;zeros(1,length(roiIndexP))];

for i=1:length(rois)
    evalc('[pointsModel,~,~]= removeRxns(pointsModel, rois{i})');
% [pointsModel,~,~]= removeRxns(pointsModel, rois{i});
end

pointsModel = addMultipleReactions(pointsModel, [rois,'Point_EX_roiPoints[roiP]_[P]'], [metsRoi,'roiPoint[roiP]','point[P]'], [metsStoich,[zeros(length(metsStoich(:,1))-2,1);-1;1]], 'lb', [roiLB.',-1000000], 'ub', [roiUB.',1000000]);



%Find solution
pointsModel = changeObjective(pointsModel,'Point_EX_Point[P]');
pointsModel.osenseStr = 'min';
pointsModelSln = optimizeWBModel(pointsModel);
% pointsModelSln.v(roiIndex)

disp(['Solution points =',num2str(pointsModelSln.f)])
disp([num2str(pointsModelSln.v(find(strcmp(pointsModel.rxns,'Point_EX_unitOfFoodChange[dP]_[P]')))),' come from diet']);
disp([num2str(pointsModelSln.v(find(strcmp(pointsModel.rxns,'Point_EX_roiPoints[roiP]_[P]')))),' come from roi']);
foodAddedIndexes=find(contains(pointsModel.rxns,'Food_Added_EX_'));
foodRemovedIndexes=find(contains(pointsModel.rxns,'Food_Removed_EX_'));
slnIndexes1=foodAddedIndexes(pointsModelSln.v(foodAddedIndexes)<0);
slnIndexes2=foodRemovedIndexes(pointsModelSln.v(foodRemovedIndexes)<0);

disp('Food items of interest are:')
T=table([pointsModel.rxns(slnIndexes1);pointsModel.rxns(slnIndexes2)],pointsModelSln.v([slnIndexes1;slnIndexes2]),'VariableNames',{'Food Rxn', 'Flux'})


%Add and remove relevant food items from diet in newDietModel
foodItemsAdd= regexprep(pointsModel.rxns(slnIndexes1),'Food_Added_EX_','Diet_EX_');
foodItemsRemove= regexprep(pointsModel.rxns(slnIndexes2),'Food_Removed_EX_','Diet_EX_');
modelOindexAdd=zeros(1,length(foodItemsAdd));
sl2IndexAdd=zeros(1,length(foodItemsAdd));
modelOindexRemove=zeros(1,length(foodItemsRemove));
sl2IndexRemove=zeros(1,length(foodItemsRemove));
for i=1:length(foodItemsAdd)
    modelOindexAdd(i)=find(contains(newDietModel.rxns,foodItemsAdd(i)));
    sl2IndexAdd(i)=find(contains(pointsModel.rxns,foodItemsAdd(i)));
end
for i=1:length(foodItemsRemove)
    modelOindexRemove(i)=find(contains(newDietModel.rxns,foodItemsRemove(i)));
    sl2IndexRemove(i)=find(contains(pointsModel.rxns,foodItemsRemove(i)));
end
newDietModel.lb(modelOindexAdd)=(pointsModelSln.v(sl2IndexAdd)+pointsModelSln.v(slnIndexes1))*1.01;
newDietModel.ub(modelOindexAdd)=(pointsModelSln.v(sl2IndexAdd)+pointsModelSln.v(slnIndexes1))*0.99;
newDietModel.lb(modelOindexRemove)=(pointsModelSln.v(sl2IndexRemove)-pointsModelSln.v(slnIndexes2))*1.01;
newDietModel.ub(modelOindexRemove)=(pointsModelSln.v(sl2IndexRemove)-pointsModelSln.v(slnIndexes2))*0.99;



if strcmp(slnType,'Quick')
    for i=1:length(rois)
        ind=find(strcmp(pointsModel.rxns,rois{i}));
        disp([rois{i},' flux = ', num2str(pointsModelSln.v(ind))])
    end
    roiFlux(i)=pointsModelSln.v(ind);
    slnRanges=pointsModelSln.v(roiIndexP);
    menuChanges=T;
    return
end


%Find new obj flux with new diet

if model.ub(objIndex)~=model.lb(objIndex)
    model_Obj = optimizeWBModel(newDietModel);
    f2=model_Obj.f;
    newDietModel=changeRxnBounds(newDietModel,obj,f2,'b'); %constrain pointsModel obj flux
else
    f2=model.ub(objIndex);
end

disp(['f1 =',num2str(f1), ' & f2=', num2str(f2)])

%Compute new min max ranges for roi with new diet
%%
for i=1:length(rois)
    disp(rois{i})
    newDietModel = changeObjective(newDietModel,rois{i});
    newDietModel.osenseStr = 'min';
    tmp=optimizeWBModel(newDietModel);
    slnMin.(['Rxn',num2str(i)]) = tmp;
    NroiFluxMin(i)=slnMin.(['Rxn',num2str(i)]).v(roiIndexO(i));
    newDietModel.osenseStr = 'max';
    tmp=optimizeWBModel(newDietModel);
    slnMax.(['Rxn',num2str(i)]) = tmp;
    NroiFluxMax(i)=slnMax.(['Rxn',num2str(i)]).v(roiIndexO(i));
    disp(['Original Diet RoI range = ', num2str(OroiFluxMin(i)), ':', num2str(OroiFluxMax(i))])
    disp(['New Diet RoI range = ', num2str(NroiFluxMin(i)), ':', num2str(NroiFluxMax(i))])
end

% slnMin
roiFlux=[slnMin.',slnMax.'];
newDietModel.ub(objIndex)=model.ub(objIndex);
newDietModel.lb(objIndex)=model.lb(objIndex);
menuChanges.itemsRemoved=[];
menuChanges.itemsAdded=[];
menuChanges.fullmenu=[];
newDietModel = changeObjective(newDietModel,obj);
newDietModel.osenseStr = objMinMax;

disp('___________________________________________________________________')


end

%TO DO:
% -Impliment itemsRemoved and itemsAdded
% -Introduce option for metabolite adjustments
% -Include intelligent default weighting for roi
% -slnMin & slnMax turn to structs or convert to flux values

%Finished:
%-introduced Quick/Detailed functionality
%-renamed variables to be more obvious
%-commented code

% figure()
% t = tiledlayout(2,2,'TileSpacing','compact');
% ax1 = nexttile;
% pie(ax1,OgMacros(1:4))
% title('Original Diet')
% ax2 = nexttile;
% x=pie(ax2,NewMacros(1:4))
% title('New Diet')
% legend(NewCategories(1:4),'Location','South')
% ax3 = nexttile;
% bar(ax3,(NewMacros(1:4)-OgMacros(1:4))./OgMacros(1:4))
% ylim([-1 14])
% ax4= nexttile;
% legend(x,NewCategories(1:4))