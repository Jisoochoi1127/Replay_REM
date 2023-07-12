%%
load('spatial_analysis.mat');
%% bin analysis

isPlaceCell= cellfun(@(c) c.IsPlaceCell, spatial_analysis.bin, 'UniformOutput', false);
isPlaceCell=cell2mat(isPlaceCell);

noplace=[];
right=[];
left=[];
bisum=[];
bi=[];

for c=1:length(isPlaceCell)
    if isPlaceCell(c,1)==0 && isPlaceCell(c,2)==0 && isPlaceCell(c,3)==0
        noplace=[noplace; c];
    elseif isPlaceCell(c,1)==1 && isPlaceCell(c,2)==0
        right=[right, c];
    elseif isPlaceCell(c,2)==1 && isPlaceCell(c,1)==0
        left=[left; c];
    elseif isPlaceCell(c,2)==1 && isPlaceCell(c,1)==1 && isPlaceCell(c,3)==0 
        bisum=[bisum; c];
    else
        bi=[bi; c];
    end
end

numPlaceFields= cellfun(@(c) c.numPlaceFields, spatial_analysis.bin, 'UniformOutput', false);
PlaceFieldPeakJointProbability= cellfun(@(c) c.PlaceFieldPeakJointProbability, spatial_analysis.bin, 'UniformOutput', false);
PlaceFieldmeanJointProbability= cellfun(@(c) c.PlaceField, spatial_analysis.bin, 'UniformOutput', false);
PlaceFieldmeanJointProbability=  cellfun(@(x) mean(nonzeros(x)), PlaceFieldmeanJointProbability, 'UniformOutput', false);
PlaceFieldStability= cellfun(@(c) c.PlaceFieldStability, spatial_analysis.bin, 'UniformOutput', false);
InFieldActivityRatio= cellfun(@(c) c.InFieldActivityRatio, spatial_analysis.bin, 'UniformOutput', false); 
SI= cellfun(@(c) c.SI, spatial_analysis.bin, 'UniformOutput', false);
sparsity= cellfun(@(c) c.sparsity, spatial_analysis.bin, 'UniformOutput', false);
CellFiringProb= cellfun(@(c) c.CellFiringProb, spatial_analysis.bin, 'UniformOutput', false);
PlaceFieldArea= cellfun(@(c) c.PlaceFieldArea, spatial_analysis.bin, 'UniformOutput', false);

v={right, left, bi, bisum, noplace};
for vi=1:length(v)
    direction =v{vi};
    PFArea_temp=[];
    if vi<4
        numPF{vi}=direction_sum(numPlaceFields, direction,vi);
        PFpeakp{vi}=direction_sum(PlaceFieldPeakJointProbability, direction,vi);
        PFmeanp{vi}=direction_sum(PlaceFieldmeanJointProbability, direction,vi);
        stability{vi}=direction_sum(PlaceFieldStability, direction,vi);
        InField{vi}=direction_sum(InFieldActivityRatio, direction,vi);
        spatial_info{vi}=direction_sum(SI, direction,vi);
        sparse{vi}=direction_sum(sparsity, direction,vi); 
        totalFiring{vi}=direction_sum(CellFiringProb, direction,vi);
        PFArea_temp=direction_sum2(PlaceFieldArea, direction,vi);
        PFArea{vi}= cellfun(@mean, PFArea_temp, 'UniformOutput', false);
    elseif vi==4
        numPF{vi-1}=[numPF{vi-1}; max(direction_sum(numPlaceFields, direction,vi-3),direction_sum(numPlaceFields, direction,vi-2)) ];
        PFpeakp{vi-1}=[PFpeakp{vi-1}; mean([direction_sum(PlaceFieldPeakJointProbability, direction,vi-3), direction_sum(PlaceFieldPeakJointProbability, direction,vi-2)],2)];
        PFmeanp{vi-1}=[PFmeanp{vi-1}; mean([direction_sum(PlaceFieldmeanJointProbability, direction,vi-3), direction_sum(PlaceFieldmeanJointProbability, direction,vi-2)],2)];
        stability{vi-1}=[stability{vi-1}; mean([direction_sum(PlaceFieldStability, direction,vi-3), direction_sum(PlaceFieldStability, direction,vi-2)],2)];
        InField{vi-1}=[InField{vi-1}; mean([direction_sum(InFieldActivityRatio, direction,vi-3), direction_sum(InFieldActivityRatio, direction,vi-2)],2)];
        spatial_info{vi-1}=[spatial_info{vi-1}; mean([direction_sum(SI, direction,vi-3), direction_sum(SI, direction,vi-2)],2)];
        sparse{vi-1}=[sparse{vi-1}; mean([direction_sum(sparsity, direction,vi-3),direction_sum(sparsity, direction,vi-2)],2)]; 
        totalFiring{vi-1}=[totalFiring{vi-1}; direction_sum(CellFiringProb, direction,vi-1)];
        PFArea_temp=direction_sum2(PlaceFieldArea, direction,vi-3);
        PFArea_temp2=direction_sum2(PlaceFieldArea, direction,vi-2);
        PFArea_temp3= cell2mat([cellfun(@mean, PFArea_temp, 'UniformOutput', false), cellfun(@mean, PFArea_temp2, 'UniformOutput', false)]);
        PFArea{vi-1}= [PFArea{vi-1}; mean(PFArea_temp3,2)];
    else
        totalFiring{vi-1}=direction_sum(CellFiringProb, direction,vi-2);
    end
end


% SA.PlaceCell_Right=right;
% SA.PlaceCell_Left=left;
% SA.PlaceCell_Bi=bi;
% SA.NonplaceCell=noplace;

SA.numPlaceFields=numPF;        
SA.PlaceFieldPeakJointProbability=PFpeakp;         
SA.PlaceFieldMeanJointProbability=PFmeanp;        
SA.PlaceFieldStability=stability;       
SA.InFieldActivityRatio=InField;        
SA.SI=spatial_info;        
SA.sparsity=sparse;                 
SA.PlaceFieldArea=cellfun(@(c) cell2mat(c), PFArea, 'UniformOutput', false);
SA.CellFiringProb=totalFiring; 
SA.classif_index=v;

clearvars -except SA
save('spatial_analysis_classif.mat', 'SA');
%%
numplace=(cellfun(@length,SA.classif_index));
numplace=[numplace(1),numplace(2),numplace(3)+numplace(4),sum(numplace(1:4)),numplace(5)];
ratioplace=100*(numplace./sum(cellfun(@length ,SA.classif_index)));

fields=fieldnames(SA);
results=[];
for f=1:length(fields)-2
    results(:,f)= [cell2mat(SA.(fields{f})(1));  cell2mat(SA.(fields{f})(2));  cell2mat(SA.(fields{f})(3))];        
end
cfp=[cell2mat(SA.CellFiringProb(1));  cell2mat(SA.CellFiringProb(2));  cell2mat(SA.CellFiringProb(3))];
%%
fields=fieldnames(SA);
meanresults=zeros(length(fields)-2,4);
for f=1:length(fields)-2
    meanresults(f,1:3)=cellfun(@mean, SA.(fields{f}));
    meanresults(f,4)= mean([cell2mat(SA.(fields{f})(1));  cell2mat(SA.(fields{f})(2));  cell2mat(SA.(fields{f})(3))]);        
end
cfp= cellfun(@mean,SA.CellFiringProb);
cfp= [cfp(1:3), mean([cell2mat(SA.CellFiringProb(1));  cell2mat(SA.CellFiringProb(2));  cell2mat(SA.CellFiringProb(3))]), cfp(4)];        

numplace=(cellfun(@length,SA.classif_index));
numplace=[numplace(1),numplace(2),numplace(3)+numplace(4),sum(numplace(1:4)),numplace(5)];
ratioplace=100*(numplace./sum(cellfun(@length ,SA.classif_index)));

%%

function out=direction_sum(variable, direction_vector, col)
    variable=variable(direction_vector, col);
    out=cell2mat(variable);
end
function out=direction_sum2(variable, direction_vector, col)
    variable=variable(direction_vector,col);
    out= cellfun(@(c) cell2mat(c), variable, 'UniformOutput', false);
end

