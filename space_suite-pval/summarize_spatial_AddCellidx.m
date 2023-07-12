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
% added by Jisoo
PlaceFieldCentroid= cellfun(@(c) c.PlaceFieldCentroid, spatial_analysis.bin, 'UniformOutput', false);

a=cell2mat(PlaceFieldStability);
Stability_Whole_cell=max(a,[],2);



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

%Added by Jisoo
SA.WholePlaceCell= [right';left;bi];
SA.PlaceCell_Right=right';
SA.PlaceCell_Left=left;
SA.PlaceCell_Bi=bi;
SA.NonplaceCell=noplace;

SA.PlaceFieldCentroid=PlaceFieldCentroid;
SA.numPlaceFields=numPF;        
SA.PlaceFieldPeakJointProbability=PFpeakp;         
SA.PlaceFieldMeanJointProbability=PFmeanp;        
SA.PlaceFieldStability=stability; 
SA.Stability_Whole_cell=Stability_Whole_cell;
SA.InFieldActivityRatio=InField;        
SA.SI=spatial_info;        
SA.sparsity=sparse;                 
SA.PlaceFieldArea=cellfun(@(c) cell2mat(c), PFArea, 'UniformOutput', false);
SA.CellFiringProb=totalFiring; 
SA.classif_index=v;

%%
mean_max_joint_R=mean(SA.PlaceFieldPeakJointProbability{1, 1});
mean_max_joint_L=mean(SA.PlaceFieldPeakJointProbability{1, 2});
mean_max_joint_B=mean(SA.PlaceFieldPeakJointProbability{1, 3});

mean_SI_R=mean(SA.SI{1, 1});
mean_SI_L=mean(SA.SI{1, 2});
mean_SI_B=mean(SA.SI{1, 3});

mean_Infield_R=mean(SA.InFieldActivityRatio{1, 1});
mean_Infield_L=mean(SA.InFieldActivityRatio{1, 2});
mean_Infield_B=mean(SA.InFieldActivityRatio{1, 3});

mean_Placefield_area_R=mean(SA.PlaceFieldArea{1, 1});
mean_Placefield_area_L=mean(SA.PlaceFieldArea{1, 2});
mean_Placefield_area_B=mean(SA.PlaceFieldArea{1, 3});

NumbPlace_R=length(SA.PlaceCell_Right);
NumbPlace_L=length(SA.PlaceCell_Left);
NumbPlace_B=length(SA.PlaceCell_Bi);

SA.mean_max_joint_R=mean_max_joint_R;
SA.mean_max_joint_L=mean_max_joint_L;
SA.mean_max_joint_B=mean_max_joint_B;

SA.mean_SI_R=mean_SI_R;
SA.mean_SI_L=mean_SI_L;
SA.mean_SI_B=mean_SI_B;

SA.mean_Infield_R=mean_Infield_R;
SA.mean_Infield_L=mean_Infield_L;
SA.mean_Infield_B=mean_Infield_B;

SA.mean_Placefield_area_R=mean_Placefield_area_R;
SA.mean_Placefield_area_L=mean_Placefield_area_L;
SA.mean_Placefield_area_B=mean_Placefield_area_B;


JointProb = cellfun(@mean,SA.PlaceFieldPeakJointProbability)
mean_JointProb=mean(JointProb);

Infield=cellfun(@mean,SA.InFieldActivityRatio)
mean_infield=mean(Infield);

stability=cellfun(@mean,SA.PlaceFieldStability)
mean_stability=mean(stability);

Firingprob=cellfun(@mean,SA.CellFiringProb)
mean_Firingprob=mean(Firingprob);

SI=cellfun(@mean,SA.SI)
mean_SI=mean(SI);

Sparsity=cellfun(@mean,SA.sparsity)
mean_sparsity=mean(Sparsity);

PlaceArea=cellfun(@mean,SA.PlaceFieldArea)
mean_PlaceArea=mean(PlaceArea);

SA.mean_JointProb=mean_JointProb;
SA.mean_infield=mean_infield;
SA.mean_stability=mean_stability;
SA.mean_Firingprob=mean_Firingprob;
SA.mean_SI=mean_SI;
SA.mean_sparsity=mean_sparsity;
SA.mean_PlaceArea=mean_PlaceArea;

%% added by Jisoo

% SA.Right_Centroid={PlaceFieldCentroid{:,1}};
% for i=1:length(SA.Right_Centroid)
%     
%     if isnan(SA.Right_Centroid{1,i})==0
%         SA.Right_Centroid{1,i}=cell2mat(SA.Right_Centroid{1,i})
%          SA.Right_Centroid{1,i}= SA.Right_Centroid(1);
% 
%  
%     else SA.Right_Centroid{1,i}=[];
%       
%     
%     
%     end
% end




clearvars -except SA
save('spatial_analysis_classif.mat', 'SA');
%%
% numplace=(cellfun(@length,SA.classif_index));
% numplace=[numplace(1),numplace(2),numplace(3)+numplace(4),sum(numplace(1:4)),numplace(5)];
% ratioplace=100*(numplace./sum(cellfun(@length ,SA.classif_index)));
% 
% fields=fieldnames(SA);
% results=[];
% for f=1:length(fields)-2
%     results(:,f)= [cell2mat(SA.(fields{f})(1));  cell2mat(SA.(fields{f})(2));  cell2mat(SA.(fields{f})(3))];        
% end
% cfp=[cell2mat(SA.CellFiringProb(1));  cell2mat(SA.CellFiringProb(2));  cell2mat(SA.CellFiringProb(3))];
% %%
% fields=fieldnames(SA);
% meanresults=zeros(length(fields)-2,4);
% for f=1:length(fields)-2
%     meanresults(f,1:3)=cellfun(@mean, SA.(fields{f}));
%     meanresults(f,4)= mean([cell2mat(SA.(fields{f})(1));  cell2mat(SA.(fields{f})(2));  cell2mat(SA.(fields{f})(3))]);        
% end
% cfp= cellfun(@mean,SA.CellFiringProb);
% cfp= [cfp(1:3), mean([cell2mat(SA.CellFiringProb(1));  cell2mat(SA.CellFiringProb(2));  cell2mat(SA.CellFiringProb(3))]), cfp(4)];        
% 
% numplace=(cellfun(@length,SA.classif_index));
% numplace=[numplace(1),numplace(2),numplace(3)+numplace(4),sum(numplace(1:4)),numplace(5)];
% ratioplace=100*(numplace./sum(cellfun(@length ,SA.classif_index)));
% 
% 
% 



function out=direction_sum(variable, direction_vector, col)
    variable=variable(direction_vector, col);
    out=cell2mat(variable);
end
function out=direction_sum2(variable, direction_vector, col)
    variable=variable(direction_vector,col);
    out= cellfun(@(c) cell2mat(c), variable, 'UniformOutput', false);
end









