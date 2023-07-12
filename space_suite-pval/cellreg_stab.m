%load cellreg

%load SAs %load spatial_analysis by hand
spatial_analysis{1}=load('spatial_analysis.mat');
SA{1}=load('spatial_analysis_classif.mat');

%%
idx=cell_registered_struct.cell_to_index_map;
sesh=[3,4,5];
keep=idx(sum(idx(:,sesh)==0,2)==0, sesh);
pc=zeros(size(keep));
for s=1: length(sesh)
    ix=keep(:,s);
    im=zeros(length(ix),1);
    si=[];        
    for c=1:4
        clasif_conc=cell2mat(SA{s}.SA.classif_index(c));
        im=im+ismember(ix,clasif_conc);
    end
    
    for lin=1:length(ix)            
        si=[si; spatial_analysis{s}.spatial_analysis.bin{ix(lin), 3}.spatial_matrix];      
    end
    pc(:,s)=im;
    sm{s}=si;    
end
c=0;
for n=1:length(pc)
    i=logical(sum(pc(n,:)));
    if i
        c=c+1;
        stab(c,1)=abs(corr2(sm{1,1}(n,:), sm{1,2}(n,:)));
        stab(c,2)=abs(corr2(sm{1,1}(n,:), sm{1,3}(n,:)));
        stab(c,3)=abs(corr2(sm{1,2}(n,:), sm{1,3}(n,:)));
    end
end
stab(isnan(stab))=0;