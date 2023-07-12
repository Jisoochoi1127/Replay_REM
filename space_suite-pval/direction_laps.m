function [dVector]= direction_laps(behav, direction)
% direction_laps returns vector =1 when animal traverses 1D environment in
% a certain direction; Takes interpolated_X vector and direction 'l'
% (left), 'r' right.

if nargin<2
    direction ='bi';
end

X_bin_vector=behav.X_bin_vector;
bin_size=X_bin_vector(2)-X_bin_vector(1);
interpolated_X=behav.interpolated_X;
starts=[];
ends=[];
laps=0;
if direction == 'r'
    idx=find(diff(interpolated_X)>0);
     
    for x =1: length(interpolated_X)-1
        if interpolated_X(x)<= X_bin_vector(1)-bin_size & interpolated_X(x+1) > X_bin_vector(1)-bin_size
         
            %Added by jisoo - just for pv1254_HATD1 session
                  % if interpolated_X(x)<=  min(interpolated_X) & interpolated_X(x+1) > min(interpolated_X)

            
            starts=[starts; x];
            laps=laps+1;
            if length(starts)>1 & (isempty(ends) | starts(end-1)>ends(end))
                starts(end)=[];
            end
        end

        if interpolated_X(x)<= X_bin_vector(end) & interpolated_X(x+1) > X_bin_vector(end) & laps>0 
            ends=[ends; x];
            if length(ends)>1 & starts(end)<ends(end-1)
                ends(end)=[];
            end
            
        end
    end
elseif direction == 'l'
    idx=find(diff(interpolated_X)<0);      
    for x =1: length(interpolated_X)-1
        if interpolated_X(x)> X_bin_vector(end) & interpolated_X(x+1) <= X_bin_vector(end)
            starts=[starts; x];
            laps=laps+1;
            if length(starts)>1 & (isempty(ends) | starts(end-1)>ends(end))
                starts(end)=[];
            end
        end

        if interpolated_X(x)>= X_bin_vector(1)-bin_size & interpolated_X(x+1) < X_bin_vector(1)-bin_size  & laps>0
            ends=[ends; x];
            if length(ends)>1 & starts(end)<ends(end-1)
                ends(end)=[];
            end
        end

    end
else
    idx=1:1:length(interpolated_X); 
    starts=1;
    ends=length(interpolated_X);

end

diVector= zeros(length(interpolated_X),1);
diVector(idx)=1;

if isempty(starts)==1 & isempty(ends)==1
    times=[starts ends];
dVector.diVector=[];
dVector.times=[];
dVector.laps=[];

elseif starts(end)>ends(end)
starts(end)=[];


times=[starts ends];
dVector.diVector=diVector;
dVector.times=times;
dVector.laps=length(starts);
end
end

