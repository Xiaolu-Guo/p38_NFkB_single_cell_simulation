function codon_struc = func_metrics(data,duration_thresh,per_frame)
codon_struc = struct();
codon_struc.peakVec = max(data,[],2);
codon_struc.time2halfMax= -ones(size(data,1),1);

for i_cell = 1:size(data,1)
    if codon_struc.peakVec(i_cell)<= 0
        codon_struc.time2halfMax(i_cell) = -1;
    else
        codon_struc.time2halfMax(i_cell) = (find(data(i_cell, :) >= codon_struc.peakVec(i_cell)/2,1,'first')-1)*1/12;
    end
end

codon_struc.totalIntegral = sum(data,2)*1/12; % units in hours

if nargin >1
    codon_struc.duration = sum(data >= duration_thresh,2) /12;
end

if nargin >2
    for i_hour = 1:floor(size(data,2)/per_frame)
        codon_struc.intigral_half_hour(:,i_hour) = sum(data(:,(i_hour-1)*per_frame+1:i_hour*per_frame),2)*1/12;
    end
else
    per_frame = 24; %6, 12, 24, 28
    for i_hour = 1:floor(size(data,2)/per_frame)
        codon_struc.intigral_half_hour(:,i_hour) = sum(data(:,(i_hour-1)*per_frame+1:i_hour*per_frame),2)*1/12;
    end
    
end

end