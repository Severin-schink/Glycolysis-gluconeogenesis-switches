function [] = MakeProteomeDataIntoCSV

loadProteomicsData;

% Set of proteins which will be saved
SetProt = {'aceA','aceB','fbp','maeB','pckA','pfkA','ppc','ppsA','pykF'};

%% create an array of structs, which contain the information of the selected proteins
for i=1:length(SetProt)
    
    indx = find(strcmp(gene,SetProt{i})); %find indexes where you have measurements of this protein
    SetProt{i}
    % general info
    protStruct(i).name = SetProt(i); %#ok<*AGROW>
    protStruct(i).charge = charge(indx(1)); %since charge is always the same for a protein, just keep it once
    protStruct(i).sequence = seq(indx(1));  %since sequence is always the same for a protein, just keep it once
    

    protStruct(i).weights = svmPred(indx); %weights-confidence, based on an SVM classifier
    protStruct(i).condition = condition(indx);  %condition
    protStruct(i).ratio = ratio(indx);          %ratio

    clear indx
    
    % assemble data based on conditions
    
    indx = find(strcmp(protStruct(i).condition,'ace')); %get the indexes of the specified condition
    protStruct(i).aceRatio = protStruct(i).ratio(indx); %asign ratio (using the indexes) 
    protStruct(i).aceWeights = protStruct(i).weights(indx); %asign weights (using the indexes)
    
    % calculate weighted Median based on weights and ratios
    [Low_ace, Med_ace, High_ace] = median_maker( protStruct(i).aceRatio,protStruct(i).aceWeights );
    
    clear indx

    indx = find(strcmp(protStruct(i).condition,'glu'));
    protStruct(i).gluRatio = protStruct(i).ratio(indx);
    protStruct(i).gluWeights = protStruct(i).weights(indx);
    
    % calculate weighted Median based on weights and ratios
    [Low_glu, Med_glu, High_glu] = median_maker( protStruct(i).gluRatio,protStruct(i).gluWeights );

    clear indx
  
    csvwrite([cell2mat( SetProt(i)), '.csv'],[Med_ace,Low_ace,High_ace,Med_glu,Low_glu,High_glu])
end



end

function [Low, Med, High] = median_maker(ratio,weights)

[sortx,IdxOrder] = sort(ratio);
sortw = weights(IdxOrder);

% Weighted low, med and high points
points = [sum(sortw)/4,sum(sortw)/2, sum(sortw)*3/4];

csumw = cumsum(sortw);

for i = 1:length(points)
    j = find(csumw<=points(i),1,'last');
    dj = csumw(j+1)-csumw(j);
    
    M(i) = sortx(j)*(1-(points(i)-csumw(j))/dj)+sortx(j+1)*(1-(csumw(j+1)-points(i))/dj);
end


Low = M(1);
Med = M(2);
High = M(3);

[Low, Med, High]
end