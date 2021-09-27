function [] = createParameters(modelObj,ResultsObj)
    paramsInitial = modelObj.Parameters;
    for i=1:length(modelObj.Parameters)
        assignin('caller',paramsInitial(i).Name,paramsInitial(i).Value)
    end
    
    for i=1:length(ResultsObj.DataNames)
        assignin('caller',ResultsObj.DataNames{i},ResultsObj.Data(:,i));
    end
end