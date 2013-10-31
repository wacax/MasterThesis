function finalStr = reconstructTrimmed(TrajectoryData)

    houses = 3
    for i = 1:houses
        
    trimName = sprintf('%s%i%s', 'House', i, 'DiscreteTrimmed');
    newVectorLength = size(TrajectoryData.(trimName), 1);
    discName = sprintf('%s%i%s', 'House', i, 'Discrete');
    newStr.(trimName) = TrajectoryData.(discName)(1:newVectorLength, 1);


    end
    
    finalStr = TrajectoryData;
    
    for i = 1:houses
        
        trimName = sprintf('%s%i%s', 'House', i, 'DiscreteTrimmed');

        finalStr.(trimName)(:,1) = newStr.(trimName);
    
    end

end