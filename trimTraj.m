function TrimmedTrajectory  = trimTraj(inputTrajectory)

houses = 3
fullConditionsNoTrimmed = {'Cont', 'Discrete', 'Benchmark', 'Random', 'DiscreteTrimmed'}

for i = 1:houses
    
    for ii = 1:length(fullConditionsNoTrimmed)
        currentCond = sprintf('%s%i%s', 'House', i, fullConditionsNoTrimmed{ii});
        trimCond = sprintf('%s%i%s', 'House', i, 'DiscreteTrimmed');
                
        lel = size(inputTrajectory.(trimCond), 1);
        
        if size(inputTrajectory.(currentCond), 1) >= size(inputTrajectory.(trimCond), 1)
            TrimmedTrajectory.(currentCond) = inputTrajectory.(currentCond)(1:lel,:);
        else
            TrimmedTrajectory.(currentCond) = inputTrajectory.(currentCond)(1:end,:);
        end
        
    end
    
end

end