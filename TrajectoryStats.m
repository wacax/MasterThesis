function [indexClosest optimalValues experimentalValues] = TrajectoryStats(VPStructure, optimalTrajectory, resolution, Times2Check)

%Define some pre-parameters
houses = 3;
ConditionsNames = {'Discrete','Cont', 'DiscreteTrimmed','Random'};

indexClosest = zeros(houses, length(ConditionsNames), length(Times2Check));
optimalValues = zeros(length(Times2Check), 8);
experimentalValues = zeros(length(ConditionsNames), 8, houses, length(Times2Check));

for iv = 1:length(Times2Check)
    val = Times2Check(iv)*1000; %value to find
    if val > 35500
        val = 35500;
    end
    f = optimalTrajectory(:,2);
    arf = abs(f-val);
    [~, idx] = min(arf);
    optimalValues(iv,:) = optimalTrajectory(idx, :);
    
    for i = 1:houses
        for ii = 1:length(ConditionsNames)
            path1 = VPStructure.(sprintf('%s%i%s', 'House', i, ConditionsNames{ii}));
            
            %Test with the other function
            newSizes1 = size(path1, 1);
            newSizes2 = size(optimalTrajectory, 1);
            
            [~, Indices1] = datasample(path1, floor(newSizes1/resolution), 1);
            
            thaIndices = sort(Indices1);
            thaDistances = zeros(newSizes2, 1);
            
            %find closest time
            val = Times2Check(iv)*1000; %time to find values on
            path1 = path1(thaIndices, :);
            f = path1(:,2);
            arf = abs(f-val);
            [~, idx] = min(arf);
            
            for iii = 1:newSizes2
                thaDistances(iii) = pdist2([path1(idx,3) path1(idx,4)], [optimalTrajectory(iii,3) optimalTrajectory(iii,4)], 'euclidean');
            end
            
            [~, coolIndex] = min(thaDistances);
            indexClosest(i, ii, iv) = coolIndex;
            experimentalValues(ii,:, i, iv) = optimalTrajectory(coolIndex, :);
            
        end
    end
end


end
