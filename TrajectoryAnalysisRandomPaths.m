function DistancesMatrix = TrajectoryAnalysisRandomPaths(VPStructure, randomStructure, resolution)
 
houses = 3;
ConditionsNames = {'Benchmark', 'Discrete','Cont', 'DiscreteTrimmed'};
comparisons = 1:105;

DistancesMatrix = zeros(houses, length(comparisons), length(ConditionsNames));

for iii = 1:length(ConditionsNames)
    for i = 1:houses
        for ii = 1:length(comparisons)
            try
            path1 = VPStructure.(sprintf('%s%i%s', 'House', i, ConditionsNames{comparisons(iii)}));
            path2 = randomStructure.(sprintf('%s%i', 'Path', ii));
            
            %Test with the other function
            newSizes1 = size(path1, 1)
            newSizes2 = size(path2, 1)
            
            [~, Indices1] = datasample(path1, floor(newSizes1/resolution), 1);
            [~, Indices2] = datasample(path2, floor(newSizes2/resolution), 1);
            
            DistancesMatrix(i, ii, iii) = DiscreteFrechetDist([path1(sort(Indices1),3), ...
                path1(sort(Indices1),4)],...
                [path2(sort(Indices2),3), ...
                path2(sort(Indices2),4)])
            end                          
            
        end
    end
end





end
