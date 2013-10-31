function [DistancesMatrix] = TrajectoryAnalysis(VPStructure, resolution)

% 
% House1 Analysis
% Distance Metrics
% f12_1 =frechet(VPStructure.House1Benchmark(:,3), VPStructure.House1Benchmark(:,4), VPStructure.House1Discrete(:,3), VPStructure.House1Discrete(:,4));
% f13_1 =frechet(VPStructure.House1Benchmark(:,3), VPStructure.House1Benchmark(:,4), VPStructure.House1Cont(:,3), VPStructure.House1Cont(:,4));
% f14_1 =frechet(VPStructure.House1Benchmark(:,3), VPStructure.House1Benchmark(:,4), VPStructure.House1DiscreteTrimmed(:,3), VPStructure.House1DiscreteTrimmed(:,4));
% f15_1 =frechet(VPStructure.House1Benchmark(:,3), VPStructure.House1Benchmark(:,4), VPStructure.House1Random(:,3), VPStructure.House1Random(:,4));
% f23_1 = 0;
% f24_1 =frechet(VPStructure.House1Discrete(:,3), VPStructure.House1Discrete(:,4), VPStructure.House1DiscreteTrimmed(:,3), VPStructure.House1DiscreteTrimmed(:,4));
% f34_1 = 0;
% f11_1 = 0;
% f22_1 = 0;
% f33_1 = 0;
% f44_1 = 0;
% f52_1 = frechet(VPStructure.House1Random(:,3), VPStructure.House1Random(:,4), VPStructure.House1Discrete(:,3), VPStructure.House1Discrete(:,4));
% f53_1 = frechet(VPStructure.House1Random(:,3), VPStructure.House1Random(:,4), VPStructure.House1Cont(:,3), VPStructure.House1Cont(:,4));
% f54_1 = frechet(VPStructure.House1Random(:,3), VPStructure.House1Random(:,4), VPStructure.House1DiscreteTrimmed(:,3), VPStructure.House1DiscreteTrimmed(:,4));
% f55_1 = 0;

%Define some pre-parameters
houses = 3
ConditionsNames = {'Benchmark', 'Discrete','Cont', 'DiscreteTrimmed','Random'};
comparisons = [1,2;1,3;1,4;1,5;2,3;5,2;5,3;5,4];

DistancesMatrix = zeros(houses, size(comparisons, 1));

for i = 1:houses
    for ii = 1:size(comparisons, 1)
        path1 = VPStructure.(sprintf('%s%i%s', 'House', i, ConditionsNames{comparisons(ii, 1)}));
        path2 = VPStructure.(sprintf('%s%i%s', 'House', i, ConditionsNames{comparisons(ii, 2)}));
        
        %Test with the other function
        newSizes1 = size(path1, 1)
        newSizes2 = size(path2, 1)
        
        [~, Indices1] = datasample(path1, floor(newSizes1/resolution), 1);
        [~, Indices2] = datasample(path2, floor(newSizes2/resolution), 1);
        
        DistancesMatrix(i, ii) = DiscreteFrechetDist([path1(sort(Indices1),3), ...
            path1(sort(Indices1),4)],...
            [path2(sort(Indices2),3), ...
            path2(sort(Indices2),4)])
    end
end


end

