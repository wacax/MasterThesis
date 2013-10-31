%Preprocess Raw2Mat
function [EV TR] = TrajectoriesRaw2Mat(VPCode)

dire = 'D:\Wacax\TU Berlin Thesis\bbciRaw\';
str = sprintf('%s%s', dire, VPCode);
cd(str)

zHouse1Random = importfile('RandomH1.csv', 'RandomH1');
zHouse1Random = zHouse1Random(1:end-3, :);
zHouse2Random = importfile('RandomH2.csv', 'RandomH2');
zHouse2Random = zHouse2Random(1:end-3, :);
zHouse3Random = importfile('RandomH3.csv', 'RandomH3');
zHouse3Random = zHouse3Random(1:end-3, :);
zHouse1Cont = importfile('House1Cont.csv', 'House1Cont');
zHouse1Cont = zHouse1Cont(1:end-3, :);
zHouse2Cont = importfile('House2Cont.csv', 'House2Cont');
zHouse2Cont = zHouse2Cont(1:end-3, :);
zHouse3Cont = importfile('House3Cont.csv', 'House3Cont');
zHouse3Cont = zHouse3Cont(1:end-3, :);
zHouse1Discrete = importfile('House1Discrete.csv', 'House1Discrete');
zHouse1Discrete = zHouse1Discrete(1:end-3, :);
zHouse2Discrete = importfile('House2Discrete.csv', 'House2Discrete');
zHouse2Discrete = zHouse2Discrete(1:end-3, :);
zHouse3Discrete = importfile('House3Discrete.csv', 'House3Discrete');
zHouse3Discrete = zHouse3Discrete(1:end-3, :);
House1Random = zHouse1Random(zHouse1Random(:,3)~=0, :);
House2Random = zHouse2Random(zHouse2Random(:,3)~=0, :);
House3Random = zHouse3Random(zHouse3Random(:,3)~=0, :);
House1Cont = zHouse1Cont(zHouse1Cont(:,3)~=0, :);
House2Cont = zHouse2Cont(zHouse2Cont(:,3)~=0, :);
House3Cont = zHouse3Cont(zHouse3Cont(:,3)~=0, :);
House1Discrete = zHouse1Discrete(zHouse1Discrete(:,3)~=0, :);
House2Discrete = zHouse2Discrete(zHouse2Discrete(:,3)~=0, :);
House3Discrete = zHouse3Discrete(zHouse3Discrete(:,3)~=0, :);

zJoy1_1 = importfile('Joy1.csv', 'Joy1');
Joy1_1 = zJoy1_1(zJoy1_1(:,3)~=0, :);
zJoy1_2 = importfile('Joy4.csv', 'Joy4');
Joy1_2 = zJoy1_2(zJoy1_2(:,3)~=0, :);
zJoy2_1 = importfile('Joy2.csv', 'Joy2');
Joy2_1 = zJoy2_1(zJoy2_1(:,3)~=0, :);
zJoy2_2 = importfile('Joy5.csv', 'Joy5');
Joy2_2 = zJoy2_2(zJoy2_2(:,3)~=0, :);
zJoy3_1 = importfile('Joy3.csv', 'Joy3');
Joy3_1 = zJoy3_1(zJoy3_1(:,3)~=0, :);
zJoy3_2 = importfile('Joy6.csv', 'Joy6');
Joy3_2 = zJoy3_2(zJoy3_2(:,3)~=0, :);


%Create Trimmed Discrete Condition
[Xrep_1 Yrep_1] = FindStatic(House1Discrete);
EmptyPointsHouse1 = intersect(Xrep_1, Yrep_1);
[Xrep_2 Yrep_2] = FindStatic(House2Discrete);
EmptyPointsHouse2 = intersect(Xrep_2, Yrep_2);
[Xrep_3 Yrep_3] = FindStatic(House3Discrete);
EmptyPointsHouse3 = intersect(Xrep_3, Yrep_3);

indicesMovement1 = ones(size(House1Discrete, 1), 1);
indicesMovement1(EmptyPointsHouse1) = 0;
House1DiscreteTrimmed = House1Discrete(indicesMovement1 == 1, :);
indicesMovement2 = ones(size(House2Discrete, 1), 1);
indicesMovement2(EmptyPointsHouse2) = 0;
House2DiscreteTrimmed = House2Discrete(indicesMovement2 == 1, :);
indicesMovement3 = ones(size(House3Discrete, 1), 1);
indicesMovement3(EmptyPointsHouse3) = 0;
House3DiscreteTrimmed = House3Discrete(indicesMovement3 == 1, :);

%Calculate distances
numSamples1 =  floor(mean([size(Joy1_1, 1) size(Joy1_2, 1)]));
numSamples2 =  floor(mean([size(Joy2_1, 1) size(Joy2_2, 1)]));
numSamples3 =  floor(mean([size(Joy3_1, 1) size(Joy3_2, 1)]));

[~, indices1_1] = datasample(Joy1_1, numSamples1, 1);
[~, indices1_2] = datasample(Joy1_2, numSamples1, 1);
[~, indices2_1] = datasample(Joy2_1, numSamples2, 1);
[~, indices2_2] = datasample(Joy2_2, numSamples2, 1);
[~, indices3_1] = datasample(Joy3_1, numSamples3, 1);
[~, indices3_2] = datasample(Joy3_2, numSamples3, 1);

House1Benchmark = (Joy1_1(sort(indices1_1), :) + Joy1_2(sort(indices1_2), :))./2;
House2Benchmark = (Joy2_1(sort(indices2_1), :) + Joy2_2(sort(indices2_2), :))./2;
House3Benchmark = (Joy3_1(sort(indices3_1), :) + Joy3_2(sort(indices3_2), :))./2;

%Events
EVHouse1Cont = zHouse1Cont(zHouse1Cont(:,3)==0, :);
EVHouse2Cont = zHouse2Cont(zHouse2Cont(:,3)==0, :);
EVHouse3Cont = zHouse3Cont(zHouse3Cont(:,3)==0, :);
EVHouse1Discrete = zHouse1Discrete(zHouse1Discrete(:,3)==0, :);
EVHouse2Discrete = zHouse2Discrete(zHouse2Discrete(:,3)==0, :);
EVHouse3Discrete = zHouse3Discrete(zHouse3Discrete(:,3)==0, :);
EVHouse1DiscreteTrimmed = House1DiscreteTrimmed(House1DiscreteTrimmed(:,3)==0, :);
EVHouse2DiscreteTrimmed = House2DiscreteTrimmed(House2DiscreteTrimmed(:,3)==0, :);
EVHouse3DiscreteTrimmed = House3DiscreteTrimmed(House3DiscreteTrimmed(:,3)==0, :);
EVHouse1Benchmark = House1Benchmark(House1Benchmark(:,3)==0, :);
EVHouse2Benchmark = House2Benchmark(House2Benchmark(:,3)==0, :);
EVHouse3Benchmark = House3Benchmark(House3Benchmark(:,3)==0, :);
EVHouse1Random = House1Random(House1Random(:,3)==0, :);
EVHouse2Random = zHouse2Random(zHouse2Random(:,3)==0, :);
EVHouse3Random = zHouse3Random(zHouse3Random(:,3)==0, :);

EV = struct('EVHouse1Cont', EVHouse1Cont,...
    'EVHouse2Cont', EVHouse2Cont, ...
    'EVHouse3Cont', EVHouse3Cont, ...
    'EVHouse1Discrete', EVHouse1Discrete, ...
    'EVHouse2Discrete', EVHouse2Discrete, ...
    'EVHouse3Discrete', EVHouse3Discrete, ...
    'EVHouse1DiscreteTrimmed', EVHouse1DiscreteTrimmed, ...
    'EVHouse2DiscreteTrimmed', EVHouse2DiscreteTrimmed, ...
    'EVHouse3DiscreteTrimmed', EVHouse3DiscreteTrimmed, ...
    'EVHouse1Benchmark', EVHouse1Benchmark, ...
    'EVHouse2Benchmark', EVHouse2Benchmark, ...
    'EVHouse3Benchmark', EVHouse3Benchmark, ...
    'EVHouse1Random', EVHouse1Random, ...
    'EVHouse2Random', EVHouse2Random, ...
    'EVHouse3Random', EVHouse3Random);

TR = struct('House1Cont', House1Cont,...
    'House2Cont', House2Cont, ...
    'House3Cont', House3Cont, ...
    'House1Discrete', House1Discrete, ...
    'House2Discrete', House2Discrete, ...
    'House3Discrete', House3Discrete, ...
    'House1DiscreteTrimmed', House1DiscreteTrimmed, ...
    'House2DiscreteTrimmed', House2DiscreteTrimmed, ...
    'House3DiscreteTrimmed', House3DiscreteTrimmed, ...
    'House1Benchmark', House1Benchmark, ...
    'House2Benchmark', House2Benchmark, ...
    'House3Benchmark', House3Benchmark, ...
    'House1Random', House1Random, ...
    'House2Random', House2Random, ...
    'House3Random', House3Random);
    
      

end