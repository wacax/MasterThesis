%Trajectory Analysis Script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all

addpath('D:\Wacax\TU Berlin Thesis\bbciMat');

importfile1('VPCodes.txt')
VPs = VPCodes(9:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import the 105 random Paths

for i = 1:105
    
    try
    str = sprintf('%i%s', i, 'Trial.txt')
    name = sprintf('%s%i', 'Path', i)
    importfile3(str)
    randomPaths.(name) = data
    end
    
end

save('randomPaths.mat', 'randomPaths')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Trim random trajectories
for i = 1:105
   try
       name = sprintf('%s%i', 'Path', i)
       data = randomPaths.(name);
       randomTrimmedPaths.(name) = data(1:1500, :);
       
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%calcuate the frechet distances between the experimental condition paths
%and the random paths

RandomFrechetTrajectory1 = TrajectoryAnalysisRandomPaths(Trajectories1, randomPaths, 200);
RandomFrechetTrajectory2 = TrajectoryAnalysisRandomPaths(Trajectories2, randomPaths, 200);
RandomFrechetTrajectory3 = TrajectoryAnalysisRandomPaths(Trajectories3, randomPaths, 200);
RandomFrechetTrajectory4 = TrajectoryAnalysisRandomPaths(Trajectories4, randomPaths, 200);
RandomFrechetTrajectory5 = TrajectoryAnalysisRandomPaths(Trajectories5, randomPaths, 200);
RandomFrechetTrajectory6 = TrajectoryAnalysisRandomPaths(Trajectories6, randomPaths, 200);
RandomFrechetTrajectory7 = TrajectoryAnalysisRandomPaths(Trajectories7, randomPaths, 200);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%better one with same time frames

Trajectories1 = reconstructTrimmed(Trajectories1);
Trajectories2 = reconstructTrimmed(Trajectories2);
Trajectories3 = reconstructTrimmed(Trajectories3);
Trajectories4 = reconstructTrimmed(Trajectories4);
Trajectories5 = reconstructTrimmed(Trajectories5);
Trajectories6 = reconstructTrimmed(Trajectories6);
Trajectories7 = reconstructTrimmed(Trajectories7);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Colapsed Distances describe condition + house
Distances1 = TrajectoryAnalysis(Trajectories1, 200);
Distances2 = TrajectoryAnalysis(Trajectories2, 200);
Distances3 = TrajectoryAnalysis(Trajectories3, 200);
Distances4 = TrajectoryAnalysis(Trajectories4, 200);
Distances5 = TrajectoryAnalysis(Trajectories5, 200);
Distances6 = TrajectoryAnalysis(Trajectories6, 200);
Distances7 = TrajectoryAnalysis(Trajectories7, 200);

ColapsedDistances = [Distances1; Distances2; Distances3; Distances4;
    Distances5; Distances6; Distances7];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
participants = 7;
uncolapsedDistances = zeros(participants, 24);
a = 0;
for i = 1:participants
    a = a+1
    b = a+2
    c = ColapsedDistances(a:b,:)
    x = c'
    uncolapsedDistances(i,:) = x(:)
end

TrajectoryPlots(uncolapsedDistances, Trajectories3, 1, 3)

LabelsBoxplot = {'BenchDisc', 'BenchCont', 'BenchDiscTrim',...
    'BenchRand', 'DiscDiscTrim', 'RandDisc', ...
    'RandCont', 'RandDiscTrim'}

figure
boxplot(ColapsedDistances, 'labels', LabelsBoxplot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% participants = 7;
% houses = 3;
% Traject = {'Trajectories1', 'Trajectories2', 'Trajectories3', 'Trajectories4'...
%     'Trajectories5', 'Trajectories6', 'Trajectories7'};
% Cont = []
% for i = 1:participants
%     
%     for ii = 1:houses
%         
%         theString = sprintf('%s%i%s', 'House', ii, 'Cont');
%         traj = Traject{i};
%         Cont = [Cont; traj.(theString)]
%     
%     end
%     
% end

Cont = [Trajectories1.House1Cont; Trajectories1.House2Cont; Trajectories1.House3Cont;... 
    Trajectories2.House1Cont; Trajectories2.House2Cont; Trajectories2.House3Cont;... 
    Trajectories3.House1Cont; Trajectories3.House2Cont; Trajectories3.House3Cont;...  
    Trajectories4.House1Cont; Trajectories4.House2Cont; Trajectories4.House3Cont;... 
    Trajectories5.House1Cont; Trajectories5.House2Cont; Trajectories5.House3Cont;... 
    Trajectories6.House1Cont; Trajectories6.House2Cont; Trajectories6.House3Cont;... 
    Trajectories7.House1Cont; Trajectories7.House2Cont; Trajectories7.House3Cont];

Disc = [Trajectories1.House1DiscreteTrimmed; Trajectories1.House2DiscreteTrimmed; Trajectories1.House3DiscreteTrimmed;... 
    Trajectories2.House1DiscreteTrimmed; Trajectories2.House2DiscreteTrimmed; Trajectories2.House3DiscreteTrimmed;... 
    Trajectories3.House1DiscreteTrimmed; Trajectories3.House2DiscreteTrimmed; Trajectories3.House3DiscreteTrimmed;...  
    Trajectories4.House1DiscreteTrimmed; Trajectories4.House2DiscreteTrimmed; Trajectories4.House3DiscreteTrimmed;... 
    Trajectories5.House1DiscreteTrimmed; Trajectories5.House2DiscreteTrimmed; Trajectories5.House3DiscreteTrimmed;... 
    Trajectories6.House1DiscreteTrimmed; Trajectories6.House2DiscreteTrimmed; Trajectories6.House3DiscreteTrimmed;... 
    Trajectories7.House1DiscreteTrimmed; Trajectories7.House2DiscreteTrimmed; Trajectories7.House3DiscreteTrimmed];

% Rand = [Trajectories1.House1Random; Trajectories1.House2Random; Trajectories1.House3Random;... 
%     Trajectories2.House1Random; Trajectories2.House2Random; Trajectories2.House3Random;... 
%     Trajectories3.House1Random; Trajectories3.House2Random; Trajectories3.House3Random;...  
%     Trajectories4.House1Random; Trajectories4.House2Random; Trajectories4.House3Random;... 
%     Trajectories5.House1Random; Trajectories5.House2Random; Trajectories5.House3Random;... 
%     Trajectories6.House1Random; Trajectories6.House2Random; Trajectories6.House3Random;... 
%     Trajectories7.House1Random; Trajectories7.House2Random; Trajectories7.House3Random];

Rand = []

for i = 1:4
    try
        newNumber = randi([1, 99])
        %newNumber = [3 4 5]
        str = sprintf('%s%i', 'Path', newNumber);
        Rand = [Rand; randomPaths.(str)];
    end
    
end

Bench = [Trajectories1.House1Benchmark; Trajectories1.House2Benchmark; Trajectories1.House3Benchmark;... 
    Trajectories2.House1Benchmark; Trajectories2.House2Benchmark; Trajectories2.House3Benchmark;... 
    Trajectories3.House1Benchmark; Trajectories3.House2Benchmark; Trajectories3.House3Benchmark;...  
    Trajectories4.House1Benchmark; Trajectories4.House2Benchmark; Trajectories4.House3Benchmark;... 
    Trajectories5.House1Benchmark; Trajectories5.House2Benchmark; Trajectories5.House3Benchmark;... 
    Trajectories6.House1Benchmark; Trajectories6.House2Benchmark; Trajectories6.House3Benchmark;... 
    Trajectories7.House1Benchmark; Trajectories7.House2Benchmark; Trajectories7.House3Benchmark];

%jet colors seemed better fit for trajectories
figure
subplot(2,2,1)
scattercloud(Cont(:,3), Cont(:,4), 20, 180, 'ko', jet(256))
title(['Continuous'])
axis equal tight
subplot(2,2,2)
scattercloud(Disc(:,3), Disc(:,4), 20, 180, 'ko', jet(256))
title(['Discrete'])
axis equal tight
subplot(2,2,3)
scattercloud(Rand(:,3), Rand(:,4), 20, 180, 'ko', jet(256))
title(['Random'])
axis equal tight
subplot(2,2,4)
scattercloud(Bench(:,3), Bench(:,4), 20, 5, 'ko', jet(256))
title(['Joystick'])
axis equal tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Frechet Distance with Trimmed Conditions to equate to trimmed continuous

%Trimming
Trajectories1Trimmed = trimTraj(Trajectories1)
Trajectories2Trimmed = trimTraj(Trajectories2)
Trajectories3Trimmed = trimTraj(Trajectories3)
Trajectories4Trimmed = trimTraj(Trajectories4)
Trajectories5Trimmed = trimTraj(Trajectories5)
Trajectories6Trimmed = trimTraj(Trajectories6)
Trajectories7Trimmed = trimTraj(Trajectories7)

%Frechet 

%Colapsed Distances describe condition + house
%NOTE: When it comes to resolution in the trajectoris analysis script a
%lower number means higher resolution as oposed to the frechet.m script

DistancesTrim1 = TrajectoryAnalysis(Trajectories1Trimmed, 50)
DistancesTrim2 = TrajectoryAnalysis(Trajectories2Trimmed, 50)
DistancesTrim3 = TrajectoryAnalysis(Trajectories3Trimmed, 50)
DistancesTrim4 = TrajectoryAnalysis(Trajectories4Trimmed, 50)
DistancesTrim5 = TrajectoryAnalysis(Trajectories5Trimmed, 50)
DistancesTrim6 = TrajectoryAnalysis(Trajectories6Trimmed, 50)
DistancesTrim7 = TrajectoryAnalysis(Trajectories7Trimmed, 50)

ColapsedTrimmedDistances = [DistancesTrim1; DistancesTrim2; DistancesTrim3; DistancesTrim4;
    DistancesTrim5; DistancesTrim6; DistancesTrim7]

LabelsBoxplot = {'BenchDisc', 'BenchCont', 'BenchDiscTrim',...
    'BenchRand', 'DiscDiscTrim', 'RandDisc', ...
    'RandCont', 'RandDiscTrim'}

figure
boxplot(ColapsedTrimmedDistances, 'labels', LabelsBoxplot)


ContTr = [Trajectories1Trimmed.House1Cont; Trajectories1Trimmed.House2Cont; Trajectories1Trimmed.House3Cont;... 
    Trajectories2Trimmed.House1Cont; Trajectories2Trimmed.House2Cont; Trajectories2Trimmed.House3Cont;... 
    Trajectories3Trimmed.House1Cont; Trajectories3Trimmed.House2Cont; Trajectories3Trimmed.House3Cont;...  
    Trajectories4Trimmed.House1Cont; Trajectories4Trimmed.House2Cont; Trajectories4Trimmed.House3Cont;... 
    Trajectories5Trimmed.House1Cont; Trajectories5Trimmed.House2Cont; Trajectories5Trimmed.House3Cont;... 
    Trajectories6Trimmed.House1Cont; Trajectories6Trimmed.House2Cont; Trajectories6Trimmed.House3Cont;... 
    Trajectories7Trimmed.House1Cont; Trajectories7Trimmed.House2Cont; Trajectories7Trimmed.House3Cont];

DiscTr = [Trajectories1Trimmed.House1DiscreteTrimmed; Trajectories1Trimmed.House2DiscreteTrimmed; Trajectories1Trimmed.House3DiscreteTrimmed;... 
    Trajectories2Trimmed.House1DiscreteTrimmed; Trajectories2Trimmed.House2DiscreteTrimmed; Trajectories2Trimmed.House3DiscreteTrimmed;... 
    Trajectories3Trimmed.House1DiscreteTrimmed; Trajectories3Trimmed.House2DiscreteTrimmed; Trajectories3Trimmed.House3DiscreteTrimmed;...  
    Trajectories4Trimmed.House1DiscreteTrimmed; Trajectories4Trimmed.House2DiscreteTrimmed; Trajectories4Trimmed.House3DiscreteTrimmed;... 
    Trajectories5Trimmed.House1DiscreteTrimmed; Trajectories5Trimmed.House2DiscreteTrimmed; Trajectories5Trimmed.House3DiscreteTrimmed;... 
    Trajectories6Trimmed.House1DiscreteTrimmed; Trajectories6Trimmed.House2DiscreteTrimmed; Trajectories6Trimmed.House3DiscreteTrimmed;... 
    Trajectories7Trimmed.House1DiscreteTrimmed; Trajectories7Trimmed.House2DiscreteTrimmed; Trajectories7Trimmed.House3DiscreteTrimmed];

% RandTr = [Trajectories1Trimmed.House1Random; Trajectories1Trimmed.House2Random; Trajectories1Trimmed.House3Random;... 
%     Trajectories2Trimmed.House1Random; Trajectories2Trimmed.House2Random; Trajectories2Trimmed.House3Random;... 
%     Trajectories3Trimmed.House1Random; Trajectories3Trimmed.House2Random; Trajectories3Trimmed.House3Random;...  
%     Trajectories4Trimmed.House1Random; Trajectories4Trimmed.House2Random; Trajectories4Trimmed.House3Random;... 
%     Trajectories5Trimmed.House1Random; Trajectories5Trimmed.House2Random; Trajectories5Trimmed.House3Random;... 
%     Trajectories6Trimmed.House1Random; Trajectories6Trimmed.House2Random; Trajectories6Trimmed.House3Random;... 
%     Trajectories7Trimmed.House1Random; Trajectories7Trimmed.House2Random; Trajectories7Trimmed.House3Random];

RandTr = []

for i = 1:8
    try
        newNumber = randi([1, 99])
        %newNumber = [3 4 5]
        str = sprintf('%s%i', 'Path', newNumber);
        RandTr = [RandTr; randomTrimmedPaths.(str)];
    end
    
end

BenchTr = [Trajectories1Trimmed.House1Benchmark; Trajectories1Trimmed.House2Benchmark; Trajectories1Trimmed.House3Benchmark;... 
    Trajectories2Trimmed.House1Benchmark; Trajectories2Trimmed.House2Benchmark; Trajectories2Trimmed.House3Benchmark;... 
    Trajectories3Trimmed.House1Benchmark; Trajectories3Trimmed.House2Benchmark; Trajectories3Trimmed.House3Benchmark;...  
    Trajectories4Trimmed.House1Benchmark; Trajectories4Trimmed.House2Benchmark; Trajectories4Trimmed.House3Benchmark;... 
    Trajectories5Trimmed.House1Benchmark; Trajectories5Trimmed.House2Benchmark; Trajectories5Trimmed.House3Benchmark;... 
    Trajectories6Trimmed.House1Benchmark; Trajectories6Trimmed.House2Benchmark; Trajectories6Trimmed.House3Benchmark;... 
    Trajectories7Trimmed.House1Benchmark; Trajectories7Trimmed.House2Benchmark; Trajectories7Trimmed.House3Benchmark];

%jet colors seemed better fit for trajectories
figure
subplot(2,2,1)
scattercloud(ContTr(:,3), ContTr(:,4), 20, 180, 'ko', jet(256))
title(['Continuous Trimmed'])
axis equal tight
subplot(2,2,2)
scattercloud(DiscTr(:,3), DiscTr(:,4), 20, 180, 'ko', jet(256))
title(['Discrete Trimmed'])
axis equal tight
subplot(2,2,3)
scattercloud(RandTr(:,3), RandTr(:,4), 20, 180, 'ko', jet(256))
title(['Random Trimmed'])
axis equal tight
subplot(2,2,4)
scattercloud(BenchTr(:,3), BenchTr(:,4), 20, 5, 'ko', jet(256))
title(['Joystick Trimmed'])
axis equal tight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%+
%Random Boxplot Data (4 last bars)

RandomComparisonsUncollapsed = [RandomFrechetTrajectory1 RandomFrechetTrajectory2 RandomFrechetTrajectory3 RandomFrechetTrajectory4 RandomFrechetTrajectory5 RandomFrechetTrajectory6 RandomFrechetTrajectory7];

BoxplotRandomData = zeros(size(RandomComparisonsUncollapsed, 1) * size(RandomComparisonsUncollapsed, 2),...
    4);
for i = 1:4
    dummy = RandomComparisonsUncollapsed(:,:,i);
    BoxplotRandomData(:,i) = dummy(:);
end

ConditionsNames = {'Benchmark - Random', 'Discrete - Random','Cont - Random', 'DiscreteTrimmed -Random'};

boxplot(BoxplotRandomData, 'labels', ConditionsNames, 'labelorientation', 'inline')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fused Boxplot
%IMPORTANT ONE

group = [repmat({'BenchDisc'}, size(ColapsedDistances, 1), 1);
    repmat({'BenchCont'}, size(ColapsedDistances, 1), 1);
    repmat({'BenchDiscTrim'}, size(ColapsedDistances, 1), 1);
    repmat({'BenchRand'}, size(BoxplotRandomData, 1), 1)];
    %repmat({'RandDisc'}, size(BoxplotRandomData, 1), 1);
    %repmat({'RandCont'}, size(BoxplotRandomData, 1), 1);
    %repmat({'RandDiscTrim'}, size(BoxplotRandomData, 1), 1)];

figure
boxplot([(ColapsedDistances(:, 1) -2.6);
    ColapsedDistances(:, 2);
    (ColapsedDistances(:, 3)-5);
    BoxplotRandomData(:, 1)], group, 'labelorientation', 'inline')
    ylabel(['Frechet Distances']) 
    title(['Frechet Distances - Full Trajectories'])
    %BoxplotRandomData(:, 2);
    %BoxplotRandomData(:, 3);
    %BoxplotRandomData(:, 4)], group, 'labelorientation', 'inline')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

RandomTrimmedFrechetTrajectory1 = TrajectoryAnalysisRandomPaths(Trajectories1Trimmed, randomPaths, 200)
RandomTrimmedFrechetTrajectory2 = TrajectoryAnalysisRandomPaths(Trajectories2Trimmed, randomPaths, 200)
RandomTrimmedFrechetTrajectory3 = TrajectoryAnalysisRandomPaths(Trajectories3Trimmed, randomPaths, 200)
RandomTrimmedFrechetTrajectory4 = TrajectoryAnalysisRandomPaths(Trajectories4Trimmed, randomPaths, 200)
RandomTrimmedFrechetTrajectory5 = TrajectoryAnalysisRandomPaths(Trajectories5Trimmed, randomPaths, 200)
RandomTrimmedFrechetTrajectory6 = TrajectoryAnalysisRandomPaths(Trajectories6Trimmed, randomPaths, 200)
RandomTrimmedFrechetTrajectory7 = TrajectoryAnalysisRandomPaths(Trajectories7Trimmed, randomPaths, 200)

%Random Boxplot Data (4 last bars)

RandomTrimmedComparisonsUncollapsed = [RandomTrimmedFrechetTrajectory1 RandomTrimmedFrechetTrajectory2 RandomTrimmedFrechetTrajectory3 RandomTrimmedFrechetTrajectory4 RandomTrimmedFrechetTrajectory5 RandomTrimmedFrechetTrajectory6 RandomTrimmedFrechetTrajectory7];

BoxplotRandomTrimmedData = zeros(size(RandomTrimmedComparisonsUncollapsed, 1) * size(RandomTrimmedComparisonsUncollapsed, 2),...
    4);
for i = 1:4
    dummy = RandomTrimmedComparisonsUncollapsed(:,:,i);
    BoxplotRandomTrimmedData(:,i) = dummy(:);
end

ConditionsNames = {'Benchmark - Random', 'Discrete - Random','Cont - Random', 'DiscreteTrimmed -Random'};

boxplot(BoxplotRandomTrimmedData, 'labels', ConditionsNames, 'labelorientation', 'inline')

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fused Boxplot TRIMMED
%SECOND IMPORTANT ONE

group = [repmat({'BenchDisc'}, size(ColapsedTrimmedDistances, 1), 1);
    repmat({'BenchCont'}, size(ColapsedTrimmedDistances, 1), 1);
    repmat({'BenchDiscTrim'}, size(ColapsedTrimmedDistances, 1), 1);
    repmat({'BenchRand'}, size(BoxplotRandomTrimmedData, 1), 1)];
    %repmat({'RandDisc'}, size(BoxplotRandomTrimmedData, 1), 1);
    %repmat({'RandCont'}, size(BoxplotRandomTrimmedData, 1), 1);
    %repmat({'RandDiscTrim'}, size(BoxplotRandomTrimmedData, 1), 1)];

figure
boxplot([(ColapsedTrimmedDistances(:, 1) - 6.54);
    (ColapsedTrimmedDistances(:, 2) - 4.2);
    ColapsedTrimmedDistances(:, 3);
    BoxplotRandomTrimmedData(:, 1)], group, 'labelorientation', 'inline')
    ylabel(['Frechet Distances']) 
    title(['Frechet Distances - Trimmed Trajectories'])
    %BoxplotRandomTrimmedData(:, 2);
    %BoxplotRandomTrimmedData(:, 3);
    %BoxplotRandomTrimmedData(:, 4)], group, 'labelorientation', 'inline')
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Differences on Speed and Trajectory

preStats1 = mazeStatistics2mat(VPCodes(1))
preStats2 = mazeStatistics2mat(VPCodes(2))
preStats3 = mazeStatistics2mat(VPCodes(3))
preStats4 = mazeStatistics2mat(VPCodes(4))
preStats5 = mazeStatistics2mat(VPCodes(5))
preStats6 = mazeStatistics2mat(VPCodes(6))
preStats7 = mazeStatistics2mat(VPCodes(7))

preStats = {preStats1 preStats2 preStats3 preStats4 preStats5...
    preStats6 preStats7}

participants = 7
comparisons = 6
differences = zeros(7, 4, comparisons, participants)
conditions = {'Cont', 'Discrete'}
houses = 3

for i = 1:participants
    
    for ii = 1:houses
        
        joy = sprintf('%s%i%s', 'House', ii, 'Joystick')
        
        for iii = 1:length(conditions)
            cond = sprintf('%s%i%s', 'House', ii, conditions{iii})
                       
            differences(i).(cond) = abs(preStats{i}.(cond) - preStats{i}.(joy));
        end
        
    end
         
end

%Boxplots of path length and time

participants = 7;
houses = 3;
regions = 7;
comparisons = 6;
differencesNames = {'Path Lenght', 'Time Elapsed', 'Velocity', 'Region'};
conditions = {'Cont', 'Discrete'};
values2plot = 4

%boxplotData = zeros(comparisons * participants, length(differencesNames));
%this is the one that is still obscure!
for i = 1:regions
    
    for iv = 1:values2plot
        
        boxplotData = []

        
        for ii = 1:participants
            
            for iii = 1:houses
                
                intString1 = sprintf('%s%i%s', 'House', iii, conditions{1})
                intString2 = sprintf('%s%i%s', 'House', iii, conditions{2})
                
                boxplotData = [boxplotData; [differences(ii).(intString1)(i, iv)...
                    differences(ii).(intString2)(i, iv)]]
            end  
                
            
        end
        
        figure
        boxplot(boxplotData, 'labels', conditions) 
        
    end
    
        
       
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ERP Amplitude and Latency Extraction
participants = 7;

amplitudesandLatenciesMatrix = zeros(6, participants);

for i = 1:participants
    [amplitudesandLatenciesMatrix(1, i) amplitudesandLatenciesMatrix(2, i) amplitudesandLatenciesMatrix(3, i) amplitudesandLatenciesMatrix(4, i) amplitudesandLatenciesMatrix(5, i) amplitudesandLatenciesMatrix(6, i)] = EEGLoad(VPCodes(i), VPs(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Final Correlations Table

%First Table
%Path length vs amplitude
%7th region represents the overall trajectory
data1 = extractRegionInfo(preStats, 7, 'Cont', 1)
data2 = extractRegionInfo(preStats, 7, 'Discrete', 1)
data3 = extractRegionInfo(preStats, 7, 'Joystick', 1)

data = [data1 data2 data3]

pValues = zeros(1, 3);
conditions = 3;
houses = 3;
amplitudes = amplitudesandLatenciesMatrix(1, :)%here must be a vector with 7 values one for each participant
amplitudesVector = repmat(amplitudes, houses, 1);
amplitudesVector = amplitudesVector(:);
Labelsscatterplot = {'Continuous', 'Discrete', 'Joystick'}
figure
for i = 1:conditions
    [R,P,RLO,RUP] = corrcoef(data(:,i), amplitudesVector)
    pValues(i) = P
    subplot(1, 3, i)
    hold on
    scatter(data(:,i), amplitudesVector)
    lsline
    title(Labelsscatterplot{i})
    xlabel('Path Length')
    ylabel('ERP Amplitude')
end

%Second Table
%Path length vs latency
data1 = extractRegionInfo(preStats, 7, 'Cont', 1);
data2 = extractRegionInfo(preStats, 7, 'Discrete', 1);
data3 = extractRegionInfo(preStats, 7, 'Joystick', 1);

data = [data1 data2 data3];

pValues = zeros(1, 3);
conditions = 3;
houses = 3;
latencies = amplitudesandLatenciesMatrix(2, :)%here must be a vector with 7 values one for each participant
latenciesVector = repmat(latencies, houses, 1);
latenciesVector = latenciesVector(:);
Labelsscatterplot = {'Continuous', 'Discrete', 'Joystick'}
figure
for i = 1:conditions
    [R,P,RLO,RUP] = corrcoef(data(:,i), latenciesVector)
    pValues(i) = P
    subplot(1, 3, i)
    hold on
    scatter(data(:,i), latenciesVector)
    lsline
    title(Labelsscatterplot{i})
    xlabel('Path Length')
    ylabel('ERP Amplitude')
end

%Third Table
%Path time elapsed vs amplitude
pValues = zeros(1, 3);
data1 = extractRegionInfo(preStats, 7, 'Cont', 2);
data2 = extractRegionInfo(preStats, 7, 'Discrete', 2);
data3 = extractRegionInfo(preStats, 7, 'Joystick', 2);

data = [data1 data2 data3];

conditions = 3;
houses = 3;
amplitudes = amplitudesandLatenciesMatrix(1, :) %here must be a vector with 7 values one for each participant
amplitudesVector = repmat(amplitudes, houses, 1);
amplitudesVector = amplitudesVector(:);
Labelsscatterplot = {'Continuous', 'Discrete', 'Joystick'}
figure
for i = 1:conditions
    [R,P,RLO,RUP] = corrcoef(data(:,i), amplitudesVector)
    pValues(i) = P
    subplot(1, 3, i)
    hold on
    scatter(data(:,i), amplitudesVector)
    lsline
    title(Labelsscatterplot{i})
    xlabel('Time Elapsed')
    ylabel('ERP Amplitude')
end

%Fourth Table
%Path time elapsed vs latency
pValues = zeros(1, 3);
data1 = extractRegionInfo(preStats, 7, 'Cont', 2);
data2 = extractRegionInfo(preStats, 7, 'Discrete', 2);
data3 = extractRegionInfo(preStats, 7, 'Joystick', 2);

data = [data1 data2 data3];

conditions = 3;
houses = 3;
latencies = amplitudesandLatenciesMatrix(2, :)%here must be a vector with 7 values one for each participant
latenciesVector = repmat(latencies, houses, 1);
latenciesVector = latenciesVector(:);
Labelsscatterplot = {'Continuous', 'Discrete', 'Joystick'}
figure
for i = 1:conditions
    [R,P,RLO,RUP] = corrcoef(data(:,i), latenciesVector)
     pValues(i) = P
     subplot(1, 3, i)
    hold on
    scatter(data(:,i), latenciesVector)
    lsline
    title(Labelsscatterplot{i})
    xlabel('Time Elapsed')
    ylabel('ERP Amplitude')
end

%Fifth Table
%Frechet vs amplitude

%PD IT CAN ALSO BE DONE WITH THE gplotmatrix(ColapsedDistances,
%amplitudesVector) COMMAND ALTHOUGH THERE IS NO LINE FIT

pValues = zeros(1, participants);
conditions = 7;
houses = 3;
amplitudes = amplitudesandLatenciesMatrix(1, :)%here must be a vector with 7 values one for each participant
amplitudesVector = repmat(amplitudes, houses, 1);
amplitudesVector = amplitudesVector(:);
Labelsscatterplot = {'Benchmark Vs. Discrete', 'Benchmark Vs. Continuous', 'Benchmark Vs. Discrete Trimmed',...
    'Benchmark Vs. Random', 'Discrete Vs. Discrete Trimmed', 'Random Vs. Discrete', ...
    'Random Vs. Continuous', 'Random Vs. Discrete Trimmed'}
figure
for i = 1:conditions
    [R,P,RLO,RUP] = corrcoef(ColapsedDistances(:,i), amplitudesVector)
    pValues(i) = P
    subplot(2, 4, i)
    hold on
    scatter(ColapsedDistances(:,i), amplitudesVector)
    lsline
    title(Labelsscatterplot{i})
    xlabel('Frechet Distance')
    ylabel('ERP Amplitude')
end

%Sixth Table
%Frechet vs latency
pValues = zeros(1, participants);
participants = 7;
houses = 3;
latencies = amplitudesandLatenciesMatrix(2, :)%here must be a vector with 7 values one for each participant
latenciesVector = repmat(latencies, houses, 1);
latenciesVector = latenciesVector(:);
Labelsscatterplot = {'Benchmark Vs. Discrete', 'Benchmark Vs. Continuous', 'Benchmark Vs. Discrete Trimmed',...
    'Benchmark Vs. Random', 'Discrete Vs. Discrete Trimmed', 'Random Vs. Discrete', ...
    'Random Vs. Continuous', 'Random Vs. Discrete Trimmed'}
figure
for i = 1:conditions
    [R,P,RLO,RUP] = corrcoef(ColapsedDistances(:,i), latenciesVector)
    pValues(i) = P
    subplot(2, 4, i)
    hold on
    scatter(ColapsedDistances(:,i), latenciesVector)
    lsline
    title(Labelsscatterplot{i})
    xlabel('Frechet Distance')
    ylabel('ERP Latency')
end

%Seventh Table
%FrechetTrimmed vs amplitude
pValues = zeros(1, participants);
conditions = 7;
houses = 3;
amplitudes = amplitudesandLatenciesMatrix(1, :)%here must be a vector with 7 values one for each participant
amplitudesVector = repmat(amplitudes, houses, 1);
amplitudesVector = amplitudesVector(:);
Labelsscatterplot = {'Benchmark Vs. Discrete', 'Benchmark Vs. Continuous', 'Benchmark Vs. Discrete Trimmed',...
    'Benchmark Vs. Random', 'Discrete Vs. Discrete Trimmed', 'Random Vs. Discrete', ...
    'Random Vs. Continuous', 'Random Vs. Discrete Trimmed'}
figure
for i = 1:conditions
    [R,P,RLO,RUP] = corrcoef(ColapsedTrimmedDistances(:,i), amplitudesVector)
    pValues(i) = P
    subplot(2, 4, i)
    hold on
    scatter(ColapsedTrimmedDistances(:,i), amplitudesVector)
    lsline
    title(Labelsscatterplot{i})
    xlabel('Frechet Distance')
    ylabel('ERP Amplitude')
end


%Eigth Table
%FrechetTrimmed vs latency
pValues = zeros(1, participants);
conditions = 7;
houses = 3;
latencies = amplitudesandLatenciesMatrix(2, :)%here must be a vector with 7 values one for each participant
latenciesVector = repmat(latencies, houses, 1);
latenciesVector = latenciesVector(:);
Labelsscatterplot = {'Benchmark Vs. Discrete', 'Benchmark Vs. Continuous', 'Benchmark Vs. Discrete Trimmed',...
    'Benchmark Vs. Random', 'Discrete Vs. Discrete Trimmed', 'Random Vs. Discrete', ...
    'Random Vs. Continuous', 'Random Vs. Discrete Trimmed'}
figure
for i = 1:conditions
    [R,P,RLO,RUP] = corrcoef(ColapsedTrimmedDistances(:,i), latenciesVector)
    pValues(i) = P
    subplot(2, 4, i)
    hold on
    scatter(ColapsedTrimmedDistances(:,i), latenciesVector)
    lsline
    title(Labelsscatterplot{i})
    xlabel('Frechet Distance')
    ylabel('ERP Latency')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EEG ANALYSIS
Participant = [VPCodes(1:5); VPCodes(7:8)];
VP = [VPs(1:5); VPs(7:8)];
eeg = EEGLoad(Participant, VP);


erp_av = {eeg.erp_av1 eeg.erp_av2 eeg.erp_av3 eeg.erp_av4 eeg.erp_av5 eeg.erp_av6 eeg.erp_av7};
erp_r_av = {eeg.erp_r_av1 eeg.erp_r_av2 eeg.erp_r_av3 eeg.erp_r_av4...
    eeg.erp_r_av5 eeg.erp_r_av6 eeg.erp_r_av7};

erp_av = proc_grandAverage(erp_av);
erp_r_av = proc_grandAverage(erp_r_av);

stdERPplots(erp_av, erp_r_av)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Summary Statistics and that's it
%Summary from length Time and Velocity
dataCont = extractRegionInfo(preStats, 1, 'Cont', 1:3);
dataDisc = extractRegionInfo(preStats, 1, 'Discrete', 1:3);
dataJoy = extractRegionInfo(preStats, 1, 'Joystick', 1:3);

%Remove outliers
cleanDataCont = []
for i = 1:3
    if findOutliers(dataCont(:,i)) ~= 0    
        cleanDataCont = [cleanDataCont dataCont(dataCont(:,i) ~= dataCont(findOutliers(dataCont(:,i)), i))]
    else
        cleanDataCont = [cleanDataCont dataCont(:,i)];
    end
    
end

cleanDataDisc = []
for i = 1:3
    if findOutliers(dataDisc(:,i)) ~= 0 
        cleanDataDisc = [cleanDataDisc dataDisc(dataDisc(:,i) ~= dataDisc(findOutliers(dataDisc(:,i)), i))]
    else
        cleanDataDisc = [cleanDataDisc dataDisc(:,i)];
    end
        
end

cleanDataJoy = []
for i = 1:3
    if findOutliers(dataJoy(:,i)) ~= 0
        cleanDataJoy = [cleanDataJoy dataJoy(dataJoy(:,i) ~= dataJoy(findOutliers(dataJoy(:,i)), i))]
    else
        cleanDataJoy = [cleanDataJoy dataJoy(:,i)];
    end
        
    
end

[DataContMeans DataContMins DataContMaxs DataContSkewds DataContKurt] = grpstats(cleanDataCont,[],{'mean','min','max', 'skewness','kurtosis'});
[DataDiscMeans DataDiscMins DataDiscMaxs DataDiscSkewds DataDiscKurt] = grpstats(cleanDataDisc,[],{'mean','min','max', 'skewness','kurtosis'});
[DataJoyMeans DataJoyMins DataJoyMaxs DataJoySkewds DataJoyKurt] = grpstats(cleanDataJoy,[],{'mean','min','max','skewness','kurtosis'});

group = [repmat({'Continuous'}, size(cleanDataCont, 1), 1);
    repmat({'Discrete'}, size(cleanDataDisc, 1), 1);
    repmat({'Joystick'}, size(cleanDataJoy, 1), 1)];


figure
boxplot([cleanDataCont(:,1); cleanDataDisc(:,1); cleanDataJoy(:,1)], group)
title(['Path Length'])
ylabel('Distance in Maze Units')
figure
hist([cleanDataCont(:,1); cleanDataDisc(:,1); cleanDataJoy(:,1)])
figure
boxplot([cleanDataCont(:,2); cleanDataDisc(:,2); cleanDataJoy(:,2)], group)
title(['Time Elapsed in Seconds'])
ylabel('Seconds')
figure
hist([cleanDataCont(:,2); cleanDataDisc(:,2); cleanDataJoy(:,2)])
figure
boxplot([cleanDataCont(:,3); cleanDataDisc(:,3); cleanDataJoy(:,3)], group)
title(['Velocity'])
ylabel('Maze Distance / Sec')
figure
hist([cleanDataCont(:,3); cleanDataDisc(:,3); cleanDataJoy(:,3)])

%Summary from Frechet Distances
%Full Trajectories
group = [repmat({'BenchDisc'}, size(ColapsedTrimmedDistances, 1), 1);
    repmat({'BenchCont'}, size(ColapsedTrimmedDistances, 1), 1);
    repmat({'BenchDiscTrim'}, size(ColapsedTrimmedDistances, 1), 1);
    repmat({'BenchRand'}, size(BoxplotRandomTrimmedData, 1), 1);
    repmat({'RandDisc'}, size(BoxplotRandomTrimmedData, 1), 1);
    repmat({'RandCont'}, size(BoxplotRandomTrimmedData, 1), 1);
    repmat({'RandDiscTrim'}, size(BoxplotRandomTrimmedData, 1), 1)];

frechData = [(ColapsedTrimmedDistances(:, 1) -2.6);
    ColapsedTrimmedDistances(:, 2);
    (ColapsedTrimmedDistances(:, 3) - 5);
    BoxplotRandomTrimmedData(:, 1);
    BoxplotRandomTrimmedData(:, 2);
    BoxplotRandomTrimmedData(:, 3);
    BoxplotRandomTrimmedData(:, 4)]

[FrechMeans FrechMins FrechMaxs FrechSkewds FrechKurt] = grpstats(frechData,group,{'mean','min','max', 'skewness','kurtosis'});

%Trimmed Trajectories
group = [repmat({'BenchDisc'}, size(ColapsedTrimmedDistances, 1), 1);
    repmat({'BenchCont'}, size(ColapsedTrimmedDistances, 1), 1);
    repmat({'BenchDiscTrim'}, size(ColapsedTrimmedDistances, 1), 1);
    repmat({'BenchRand'}, size(BoxplotRandomTrimmedData, 1), 1);
    repmat({'RandDisc'}, size(BoxplotRandomTrimmedData, 1), 1);
    repmat({'RandCont'}, size(BoxplotRandomTrimmedData, 1), 1);
    repmat({'RandDiscTrim'}, size(BoxplotRandomTrimmedData, 1), 1)];

TrimFrechData = [(ColapsedTrimmedDistances(:, 1) -6.54);
    (ColapsedTrimmedDistances(:, 2) - 4.2);
    ColapsedTrimmedDistances(:, 3);
    BoxplotRandomTrimmedData(:, 1);
    BoxplotRandomTrimmedData(:, 2);
    BoxplotRandomTrimmedData(:, 3);
    BoxplotRandomTrimmedData(:, 4)];

[FrechTrimMeans FrechTrimMins FrechTrimMaxs FrechTrimSkewds FrechTrimKurt] = grpstats(TrimFrechData, group,{'mean','min','max', 'skewness','kurtosis'});

%Extract Data From Region 1
dataContReg1 = extractRegionInfo(preStats, 1, 'Cont', 1:3);
dataDiscReg1 = extractRegionInfo(preStats, 1, 'Discrete', 1:3);
dataJoyReg1 = extractRegionInfo(preStats, 1, 'Joystick', 1:3);

%Remove outliers
cleanDataCont = []
for i = 1:3
    if findOutliers(dataContReg1(:,i)) ~= 0 
        cleanDataCont = [cleanDataCont dataContReg1(dataContReg1(:,i) ~= dataContReg1(findOutliers(dataContReg1(:,i)), i))]
    else
        cleanDataCont = [cleanDataCont dataContReg1(:,i)];
    end
end

cleanDataDisc = []
for i = 1:3
    if findOutliers(dataDiscReg1(:,i)) ~= 0 
        cleanDataDisc = [cleanDataDisc dataDiscReg1(dataDiscReg1(:,i) ~= dataDiscReg1(findOutliers(dataDiscReg1(:,i)), i))]
    else
        cleanDataDisc = [cleanDataDisc dataDiscReg1(:,i)];
    end
end

cleanDataJoy = []
for i = 1:3
    if findOutliers(dataJoyReg1(:,i)) ~= 0
        cleanDataJoy = [cleanDataJoy dataJoyReg1(dataJoyReg1(:,i) ~= dataJoyReg1(findOutliers(dataJoyReg1(:,i)), i))]
    else
        cleanDataJoy = [cleanDataJoy dataJoyReg1(:,i)];
    end
end

[DataContMeans DataContMins DataContMaxs DataContSkewds DataContKurt] = grpstats(cleanDataCont,[],{'mean','min','max', 'skewness','kurtosis'});
[DataDiscMeans DataDiscMins DataDiscMaxs DataDiscSkewds DataDiscKurt] = grpstats(cleanDataDisc,[],{'mean','min','max', 'skewness','kurtosis'});
[DataJoyMeans DataJoyMins DataJoyMaxs DataJoySkewds DataJoyKurt] = grpstats(cleanDataJoy,[],{'mean','min','max','skewness','kurtosis'});

group = [repmat({'Continuous'}, size(cleanDataCont, 1), 1);
    repmat({'Discrete'}, size(cleanDataDisc, 1), 1);
    repmat({'Joystick'}, size(cleanDataJoy, 1), 1)];

group2 = [repmat({'Continuous'}, size(cleanDataCont, 1), 1);
    repmat({'Discrete'}, size(cleanDataDisc, 1), 1)];

figure
boxplot([cleanDataCont(:,1); cleanDataDisc(:,1); cleanDataJoy(:,1)], group)
title(['Path Length'])
ylabel('Distance in Maze Units')
figure
boxplot([cleanDataCont(:,1); cleanDataDisc(:,1)], group2)
title(['Path Length'])
ylabel('Distance in Maze Units')
figure
hist([cleanDataCont(:,1); cleanDataDisc(:,1); cleanDataJoy(:,1)])
figure
boxplot([cleanDataCont(:,2); cleanDataDisc(:,2); cleanDataJoy(:,2)], group)
title(['Time Elapsed in Seconds'])
ylabel('Seconds')
figure
hist([cleanDataCont(:,2); cleanDataDisc(:,2); cleanDataJoy(:,2)])
figure
boxplot([cleanDataCont(:,3); cleanDataDisc(:,3); cleanDataJoy(:,3)], group)
title(['Velocity'])
ylabel('Maze Distance / Sec')
figure
hist([cleanDataCont(:,3); cleanDataDisc(:,3); cleanDataJoy(:,3)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Benchmark Average

Traject = {Trajectories1, Trajectories2, Trajectories3, Trajectories4, Trajectories5, Trajectories6, Trajectories7};

AveragePath = [];
for i = 1:length(Traject)
	for ii = 1:3
		string = sprintf('%s%i%s', 'House', ii, 'Benchmark');
		AveragePath{i, ii} = Traject{i}.(string);

	end
end

AveragePath = reshape(AveragePath, 21, 1);

AverageSize = 1519;
NewAveragePath = zeros(1519, 8, 21);

for i = 1:21
    CurrentAveragePath = AveragePath{i};
    [~, Indices] = datasample(CurrentAveragePath, AverageSize, 1);
    NewAveragePath(:,:,i) = CurrentAveragePath(sort(Indices), :);
end
    
NewAveragePath = mean(NewAveragePath, 3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%That other Data
[indexClosest1 optimalValues1 experimentalValues1] = TrajectoryStats(Trajectories1, NewAveragePath, 1, [10 30 37 120 400])
[indexClosest2 optimalValues2 experimentalValues2] = TrajectoryStats(Trajectories2, NewAveragePath, 1, [10 30 37 120 400])
[indexClosest3 optimalValues3 experimentalValues3] = TrajectoryStats(Trajectories3, NewAveragePath, 1, [10 30 37 120 400])
[indexClosest4 optimalValues4 experimentalValues4] = TrajectoryStats(Trajectories4, NewAveragePath, 1, [10 30 37 120 400])
[indexClosest5 optimalValues5 experimentalValues5] = TrajectoryStats(Trajectories5, NewAveragePath, 1, [10 30 37 120 400])
[indexClosest6 optimalValues6 experimentalValues6] = TrajectoryStats(Trajectories6, NewAveragePath, 1, [10 30 37 120 400])
[indexClosest7 optimalValues7 experimentalValues7] = TrajectoryStats(Trajectories7, NewAveragePath, 1, [10 30 37 120 400])

%%%%%%%%%%
%Extra Trajectory/Time
    
houses = 3
for i = 1:houses
    
    Trim = sprintf('%s%i%s', 'House', i, 'DiscreteTrimmed');
    Disc = sprintf('%s%i%s', 'House', i, 'Discrete');
    
    lel = size(Trajectories7.(Trim), 1);

    Trajectories7.(Trim)(:,2) = Trajectories7.(Disc)(1:lel,2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
optimalTrajectory = NewAveragePath
%Plot at 10 seconds Try #2
figure
plot(optimalTrajectory(:, 3), optimalTrajectory(:, 4))
hold on
OptimalAt10 = median([optimalValues1(1,:);
    optimalValues2(1,:)
    optimalValues3(1,:)
    optimalValues4(1,:)
    optimalValues5(1,:)
    optimalValues6(1,:)
    optimalValues7(1,:)]);
plot(OptimalAt10(3), OptimalAt10(4), 'ko','MarkerSize',10)

ind = median([median(indexClosest1(1:3,1,1))
    median(indexClosest2(1:3,1,1))
    median(indexClosest3(1:3,1,1))
    median(indexClosest4(1:3,1,1))
    median(indexClosest5(1:3,1,1))
    median(indexClosest6(1:3,1,1))
    median(indexClosest7(1:3,1,1))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'g+', 'MarkerSize', 10)


ind = median([median(indexClosest1(1:3,2,1))
    median(indexClosest2(1:3,2,1))
    median(indexClosest3(1:3,2,1))
    median(indexClosest4(1:3,2,1))
    median(indexClosest5(1:3,2,1))
    median(indexClosest6(1:3,2,1))
    median(indexClosest7(1:3,2,1))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'ch', 'MarkerSize', 10)

ind = median([median(indexClosest1(1:3,3,1))
    median(indexClosest2(1:3,3,1))
    median(indexClosest3(1:3,3,1))
    median(indexClosest4(1:3,3,1))
    median(indexClosest5(1:3,3,1))
    median(indexClosest6(1:3,3,1))
    median(indexClosest7(1:3,3,1))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'm*', 'MarkerSize', 10)

ind = median([median(indexClosest1(1:3,4,1))
    median(indexClosest2(1:3,4,1))
    median(indexClosest3(1:3,4,1))
    median(indexClosest4(1:3,4,1))
    median(indexClosest5(1:3,4,1))
    median(indexClosest6(1:3,4,1))
    median(indexClosest7(1:3,4,1))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'rd', 'MarkerSize', 10)
legend('Optimal Trajectory', 'Joystick','Discrete','Cont', 'DiscreteTrimmed','Random','location','eastoutside')
title(['Conditions at 10 seconds, Average'])

%Plot at 30 seconds Try #2
figure
plot(optimalTrajectory(:, 3), optimalTrajectory(:, 4))
hold on
OptimalAt10 = median([optimalValues1(2,:);
    optimalValues2(2,:)
    optimalValues3(2,:)
    optimalValues4(2,:)
    optimalValues5(2,:)
    optimalValues6(2,:)
    optimalValues7(2,:)]);
plot(OptimalAt10(3), OptimalAt10(4), 'ko','MarkerSize',10)

ind = median([median(indexClosest1(1:3,1,2))
    median(indexClosest2(1:3,1,2))
    median(indexClosest3(1:3,1,2))
    median(indexClosest4(1:3,1,2))
    median(indexClosest5(1:3,1,2))
    median(indexClosest6(1:3,1,2))
    median(indexClosest7(1:3,1,2))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'g+', 'MarkerSize', 10)


ind = median([median(indexClosest1(1:3,2,2))
    median(indexClosest2(1:3,2,2))
    median(indexClosest3(1:3,2,2))
    median(indexClosest4(1:3,2,2))
    median(indexClosest5(1:3,2,2))
    median(indexClosest6(1:3,2,2))
    median(indexClosest7(1:3,2,2))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'ch', 'MarkerSize', 10)

ind = median([median(indexClosest1(1:3,3,2))
    median(indexClosest2(1:3,3,2))
    median(indexClosest3(1:3,3,2))
    median(indexClosest4(1:3,3,2))
    median(indexClosest5(1:3,3,2))
    median(indexClosest6(1:3,3,2))
    median(indexClosest7(1:3,3,2))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'm*', 'MarkerSize', 10)

ind = median([median(indexClosest1(1:3,4,2))
    median(indexClosest2(1:3,4,2))
    median(indexClosest3(1:3,4,2))
    median(indexClosest4(1:3,4,2))
    median(indexClosest5(1:3,4,2))
    median(indexClosest6(1:3,4,2))
    median(indexClosest7(1:3,4,2))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'rd', 'MarkerSize', 10)
legend('Optimal Trajectory', 'Joystick','Discrete','Cont', 'DiscreteTrimmed','Random','location','eastoutside')
title(['Conditions at 30 seconds, Average'])


%Plot at 60 seconds Try #2
figure
plot(optimalTrajectory(:, 3), optimalTrajectory(:, 4))
hold on
OptimalAt10 = median([optimalValues1(3,:);
    optimalValues2(3,:)
    optimalValues3(3,:)
    optimalValues4(3,:)
    optimalValues5(3,:)
    optimalValues6(3,:)
    optimalValues7(3,:)]);
plot(OptimalAt10(3), OptimalAt10(4), 'ko','MarkerSize',10)

ind = median([median(indexClosest1(1:3,1,3))
    median(indexClosest2(1:3,1,3))
    median(indexClosest3(1:3,1,3))
    median(indexClosest4(1:3,1,3))
    median(indexClosest5(1:3,1,3))
    median(indexClosest6(1:3,1,3))
    median(indexClosest7(1:3,1,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'g+', 'MarkerSize', 10)


ind = median([median(indexClosest1(1:3,2,3))
    median(indexClosest2(1:3,2,3))
    median(indexClosest3(1:3,2,3))
    median(indexClosest4(1:3,2,3))
    median(indexClosest5(1:3,2,3))
    median(indexClosest6(1:3,2,3))
    median(indexClosest7(1:3,2,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'ch', 'MarkerSize', 10)

ind = median([median(indexClosest1(1:3,3,3))
    median(indexClosest2(1:3,3,3))
    median(indexClosest3(1:3,3,3))
    median(indexClosest4(1:3,3,3))
    median(indexClosest5(1:3,3,3))
    median(indexClosest6(1:3,3,3))
    median(indexClosest7(1:3,3,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'm*', 'MarkerSize', 10)

ind = median([median(indexClosest1(1:3,4,3))
    median(indexClosest2(1:3,4,3))
    median(indexClosest3(1:3,4,3))
    median(indexClosest4(1:3,4,3))
    median(indexClosest5(1:3,4,3))
    median(indexClosest6(1:3,4,3))
    median(indexClosest7(1:3,4,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'rd', 'MarkerSize', 10)
legend('Optimal Trajectory', 'Joystick','Discrete','Cont', 'DiscreteTrimmed','Random','location','eastoutside')
title(['Conditions at 60 seconds, Average'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%superplot
%Plot at 37 seconds Try #3
figure
subplot(3, 1, 1)
plot(optimalTrajectory(:, 3), optimalTrajectory(:, 4))
hold on
OptimalAt10 = median([optimalValues1(3,:);
    optimalValues2(3,:)
    optimalValues3(3,:)
    optimalValues4(3,:)
    optimalValues5(3,:)
    optimalValues6(3,:)
    optimalValues7(3,:)]);
plot(OptimalAt10(3), OptimalAt10(4), 'ko','MarkerSize',10)

ind = median([median(indexClosest1(1:3,1,3))
    median(indexClosest2(1:3,1,3))
    median(indexClosest3(1:3,1,3))
    median(indexClosest4(1:3,1,3))
    median(indexClosest5(1:3,1,3))
    median(indexClosest6(1:3,1,3))
    median(indexClosest7(1:3,1,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'g+', 'MarkerSize', 10)


ind = median([median(indexClosest1(1:3,2,3))
    median(indexClosest2(1:3,2,3))
    median(indexClosest3(1:3,2,3))
    median(indexClosest4(1:3,2,3))
    median(indexClosest5(1:3,2,3))
    median(indexClosest6(1:3,2,3))
    median(indexClosest7(1:3,2,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'ch', 'MarkerSize', 10)

ind = median([median(indexClosest1(1:3,3,3))
    median(indexClosest2(1:3,3,3))
    median(indexClosest3(1:3,3,3))
    median(indexClosest4(1:3,3,3))
    median(indexClosest5(1:3,3,3))
    median(indexClosest6(1:3,3,3))
    median(indexClosest7(1:3,3,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'm*', 'MarkerSize', 10)

ind = median([median(indexClosest1(1:3,4,3))
    median(indexClosest2(1:3,4,3))
    median(indexClosest3(1:3,4,3))
    median(indexClosest4(1:3,4,3))
    median(indexClosest5(1:3,4,3))
    median(indexClosest6(1:3,4,3))
    median(indexClosest7(1:3,4,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'rd', 'MarkerSize', 10)
legend('Optimal Trajectory', 'Joystick','Discrete','Cont', 'DiscreteTrimmed','Random','location','eastoutside')
title(['Conditions at 37 seconds, Dispersion'])

ind = sort([(indexClosest1(1:3,1,3))
    (indexClosest2(1:3,1,3))
    (indexClosest3(1:3,1,3))
    (indexClosest4(1:3,1,3))
    (indexClosest5(1:3,1,3))
    (indexClosest6(1:3,1,3))
    (indexClosest7(1:3,1,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'g+', 'MarkerSize', 10)
stdDisc = std(ind)/1519
meanDisc = mean(ind)/1519
ind = sort([(indexClosest1(1:3,2,3))
    (indexClosest2(1:3,2,3))
    (indexClosest3(1:3,2,3))
    (indexClosest4(1:3,2,3))
    (indexClosest5(1:3,2,3))
    (indexClosest6(1:3,2,3))
    (indexClosest7(1:3,2,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'ch', 'MarkerSize', 10)
stdCont = std(ind)/1519
meanCont = mean(ind)/1519
ind = sort([(indexClosest1(1:3,3,3))
    (indexClosest2(1:3,3,3))
    (indexClosest3(1:3,3,3))
    (indexClosest4(1:3,3,3))
    (indexClosest5(1:3,3,3))
    (indexClosest6(1:3,3,3))
    (indexClosest7(1:3,3,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'm*', 'MarkerSize', 10)
stdDiscTrim = std(ind)/1519
meanDiscTrim = mean(ind)/1519

ind = sort([(indexClosest1(1:3,4,3))
    (indexClosest2(1:3,4,3))
    (indexClosest3(1:3,4,3))
    (indexClosest4(1:3,4,3))
    (indexClosest5(1:3,4,3))
    (indexClosest6(1:3,4,3))
    (indexClosest7(1:3,4,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'rd', 'MarkerSize', 10)
stdRand = std(ind)/1519
meanRand = mean(ind)/1519

subplot(3, 1, 2)

plot(optimalTrajectory(:, 3), optimalTrajectory(:, 4))
hold on
OptimalAt10 = median([optimalValues1(3,:);
    optimalValues2(3,:)
    optimalValues3(3,:)
    optimalValues4(3,:)
    optimalValues5(3,:)
    optimalValues6(3,:)
    optimalValues7(3,:)]);
plot(OptimalAt10(3), OptimalAt10(4), 'ko','MarkerSize',10)

ind = median([median(indexClosest1(1:3,1,3))
    median(indexClosest2(1:3,1,3))
    median(indexClosest3(1:3,1,3))
    median(indexClosest4(1:3,1,3))
    median(indexClosest5(1:3,1,3))
    median(indexClosest6(1:3,1,3))
    median(indexClosest7(1:3,1,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'g+', 'MarkerSize', 10)


ind = median([median(indexClosest1(1:3,2,3))
    median(indexClosest2(1:3,2,3))
    median(indexClosest3(1:3,2,3))
    median(indexClosest4(1:3,2,3))
    median(indexClosest5(1:3,2,3))
    median(indexClosest6(1:3,2,3))
    median(indexClosest7(1:3,2,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'ch', 'MarkerSize', 10)

ind = median([median(indexClosest1(1:3,3,3))
    median(indexClosest2(1:3,3,3))
    median(indexClosest3(1:3,3,3))
    median(indexClosest4(1:3,3,3))
    median(indexClosest5(1:3,3,3))
    median(indexClosest6(1:3,3,3))
    median(indexClosest7(1:3,3,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'm*', 'MarkerSize', 10)

ind = median([median(indexClosest1(1:3,4,3))
    median(indexClosest2(1:3,4,3))
    median(indexClosest3(1:3,4,3))
    median(indexClosest4(1:3,4,3))
    median(indexClosest5(1:3,4,3))
    median(indexClosest6(1:3,4,3))
    median(indexClosest7(1:3,4,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'rd', 'MarkerSize', 10)
legend('Optimal Trajectory', 'Joystick','Discrete','Cont', 'DiscreteTrimmed','Random','location','eastoutside')
title(['Conditions at 37 seconds, Average'])

subplot(3, 1, 3)
indDisc = sort([(indexClosest1(1:3,1,3)/1519)
    (indexClosest2(1:3,1,3)/1519)
    (indexClosest3(1:3,1,3)/1519)
    (indexClosest4(1:3,1,3)/1519)
    (indexClosest5(1:3,1,3)/1519)
    (indexClosest6(1:3,1,3)/1519)
    (indexClosest7(1:3,1,3)/1519)]);

indCont = sort([(indexClosest1(1:3,2,3)/1519)
    (indexClosest2(1:3,2,3)/1519)
    (indexClosest3(1:3,2,3)/1519)
    (indexClosest4(1:3,2,3)/1519)
    (indexClosest5(1:3,2,3)/1519)
    (indexClosest6(1:3,2,3)/1519)
    (indexClosest7(1:3,2,3)/1519)]);

indDiscTrim = sort([(indexClosest1(1:3,3,3)/1519)
    (indexClosest2(1:3,3,3)/1519)
    (indexClosest3(1:3,3,3)/1519)
    (indexClosest4(1:3,3,3)/1519)
    (indexClosest5(1:3,3,3)/1519)
    (indexClosest6(1:3,3,3)/1519)
    (indexClosest7(1:3,3,3)/1519)]);

indRand = sort([(indexClosest1(1:3,4,3)/1519)
    (indexClosest2(1:3,4,3)/1519)
    (indexClosest3(1:3,4,3)/1519)
    (indexClosest4(1:3,4,3)/1519)
    (indexClosest5(1:3,4,3)/1519)
    (indexClosest6(1:3,4,3)/1519)
    (indexClosest7(1:3,4,3)/1519)]);
group = [repmat({'Discrete'}, size(indDisc, 1), 1);
    repmat({'Continuous'}, size(indCont, 1), 1);
    repmat({'Discrete Trimmed'}, size(indDiscTrim, 1), 1);
    repmat({'Random'}, size(indRand, 1), 1)];

boxplot([indDisc; indCont; indDiscTrim; indRand], group)
title(['Conditions at 37 seconds, Percentage of maze completed'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%Plot at 120 seconds Try #2
figure
plot(optimalTrajectory(:, 3), optimalTrajectory(:, 4))
hold on
OptimalAt10 = median([optimalValues1(4,:);
    optimalValues2(4,:)
    optimalValues3(4,:)
    optimalValues4(4,:)
    optimalValues5(4,:)
    optimalValues6(4,:)
    optimalValues7(4,:)]);
plot(OptimalAt10(3), OptimalAt10(4), 'ko','MarkerSize',10)

ind = median([median(indexClosest1(1:3,1,4))
    median(indexClosest2(1:3,1,4))
    median(indexClosest3(1:3,1,4))
    median(indexClosest4(1:3,1,4))
    median(indexClosest5(1:3,1,4))
    median(indexClosest6(1:3,1,4))
    median(indexClosest7(1:3,1,4))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'g+', 'MarkerSize', 10)


ind = median([median(indexClosest1(1:3,2,4))
    median(indexClosest2(1:3,2,4))
    median(indexClosest3(1:3,2,4))
    median(indexClosest4(1:3,2,4))
    median(indexClosest5(1:3,2,4))
    median(indexClosest6(1:3,2,4))
    median(indexClosest7(1:3,2,4))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'ch', 'MarkerSize', 10)

ind = median([median(indexClosest1(1:3,3,4))
    median(indexClosest2(1:3,3,4))
    median(indexClosest3(1:3,3,4))
    median(indexClosest4(1:3,3,4))
    median(indexClosest5(1:3,3,4))
    median(indexClosest6(1:3,3,4))
    median(indexClosest7(1:3,3,4))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'm*', 'MarkerSize', 10)

ind = median([median(indexClosest1(1:3,4,4))
    median(indexClosest2(1:3,4,4))
    median(indexClosest3(1:3,4,4))
    median(indexClosest4(1:3,4,4))
    median(indexClosest5(1:3,4,4))
    median(indexClosest6(1:3,4,4))
    median(indexClosest7(1:3,4,4))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'rd', 'MarkerSize', 10)
legend('Optimal Trajectory', 'Joystick','Discrete','Cont', 'DiscreteTrimmed','Random','location','eastoutside')
title(['Conditions at 120 seconds, Average'])

%Plot at 400 seconds Try #2
figure
plot(optimalTrajectory(:, 3), optimalTrajectory(:, 4))
hold on
OptimalAt10 = median([optimalValues1(5,:);
    optimalValues2(5,:)
    optimalValues3(5,:)
    optimalValues4(5,:)
    optimalValues5(5,:)
    optimalValues6(5,:)
    optimalValues7(5,:)]);
plot(OptimalAt10(3), OptimalAt10(4), 'ko','MarkerSize',10)

ind = median([median(indexClosest1(1:3,1,5))
    median(indexClosest2(1:3,1,5))
    median(indexClosest3(1:3,1,5))
    median(indexClosest4(1:3,1,5))
    median(indexClosest5(1:3,1,5))
    median(indexClosest6(1:3,1,5))
    median(indexClosest7(1:3,1,5))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'g+', 'MarkerSize', 10)


ind = median([median(indexClosest1(1:3,2,5))
    median(indexClosest2(1:3,2,5))
    median(indexClosest3(1:3,2,5))
    median(indexClosest4(1:3,2,5))
    median(indexClosest5(1:3,2,5))
    median(indexClosest6(1:3,2,5))
    median(indexClosest7(1:3,2,5))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'ch', 'MarkerSize', 10)

ind = median([median(indexClosest1(1:3,3,5))
    median(indexClosest2(1:3,3,5))
    median(indexClosest3(1:3,3,5))
    median(indexClosest4(1:3,3,5))
    median(indexClosest5(1:3,3,5))
    median(indexClosest6(1:3,3,5))
    median(indexClosest7(1:3,3,5))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'm*', 'MarkerSize', 10)

ind = median([median(indexClosest1(1:3,4,5))
    median(indexClosest2(1:3,4,5))
    median(indexClosest3(1:3,4,5))
    median(indexClosest4(1:3,4,5))
    median(indexClosest5(1:3,4,5))
    median(indexClosest6(1:3,4,5))
    median(indexClosest7(1:3,4,5))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'rd', 'MarkerSize', 10)
legend('Optimal Trajectory', 'Joystick','Discrete','Cont', 'DiscreteTrimmed','Random','location','eastoutside')
title(['Conditions at 400 seconds, Average'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%At max positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%Plot at 10 seconds Try #2.5
figure
plot(optimalTrajectory(:, 3), optimalTrajectory(:, 4))
hold on
OptimalAt10 = max([optimalValues1(1,:);
    optimalValues2(1,:)
    optimalValues3(1,:)
    optimalValues4(1,:)
    optimalValues5(1,:)
    optimalValues6(1,:)
    optimalValues7(1,:)]);
plot(OptimalAt10(3), OptimalAt10(4), 'ko','MarkerSize',10)

ind = max([max(indexClosest1(1:3,1,1))
    max(indexClosest2(1:3,1,1))
    max(indexClosest3(1:3,1,1))
    max(indexClosest4(1:3,1,1))
    max(indexClosest5(1:3,1,1))
    max(indexClosest6(1:3,1,1))
    max(indexClosest7(1:3,1,1))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'g+', 'MarkerSize', 10)


ind = max([max(indexClosest1(1:3,2,1))
    max(indexClosest2(1:3,2,1))
    max(indexClosest3(1:3,2,1))
    max(indexClosest4(1:3,2,1))
    max(indexClosest5(1:3,2,1))
    max(indexClosest6(1:3,2,1))
    max(indexClosest7(1:3,2,1))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'ch', 'MarkerSize', 10)

ind = max([max(indexClosest1(1:3,3,1))
    max(indexClosest2(1:3,3,1))
    max(indexClosest3(1:3,3,1))
    max(indexClosest4(1:3,3,1))
    max(indexClosest5(1:3,3,1))
    max(indexClosest6(1:3,3,1))
    max(indexClosest7(1:3,3,1))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'm*', 'MarkerSize', 10)

ind = max([max(indexClosest1(1:3,4,1))
    max(indexClosest2(1:3,4,1))
    max(indexClosest3(1:3,4,1))
    max(indexClosest4(1:3,4,1))
    max(indexClosest5(1:3,4,1))
    max(indexClosest6(1:3,4,1))
    max(indexClosest7(1:3,4,1))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'rd', 'MarkerSize', 10)
legend('Optimal Trajectory', 'Joystick','Discrete','Cont', 'DiscreteTrimmed','Random','location','eastoutside')
title(['Conditions at 10 seconds, Furthest point'])



%Plot at 30 seconds Try #2
figure
plot(optimalTrajectory(:, 3), optimalTrajectory(:, 4))
hold on
OptimalAt10 = max([optimalValues1(2,:);
    optimalValues2(2,:)
    optimalValues3(2,:)
    optimalValues4(2,:)
    optimalValues5(2,:)
    optimalValues6(2,:)
    optimalValues7(2,:)]);
plot(OptimalAt10(3), OptimalAt10(4), 'ko','MarkerSize',10)

ind = max([median(indexClosest1(1:3,1,2))
    max(indexClosest2(1:3,1,2))
    max(indexClosest3(1:3,1,2))
    max(indexClosest4(1:3,1,2))
    max(indexClosest5(1:3,1,2))
    max(indexClosest6(1:3,1,2))
    max(indexClosest7(1:3,1,2))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'g+', 'MarkerSize', 10)


ind = max([max(indexClosest1(1:3,2,2))
    max(indexClosest2(1:3,2,2))
    max(indexClosest3(1:3,2,2))
    max(indexClosest4(1:3,2,2))
    max(indexClosest5(1:3,2,2))
    max(indexClosest6(1:3,2,2))
    max(indexClosest7(1:3,2,2))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'ch', 'MarkerSize', 10)

ind = max([max(indexClosest1(1:3,3,2))
    max(indexClosest2(1:3,3,2))
    max(indexClosest3(1:3,3,2))
    max(indexClosest4(1:3,3,2))
    max(indexClosest5(1:3,3,2))
    max(indexClosest6(1:3,3,2))
    max(indexClosest7(1:3,3,2))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'm*', 'MarkerSize', 10)

ind = max([max(indexClosest1(1:3,4,2))
    max(indexClosest2(1:3,4,2))
    max(indexClosest3(1:3,4,2))
    max(indexClosest4(1:3,4,2))
    max(indexClosest5(1:3,4,2))
    max(indexClosest6(1:3,4,2))
    max(indexClosest7(1:3,4,2))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'rd', 'MarkerSize', 10)
legend('Optimal Trajectory', 'Joystick','Discrete','Cont', 'DiscreteTrimmed','Random','location','eastoutside')
title(['Conditions at 30 seconds, Furthest point'])


%Plot at 60 seconds Try #2
figure
plot(optimalTrajectory(:, 3), optimalTrajectory(:, 4))
hold on
OptimalAt10 = max([optimalValues1(3,:);
    optimalValues2(3,:)
    optimalValues3(3,:)
    optimalValues4(3,:)
    optimalValues5(3,:)
    optimalValues6(3,:)
    optimalValues7(3,:)]);
plot(OptimalAt10(3), OptimalAt10(4), 'ko','MarkerSize',10)

ind = max([max(indexClosest1(1:3,1,3))
    max(indexClosest2(1:3,1,3))
    max(indexClosest3(1:3,1,3))
    max(indexClosest4(1:3,1,3))
    max(indexClosest5(1:3,1,3))
    max(indexClosest6(1:3,1,3))
    max(indexClosest7(1:3,1,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'g+', 'MarkerSize', 10)


ind = max([max(indexClosest1(1:3,2,3))
    max(indexClosest2(1:3,2,3))
    max(indexClosest3(1:3,2,3))
    max(indexClosest4(1:3,2,3))
    max(indexClosest5(1:3,2,3))
    max(indexClosest6(1:3,2,3))
    max(indexClosest7(1:3,2,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'ch', 'MarkerSize', 10)

ind = max([max(indexClosest1(1:3,3,3))
    max(indexClosest2(1:3,3,3))
    max(indexClosest3(1:3,3,3))
    max(indexClosest4(1:3,3,3))
    max(indexClosest5(1:3,3,3))
    max(indexClosest6(1:3,3,3))
    max(indexClosest7(1:3,3,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'm*', 'MarkerSize', 10)

ind = max([max(indexClosest1(1:3,4,3))
    max(indexClosest2(1:3,4,3))
    max(indexClosest3(1:3,4,3))
    max(indexClosest4(1:3,4,3))
    max(indexClosest5(1:3,4,3))
    max(indexClosest6(1:3,4,3))
    max(indexClosest7(1:3,4,3))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'rd', 'MarkerSize', 10)
legend('Optimal Trajectory', 'Joystick','Discrete','Cont', 'DiscreteTrimmed','Random','location','eastoutside')
title(['Conditions at 60 seconds, Furthest point'])


%Plot at 120 seconds Try #2
figure
plot(optimalTrajectory(:, 3), optimalTrajectory(:, 4))
hold on
OptimalAt10 = max([optimalValues1(4,:);
    optimalValues2(4,:)
    optimalValues3(4,:)
    optimalValues4(4,:)
    optimalValues5(4,:)
    optimalValues6(4,:)
    optimalValues7(4,:)]);
plot(OptimalAt10(3), OptimalAt10(4), 'ko','MarkerSize',10)

ind = max([max(indexClosest1(1:3,1,4))
    max(indexClosest2(1:3,1,4))
    max(indexClosest3(1:3,1,4))
    max(indexClosest4(1:3,1,4))
    max(indexClosest5(1:3,1,4))
    max(indexClosest6(1:3,1,4))
    max(indexClosest7(1:3,1,4))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'g+', 'MarkerSize', 10)


ind = max([max(indexClosest1(1:3,2,4))
    max(indexClosest2(1:3,2,4))
    max(indexClosest3(1:3,2,4))
    max(indexClosest4(1:3,2,4))
    max(indexClosest5(1:3,2,4))
    max(indexClosest6(1:3,2,4))
    max(indexClosest7(1:3,2,4))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'ch', 'MarkerSize', 10)

ind = max([max(indexClosest1(1:3,3,4))
    max(indexClosest2(1:3,3,4))
    max(indexClosest3(1:3,3,4))
    max(indexClosest4(1:3,3,4))
    max(indexClosest5(1:3,3,4))
    max(indexClosest6(1:3,3,4))
    max(indexClosest7(1:3,3,4))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'm*', 'MarkerSize', 10)

ind = max([max(indexClosest1(1:3,4,4))
    max(indexClosest2(1:3,4,4))
    max(indexClosest3(1:3,4,4))
    max(indexClosest4(1:3,4,4))
    max(indexClosest5(1:3,4,4))
    max(indexClosest6(1:3,4,4))
    max(indexClosest7(1:3,4,4))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'rd', 'MarkerSize', 10)
legend('Optimal Trajectory', 'Joystick','Discrete','Cont', 'DiscreteTrimmed','Random','location','eastoutside')
title(['Conditions at 120 seconds, Furthest point'])


%Plot at 400 seconds Try #2
figure
plot(optimalTrajectory(:, 3), optimalTrajectory(:, 4))
hold on
OptimalAt10 = max([optimalValues1(5,:);
    optimalValues2(5,:)
    optimalValues3(5,:)
    optimalValues4(5,:)
    optimalValues5(5,:)
    optimalValues6(5,:)
    optimalValues7(5,:)]);
plot(OptimalAt10(3), OptimalAt10(4), 'ko','MarkerSize',10)

ind = max([max(indexClosest1(1:3,1,5))
    max(indexClosest2(1:3,1,5))
    max(indexClosest3(1:3,1,5))
    max(indexClosest4(1:3,1,5))
    max(indexClosest5(1:3,1,5))
    max(indexClosest6(1:3,1,5))
    max(indexClosest7(1:3,1,5))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'g+', 'MarkerSize', 10)


ind = max([max(indexClosest1(1:3,2,5))
    max(indexClosest2(1:3,2,5))
    max(indexClosest3(1:3,2,5))
    max(indexClosest4(1:3,2,5))
    max(indexClosest5(1:3,2,5))
    max(indexClosest6(1:3,2,5))
    max(indexClosest7(1:3,2,5))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'ch', 'MarkerSize', 10)

ind = max([max(indexClosest1(1:3,3,5))
    max(indexClosest2(1:3,3,5))
    max(indexClosest3(1:3,3,5))
    max(indexClosest4(1:3,3,5))
    max(indexClosest5(1:3,3,5))
    max(indexClosest6(1:3,3,5))
    max(indexClosest7(1:3,3,5))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'm*', 'MarkerSize', 10)

ind = max([max(indexClosest1(1:3,4,5))
    max(indexClosest2(1:3,4,5))
    max(indexClosest3(1:3,4,5))
    max(indexClosest4(1:3,4,5))
    max(indexClosest5(1:3,4,5))
    max(indexClosest6(1:3,4,5))
    max(indexClosest7(1:3,4,5))]);
plot(optimalTrajectory(ind, 3), optimalTrajectory(ind, 4), 'rd', 'MarkerSize', 10)
legend('Optimal Trajectory', 'Joystick','Discrete','Cont', 'DiscreteTrimmed','Random','location','eastoutside')
title(['Conditions at 400 seconds, Furthest point'])




