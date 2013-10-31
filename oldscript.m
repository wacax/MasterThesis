function [TrajectoryStructure ResultsStructure EventsStructure] = TrajetoryAnalysis(VPStructure)


%Import Trajectories
% 
% House1Random = importfile('RandomH1.csv', 'RandomH1');
% House2Random = importfile('RandomH2.csv', 'RandomH2');
% House3Random = importfile('RandomH3.csv', 'RandomH3');
% House1Random = House1Random(House1Random(:,3)~=0, :);
% House2Random = House2Random(House2Random(:,3)~=0, :);
% House3Random = House3Random(House3Random(:,3)~=0, :);
% House1Cont = importfile('House1Cont.csv', 'House1Cont');
% House2Cont = importfile('House2Cont.csv', 'House2Cont');
% House3Cont = importfile('House3Cont.csv', 'House3Cont');
% House1Cont = House1Cont(House1Cont(:,3)~=0, :);
% House2Cont = House2Cont(House2Cont(:,3)~=0, :);
% House3Cont = House3Cont(House2Cont(:,3)~=0, :);
% House1Discrete = importfile('House1Discrete.csv', 'House1Discrete');
% House2Discrete = importfile('House2Discrete.csv', 'House2Discrete');
% House3Discrete = importfile('House3Discrete.csv', 'House3Discrete');
% House1Discrete = House1Discrete(House1Discrete(:,3)~=0, :);
% House2Discrete = House2Discrete(House2Discrete(:,3)~=0, :);
% House3Discrete = House3Discrete(House3Discrete(:,3)~=0, :);
% 
% Joy1_1 = importfile('Joy1.csv', 'Joy1');
% Joy1_1 = Joy1_1(Joy1_1(:,3)~=0, :);
% Joy1_2 = importfile('Joy4.csv', 'Joy4');
% Joy1_2 = Joy1_2(Joy1_2(:,3)~=0, :);
% Joy2_1 = importfile('Joy2.csv', 'Joy2');
% Joy2_1 = Joy2_1(Joy2_1(:,3)~=0, :);
% Joy2_2 = importfile('Joy5.csv', 'Joy5');
% Joy2_2 = Joy2_2(Joy2_2(:,3)~=0, :);
% Joy3_1 = importfile('Joy3.csv', 'Joy3');
% Joy3_1 = Joy3_1(Joy3_1(:,3)~=0, :);
% Joy3_2 = importfile('Joy6.csv', 'Joy6');
% Joy3_2 = Joy3_2(Joy3_2(:,3)~=0, :);

House1Random = importfile(RandomH1.csv, 'RandomH1');
House2Random = importfile(RandomH2.csv, 'RandomH2');
House3Random = importfile(RandomH3.csv, 'RandomH3');
House1Random = House1Random(House1Random(:,3)~=0, :);
House2Random = House2Random(House2Random(:,3)~=0, :);
House3Random = House3Random(House3Random(:,3)~=0, :);
House1Cont = importfile(House1Cont.csv, 'House1Cont');
House2Cont = importfile(House2Cont.csv, 'House2Cont');
House3Cont = importfile(House3Cont.csv, 'House3Cont');
House1Cont = House1Cont(House1Cont(:,3)~=0, :);
House2Cont = House2Cont(House2Cont(:,3)~=0, :);
House3Cont = House3Cont(House2Cont(:,3)~=0, :);
House1Discrete = importfile(House1Discrete.csv, 'House1Discrete');
House2Discrete = importfile(House2Discrete.csv, 'House2Discrete');
House3Discrete = importfile(House3Discrete.csv, 'House3Discrete');
House1Discrete = House1Discrete(House1Discrete(:,3)~=0, :);
House2Discrete = House2Discrete(House2Discrete(:,3)~=0, :);
House3Discrete = House3Discrete(House3Discrete(:,3)~=0, :);

Joy1_1 = importfile(Joy1.csv, 'Joy1');
Joy1_1 = Joy1_1(Joy1_1(:,3)~=0, :);
Joy1_2 = importfile(Joy4.csv, 'Joy4');
Joy1_2 = Joy1_2(Joy1_2(:,3)~=0, :);
Joy2_1 = importfile(Joy2.csv, 'Joy2');
Joy2_1 = Joy2_1(Joy2_1(:,3)~=0, :);
Joy2_2 = importfile(Joy5.csv, 'Joy5');
Joy2_2 = Joy2_2(Joy2_2(:,3)~=0, :);
Joy3_1 = importfile(Joy3.csv, 'Joy3');
Joy3_1 = Joy3_1(Joy3_1(:,3)~=0, :);
Joy3_2 = importfile(Joy6.csv, 'Joy6');
Joy3_2 = Joy3_2(Joy3_2(:,3)~=0, :);

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

%Joy1_1(round(linspace(1,2700,2400)),4)

[~, indices1_1] = datasample(Joy1_1, numSamples1, 1);
[~, indices1_2] = datasample(Joy1_2, numSamples1, 1);
[~, indices2_1] = datasample(Joy2_1, numSamples2, 1);
[~, indices2_2] = datasample(Joy2_2, numSamples2, 1);
[~, indices3_1] = datasample(Joy3_1, numSamples3, 1);
[~, indices3_2] = datasample(Joy3_2, numSamples3, 1);

House1Benchmark = (Joy1_1(sort(indices1_1), :) + Joy1_2(sort(indices1_2), :))./2;
House2Benchmark = (Joy2_1(sort(indices2_1), :) + Joy2_2(sort(indices2_2), :))./2;
House3Benchmark = (Joy3_1(sort(indices3_1), :) + Joy3_2(sort(indices3_2), :))./2;

%House1 Analysis
%Distance Metrics
f12_1 =frechet(House1Benchmark(:,3), House1Benchmark(:,4), House1Discrete(:,3), House1Discrete(:,4));
f13_1 =frechet(House1Benchmark(:,3), House1Benchmark(:,4), House1Cont(:,3), House1Cont(:,4));
f14_1 =frechet(House1Benchmark(:,3), House1Benchmark(:,4), House1DiscreteTrimmed(:,3), House1DiscreteTrimmed(:,4));
f15_1 =frechet(House1Benchmark(:,3), House1Benchmark(:,4), House1Random(:,3), House1Random(:,4));
f23_1 = 0;
%f23_1 =frechet(House1Discrete(:,3), House1Discrete(:,4), House1Cont(:,3), House1Cont(:,4));
f24_1 =frechet(House1Discrete(:,3), House1Discrete(:,4), House1DiscreteTrimmed(:,3), House1DiscreteTrimmed(:,4));
f34_1 = 0;
%f34_1 =frechet(House1Cont(:,3), House1Cont(:,4), House1DiscreteTrimmed(:,3), House1DiscreteTrimmed(:,4));
f11_1 = 0;
f22_1 = 0;
f33_1 = 0;
f44_1 = 0;
%f11_1 =frechet(House1Benchmark(:,3), House1Benchmark(:,4) ,House1Benchmark(:,3), House1Benchmark(:,4));
%f22_1 =frechet(House1Discrete(:,3), House1Discrete(:,4), House1Discrete(:,3), House1Discrete(:,4));
%f33_1 =frechet(House1Cont(:,3), House1Cont(:,4), House1Cont(:,3), House1Cont(:,4));
%f44_1 =frechet(House1DiscreteTrimmed(:,3), House1DiscreteTrimmed(:,4), House1DiscreteTrimmed(:,3), House1DiscreteTrimmed(:,4));
f52_1 = frechet(House1Random(:,3), House1Random(:,4), House1Discrete(:,3), House1Discrete(:,4));
f53_1 = frechet(House1Random(:,3), House1Random(:,4), House1Cont(:,3), House1Cont(:,4));
f54_1 = frechet(House1Random(:,3), House1Random(:,4), House1DiscreteTrimmed(:,3), House1DiscreteTrimmed(:,4));
f55_1 = 0;
%f55_1 = frechet(House1Random(:,3), House1Random(:,4), House1Random(:,3), House1Random(:,4));


figure;
subplot(2,1,1)
hold on
plot(House1Benchmark(:,3) , House1Benchmark(:,4), 'r', 'linewidth', 2)
plot(House1Discrete(:,3), House1Discrete(:,4),'g','linewidth',2)
plot(House1Cont(:,3), House1Cont(:,4),'b','linewidth',2)
plot(House1DiscreteTrimmed(:,3), House1DiscreteTrimmed(:,4),'g','linewidth',2)
plot(House1Random(:,3), House1Random(:,4),'c','linewidth',2)
legend('Benchmark','Discrete','Continuous', 'Discrete Trimmed', 'Random', 'location','eastoutside')
xlabel('X')
ylabel('Y')
axis equal tight
box on
title(['four space curves to compare'])
legend

subplot(2,1,2)
%imagesc([[f11_1,f12_1,f13_1];[f12_1,f22_1,f23_1];[f13_1,f23_1,f33_1]])
imagesc([[f11_1,f12_1,f13_1, f14_1, f15_1];[f12_1,f22_1,f23_1, f24_1, f52_1];[f13_1,f23_1,f33_1, f43_1, f53_1]; [f14_1,f24_1,f34_1, f44_1, f54_1]; [f15_1,f52_1,f53_1, f54_1, f55_1]])
xlabel('curve')
ylabel('curve')
cb1=colorbar('peer',gca);
set(get(cb1,'Ylabel'),'String','Frechet Distance')
axis equal tight

%House2 Analysis
%Distance Metrics
f12_2=frechet(House2Benchmark(:,3), House2Benchmark(:,4), House2Discrete(:,3), House2Discrete(:,4));
f13_2=frechet(House2Benchmark(:,3), House2Benchmark(:,4), House2Cont(:,3), House2Cont(:,4));
f14_2 =frechet(House2Benchmark(:,3), House2Benchmark(:,4), House2DiscreteTrimmed(:,3), House2DiscreteTrimmed(:,4));
f15_2 =frechet(House2Benchmark(:,3), House2Benchmark(:,4), House2Random(:,3), House2Random(:,4));
f23_2 = 0;
%f23_2=frechet(House2Discrete(:,3), House2Discrete(:,4), House2Cont(:,3), House2Cont(:,4));
f24_2 =frechet(House2Discrete(:,3), House2Discrete(:,4), House2DiscreteTrimmed(:,3), House2DiscreteTrimmed(:,4));
f34_2 = 0;
%f34_2 =frechet(House2Cont(:,3), House2Cont(:,4), House2DiscreteTrimmed(:,3), House2DiscreteTrimmed(:,4));
f11_2 = 0;
f22_2 = 0;
f33_2 = 0;
f44_2 = 0;
%f11_2=frechet(House2Benchmark(:,3), House2Benchmark(:,4) ,House2Benchmark(:,3), House2Benchmark(:,4));
%f22_2=frechet(House2Discrete(:,3), House2Discrete(:,4), House2Discrete(:,3), House2Discrete(:,4));
%f33_2=frechet(House2Cont(:,3), House2Cont(:,4), House2Cont(:,3), House2Cont(:,4));
%f44_2 =frechet(House1DiscreteTrimmed(:,3), House2DiscreteTrimmed(:,4), House2DiscreteTrimmed(:,3), House2DiscreteTrimmed(:,4));
f52_2 = frechet(House2Random(:,3), House2Random(:,4), House2Discrete(:,3), House2Discrete(:,4));
f53_2 = frechet(House2Random(:,3), House2Random(:,4), House2Cont(:,3), House2Cont(:,4));
f54_2 = frechet(House2Random(:,3), House2Random(:,4), House2DiscreteTrimmed(:,3), House2DiscreteTrimmed(:,4));
f55_2 = 0;
%f55_2 = frechet(House2Random(:,3), House2Random(:,4), House2Random(:,3), House2Random(:,4));


figure;
subplot(2,1,1)
hold on
plot(House2Benchmark(:,3) , House2Benchmark(:,4), 'r', 'linewidth', 2)
plot(House2Discrete(:,3), House2Discrete(:,4),'g','linewidth',2)
plot(House2Cont(:,3), House2Cont(:,4),'b','linewidth',2)
plot(House2DiscreteTrimmed(:,3), House2DiscreteTrimmed(:,4),'g','linewidth',2)
plot(House2Random(:,3), House2Random(:,4),'c','linewidth',2)
legend('Benchmark','Discrete','Continuous', 'Discrete Trimmed', 'Random', 'location','eastoutside')
xlabel('X')
ylabel('Y')
axis equal tight
box on
title(['four space curves to compare'])
legend

subplot(2,1,2)
%imagesc([[f11_2,f12_2,f13_2];[f12_2,f22_2,f23_2];[f13_2,f23_2,f33_2]])
imagesc([[f11_2,f12_2,f13_2, f14_2, f15_2];[f12_2,f22_2,f23_2, f24_2, f52_2];[f13_2,f23_2,f33_2, f43_2, f53_2]; [f14_2,f24_2,f34_2, f44_2, f54_2]; [f15_2,f52_2,f53_2, f54_2, f55_2]])
xlabel('curve')
ylabel('curve')
cb1=colorbar('peer',gca);
set(get(cb1,'Ylabel'),'String','Frechet Distance')
axis equal tight

%House3 Analysis
%Distance Metrics
f12_3=frechet(House3Benchmark(:,3), House3Benchmark(:,4), House3Discrete(:,3), House3Discrete(:,4));
f13_3=frechet(House3Benchmark(:,3), House3Benchmark(:,4), House3Cont(:,3), House3Cont(:,4));
f14_3=frechet(House3Benchmark(:,3), House3Benchmark(:,4), House3DiscreteTrimmed(:,3), House3DiscreteTrimmed(:,4));
f15_3 =frechet(House2Benchmark(:,3), House2Benchmark(:,4), House2Random(:,3), House2Random(:,4));
f23_3 = 0;
%f23_3=frechet(House3Discrete(:,3), House3Discrete(:,4), House3Cont(:,3), House3Cont(:,4));
f24_3 =frechet(House3Discrete(:,3), House3Discrete(:,4), House3DiscreteTrimmed(:,3), House3DiscreteTrimmed(:,4));
f34_3 = 0;
%f34_3 =frechet(House3Cont(:,3), House3Cont(:,4), House3DiscreteTrimmed(:,3), House3DiscreteTrimmed(:,4));
f11_3 = 0;
f22_3 = 0;
f33_3 = 0;
f44_3 = 0;
%f11_3=frechet(House3Benchmark(:,3), House3Benchmark(:,4) ,House3Benchmark(:,3), House3Benchmark(:,4));
%f22_3=frechet(House3Discrete(:,3), House3Discrete(:,4), House3Discrete(:,3), House3Discrete(:,4));
%f33_3=frechet(House3Cont(:,3), House3Cont(:,4), House3Cont(:,3), House3Cont(:,4));
%f44_3 =frechet(House3DiscreteTrimmed(:,3), House3DiscreteTrimmed(:,4), House3DiscreteTrimmed(:,3), House3DiscreteTrimmed(:,4));
f52_3 = frechet(House3Random(:,3), House3Random(:,4), House3Discrete(:,3), House3Discrete(:,4));
f53_3 = frechet(House3Random(:,3), House3Random(:,4), House3Cont(:,3), House3Cont(:,4));
f54_3 = frechet(House3Random(:,3), House3Random(:,4), House3DiscreteTrimmed(:,3), House3DiscreteTrimmed(:,4));
f55_3 = 0;
%f55_3 = frechet(House3Random(:,3), House3Random(:,4), House3Random(:,3), House3Random(:,4));

figure;
subplot(2,1,1)
hold on
plot(House3Benchmark(:,3) , House3Benchmark(:,4), 'r', 'linewidth', 2)
plot(House3Discrete(:,3), House3Discrete(:,4),'g','linewidth',2)
plot(House3Cont(:,3), House3Cont(:,4),'b','linewidth',2)
plot(House3DiscreteTrimmed(:,3), House3DiscreteTrimmed(:,4),'g','linewidth',2)
plot(House3Random(:,3), House3Random(:,4),'c','linewidth',2)
legend('Benchmark','Discrete','Continuous', 'Discrete Trimmed', 'Random', 'location','eastoutside')
xlabel('X')
ylabel('Y')
axis equal tight
box on
title(['four space curves to compare'])
legend

subplot(2,1,2)
%imagesc([[f11_3,f12_3,f13_3];[f12_3,f22_3,f23_3];[f13_3,f23_3,f33_3]])
imagesc([[f11_3,f12_3,f13_3, f14_3, f15_3];[f12_3,f22_3,f23_3, f24_3, f52_3];[f13_3,f23_3,f33_3, f43_3, f53_3]; [f14_3,f24_3,f34_3, f44_3, f54_3]; [f15_3,f52_3,f53_3, f54_3, f55_3]])
xlabel('curve')
ylabel('curve')
cb1=colorbar('peer',gca);
set(get(cb1,'Ylabel'),'String','Frechet Distance')
axis equal tight

%Trajectories
TrajectoryStructure.House1Cont = House1Cont;
TrajectoryStructure.House2Cont = House2Cont;
TrajectoryStructure.House3Cont = House3Cont;
TrajectoryStructure.House1Discrete = House1Discrete;
TrajectoryStructure.House2Discrete = House2Discrete;
TrajectoryStructure.House3Discrete = House3Discrete;
TrajectoryStructure.House1DiscreteTrimmed = House1DiscreteTrimmed;
TrajectoryStructure.House2DiscreteTrimmed = House2DiscreteTrimmed;
TrajectoryStructure.House3DiscreteTrimmed = House3DiscreteTrimmed;
TrajectoryStructure.House1Benchmark = House1Benchmark;
TrajectoryStructure.House2Benchmark = House2Benchmark;
TrajectoryStructure.House3Benchmark = House3Benchmark;
TrajectoryStructure.House1Random = House1Random;
TrajectoryStructure.House2Random = House2Random;
TrajectoryStructure.House3Random = House3Random;

%Distances
ResultsStructure.DistanceH1BenchmarkDiscrete = f12_1;
ResultsStructure.DistanceH2BenchmarkDiscrete = f12_2;
ResultsStructure.DistanceH3BenchmarkDiscrete = f12_3;
ResultsStructure.DistanceH1BenchmarkContinuous = f13_1;
ResultsStructure.DistanceH2BenchmarkContinuous = f13_2;
ResultsStructure.DistanceH3BenchmarkContinuous = f13_3;
ResultsStructure.DistanceH1BenchmarkDiscreteTrimmed = f14_1;
ResultsStructure.DistanceH2BenchmarkDiscreteTrimmed = f14_2;
ResultsStructure.DistanceH3BenchmarkDiscreteTrimmed = f14_3;
ResultsStructure.DistanceH1DiscreteDiscreteTrimmed =f24_1;
ResultsStructure.DistanceH2DiscreteDiscreteTrimmed =f24_2;
ResultsStructure.DistanceH3DiscreteDiscreteTrimmed =f24_3;
ResultsStructure.DistanceH1BenchmarkRandom =f15_1;
ResultsStructure.DistanceH2BenchmarkRandom =f15_2;
ResultsStructure.DistanceH3BenchmarkRandom =f15_3;
ResultsStructure.DistanceH1RandomDiscrete =f52_1;
ResultsStructure.DistanceH2RandomDiscrete =f52_2;
ResultsStructure.DistanceH3RandomDiscrete =f52_3;
ResultsStructure.DistanceH1RandomContinuous =f53_1;
ResultsStructure.DistanceH2RandomContinuous =f53_2;
ResultsStructure.DistanceH3RandomContinuous =f53_3;
ResultsStructure.DistanceH1RandomDiscreteTrimmed =f54_1;
ResultsStructure.DistanceH2RandomDiscreteTrimmed =f54_2;
ResultsStructure.DistanceH3RandomDiscreteTrimmed =f54_3;

%Events
EventsStructure.House1Cont = House1Cont(House1Cont(:,3)==0, :);
EventsStructure.House2Cont = House2Cont(House2Cont(:,3)==0, :);
EventsStructure.House3Cont = House3Cont(House3Cont(:,3)==0, :);
EventsStructure.House1Discrete = House1Discrete(House1Discrete(:,3)==0, :);
EventsStructure.House2Discrete = House2Discrete(House2Discrete(:,3)==0, :);
EventsStructure.House3Discrete = House3Discrete(House3Discrete(:,3)==0, :);
EventsStructure.House1DiscreteTrimmed = House1DiscreteTrimmed(House1DiscreteTrimmed(:,3)==0, :);
EventsStructure.House2DiscreteTrimmed = House2DiscreteTrimmed(House2DiscreteTrimmed(:,3)==0, :);
EventsStructure.House3DiscreteTrimmed = House3DiscreteTrimmed(House3DiscreteTrimmed(:,3)==0, :);
EventsStructure.House1Benchmark = House1Benchmark(House1Benchmark(:,3)==0, :);
EventsStructure.House2Benchmark = House2Benchmark(House2Benchmark(:,3)==0, :);
EventsStructure.House3Benchmark = House3Benchmark(House3Benchmark(:,3)==0, :);
EventsStructure.House1Random = House1Random(House1Random(:,3)==0, :);
EventsStructure.House2Random = House2Random(House2Random(:,3)==0, :);
EventsStructure.House3Random = House3Random(House3Random(:,3)==0, :);

end