function f = TrajectoryPlots(DistancesMatrix, TrajectoriesMatrix, HouseNumber, Participant)

% figure;
% subplot(2,1,1)
% hold on
% plot(TrajectoriesMatrix.House1Benchmark(:,3) , TrajectoriesMatrix.House1Benchmark(:,4), 'r', 'linewidth', 2)
% plot(TrajectoriesMatrix.House1Discrete(:,3), TrajectoriesMatrix.House1Discrete(:,4),'g','linewidth',2)
% plot(TrajectoriesMatrix.House1Cont(:,3), TrajectoriesMatrix.House1Cont(:,4),'b','linewidth',2)
% plot(TrajectoriesMatrix.House1DiscreteTrimmed(:,3), TrajectoriesMatrix.House1DiscreteTrimmed(:,4),'g','linewidth',2)
% plot(TrajectoriesMatrix.House1Random(:,3), TrajectoriesMatrix.House1Random(:,4),'c','linewidth',2)
% legend('Benchmark','Discrete','Continuous', 'Discrete Trimmed', 'Random', 'location','eastoutside')
% xlabel('X')
% ylabel('Y')
% axis equal tight
% box on
% title(['five space curves to compare'])
% legend

%HERE COMES THE QUESTION OF HOW TO PUT THE NUMBER INTO THE STRING

if HouseNumber == 1
figure;
subplot(2,1,1)
hold on
plot(TrajectoriesMatrix.House1Benchmark(:,3) , TrajectoriesMatrix.House1Benchmark(:,4), 'r', 'linewidth', 2)
plot(TrajectoriesMatrix.House1Discrete(:,3), TrajectoriesMatrix.House1Discrete(:,4),'g','linewidth',2)
plot(TrajectoriesMatrix.House1Cont(:,3), TrajectoriesMatrix.House1Cont(:,4),'b','linewidth',2)
plot(TrajectoriesMatrix.House1DiscreteTrimmed(:,3), TrajectoriesMatrix.House1DiscreteTrimmed(:,4),'g','linewidth',2)
plot(TrajectoriesMatrix.House1Random(:,3), TrajectoriesMatrix.House1Random(:,4),'c','linewidth',2)
legend('Benchmark','Discrete','Continuous', 'Discrete Trimmed', 'Random', 'location','eastoutside')
xlabel('X')
ylabel('Y')
axis equal tight
box on
title(['five space curves to compare'])
legend
elseif HouseNumber == 2
figure;
subplot(2,1,1)
hold on
plot(TrajectoriesMatrix.House2Benchmark(:,3) , TrajectoriesMatrix.House2Benchmark(:,4), 'r', 'linewidth', 2)
plot(TrajectoriesMatrix.House2Discrete(:,3), TrajectoriesMatrix.House2Discrete(:,4),'g','linewidth',2)
plot(TrajectoriesMatrix.House2Cont(:,3), TrajectoriesMatrix.House2Cont(:,4),'b','linewidth',2)
plot(TrajectoriesMatrix.House2DiscreteTrimmed(:,3), TrajectoriesMatrix.House2DiscreteTrimmed(:,4),'g','linewidth',2)
plot(TrajectoriesMatrix.House2Random(:,3), TrajectoriesMatrix.House2Random(:,4),'c','linewidth',2)
legend('Benchmark','Discrete','Continuous', 'Discrete Trimmed', 'Random', 'location','eastoutside')
xlabel('X')
ylabel('Y')
axis equal tight
box on
title(['five space curves to compare'])
legend  
elseif HouseNumber == 3
figure;
subplot(2,1,1)
hold on
plot(TrajectoriesMatrix.House3Benchmark(:,3) , TrajectoriesMatrix.House3Benchmark(:,4), 'r', 'linewidth', 2)
plot(TrajectoriesMatrix.House3Discrete(:,3), TrajectoriesMatrix.House3Discrete(:,4),'g','linewidth',2)
plot(TrajectoriesMatrix.House3Cont(:,3), TrajectoriesMatrix.House3Cont(:,4),'b','linewidth',2)
plot(TrajectoriesMatrix.House3DiscreteTrimmed(:,3), TrajectoriesMatrix.House3DiscreteTrimmed(:,4),'g','linewidth',2)
plot(TrajectoriesMatrix.House3Random(:,3), TrajectoriesMatrix.House3Random(:,4),'c','linewidth',2)
legend('Benchmark','Discrete','Continuous', 'Discrete Trimmed', 'Random', 'location','eastoutside')
xlabel('X')
ylabel('Y')
axis equal tight
box on
title(['five space curves to compare'])
legend    
end

% subplot(2,1,2)
% imagesc([[DistancesMatrix.f11_1, DistancesMatrix.f12_1, DistancesMatrix.f13_1, DistancesMatrix.f14_1, DistancesMatrix.f15_1];[DistancesMatrix.f12_1, DistancesMatrix.f22_1, DistancesMatrix.f23_1, DistancesMatrix.f24_1, DistancesMatrix.f52_1];[DistancesMatrix.f13_1, DistancesMatrix.f23_1, DistancesMatrix.f33_1, DistancesMatrix.f43_1, DistancesMatrix.f53_1]; [DistancesMatrix.f14_1, DistancesMatrix.f24_1, DistancesMatrix.f34_1, DistancesMatrix.f44_1, DistancesMatrix.f54_1]; [DistancesMatrix.f15_1, DistancesMatrix.f52_1, DistancesMatrix.f53_1, DistancesMatrix.f54_1, DistancesMatrix.f55_1]])
% xlabel('curve')
% ylabel('curve')
% cb1=colorbar('peer',gca);
% set(get(cb1,'Ylabel'),'String','Frechet Distance')
% axis equal tight

HouseNumber = HouseNumber-1

subplot(2,1,2)
imagesc([[0, DistancesMatrix(Participant, 1+HouseNumber), DistancesMatrix(Participant, 4+HouseNumber),...
    DistancesMatrix(Participant, 7+HouseNumber), DistancesMatrix(Participant, 13+HouseNumber)];...
    [DistancesMatrix(Participant, 1+HouseNumber), 0, 0, DistancesMatrix(Participant, 10+HouseNumber), DistancesMatrix(Participant, 16+HouseNumber)];...
    [DistancesMatrix(Participant, 4+HouseNumber), 0, 0, 0, DistancesMatrix(Participant, 19+HouseNumber)];...
    [DistancesMatrix(Participant, 7+HouseNumber), DistancesMatrix(Participant, 10+HouseNumber),...
    0, 0, DistancesMatrix(Participant, 22+HouseNumber)]; [DistancesMatrix(Participant, 13+HouseNumber),...
    DistancesMatrix(Participant, 16+HouseNumber), DistancesMatrix(Participant, 19+HouseNumber), DistancesMatrix(Participant, 22+HouseNumber), 0]])

xlabel('curve')
ylabel('curve')
cb1=colorbar('peer',gca);
set(get(cb1,'Ylabel'),'String','Frechet Distance')
axis equal tight

end