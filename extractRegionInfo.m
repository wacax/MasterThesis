function CollapsedLengthTimeandSpeed = extractRegionInfo(statFiles,...
    Region, Condition, dataType)

%First column is path length, second is time and the third is speed
houses = 3
participants = 7
CollapsedLengthTimeandSpeed = [];

for ii = 1:participants
    for iii = 1:houses
        intString1 = sprintf('%s%i%s', 'House', iii, char(Condition))
        structure = statFiles{ii}.(intString1)
        CollapsedLengthTimeandSpeed = [CollapsedLengthTimeandSpeed
            structure(Region, dataType)]
    end
end
end