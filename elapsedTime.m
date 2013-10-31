function timeElapsedDiscTrim  = elapsedTime(inputTime)

%find out elapsed time in the discrete condition

indicesBreak = []

for iii = 1:length(inputTime)
    if iii == size(inputTime, 1)
        break
    end
    if inputTime(iii+1) - inputTime(iii) > 300
        indicesBreak = [indicesBreak  iii];
    end
end

differenceTotal = 0

for b = 1:length(indicesBreak)
   
    try
    dif = inputTime(indicesBreak(b + 1)) - inputTime(indicesBreak(b))
    differenceTotal = differenceTotal + dif
    end

end

timeElapsedDiscTrim =  max(inputTime) - differenceTotal
end