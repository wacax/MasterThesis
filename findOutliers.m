function OutlierIndices = findOutliers(data) 

Z = zscore(data);
OutlierIndices = find(abs(Z)>4);

end
