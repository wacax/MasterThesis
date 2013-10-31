function [Xrep Yrep] = FindStatic(file)

Xrep = [];
for i = 1:size(file, 1)
    if i == size(file, 1)
        break
    end
    if file(i, 3) == file(i+1,3)
        Xrep = [Xrep  i];
    end
end

Yrep = [];
for i = 1:size(file, 1)
    if i == size(file, 1)
        break
    end
    if file(i, 4) == file(i+1, 4)
        Yrep = [Yrep  i];
    end
end