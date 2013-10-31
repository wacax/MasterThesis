function struct = mazeStatistics2mat(VPCode)
dire = 'D:\Wacax\TU Berlin Thesis\bbciRaw\';
str = sprintf('%s%s', dire, char(VPCode));
cd(str)

Houses = 3;
Conditions = {'Cont', 'Discrete', 'Joystick'};

for i = 1:Houses
    for ii = 1:length(Conditions)
        condition = sprintf('%s%i%s', 'House', i, Conditions{ii})
        name2import = sprintf('%s%i%s%s', 'House', i, Conditions{ii}, 'Measurements.txt')
        new = importfile2(name2import)
        struct.(condition) = new.data
    end
end

end
