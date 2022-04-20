function [filename] = chooseDataFile(num)
%% choose a network to load from networks folder
clc
dirName = [pwd,'\real_data\*.mat'];
Allinpnames = dir(dirName);

if isempty(num)
    disp(sprintf('\nChoose real data file:'))
    for i=1:length(Allinpnames)
        disp([num2str(i),'. ', Allinpnames(i).name])
    end
    x = input(sprintf('\nEnter File Number: '));
else
    x = num;
end
filename=['real_data\',Allinpnames(x).name];    
end

