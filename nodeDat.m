function [StressData, count] = nodeDat(filePath, numNodes)
%% Columns of StressData have the structure:
% StressData = [Node, SX, SY, SZ, SXY, SYZ, SXZ]
    numOfResults = 7;
     if ~exist('filePath', 'var')
        filePath = 'C:\Users\z5020362\Desktop\Uni\OneDrive_2_9-6-2017\test1_files\dp0\SYS\MECH\nodalSolution.txt';
    else
        filePath = [filePath  '\nodalSolution.txt'];
    end
    
    datafile = fopen(filePath);

    trashdata = 'a';
    startelements = '    NODE';

    while ~strncmpi(trashdata,startelements, length(startelements))
        trashdata = fgetl(datafile);
    end
    StressData = nan(numOfResults,numNodes);
    count = 1;
    for i = 1:numNodes
        linetest = strtrim(fgetl(datafile));
        
        if isempty(linetest)
            break
        end
        
        linetest = str2double(strsplit(linetest))';
        
        StressData(:,linetest(1)) = linetest;
        count = count+1;
    end
    fclose(datafile);
end