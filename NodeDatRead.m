function [nodes] = NodeDatRead(fname, StressData, numNodes)
    clear nodes
    if ~exist('fname', 'var')
        fname = 'C:\Users\z5020362\Desktop\Uni\OneDrive_2_9-6-2017\test_files\dp0\SYS\MECH\ds.dat';
    else
        fname = [fname  '\ds.dat'];
    end
    
    nodeNums = StressData(1,:);
    xstress = StressData(2,:);
    ystress = StressData(3,:);
    zstress = StressData(4,:);
    xystress = StressData(5,:);
    yzstress = StressData(6,:);
    xzstress = StressData(7,:);
    
    
    datafile = fopen(fname);

    trashdata = 'a';
    startelements = '/com,*********** Nodes ';

    while ~strncmpi(trashdata,startelements, 23)
        trashdata = fgetl(datafile);
    end
    
    nid = 'a';
    
    while ~strncmpi(nid(1), '(', 1)
        nid = strtrim(fgetl(datafile));
    end
    
    nodes(1,numNodes) = Node();
    
    linetest = fgetl(datafile);
    linetest = strsplit(linetest);
    linetest = str2double(linetest);
    
    i = 1;
    while linetest(1) ~= -1
        if i == linetest(2)
            
            nodes(i) = Node(linetest(2), linetest(3), linetest(4),linetest(5),...
                       xstress(i), ystress(i),zstress(i),xystress(i),...
                       yzstress(i),xzstress(i));
            linetest = fgetl(datafile);
            linetest = strsplit(linetest);
            linetest = str2double(linetest);
        else
            nodes(i) = Node(i, NaN, NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN);
        end
    
    i = i+1;
    end
    fclose(datafile); 
end

