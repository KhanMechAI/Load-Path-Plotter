function [Mesh, numNodes,PartArr] = datread(fpath, nodes)

    Mesh =[];
%     numConnectingPoints = 2;
    
    if ~exist('fpath', 'var')
        fname = 'C:\Users\z5020362\Desktop\Uni\OneDrive_2_9-6-2017\test_files\dp0\SYS\MECH\ds.dat';
    else
        fname = [fpath  '\ds.dat'];
    end
    
    datafile = fopen(fname);

    trashdata = 'a';
    startelements = '/com,*********** Elements';

    while ~strncmpi(trashdata,startelements, 25)
        trashdata = fgetl(datafile);
    end

    elid = 'a';
    numNodes = 0;
    while ~strncmpi(elid, 'eblock', 5)
        elid = strtrim(fgetl(datafile));
        if strncmpi(elid, 'et', 2)
            ElTypeCheck = strsplit(elid, ',');
            type = char(ElTypeCheck(3));
            switch char(ElTypeCheck(3))
                case {'183'}
                    numNodes = 8;
                    skipLine = 0;
                    xver =@(seedel)  [seedel.nodes(1).xCoordinate,...
                                      seedel.nodes(5).xCoordinate,...
                                      seedel.nodes(2).xCoordinate,...
                                      seedel.nodes(6).xCoordinate,...
                                      seedel.nodes(3).xCoordinate,...
                                      seedel.nodes(7).xCoordinate,...
                                      seedel.nodes(4).xCoordinate,...
                                      seedel.nodes(8).xCoordinate,...
                                      seedel.nodes(1).xCoordinate];
                    yver =@(seedel)  [seedel.nodes(1).yCoordinate,...
                                      seedel.nodes(5).yCoordinate,...
                                      seedel.nodes(2).yCoordinate,...
                                      seedel.nodes(6).yCoordinate,...
                                      seedel.nodes(3).yCoordinate,...
                                      seedel.nodes(7).yCoordinate,...
                                      seedel.nodes(4).yCoordinate,...
                                      seedel.nodes(8).yCoordinate,...
                                      seedel.nodes(1).yCoordinate];
                    numEdges = 4;
                case {'182'}
                    numNodes = 4;
                    skipLine = 0;
                    xver =@(seedel)  [seedel.nodes(:).xCoordinate,...
                                        seedel.nodes(1).xCoordinate];
                    yver =@(seedel)  [seedel.nodes(:).yCoordinate,...
                                        seedel.nodes(1).yCoordinate];
                case{'185'}
                    numNodes = 8;
                    skipLine = 0;
                    numEdges = 12;
                case{'186'}
                    numNodes = 8;
                    skipLine = 1;
                    numEdges = 12;
            end
        end
    end




% %     Code relating to building a mesh is now unused. It is left in case
% future developements require full mesh visualisation.
% %     Mesh = nan(numConnectingPoints,numEdges*numel);

    temp = strsplit(elid,',');
    numel = str2double(temp(end));
    
    linetest = fgetl(datafile);
    linetest = fgetl(datafile);
    start = 1;
    numParts = 1;
%     eltype = char(ElTypeCheck(3));
    PartArr(numParts).elements(numel) = Element1();
    PartArr(numParts).range = [];
    PartArr(numParts).span = numel;
    counter=1;
    max_rad = 0;
    faces = faceDefinition(type);
    while ~strcmpi(linetest, '-1')
        
        nums = strsplit(linetest);
        nums = str2double(nums(end-numNodes:end));
        
        if start
            start = false;
            stidx = nums(1);
            PartArr(numParts).range(1) = stidx;
        end
        nodes_nums = nums(2:end);
        element_nodes = nodes(nodes_nums);
        PartArr(numParts).elements(counter) = Element1(nums(1), element_nodes,1);
        Faces.A= [];  
        Faces.E1= [];
        Faces.E2= [];
        Faces.ECross= [];
        for k = 1:size(faces,2)
            face = faces(:,k);
            [~,idx_dup,~] = unique(nodes_nums(face));
            if length(idx_dup) < 3
                continue
            end
            [a,e1,e2,ecross] = faceTriangulation(element_nodes(face(idx_dup)));
            
            try
                Faces.A(:,:,k)= a;
                Faces.E1(:,:,k)= e1;
                Faces.E2(:,:,k)= e2;
                Faces.ECross(:,:,k)= ecross;
                
            catch
            end
            
        end
        PartArr(numParts).elements(counter).Faces = Faces;
        PartArr(numParts).elements(counter).part_num = numParts;
        sphere_influence_tracker = PartArr(numParts).elements(counter).sphere_radius;
        if sphere_influence_tracker >= max_rad
            max_rad = sphere_influence_tracker;
        end
%         Edges = elEdges(elarray(counter),eltype);
%         Mesh(:, (nums(1)-1)*numEdges+1:nums(1)*numEdges) = Edges;
        linetest = fgetl(datafile);
        counter =counter+1;
        
        if skipLine
            linetest = fgetl(datafile);
        end
        
        if strcmpi(linetest, '-1')
            endidx = nums(1);
            PartArr(numParts).range(2) = endidx;
            linetest = fgetl(datafile);
            linetest = fgetl(datafile);
            PartArr(numParts).span = PartArr(numParts).range(2) - PartArr(numParts).range(1);
            if strncmpi(linetest, '/com,*********** Elements for Body',34)
                numParts = numParts+1;
                temp = {};
                while length(temp) <9
                    if strncmpi(linetest, 'et', 2)
                        [skipLine, numNodes,numElements, type] = caseCheck(linetest);
                        faces = faceDefinition(type);
                    end
                    linetest = fgetl(datafile);
                    temp = strsplit(linetest);
                end
                start = true; 
                counter = 1;
            else
                linetest = '-1';
            end            
        end
    end
    PartArr(1).maxRadius = max_rad;
    fclose(datafile); 
end

function [faces] = faceDefinition(element_type)

switch lower(element_type)
    case {'168', '187'}
        I = 1;
        J = 2;
        K = 3;
        L = 4;
        
        faces = [[J;I;K],[I;J;L],[J;K;L],[K;I;L]];
    case {'164','185','186'}
        I = 1;
        J = 2;
        K = 3;
        L = 4;
        M = 5;
        N = 6;
        O = 7;
        P = 8;
        
        
        
        faces = [[J;I;L;K],[I;J;N;M],[J;K;O;N],[K;L;P;O],[L;I;M;P],[M;N;O;P]];
        
        
    case {'182', '181','183'}
        I = 1;
        J = 2;
        K = 3;
        L = 4;
        
        faces = [J;I;L;K];        
end
end

function [Edges] = elEdges(seedel,eltype)

    switch eltype
        case {'183'}
            I = 1;
            J = 3;
            K = 5;
            L = 7;
            E1 = [I; J];
            E2 = [J; K];
            E3 = [K; L];
            E4 = [I; L];

            EdgeIDX = [E1,E2,E3,E4];
        case{'185', '186'}
            I = 1;
            J = 2;
            K = 3;
            L = 4;
            M = 5;
            N = 6;
            O = 7;
            P = 8;

            E1 = [I; J];
            E2 = [J; K];
            E3 = [K; L];
            E4 = [I; L];
            E5 = [M; N];
            E6 = [N; O];
            E7 = [O; P];
            E8 = [M; P];
            E9 = [I; M];
            E10 = [J; N];
            E11 = [K; O];
            E12 = [L; P];

            EdgeIDX = [E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12];
    end
    
    Edges = seedel.nodenums(EdgeIDX);
end
function [skipLine, numNodes, numElements, type] = caseCheck(linetest)
    ElTypeCheck = strsplit(linetest, ',');
    switch char(ElTypeCheck(3))
        case {'183'}
            numNodes = 8;
            skipLine = 0;
        case {'182'}
            numNodes = 4;
            skipLine = 0;
        case{'185'}
            numNodes = 8;
            skipLine = 0;
        case{'186'}
            numNodes = 8;
            skipLine = 1;
        otherwise
            numNodes = 8;
            skipLine = 0;
    end
    numElements = str2double(ElTypeCheck(end));
    type = char(ElTypeCheck(3));
end


