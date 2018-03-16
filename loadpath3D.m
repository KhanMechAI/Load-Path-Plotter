function loadpath3D(dim, sim_dir, seed_dir, save_dir, model_name,path_dir,...
                    overlay, parallel, newPDF, recompute, step_size, path_length)    
%% ********************  Seed Coordinates  ******************************
    F = findall(0,'type','figure','tag','TMWWaitbar'); %%%Unsupress this
%     if wait bar wont close
    delete(F);
    Seed = importdata(seed_dir, ',');
    xseed = Seed(:,1);
    yseed = Seed(:,2);
    zseed = Seed(:,3);
    numSeeds = length(xseed);
%% HouseKeeping
    
    wb = waitbar(0,'1','Name','Computing Load Paths',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
    setappdata(wb,'canceling',0)
    data_read_time = 7;
    plot_time = 3;
    print_time = 10;
    path_time = 80;
    total_time = data_read_time + plot_time + print_time + path_time;
    current_time = 0;
    warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
    warning('off','MATLAB:MKDIR:DirectoryExists')

    %% Scan Data
    
    model_name = [model_name, ' - ',upper(path_dir), ' Path'];
    
    model_data_name = regexprep(model_name, ' ', '_');
    
    system(['taskkill /fi "WINDOWTITLE eq ', model_name,'.pdf"']);

    nodei = [sim_dir '\' 'nodeInfo.txt'];

    numNodes = importdata(nodei);
    numNodes = numNodes(2);

    
    %% ******************  Populate Nodes and Elements  *********************
    
    if ~exist([save_dir, '\Path Data\data_',model_data_name,'.mat'], 'file') || recompute
        fprintf('New model or user nominated to recompute data. Starting now.\n')
        waitbar(current_time/total_time,wb,sprintf('Computing Initial Data'))
        
        [StressData, numNodes] = nodeDat(sim_dir, numNodes);
        current_time = current_time + data_read_time/3;
        waitbar(current_time/total_time,wb,sprintf('Computing Initial Data'))
        fprintf('Nodal information complete. Starting stress population.\n')
        
        [nodes] = NodeDatRead(sim_dir, StressData, numNodes);
        current_time = current_time + data_read_time/3;
        waitbar(current_time/total_time,wb,sprintf('Computing Initial Data'))
        fprintf('Nodal stresses populated. Element generation beginning.\n')
        
        [MeshNode, nodePerEl, PartArr] = datread(sim_dir, nodes); 
        current_time = current_time + data_read_time/3;
        
        fprintf('Elements constructed, directories being created and data being saved.\n')
        mkdir(save_dir, '\Path Data')
        save([save_dir,'\Path Data\data_',model_data_name,'.mat'],'PartArr','nodes', 'nodePerEl', 'MeshNode');
        
    else
        fprintf('Previous model detected, loading data.\n')
        waitbar(current_time/total_time,wb,sprintf('Loading Data'))
        load([save_dir, '\Path Data\data_',model_data_name,'.mat']);
        current_time = current_time + data_read_time;
        
        fprintf('Data loaded. Starting path computation.\n')
        
    end
    if getappdata(wb,'canceling')
        delete(wb)
        return
    end
    
    waitbar(current_time/total_time,wb,sprintf('Starting Paths'))
    %% ****************  Load Path Generation  ******************************
    
    
    
    Paths(numSeeds).X.forward = [];
    Paths(numSeeds).Y.forward = [];
    Paths(numSeeds).Z.forward = [];
    Paths(numSeeds).I.forward = [];
    Paths(numSeeds).X.total = [];
    Paths(numSeeds).Y.total = [];
    Paths(numSeeds).Z.total = [];
    Paths(numSeeds).I.total = [];
        
    switch parallel
        case 1
            workers = 4;
            switch dim
                case {'3D'}
                    parfor (i = 1:numSeeds, workers)
                        fprintf('Starting path %i\n',i)
                        warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
                        [x, y, z, intense] = rungekuttaNatInter3D(xseed(i),...
                            yseed(i),zseed(i), PartArr, path_dir, nodePerEl,path_length,false,step_size,wb);
                        if isempty(x)
                            fprintf('Path %i unsuccessful\n',i)
                            continue
                        end
                        Paths(i).X.forward = x;
                        Paths(i).Y.forward = y;
                        Paths(i).Z.forward = z;
                        Paths(i).I.forward = intense;

                        [x, y, z, intense ] = rungekuttaNatInter3D(xseed(i), yseed(i),...
                            zseed(i), PartArr, path_dir, nodePerEl,path_length, true,step_size, wb);
                        Paths(i).X.total = [fliplr(x), Paths(i).X.forward];
                        Paths(i).Y.total = [fliplr(y), Paths(i).Y.forward];
                        Paths(i).Z.total = [fliplr(z), Paths(i).Z.forward];
                        Paths(i).I.total = [fliplr(intense), Paths(i).I.forward];
                        fprintf('Path %i done\n',i)
                    end
                    current_time = current_time +80;   
                case {'2D'}
                    parfor (i = 1:numSeeds, workers)
                        fprintf('Starting path %i\n',i)
                        warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
                        [x, y, intense] = rungekuttaNatInter2D(xseed(i),...
                            yseed(i), PartArr, path_dir, nodePerEl,path_length,false,step_size,wb);
                        if isempty(x)
                            fprintf('Path %i unsuccessful\n',i)
                            continue
                        end
                        Paths(i).X.forward = x;
                        Paths(i).Y.forward = y;

                        Paths(i).I.forward = intense;

                        [x, y, intense ] = rungekuttaNatInter2D(xseed(i), yseed(i),...
                             PartArr, path_dir, nodePerEl,path_length, true,step_size, wb);
                        Paths(i).X.total = [fliplr(x), Paths(i).X.forward];
                        Paths(i).Y.total = [fliplr(y), Paths(i).Y.forward];

                        Paths(i).I.total = [fliplr(intense), Paths(i).I.forward];
                        fprintf('Path %i done\n',i)
                    end
                    current_time = current_time +80;  
            end
        case 0
            switch dim
                case {'3D'}

                    for i = 1:numSeeds
                        fprintf('Starting path %i\n',i)
                        if getappdata(wb,'canceling')
                            delete(wb)
                            return
                        end

                        waitbar(current_time/total_time,wb,sprintf('Seed %i of %i Computing', i, numSeeds))
                        warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');

                        [x, y, z, intense] = rungekuttaNatInter3D(xseed(i),...
                            yseed(i),zseed(i), PartArr, path_dir, nodePerEl,path_length,false,step_size, wb);
                        if isempty(x)
                            fprintf('Path %i unsuccessful\n',i)
                            continue
                        end
                        Paths(i).X.forward = x;
                        Paths(i).Y.forward = y;
                        Paths(i).Z.forward = z;
                        Paths(i).I.forward = intense;
                        current_time = current_time + 1/numSeeds *80/2;

                        [x, y, z, intense ] = rungekuttaNatInter3D(xseed(i), yseed(i),...
                            zseed(i), PartArr, path_dir, nodePerEl,path_length,true,step_size, wb);

                        Paths(i).X.total = [fliplr(x), Paths(i).X.forward];
                        Paths(i).Y.total = [fliplr(y), Paths(i).Y.forward];
                        Paths(i).Z.total = [fliplr(z), Paths(i).Z.forward];
                        Paths(i).I.total = [fliplr(intense), Paths(i).I.forward];
                        current_time = current_time + 1/numSeeds *80/2;
                        fprintf('Path %i done\n',i)
                    end
                case {'2D'}
                    
                    for i = 1:numSeeds
                        fprintf('Starting path %i\n',i)
                        if getappdata(wb,'canceling')
                            delete(wb)
                            return
                        end

                        waitbar(current_time/total_time,wb,sprintf('Seed %i of %i Computing', i, numSeeds))
                        warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');

                        [x, y, intense] = rungekuttaNatInter2D(xseed(i),...
                            yseed(i), PartArr, path_dir, nodePerEl,path_length,false,step_size, wb);
                        if isempty(x)
                            fprintf('Path %i unsuccessful\n',i)
                            continue
                        end
                        Paths(i).X.forward = x;
                        Paths(i).Y.forward = y;

                        Paths(i).I.forward = intense;
                        current_time = current_time + 1/numSeeds *80/2;

                        [x, y, intense ] = rungekuttaNatInter2D(xseed(i), yseed(i),...
                            PartArr, path_dir, nodePerEl,path_length,true,step_size, wb);

                        Paths(i).X.total = [fliplr(x), Paths(i).X.forward];
                        Paths(i).Y.total = [fliplr(y), Paths(i).Y.forward];

                        Paths(i).I.total = [fliplr(intense), Paths(i).I.forward];
                        current_time = current_time + 1/numSeeds *80/2;
                        fprintf('Path %i done\n',i)
                    end
            end
    end        
    fprintf('All seeds tested and appropriate paths computed. Saving path data.\n')
    %% ****************  Plotting and Printing  ******************************
    if getappdata(wb,'canceling')
        delete(wb)
        return
    end
    waitbar(current_time/total_time,wb,sprintf('Paths Finished Paths, Saving and Plotting now...\n'))

    save([save_dir,'\Path Data\pathdata_',model_data_name,'.mat'], 'Paths');
    fig = figure;
    fprintf('Plotting Paths\n')
    
    switch dim
        case {'3D'}
            modelPlot3D([Paths(:).X],[Paths(:).Y],[Paths(:).Z],[Paths(:).I],PartArr)
        case {'2D'}
            modelPlot2D([Paths(:).X],[Paths(:).Y],[Paths(:).I],PartArr)
    end
    
    
    if getappdata(wb,'canceling')
        delete(wb)
        return
    end
    fprintf('Printing to PDF\n')
    current_time = current_time + plot_time;
    waitbar(current_time/total_time,wb,sprintf('Printing PDF'))
    
    mkdir(save_dir,'\Path Plots')
    if newPDF
        dt = datestr(now,'HH.MM.SS_dd/mm/yy');
        dateAppenedFN = [save_dir,'\Path Plots\',model_name,'_', dt, '.pdf'];
    else
        dateAppenedFN = [save_dir,'\Path Plots\',model_name, '.pdf'];
    end
    

    switch overlay
        case 1
            gcf        
            print(fig,dateAppenedFN, '-dpdf','-r2000', '-fillpage');
            open(dateAppenedFN);
        case 0
            hold off
    end
    delete(wb)
    fprintf('Load paths complete.\n')
end
function [] = modelPlot3D(x_paths,y_paths,z_paths,Intensity,PartArr)
    Alpha = 0.1;
    Buffer = 0.35;
    wireFrame(PartArr,Alpha, Buffer);
    seedLength = size(x_paths(:),1);
    maxInt = max(max([Intensity.total]));
    minInt = min(min([Intensity.total]));
    if isempty(maxInt)
        disp('No Successful Paths')
        return
    end
    for k = 1:seedLength
        if isempty(x_paths(k).total)
            fprintf('Path %i unsuccessful', k)
            continue
        end
        cd = colormap('parula');
        cd = interp1(linspace(minInt,maxInt,length(cd)),cd,Intensity(k).total);
        cd = uint8(cd'*255);
        cd(4,:) = 255;
        paths = line(x_paths(k).total,y_paths(k).total,z_paths(k).total);
        set(paths.Edge,'ColorBinding','interpolated', 'ColorData',cd)
    end
    colorbar;
end
function [] = modelPlot2D(x_paths,y_paths,Intensity,PartArr)
    Alpha = 0.1;
    Buffer = 0.35;
    wireFrame2D(PartArr,Alpha, Buffer);
    seedLength = size(x_paths(:),1);
    maxInt = max(max([Intensity.total]));
    minInt = min(min([Intensity.total]));
    if isempty(maxInt)
        disp('No Successful Paths')
        return
    end
    for k = 1:seedLength
        if isempty(x_paths(k).total)
            fprintf('Path %i unsuccessful', k)
            continue
        end
        cd = colormap('parula');
        cd = interp1(linspace(minInt,maxInt,length(cd)),cd,Intensity(k).total);
        cd = uint8(cd'*255);
        cd(4,:) = 255;
        paths = line(x_paths(k).total,y_paths(k).total);
        set(paths.Edge,'ColorBinding','interpolated', 'ColorData',cd)
    end
    colorbar;
end