function [x_path, y_path,z_path, intensity] =  rungekuttaNatInter3D(xseed,yseed,zseed, PartArr,...
                             pathDir,~,maxPathLength, ReversePath, step_size, wb)
    
                         
    
    %x and y are seed points
    projectionMultiplier = 2;
    p0 = [xseed; yseed; zseed];


    Vx =@(stress, shearxy,shearxz) [stress; shearxy; shearxz];
    Vy =@(stress, shearxy, shearyz) [shearxy; stress; shearyz];
    Vz =@(stress, shearxz, shearyz) [shearxz; shearyz; stress];
    
    switch lower(pathDir)
        case 'x'
            V = Vx;
        case 'y'
            V = Vy;
        case 'z'
            V = Vz;
    end
    
    [in, Element] = globalFind(PartArr, PartArr(1).elements(1), p0);
    
    if in
        [F, Fs1, Fs2] = setInterpFunc(Element, pathDir);
    else
        fprintf('Seed Point (%f, %f, %f) not in solution domain\n', xseed, yseed, zseed);
        x_path = [];
        y_path = [];
        z_path =[];
        intensity = [];
        return
    end
    
    p = NaN(3,maxPathLength,'double');
    intensity = NaN(1,maxPathLength,'double');
    w = 1; 
    element_change = false;
    
    if ReversePath
        F=@(x,y, z) -F(x,y, z);
        Fs1=@(x,y,z) -Fs1(x,y,z);
        Fs2=@(x,y,z) -Fs2(x,y,z);
    end

    while w <= maxPathLength && in ~= false
        %Terminate program if cancel button is pressed
        if getappdata(wb,'canceling')
            delete(wb)
            return
        end
        %Interpolate stress initially, and get the relative normalisation
        %value.
        if ~element_change
            p(:,w) = p0;
            stress = F(p0(1), p0(2), p0(3));
            shear1 = Fs1(p0(1), p0(2), p0(3));
            shear2 = Fs2(p0(1), p0(2), p0(3));
            intensity(w) = norm([stress shear1 shear2]);
        else
            [F, Fs1, Fs2] = setInterpFunc(Element,pathDir);
            if ReversePath
                F=@(x,y, z) -F(x,y, z);
                Fs1=@(x,y,z) -Fs1(x,y,z);
                Fs2=@(x,y,z) -Fs2(x,y,z);
            end
            p(:,w) = p0;
            stress = F(p0(1,1), p0(2,1), p0(3,1));
            shear1 = Fs1(p0(1,1), p0(2,1), p0(3,1));
            shear2 = Fs2(p0(1,1), p0(2,1), p0(3,1));
            intensity(w) = norm([stress shear1 shear2]);
        end
        w=w+1;
        
        
        %Find dp1 at first interpolation.
        %Normalise the poining vector relative to the initial test point.
        %Calculate new points:
        
        
        dp1 = V(stress, shear1, shear2)*step_size/intensity(w-1);
        p1 = p0 + dp1;
        
        stress = F(p1(1), p1(2), p1(3));
        shear1 = Fs1(p1(1), p1(2), p1(3));
        shear2 = Fs2(p1(1), p1(2), p1(3));
        dp2 =  V(stress, shear1, shear2)*step_size/intensity(w-1);
        p2 = p0 + 0.5*dp2;

        stress = F(p2(1), p2(2), p2(3));
        shear1 = Fs1(p1(1), p1(2), p1(3));
        shear2 = Fs2(p1(1), p1(2), p1(3));
        dp3 =  V(stress, shear1, shear2)*step_size/intensity(w-1);
        p3 = p0 + 0.5*dp2;

        stress = F(p3(1), p3(2), p3(3));
        shear1 = Fs1(p1(1), p1(2), p1(3));
        shear2 = Fs2(p1(1), p1(2), p1(3));
        dp4 =  V(stress, shear1, shear2)*step_size/intensity(w-1);

        p0 =p0 + 1/6 * (dp1 + 2*dp2 + 2*dp3 +dp4);
        
        [in, new_Element] = localFind(p0',Element);
        if new_Element.ElementNo ~= Element.ElementNo
            element_change = true;
        else
            element_change = false;
        end
        
        Element = new_Element;
            
        if ~in
            extension = 1;
            while ~in && extension < projectionMultiplier+1
                R = (p0 - p(:,w-1))*extension*2 + p0;
                [in, return_Element] = globalFind(PartArr,Element, R);
                extension = extension+1;
            end
            if in
                p0 = R;
                Element = return_Element;
            end
        end

    end
    nancols = ~isnan(p(1,:));
    x_path = p(1,nancols);    
    y_path = p(2,nancols);
    z_path = p(3,nancols);
    intensity = intensity(nancols);
end
function [F, Fs1, Fs2] = setInterpFunc(Element, pathDir)
    
    surr_elements = surrrounding_elements(Element, Element);
    nodes = unique([surr_elements(:).nodes]);
    coordx = [nodes(:).xCoordinate]';
    coordy = [nodes(:).yCoordinate]';
    coordz = [nodes(:).zCoordinate]';
    

    if strcmpi(pathDir,'x')
        F = scatteredInterpolant(coordx(:), coordy(:), coordz(:), [nodes(:).xStress]', 'natural');
        Fs1 =  scatteredInterpolant(coordx, coordy, coordz, [nodes(:).xyStress]', 'natural');
        Fs2 = scatteredInterpolant(coordx, coordy, coordz, [nodes(:).xzStress]', 'natural');
    elseif strcmpi(pathDir,'y') 
        F = scatteredInterpolant(coordx(:), coordy(:), coordz(:), [nodes(:).yStress]', 'natural');
        Fs1 =  scatteredInterpolant(coordx, coordy, coordz, [nodes(:).xyStress]', 'natural');
        Fs2 = scatteredInterpolant(coordx, coordy, coordz, [nodes(:).yzStress]', 'natural');
    elseif strcmpi(pathDir,'z')
        F = scatteredInterpolant(coordx(:), coordy(:), coordz(:), [nodes(:).zStress]', 'natural');
        Fs1 = scatteredInterpolant(coordx, coordy, coordz, [nodes(:).yzStress]', 'natural');
        Fs2 = scatteredInterpolant(coordx, coordy, coordz, [nodes(:).xzStress]', 'natural');
    end
end