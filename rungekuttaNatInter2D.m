function [x_path, y_path, intensity] =  rungekuttaNatInter2D(xseed,yseed, PartArr,...
                             pathDir,~,maxPathLength, ReversePath, step_size, wb)
    
                         
    
    %x and y are seed points
    projectionMultiplier = 2;
    p0 = [xseed; yseed];


    Vx =@(stress, shearxy) [stress; shearxy];
    Vy =@(stress, shearxy) [shearxy; stress];

    
    switch lower(pathDir)
        case 'x'
            V = Vx;
        case 'y'
            V = Vy;
        case 'z'
            V = Vz;
    end
    
    [in, Element] = globalFind2D(PartArr, PartArr(1).elements(1), p0);
    
    if in
        [F, Fs1] = setInterpFunc(Element, pathDir);
    else
        fprintf('Seed Point (%f, %f, %f) not in solution domain\n', xseed, yseed);
        x_path = [];
        y_path = [];

        intensity = [];
        return
    end
    
    p = NaN(2,maxPathLength,'double');
    intensity = NaN(1,maxPathLength,'double');
    w = 1; 
    element_change = false;
    
    if ReversePath
        F=@(x,y) -F(x,y);
        Fs1=@(x,y) -Fs1(x,y);

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
            stress = F(p0(1), p0(2));
            shear1 = Fs1(p0(1), p0(2));

            intensity(w) = norm([stress shear1]);
        else
            [F, Fs1] = setInterpFunc(Element,pathDir);
            if ReversePath
                F=@(x,y) -F(x,y);
                Fs1=@(x,y) -Fs1(x,y);

            end
            p(:,w) = p0;
            stress = F(p0(1,1), p0(2,1));
            shear1 = Fs1(p0(1,1), p0(2,1));

            intensity(w) = norm([stress shear1]);
        end
        w=w+1;
        
        
        %Find dp1 at first interpolation.
        %Normalise the poining vector relative to the initial test point.
        %Calculate new points:
        
        
        dp1 = V(stress, shear1)*step_size/intensity(w-1);
        p1 = p0 + dp1;
        
        stress = F(p1(1), p1(2));
        shear1 = Fs1(p1(1), p1(2));

        dp2 =  V(stress, shear1)*step_size/intensity(w-1);
        p2 = p0 + 0.5*dp2;

        stress = F(p2(1), p2(2));
        shear1 = Fs1(p1(1), p1(2));

        dp3 =  V(stress, shear1)*step_size/intensity(w-1);
        p3 = p0 + 0.5*dp2;

        stress = F(p3(1), p3(2));
        shear1 = Fs1(p1(1), p1(2));

        dp4 =  V(stress, shear1)*step_size/intensity(w-1);

        p0 =p0 + 1/6 * (dp1 + 2*dp2 + 2*dp3 +dp4);
        
        [in, new_Element] = localFind2D(p0',Element);
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
                [in, return_Element] = globalFind2D(PartArr,Element, R);
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

    intensity = intensity(nancols);
end
function [F, Fs1] = setInterpFunc(Element, pathDir)
    
    surr_elements = surrrounding_elements(Element, Element);
    nodes = unique([surr_elements(:).nodes]);
    coordx = [nodes(:).xCoordinate]';
    coordy = [nodes(:).yCoordinate]';
    

    if strcmpi(pathDir,'x')
        F = scatteredInterpolant(coordx, coordy, [nodes(:).xStress]', 'natural');
        Fs1 =  scatteredInterpolant(coordx, coordy, [nodes(:).xyStress]', 'natural');
    elseif strcmpi(pathDir,'y') 
        F = scatteredInterpolant(coordx, coordy, [nodes(:).yStress]', 'natural');
        Fs1 =  scatteredInterpolant(coordx, coordy, [nodes(:).xyStress]', 'natural');
    end
end