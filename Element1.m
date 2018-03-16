classdef Element1 < FEM 
    
    properties (Constant, Access = private)
        gaussWeightW = [.3478546451 .3478546451 .6521451549 .6521451549]
        gaussPointsS = [.8611363116 -.8611363116 .3399810436 -.3399810436]
        gaussCoords = [0.861136311600000,0.861136311600000;0.861136311600000,-0.861136311600000;0.861136311600000,0.339981043600000;0.861136311600000,-0.339981043600000;-0.861136311600000,0.861136311600000;-0.861136311600000,-0.861136311600000;-0.861136311600000,0.339981043600000;-0.861136311600000,-0.339981043600000;0.339981043600000,0.861136311600000;0.339981043600000,-0.861136311600000;0.339981043600000,0.339981043600000;0.339981043600000,-0.339981043600000;-0.339981043600000,0.861136311600000;-0.339981043600000,-0.861136311600000;-0.339981043600000,0.339981043600000;-0.339981043600000,-0.339981043600000]
    end
    
    properties (SetAccess = public, GetAccess = public)
        ElementNo = 0;
        nodes;
        importeddata = 0;
        nodenums;
        Faces;
        sphere_radius;
        part_num;
        part_idx;
        verticies;
        
        %Likely Redundant below
        topFace;
        bottomFace;
        rightFace;
        leftFace;
        frontFace;
        backFace;
        
        
        
        topFaceCross;
        bottomFaceCross;
        rightFaceCross;
        leftFaceCross;
        frontFaceCross;
        backFaceCross;
        FaceCross;
        
        topFacePt;
        bottomFacePt;
        rightFacePt;
        leftFacePt;
        frontFacePt;
        backFacePt;
        FacePoint;
        
    end
    
    properties
        elementtype = '4node';
    end
    
    properties (Access = private)
        xi;
        eta;
        D;
        xmap;
        ymap;
    end
    
    properties (Dependent=true)
        del_Ni_xi;
        del_Ni_eta;
        connectivity;        
        assembly;
    end
    
    properties (Dependent = true)
        Jacobian;
        Bmatrix;
        Stiffnessmatrix;
        NodeDisp;
        eStress;
    end
    
    
    methods
        %%Initialise Element%%
        function obj = Element1(ElementNumber, xNodes, Import)

            if ~exist('ElementNumber', 'var')
                ElementNumber = NaN;
            end
            if ~exist('xNodes', 'var')
%                 xNodes = Node();
            end
            if ~exist('Import', 'var')
                obj.importeddata = 0;
            else
                obj.importeddata = 1;                
            end
            obj.ElementNo = ElementNumber;
            
            if exist('Import', 'var')
                                                              
                obj.nodes = xNodes;
                obj.nodenums = [xNodes(:).NodeNum];
                obj.verticies = [xNodes(:).Coordinates]';
                
                obj.sphere_radius = max(pdist([[xNodes(:).xCoordinate];...
                           [xNodes(:).yCoordinate];...
                           [xNodes(:).zCoordinate]]'))/2;
            end
        end
        
        function obj = set.ElementNo(obj, num)
            obj.ElementNo = num;
        end
        function obj = set.nodes(obj, Nodes)
            obj.nodes = Nodes;
            for k=1:length(Nodes)
                Nodes(k).inelement = obj;
            end
        end
        
        %%Unused assignment at this stage
        function AsgNodes(obj, node)
            if nargin~=0
                if obj.ElementNo ~= 0
                    obj.nodes = node;
                else
                    error('Not a node object')
                end
            else
                error('Invaid input')
            end
        end
        
%         function ret = get.xmap(obj)
%             a = zeros(16,1);
%             for i = 1:16
%                 obj.xi = obj.gaussCoords(i,1);
%                 obj.eta = obj.gaussCoords(i,2);
%                 xcoor = [obj.nodes.xCoordinate];
%                 a(i) = dot(xcoor, obj.Shapefunctions);
%             end
%             ret = a;
%         end
%         
%         function ret = get.ymap(obj)
%             a = zeros(16,1);
%             for i = 1:16
%                 obj.xi = obj.gaussCoords(i,1);
%                 obj.eta = obj.gaussCoords(i,2);
%                 ycoor = [obj.nodes.yCoordinate];
%                 a(i) = dot(ycoor, obj.Shapefunctions);
%             end
%             ret = a;
%         end
        
        function ret = IrregBiLinInterp(obj, x, y)

            %Interpolation data
            x1 = obj.nodes(1).xCoordinate;
            x2 = obj.nodes(2).xCoordinate;
            x21 = x2 - x1;
            x31 = obj.nodes(3).xCoordinate - obj.nodes(1).xCoordinate;
            x42 = obj.nodes(4).xCoordinate - obj.nodes(2).xCoordinate;
            
            y1 = obj.nodes(1).yCoordinate;
            y2 = obj.nodes(2).yCoordinate;
            y21 = obj.nodes(2).yCoordinate - obj.nodes(1).yCoordinate;
            y31 = obj.nodes(3).yCoordinate - obj.nodes(1).yCoordinate;
            y42 = obj.nodes(4).yCoordinate - obj.nodes(2).yCoordinate;
            
            A = x31*y42 - y31*x42;
            
            B = y*(x42 - x31) - x*(y42 - y31) + x31*y2 - y31*x2 + x1*y42...
                -y1*x42;
            C = y * x21 - x*y21 + x1*y2 -x2*y1;
            
            t1 = single((-B + sqrt(B^2 - 4* A * C))/(2*C));
            t2 = single((-B - sqrt(B^2 - 4* A * C))/(2*C));
            
            
            if t1 <= 1 && t1 >= 0
                t = t1;
            elseif t2 <= 1 && t2 >= 0
                t = t2;
            else 
                disp('didntwork')
                return
            end
            t = double(t);
            
            s = (y - y1 - y31*t)/(y2 + y42*t - y1 - y31*t);
            
            Px = @(obj, s,t) obj.nodes(1).xStress *(1-s)*(1-t) + ...
                            obj.nodes(2).xStress * s * (1-t) + ...
                            obj.nodes(3).xStress * (1-s) * t + ...
                            obj.nodes(4).xStress * s * t;   
            Py = @(obj, s,t) obj.nodes(1).yStress *(1-s)*(1-t) + ...
                            obj.nodes(2).yStress * s * (1-t) + ...
                            obj.nodes(3).yStress * (1-s) * t + ...
                            obj.nodes(4).yStress * s * t;
            Pxy = @(obj, s,t) obj.nodes(1).xyStress *(1-s)*(1-t) + ...
                            obj.nodes(2).xyStress * s * (1-t) + ...
                            obj.nodes(3).xyStress * (1-s) * t + ...
                            obj.nodes(4).xyStress * s * t;            
            
            ret = [Px(obj,s,t); Py(obj,s,t);Pxy(obj,s,t)];
            
        end
        
        
%         function InterpMat(obj)
%             row = @(x, y) [1 x y x*y];
%             A = [row(obj.nodes(1).xCoordinate, obj.nodes(1).yCoordinate);...
%                  row(obj.nodes(2).xCoordinate, obj.nodes(2).yCoordinate);...
%                  row(obj.nodes(3).xCoordinate, obj.nodes(3).yCoordinate);...
%                  row(obj.nodes(4).xCoordinate, obj.nodes(4).yCoordinate)];
%             sigma = @(node) [node.xStress node.yStress node.xyStress];
%             sigmat = [sigma(obj.nodes(1)); sigma(obj.nodes(2));...
%                       sigma(obj.nodes(3)); sigma(obj.nodes(4))];
%             obj.InterpMatrix = A\sigmat;
%         end
            
        function sigma = EInterpolate(obj, x,y)
            row = @(x, y) [1 x y x*y];
            sigma = row(x,y)*obj.InterpMatrix;
        end
        
        function ret = Interp(obj,x,y)
            row = @(x, y) [1 x y x*y];
            column =@(x,y) [1 x y x*y]';
            A =@ (x,y) column(x,y)*row(x,y) .*column(x,y);
            b = @(x,y, stress) column(x,y) * [stress(1) stress(2) stress(3)];
            AMat = zeros(4,4);
            bMat = zeros(4,3);
            if ~obj.importeddata
                for i=1:length(obj.xmap)
                    AMat = AMat + A(obj.xmap(i), obj.ymap(i));
                    bMat = bMat + b(obj.xmap(i), obj.ymap(i),...
                                            obj.nodes(i).avStress);
                end
            else
            end
            a = AMat\bMat;
            sigma_xy = row(x,y) * a;
            ret = sigma_xy';
        end

        function ret = Extrap(obj, prop, x, y)
            switch prop
                case 'xStress'
                    ret = [obj.nodes(1).xStress(obj.ElementNo) ...
                           obj.nodes(2).xStress(obj.ElementNo) ...
                           obj.nodes(3).xStress(obj.ElementNo) ...
                           obj.nodes(4).xStress(obj.ElementNo)] ...
                          * obj.Shapefunctions';
                case 'yStress'
                    ret = [obj.nodes(1).yStress(obj.ElementNo) ...
                           obj.nodes(2).yStress(obj.ElementNo) ...
                           obj.nodes(3).yStress(obj.ElementNo) ...
                           obj.nodes(4).yStress(obj.ElementNo)] ...
                          * obj.Shapefunctions(x,y)';
                case 'Stress'
                    ret1 = Extrap(obj, 'xStress', x, y);
                    ret2 = Extrap(obj, 'yStress', x, y);
                    ret = [ret1 ret2];
            end   
        end
        
        function ret = Coords(obj)
            ret = [[obj.nodes.xCoordinate]' [obj.nodes.yCoordinate]'];
        end
        
        %Access node numbers of Node objects assigned to element
        function Nodes(obj)
            ShowNum(obj.nodes)
        end
        
%         function [ret] = get.nodenums(obj)
%             [ret] = [obj.nodenums];
%         end
        
        function ret = get.del_Ni_xi(obj)
            ret = [obj.del_N1_xi(obj.xi, obj.eta) obj.del_N2_xi(obj.xi, obj.eta) ...
                   obj.del_N3_xi(obj.xi, obj.eta) obj.del_N4_xi(obj.xi, obj.eta)]';
        end
        
        function ret = get.del_Ni_eta(obj)
            ret = [obj.del_N1_eta(obj.xi, obj.eta) obj.del_N2_eta(obj.xi, obj.eta) ...
            obj.del_N3_eta(obj.xi, obj.eta) obj.del_N4_eta(obj.xi, obj.eta)]';
        end
        
%         function ret = get.Shapefunctions(obj)
%             ret = [obj.N1(obj.xi, obj.eta) obj.N2(obj.xi, obj.eta) ...
%                    obj.N3(obj.xi, obj.eta) obj.N4(obj.xi, obj.eta)];
%         end
        
        %Jacobian property calculation
        function ret = get.Jacobian(obj) 
            ret = [obj.del_Ni_xi obj.del_Ni_eta]' * Coords(obj);
        end
        
        function ret = get.Bmatrix(obj)
            Nx = [1 0] * obj.Jacobian * [obj.del_Ni_xi'; obj.del_Ni_eta'];
            Ny = [0 1] * obj.Jacobian * [obj.del_Ni_xi'; obj.del_Ni_eta'];
            B = [Nx(1) 0 Nx(2) 0 Nx(3) 0 Nx(4) 0; ...
                 0 Ny(1) 0 Ny(2) 0 Ny(3) 0 Ny(4); ...
                 Ny(1) Nx(1) Ny(2) Nx(2) Ny(3) Nx(3) Ny(4) Nx(4)];
            ret = B;
        end
        
%         function obj = set.D(obj,~)
%             obj.D = obj.E/(1-obj.nu^2).*[1 obj.nu 0 ; obj.nu 1 0; 0 0 (1-obj.nu)/2];
%         end
        
        function ret = get.D(obj)
            obj.D = obj.E/(1-obj.nu^2).*[1 obj.nu 0 ; obj.nu 1 0; 0 0 (1-obj.nu)/2];
            ret = obj.D;
        end
        
        function ret = get.Stiffnessmatrix(obj)
            K = 0;
            for i =1:4
                for j =1:4
                    obj.xi = obj.gaussPointsS(i);
                    obj.eta = obj.gaussPointsS(j);
                    K = K + obj.gaussWeightW(i) * obj.gaussWeightW(j) .* ...
                    (obj.Bmatrix' * obj.D * obj.Bmatrix) .* det(obj.Jacobian);
                end
            end
            ret = K;
        end
        
        function ret = get.assembly(obj)
            assem = @(a) [(2*a(1) -1) 2*a(1) (2*a(2) -1) 2*a(2)  (2*a(3) -1) 2*a(3)  (2*a(4) -1) 2*a(4)];
            ret = assem(obj.nodenums);
        end
        
        function ret = get.NodeDisp(obj)
            ret = [obj.nodes(1).xDisplacement obj.nodes(1).yDisplacement ...
                obj.nodes(2).xDisplacement obj.nodes(2).yDisplacement ...
                obj.nodes(3).xDisplacement obj.nodes(3).yDisplacement ...
                obj.nodes(4).xDisplacement obj.nodes(4).yDisplacement];
        end
        
        function ret = get.eStress(obj)
            ret = obj.D * obj.Bmatrix * obj.NodeDisp';
        end
        
        function nStress(obj)
            for i = 1:4
                 a = obj.D * obj.Bmatrix(:,[(2*i -1) 2*i])...
                * [obj.nodes(i).xDisplacement; obj.nodes(i).yDisplacement];
                obj.nodes(i).tStress(:,end+1) = a;
            end       
        end
        function ElPlot(obj)
            nodePerEl = length([obj.nodes(:)]);
            if nodePerEl == 8
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
                else
                xver =@(seedel)  [seedel.nodes(:).xCoordinate,...
                                            seedel.nodes(1).xCoordinate];
                yver =@(seedel)  [seedel.nodes(:).yCoordinate,...    
                                            seedel.nodes(1).yCoordinate];
            end

            XVert = xver(obj);
            YVert = yver(obj);
            hold on 
            plot(XVert, YVert, 'k')
        end
        function EPlot3(obj)
            x = [obj.nodes(:).xCoordinate];
            y = [obj.nodes(:).yCoordinate];
            z = [obj.nodes(:).zCoordinate];
            plot3(x(:),y(:),z(:), 'o')
%             a = gca;
%             a.DataAspectRatio = [1 1 1];
        end
    end 
end
 