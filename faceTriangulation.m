function [A, E1, E2, ECross] = faceTriangulation(verticies)
    X_vert = [verticies(:).xCoordinate]';
    Y_vert = [verticies(:).yCoordinate]';
    Z_vert = [verticies(:).zCoordinate]';
    vert = [X_vert,Y_vert,Z_vert];

    plane = optimalPlane(vert);
    
    edge_midpoints = [vert(end,:);vert((1:end-1),:)];
    points = [(0.5*(edge_midpoints - vert) + vert);vert];
    
    x = points(:,1);
    y = points(:,2);
    z = points(:,3);
    
    switch plane
        case 1
            tri = delaunayTriangulation(x,y);
            ti = tri.ConnectivityList;
            pt = tri.Points;

            A = [pt(ti(:,1),:), z(ti(:,1))];
            B = [pt(ti(:,2),:), z(ti(:,2))];
            C = [pt(ti(:,3),:), z(ti(:,3))];
        case 2
            tri = delaunayTriangulation(x,z);
            ti = tri.ConnectivityList;
            pt = tri.Points;      
            
            A = [pt(ti(:,1),1), y(ti(:,1)), pt(ti(:,1),2)];
            B = [pt(ti(:,2),1), y(ti(:,2)), pt(ti(:,2),2)];
            C = [pt(ti(:,3),1), y(ti(:,3)), pt(ti(:,3),2)];
        case 3
            tri = delaunayTriangulation(y,z);
            ti = tri.ConnectivityList;
            pt = tri.Points;

            A = [x(ti(:,1)), pt(ti(:,1),:)];
            B = [x(ti(:,2)), pt(ti(:,2),:)];
            C = [x(ti(:,3)), pt(ti(:,3),:)];
    end
    
    E1 = C - A;
    E2 = B - A;
    ECross = cross(E2,E1,2);
%     figure (3)
%     triplot(tri)
%     hold on
%     
%     vxlabels = arrayfun(@(n) {sprintf('P%d', n)}, (1:length(pt(:,1)))');
%     Hpl = text(pt(:,1)', pt(:,2)', vxlabels, 'FontWeight', 'bold', 'HorizontalAlignment',...
%        'center', 'BackgroundColor', 'none');
%     ic = incenter(tri);
%     numtri = size(tri,1);
%     trilabels = arrayfun(@(x) {sprintf('T%d', x)}, (1:numtri)');
%     Htl = text(ic(:,1), ic(:,2), trilabels, 'FontWeight', 'bold', ...
%        'HorizontalAlignment', 'center', 'Color', 'blue');
%     hold off
%     If the cross product close to zero, then the points a colinear
    buffer = 1e5;
    a = ~all((ECross < eps*buffer & ECross > -eps*buffer), 2);
    ECross = ECross(a,:);
    E1 = E1(a,:);
    E2 = E2(a,:);
    A = A(a,:);

    
end
    function [plane] = optimalPlane(vert)
        vert = [vert;vert(1,:)];
        sq_dist = diff(vert).^2;
        [~, plane] = max(sum([sum(sq_dist(:,[1,2]),2), sum(sq_dist(:,[1,3]),2), sum(sq_dist(:,[2,3]),2)].^(1/2),1));
    end
        
        