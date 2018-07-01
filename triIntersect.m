function [in] = triIntersect(Faces, O)
%% Triangle Intersection
%Moller Trombore Implementation of ray casting. See their paper.
    in =0;
    buffer = 1e5;
    A = Faces.A;
    ECross = Faces.ECross;
    E1 = Faces.E1;
    E2 = Faces.E2;
    T = O-A;
    
    numFaces = size(A);
    D = ones(numFaces)*-1;
    Det = dot(D,ECross,2);
    
    a = (Det > eps*buffer | Det < -eps*buffer);
    
    if all(a==0)
        return
    end
    
    invDet = 1./Det;
    
    t = dot(ECross,-T,2).*invDet;
    a = a.* (t>=0);
    
    if all(a==0)
        return 
    end
    
    P = cross(D,E2,2);

    u = dot(T,P,2) .* invDet;
    
    a =a.* (u >= 0 & u <= 1);
    
    if all(a==0)
        return 
    end

    Q = cross(T,E1,2);
    
    v = dot(D,Q,2).*invDet;
    
    a = a.* (v >= 0 & u+v <= 1);
    
    if all(a==0)
        return 
    end
    
    hit_count = sum(sum(sum(a)));
    
    if hit_count == 1 || hit_count ==0
        return
    else
        in=1;
        return
    end
%     trisurf(ti,x,y,z)
%     plot3(O(1),O(2),O(3),'*')
%     plot3(D(1),D(2),D(3),'*')
%     X = [O(1) D(1)];
%     Y = [O(2) D(2)];
%     Z = [O(3) D(3)];
%     plot3(X,Y,Z)
    
end
    
    
    

