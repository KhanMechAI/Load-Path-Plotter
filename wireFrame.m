function [] = wireFrame(PartArr,Alpha, Buffer)
    hold on
    warning('off','MATLAB:delaunayTriangulation:DupPtsWarnId');
    mmMat = ones(3,2);
    for k = 1:size(PartArr(:),1)
        tempNodearr = [PartArr(k).elements.nodes];
        x = [tempNodearr.xCoordinate];
        y = [tempNodearr.yCoordinate];
        z = [tempNodearr.zCoordinate];
        mmMat = newMinMax(x,y,z,mmMat);
        C = [x;y;z];
        DT = delaunayTriangulation(C');
        K = convexHull(DT);
        wf = trimesh(K,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3));
        wf.FaceAlpha = 0;
        wf.EdgeColor = rand(1,3);
        wf.EdgeAlpha = Alpha;
    end
    a = gca;
    a.DataAspectRatio = [1 1 1];
    a.XLabel.String = 'X';
    a.YLabel.String = 'Y';
    a.ZLabel.String = 'Z';
    alim =max((mmMat(:,2) - mmMat(:,1)))*Buffer;
    lims =[mmMat(:,1)-alim,mmMat(:,2)+alim];
    a.XLim = lims(1,:);
    a.YLim = lims(2,:);
    a.ZLim = mmMat(3,:);
end

function [mmMat] = newMinMax(x,y,z,mmMat)
    tempMinMaxVal =@(x,y,z) [[min(x);min(y);min(z)],[max(x);max(y);max(z)]];
    tmmV = tempMinMaxVal(x,y,z);
    a = mmMat(:,1) > tmmV(:,1);
    b = mmMat(:,2) < tmmV(:,2);
    mmMat(a,1) = tmmV(a,1);
    mmMat(b,2) = tmmV(b,2);
end
