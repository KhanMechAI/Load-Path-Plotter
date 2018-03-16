function [] = wireFrame2D(PartArr,Alpha, Buffer)
    hold on
    warning('off','MATLAB:delaunayTriangulation:DupPtsWarnId');
    mmMat = ones(2,2);
    for k = 1:size(PartArr(:),1)
        tempNodearr = [PartArr(k).elements.nodes];
        x = [tempNodearr.xCoordinate]';
        y = [tempNodearr.yCoordinate]';
        mmMat = newMinMax(x,y,mmMat);

        DT = delaunayTriangulation(x,y);

        wf = triplot(DT);

        wf.Color = [rand(1,3), Alpha];
    end
    a = gca;
    a.DataAspectRatio = [1 1 1];
    a.XLabel.String = 'X';
    a.YLabel.String = 'Y';

    alim =max((mmMat(:,2) - mmMat(:,1)))*Buffer;
    lims =[mmMat(:,1)-alim,mmMat(:,2)+alim];
    a.XLim = lims(1,:);
    a.YLim = lims(2,:);

end

function [mmMat] = newMinMax(x,y,mmMat)
    tempMinMaxVal =@(x,y) [[min(x);min(y)],[max(x);max(y)]];
    tmmV = tempMinMaxVal(x,y);
    a = mmMat(:,1) > tmmV(:,1);
    b = mmMat(:,2) < tmmV(:,2);
    mmMat(a,1) = tmmV(a,1);
    mmMat(b,2) = tmmV(b,2);
end
