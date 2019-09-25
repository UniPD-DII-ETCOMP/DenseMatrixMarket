function [xmin,xmax,ymin,ymax,zmin,zmax] = fun_domain_limits(NN)
%%
    xmin=min(NN(1,:));
    xmax=max(NN(1,:));
    if xmin==xmax
        xmin=xmin-0.1;
        xmax=xmax+0.1;
    end
    if xmin<0
       xmin=xmin*1.1;
    else
        xmin=xmin*0.9;
    end
    if xmax>0
       xmax=xmax*1.1;
    else
        xmax=xmax*0.9;
    end

    ymin=min(NN(2,:));
    ymax=max(NN(2,:));
    if ymin==ymax
        ymin=ymin-0.1;
        ymax=ymax+0.1;
    end
    if ymin<0
       ymin=ymin*1.1;
    else
        ymin=ymin*0.9;
    end
    if ymax>0
       ymax=ymax*1.1;
    else
        ymax=ymax*0.9;
    end

    zmin=min(NN(3,:));
    zmax=max(NN(3,:));
    if zmin==zmax
        zmin=zmin-0.1;
        zmax=zmax+0.1;
    end
    if zmin<0
       zmin=zmin*1.1;
    else
        zmin=zmin*0.9;
    end
    if zmax>0
       zmax=zmax*1.1;
    else
        zmax=zmax*0.9;
    end
%%
end

