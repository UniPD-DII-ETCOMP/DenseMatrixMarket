function tree=quicktree(edg)

istree=zeros(size(edg,1),1);

index = 1;
iedgex = edg;
edgnum = size(edg,1);

while index < edgnum
    itest = iedgex(index,2);
    iedgex(index,2) = iedgex(index,1);
    istree(index) = 1;
    for iedge = index+1:edgnum
        if iedgex(iedge,1)==itest
            iedgex(iedge,1)=iedgex(index,1);
        end
        if iedgex(iedge,2)==itest
            iedgex(iedge,2)=iedgex(index,1);
        end
    end
    while iedgex(index+1,1)==iedgex(index+1,2)
        index=index+1;
        if (index==edgnum)
            break
        end
    end
    index=index+1;
end

tree=find(istree);

end

