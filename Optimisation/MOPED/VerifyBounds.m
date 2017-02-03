function x=VerifyBounds(x,bounds)
for indi=1:size(x,1)
    for j=1:size(x,2) %   tansig    y=2./(1+exp(-.05.*x))-1
        ppc=rand(1,1); %   tansig    y=2./(1+exp(-.05.*x))-1
        if x(indi,j)<bounds(j,1)
            if ppc<.8
                x(indi,j)=bounds(j,1);
            else
                dt=abs(bounds(j,1)-bounds(j,2));
                dl=abs(x(indi,j)-bounds(j,1));
                sl=2./(1+exp(-.01*dl))-1;
                x(indi,j)=bounds(j,1)+sl*dt;
            end
        end
        if x(indi,j)>bounds(j,2)
            if ppc<.8
                x(indi,j)=bounds(j,2);
            else
                dt=abs(bounds(j,1)-bounds(j,2));
                dl=abs(x(indi,j)-bounds(j,2));
                sl=2./(1+exp(-.01*dl))-1;
                x(indi,j)=bounds(j,2)-sl*dt;
            end
        end
        if isnan(x(indi,j))
            x(indi,j)=random(1,1)*abs(bounds(j,1)-bounds(j,2))+bounds(j,1);
        end
    end
end
