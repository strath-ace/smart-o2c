function [a,v] = alpha_clip2(x,a,v,lb,hb)

% Clips alpha to the max admissible value, so that x+alpha*v is within lb and hb
%
% inputs
%	x	:	vector of coordinates of starting point
%	a	:	initial value of alpha
%	v	:	vector of components of velocity vector
%	lb	:	vector of min admissible value for each coordinate
%	hb	:	vector of max admissible value for each coordinate
% outputs
%   a   :   clipped alpha
%   v   :   clipped vector, if case be

% prevents stagnation if point is on the border and velocity pushes it out.
v = v.*~((v<0).*(x==lb)+(v>0).*(x==hb));                              				% this reads: where x==lb or x==hb, v shall be 0 wrong!!! v shall be 0 if it pushes x out of the domain

% alpha clipping
%lbc = ((lb-x)./v).*(((lb-x)./v)>=0)+(((lb-x)./v)<0);                    % for each element, compute the alpha that would be needed to get to lb: keep positive alphas and set 1 negative ones
%a2 = lbc.*(lbc<=1)+(lbc>1);                                             % clip positive alphas to 1
%hbc = ((hb-x)./v).*(((hb-x)./v)>=0)+(((hb-x)./v)<0);                    % as above, with hb
%a3 = hbc.*(hbc<=1)+(hbc>1);                     						% clip as above

%a = min([a a2 a3]);                                                     % final val

c1 = (x+a*v)<lb;     % lower bound overflow
c2 = (x+a*v)>hb;     % upper bound overflow

if any(c1~=0)
    
    % if we're here we have a violation of the lower bound
    
    id = 1:length(x);
    id = id(c1~=0);
    
    %keyboard
    
    for i=id
       
        atemp = (lb(i)-x(i))./v(i);
        a = min([a;atemp]);
        
    end    
end

if any(c2~=0)
    
    % if we're here we have a violation of the lower bound
    
    id = 1:length(x);
    id = id(c2~=0);
    
    %keyboard
    
    for i=id
       
        atemp = (hb(i)-x(i))./v(i);
        a = min([a;atemp]);
        
    end
    
end

end



