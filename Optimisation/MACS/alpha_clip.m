function a = alpha_clip(x,a,v,lb,hb)

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

	lbc = ((lb-x)./v).*(((lb-x)./v)>=0)+(((lb-x)./v)<0);                    % for each element, compute the alpha that wuold be needed to get to lb: keep positive alphas and set negative ones to 1
	a2 = lbc.*(lbc<=1)+(lbc>1);                                             % clip positive alphas to 1
	hbc = ((hb-x)./v).*(((hb-x)./v)>=0)+(((hb-x)./v)<0);                    % as above, with hb
	a3 = hbc.*(hbc<=1)+(hbc>1);                                             % clip as above
	
	a = min([a a2 a3]);                                                     % final value of alpha

end



