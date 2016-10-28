function [indicator] = kriging_EIdom(y, mse, ymin)

% % Note that it takes as many x points as objectives
% if (size(x,1)~=2)
%     error('EIdom is only defined for 2 objectives and 2 points need to be passed, one per each objective');
% end

% y = [];
% mse = [];
% for obj = 1:size(x,1)
%     [y_aux, mse_aux] = surrogate.predictor(x(obj,:), surrogate.model);
%     y(1,obj) = y_aux(1,obj);
%     mse(1,obj) = mse_aux(1,obj);
% end

s = sqrt(abs(mse));

    M0 = size(ymin,1);
    y1 = repmat(y,M0,1);
    s1 = repmat(s,M0,1);
    u = (ymin-y1)./s1;
    PHI = 0.5*erfc(-u./sqrt(2));
    phi = exp(-0.5*u.^2)./(sqrt(2*pi).*s1);
    PHI1 = PHI(1,1);
    PHIs = 0;
    E1 = y(1)*PHI(1,1) - s(1)*phi(1,1);
    E2 = y(2)*PHI(1,2) + s(2)*phi(1,2);
    E1s = 0;
    E2s = 0;
    for i = 1:M0-1
        
        PHIi = PHI(i,1);
        PHIii = PHI(i+1,1);
        PHI2ii = PHI(i+1,2);
        PHIs = PHIs + (PHIii - PHIi)*PHI2ii;
        
        E1i = y(1)*PHI(i,1) - s(1)*phi(i,1);
        E1ii = y(1)*PHI(i+1,1) - s(1)*phi(i+1,1);
        E2i = y(2)*PHI(i,2) - s(2)*phi(i,2);
        E2ii = y(2)*PHI(i+1,2) - s(2)*phi(i+1,2);
        E1s = E1s + (E1ii - E1i)*PHI2ii;
        E2s = E2s + (E2i - E2ii)*PHIi;
        
    end
    PHI1M0 = PHI(end,1);
    PHI2M0 = PHI(end,2);
    P = PHI1 + PHIs + (1-PHI1M0)*PHI2M0;
    E1M0 = y(1)*PHI(end,1) + s(1)*phi(end,1);
    E2M0 = y(2)*PHI(end,2) - s(2)*phi(end,2);
    yE(1) = (E1 + E1s + E1M0*PHI2M0)/P;
    yE(2) = (E2*PHI1 + E2s + E2M0)/P;
    dist = zeros(1,M0);
    for i = 1:M0
        dist(i) = norm(yE - ymin(i,:));
    end
    [~,ind] = min(dist);
    fE = ymin(ind,:);
    E = P*norm(yE - fE);
    
    indicator = E;
    indicator (P==0) = 0;
    
return
