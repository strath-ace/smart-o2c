% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function update = kriging_update(y, Y, mse, threshold, method, subobj, minmax)

if nargin < 7
    minmax = 0;
end

y = y(subobj);
Y = Y(:,subobj);
mse = mse(subobj);
s = sqrt(mse);
s(s <= 1e-4) = 0;
threshold = threshold(subobj);

switch lower(method)
    
    case 'm'
        % mse (classic method)
        V = s;
        
    case 'p'
        % Probability of improvement
        if minmax == -1
            ybest = -max(Y,[],1);
            y = -y;
        else
            ybest = min(Y,[],1);
        end
        I = ybest - y;
        if any(s == 0)
            P = zeros(size(I));
        else
            u = I./s;
            P = 0.5*erfc(-u./sqrt(2));
        end
        V = P;
        
    case 'e'
        % Expected improvement
        if minmax == -1
            ybest = -max(Y,[],1);
            y = -y;
        else
            ybest = min(Y,[],1);
        end
        I = ybest - y;
        if all(I <= 0) || any(s == 0)
            E = zeros(size(I));
        else
            u = I./s;
            PHI = 0.5*erfc(-u./sqrt(2));
            phi = exp(-0.5*u.^2)./(sqrt(2*pi).*s);
            E = s.*(u.*PHI + phi);
        end
        V = E;
        
    case 'pmoo1'
        % MOO Probability of improvement, front augmentation
        if any(s == 0)
            P = zeros(1,numel(subobj));
        else
            M0 = size(Y,1);
            y1 = repmat(y,M0,1);
            s1 = repmat(s,M0,1);
            u = (Y-y1)./s1;
            PHI = 0.5*erfc(-u./sqrt(2));
            PHI1 = PHI(1,1);
            PHIs = 0;
            for i = 1:M0-1
                PHIi = PHI(i,1);
                PHIii = PHI(i+1,1);
                PHI2i = PHI(i,2);
                PHIs = PHIs + (PHIii - PHIi)*PHI2i;
            end
            PHI1M0 = PHI(end,1);
            PHI2M0 = PHI(end,2);
            P = PHI1 + PHIs + (1-PHI1M0)*PHI2M0;
        end
        V = P(ones(1,numel(subobj)));
        
    case 'pmoo2'
        % MOO Probability of improvement, front domination
        if any(s == 0)
            P = zeros(1,numel(subobj));
        else
            M0 = size(Y,1);
            y1 = repmat(y,M0,1);
            s1 = repmat(s,M0,1);
            u = (Y-y1)./s1;
            PHI = 0.5*erfc(-u./sqrt(2));
            PHI1 = PHI(1,1);
            PHIs = 0;
            for i = 1:M0-1
                PHIi = PHI(i,1);
                PHIii = PHI(i+1,1);
                PHI2ii = PHI(i+1,2);
                PHIs = PHIs + (PHIii - PHIi)*PHI2ii;
            end
            PHI1M0 = PHI(end,1);
            PHI2M0 = PHI(end,2);
            P = PHI1 + PHIs + (1-PHI1M0)*PHI2M0;
        end
        V = P(ones(1,numel(subobj)));
        
    case 'emoo1'
        % MOO Expected improvement, front augmentation
        if any(s == 0)
            E = zeros(1,numel(subobj));
        else
            M0 = size(Y,1);
            y1 = repmat(y,M0,1);
            s1 = repmat(s,M0,1);
            u = (Y-y1)./s1;
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
                PHI2i = PHI(i,2);
                PHIs = PHIs + (PHIii - PHIi)*PHI2i;
                
                E1i = y(1)*PHI(i,1) - s(1)*phi(i,1);
                E1ii = y(1)*PHI(i+1,1) - s(1)*phi(i+1,1);
                E2i = y(2)*PHI(i,2) - s(2)*phi(i,2);
                E2ii = y(2)*PHI(i+1,2) - s(2)*phi(i+1,2);
                E1s = E1s + (E1ii - E1i)*PHI2i;
                E2s = E2s + (E2i - E2ii)*PHIii;
                
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
                dist(i) = norm(yE - Y(i,:));
            end
            [~,ind] = min(dist);
            fE = Y(ind,:);
            E = P*norm(yE - fE);
        end
        V = E(ones(1,numel(subobj)));
        
    case 'emoo2'
        % MOO Expected improvement, front domination
        if any(s == 0)
            E = zeros(1,numel(subobj));
        else
            M0 = size(Y,1);
            y1 = repmat(y,M0,1);
            s1 = repmat(s,M0,1);
            u = (Y-y1)./s1;
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
                dist(i) = norm(yE - Y(i,:));
            end
            [~,ind] = min(dist);
            fE = Y(ind,:);
            E = P*norm(yE - fE);
        end
        V = E(ones(1,numel(subobj)));
        
end

update = any( V > threshold );

end
