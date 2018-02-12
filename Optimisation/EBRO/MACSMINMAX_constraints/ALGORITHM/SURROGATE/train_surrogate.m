% function surrogate = train_surrogate(x, y, surrogate, minmax, moo)
% %
% %   surrogate = train_surrogate(x, y, surrogate)
% %
% %   INPUT
% %        x: samples in the parameter space
% %        y: corresponding values of the true model
% %        surrogate: structure containing the surrogate model
% %
% %   OUTPUT
% %        surrogate: surrogate model
% %

% % Change Log
% % Simone Alicino, 2012
% % Revised by Massimiliano Vasile July 2012: simplified call, revised update
% % of the surrogate

% if nargin < 5
%     moo = 0;
%     if nargin < 4
%         minmax = 0;
%     end
% end

% % Remove double (to certain precision) points, if any
% x0 = roundn(x,-6);
% [~,ix] = unique(x0,'rows');
% X = x(ix,:);
% Y = y(ix,:);
% [~,i] = sort(Y(:,1));
% Y = Y(i,:);
% X = X(i,:);

% % Reduce size of surrogate if it exceeds maximum size
% [lY,nY] = size(Y);
% if lY > surrogate.maxsize
%     if moo == 0
%         mY = mean(Y);
%         sel0 = zeros(lY,nY);
%         for k = 1:nY
%             if minmax == -1
%                 sel0(:,k) = (Y(:,k) >= mY(k));
%             else
%                 sel0(:,k) = (Y(:,k) <= mY(k));
%             end
%         end
%         sel1 = sum(sel0,2);
%         sel = find(sel1 == nY);
%     else
%         dom = dominance(Y,0);
%         Y = Y(dom == 0,:);
%         X = X(dom == 0,:);
%         m = size(Y,1);
%         if m <= surrogate.maxsize
%             sel = 1:m;
%         else
%             sel = archresize(Y, m, surrogate.maxsize);
%         end
%     end
%     if length(sel) > 1
%         Y = Y(sel,:);
%         X = X(sel,:);
%         if length(sel) > surrogate.maxsize
%             [~, iy] = sort(Y(:,1));
%             Y = Y(iy(1:surrogate.maxsize),:);
%             X = X(iy(1:surrogate.maxsize),:);
%         end
%     end
% end

% % Fit surrogate model
% if numel(ix) > 1
    
%     trainfunc = surrogate.trainfunc;
%     model = trainfunc(X,Y,surrogate);
%     surrogate.model = model;
    
% end

% surrogate.x = X;
% surrogate.y = Y;

% end
