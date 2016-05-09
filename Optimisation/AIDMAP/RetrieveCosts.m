function [Costs] = RetrieveCosts(Solutions,ListNodes)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(Solutions)
    for j = 2:length(Solutions{i})
        currnode = Solutions{i}(j);
        Costs{1,i}(j-1) = ListNodes.(char(currnode)).length;
    end
end

end

