function [J]=mexspaceartCassiniNOdsmA_MO(ya)

ll=[-1000 30   100  30   400   1000
     0    400  470  400  2000  6000];
 
y=(ll(2,:)-ll(1,:)).*ya+ll(1,:);

J=[mexspaceartCassiniNOdsm(y) sum(y(2:end))];


