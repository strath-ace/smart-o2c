function [] = SaveTrajectorySolution(Sequence,ListNodes,filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(filename,'w');
fprintf(fid,'---------------------------------------------------------------------------\n');
fprintf(fid,'\tAsteroid\tLambert T_dep\tLambert ToF\t\tLambert T_arr\t\tdV\n');
fprintf(fid,'---------------------------------------------------------------------------\n');
formatspec = '%12s\t%12s\t%12.4f\t%12s\t%8.2f\n';
%formatspec = '%s';
Asteroids = Sequence';
dV_tot = 0;
for i = 2:length(Asteroids)
%Find asteroid names
asteroidsplit = strsplit(Asteroids{i},'____');
asteroidnamesplit = strsplit(char(asteroidsplit(2)),'___');
asteoroidname = asteroidnamesplit{1};
AsteroidsNames= char(asteoroidname);

%Departure Date Lambert Arc
lambertdepdate2000 = ListNodes.(Sequence{i}).attributes.t_dep;
lambertdepdate = datestr(mjd20002date(lambertdepdate2000),'yyyy/mm/dd');

%ToF Lambert Arc [Days]
lamberttof = ListNodes.(Sequence{i}).attributes.tof;

%Arrival Date at Asteroid Node
lambertarrdate2000 = ListNodes.(Sequence{i}).attributes.t_arr;
lambertarrdate = datestr(mjd20002date(lambertarrdate2000),'yyyy/mm/dd');

%dV [km/s]
dV = ListNodes.(Sequence{i}).attributes.dV_tot;
dV_tot = dV_tot+dV;
fprintf(fid,formatspec,AsteroidsNames,lambertdepdate,lamberttof,lambertarrdate,dV);
end

fprintf(fid,'---------------------------------------------------------------------------\n');
fprintf(fid,'%12s\t%12s\t%12.4f\t%12s\t\t\t\t%8.2f\n','Total:',[],[],[],dV_tot);
% fprintf(fid,formatspec,char(AsteroidsNames),lambertdepdate,lamberttof,lamberttof,lambertarrdate,dV)

end

