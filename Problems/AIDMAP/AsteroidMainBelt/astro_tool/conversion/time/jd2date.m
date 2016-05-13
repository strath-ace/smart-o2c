function date = jd2date(jd)

%JD2DATE Gregorian calendar date from Julian day number.
%
% DATE = JD2DATE(JD) returns the Gregorian calendar date (year, month, day,
% hour, minute, and second) corresponding to the Julian day number JD.
%
% Note: jd must be a non-negative real. This means that the function is
%       valid for the whole range of dates since 12:00 noon 24 November
%       -4713, Gregorian calendar. (This bound is due to the function
%       FRACDAY2HMS that does not work for negative inputs)
%
% INPUT
%   jd          Date in Julian Day. The JD (Julian day) count is from 0 at
%               12:00 noon, 1 January -4712 (4713 BC), Julian proleptic
%               calendar. The corresponding date in Gregorian calendar is
%               12:00 noon, 24 November -4713. It must be a non-negative
%               real.
%
% OUTPUT
%	date        Date in the Gregorian calendar, as a 6-element vector
%               [year, month, day, hour, minute, second]. For dates before
%               1582, the resulting date components are valid only in the
%               Gregorian proleptic calendar. This is based on the
%               Gregorian calendar but extended to cover dates before its
%               introduction.
%
% REFERENCES
%   Formula from http://en.wikipedia.org/wiki/Julian_day
%   (last visited 16/02/2008)
%   Compared to http://pdc.ro.nu/mjd.cgi for a few dates, the same results
%   were found
%
% See also DATE2JD.
%
% FUNCTIONS CALLED
%   FRACDAY2HMS
%
% - Nicolas Croisard - 16/02/2008
% - Revised by Matteo Ceriotti - 21/02/2008 - Validated with:
%   - Fliegel, Van Flandern "A machine Algorithm for Processing Calendar
%     Dates", Communications of the ACM, 1968. Also on Wertz, "Space
%     Mission Analysis and Design".
%   - A revised version of the algorithm on Vallado, "Fundamentals of
%     Astrodynamics and Applications", third edition, for dates from year
%     1900 to year 2100.
%
% ------------------------- - SpaceART Toolbox - --------------------------

% Check the input
if nargin ~= 1 || numel(jd) ~= 1
    error('JD2DATE:incorrectInput','The input should be a single real');
end
% if jd < 0
%     error('JD2DATE:jdLessThanZero','The input jd value cannot be negative');
% end

% Adding 0.5 to JD and taking FLOOR ensures that the date is correct.
j = floor(jd+0.5) + 32044;
g = floor(j/146097);
dg = mod(j,146097);
c = floor((floor(dg/36524)+1) * 3/4);
dc = dg - c*36524;
b = floor(dc/1461);
db = mod(dc,1461);
a = floor((floor(db/365)+1) * 3/4);
da = db - a*365;
y = g*400 + c*100 + b*4 + a;
m = floor((da*5 + 308)/153) - 2;
d = da - floor((m+4)*153/5) + 122;

% Year, Month and Day
Y = y-4800 + floor((m+2)/12);
M = mod((m+2),12) + 1;
D = floor(d+1);

% Hour, Minute and Second
[hrs, mn, sec] = fracday2hms(mod(jd+0.5,floor(jd+0.5)));

% Prepare output
date = [Y, M, D, hrs, mn, sec];


return