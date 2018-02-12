function [hrs, mn, sec] = fracday2hms(fracDay)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
% FRACDAY2HMS Convert a fraction of day into hours, minutes, and seconds.
% 
% [HRS, MN, SEC] = FRACDAY2HMS(FRACDAY) converts the fraction of day to
% hours, minutes, and seconds.
% 
% INPUT
%   fracDay    A single real greater or equal to 0 and strictly lower than
%              1.
% 
% OUTPUT
%   hrs        Number of hours as integer greater or equal to 0 and lower
%              or equal to 23.
%   mn         Number of minutes as integer greater or equal to 0 and
%              lower or equal to 59.
%   sec        Number of seconds as a real greater or equal to 0 and
%              strictly lower than 60.
% 
% See also HMS2FRACDAY.
% 
% FUNCTIONS CALLED
%  None
% 
% - Nicolas Croisard - 16/02/2008
% - Revised by Matteo Ceriotti - 21/02/2008
% 
% ------------------------- - SpaceART Toolbox - --------------------------


% Check the input
if nargin ~= 1 || numel(fracDay) ~= 1
    error('FRACDAY2HMS:incorrectInput', ...
          'The input should a single element');
end
if fracDay<0 || fracDay>=1
    error('FRACDAY2HMS:incorrectInput', ...
          ['The input should be real greater or equal to 0 , '...
           'and strictly lower than 1']);
end

temp = fracDay*24;
hrs = fix(temp);
mn = fix((temp-hrs)*60);
sec = (temp-hrs-mn/60)*3600;


return
