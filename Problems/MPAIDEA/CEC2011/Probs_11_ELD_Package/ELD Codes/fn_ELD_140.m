<<<<<<< HEAD
% Downloaded from:
% http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2fepnsugan%2fPublicSite%2fShared%20Documents%2fCEC%202011%2d%20RWP&FolderCTID=&View=%7bDAF31868%2d97D8%2d4779%2dAE49%2d9CEC4DC3F310%7d

=======
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
>>>>>>> 56eeb4b328a0319b2f58a2e2413248a83fcc168b
%% ECONOMIC LOAD DISPATCH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Input  => Population of Row vectors (generation units' generations)
%%%  Output => Each is a column vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Total_Cost Cost Total_Penalty] = fn_ELD_140(Input_Population,Display)
%% DATA REQUIRED
[Pop_Size No_of_Units] = size(Input_Population);
Power_Demand = 49342; % in MW 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % ============= 140 unit system data ==========
%% [Pmin       Pmax          c              b               a   ];
Data1= ...
    [71        119        1220.645        61.242        0.032888;
    120        189        1315.118        41.095        0.00828;
    125        190        874.288         46.31         0.003849;
    125        190        874.288         46.31         0.003849;
    90         190        1976.469        54.242        0.042468;
    90         190        1338.087        61.215        0.014992;
    280        490        1818.299        11.791        0.007039;
    280        490        1133.978        15.055        0.003079;
    260        496        1320.636        13.226        0.005063;
    260        496        1320.636        13.226        0.005063;
    260        496        1320.636        13.226        0.005063;
    260        496        1106.539        14.498        0.003552;
    260        506        1176.504        14.651        0.003901;
    260        509        1176.504        14.651        0.003901;
    260        506        1176.504        14.651        0.003901;
    260        505        1176.504        14.651        0.003901;
    260        506        1017.406        15.669        0.002393;
    260        506        1017.406        15.669        0.002393;
    260        505        1229.131        14.656        0.003684;
    260        505        1229.131        14.656        0.003684;
    260        505        1229.131        14.656        0.003684;
    260        505        1229.131        14.656        0.003684;
    260        505        1267.894        14.378        0.004004;
    260        505        1229.131        14.656        0.003684;
    280        537        975.926         16.261        0.001619;
    280        537        1532.093        13.362        0.005093;
    280        549        641.989         17.203        0.000993;
    280        549        641.989         17.203        0.000993;
    260        501        911.533         15.274        0.002473;
    260        501        910.533         15.212        0.002547;
    260        506        1074.81         15.033        0.003542;
    260        506        1074.81         15.033        0.003542;
    260        506        1074.81         15.033        0.003542;
    260        506        1074.81         15.033        0.003542;
    260        500        1278.46         13.992        0.003132;
    260        500        861.742         15.679        0.001323;
    120        241        408.834         16.542        0.00295;
    120        241        408.834         16.542        0.00295;
    423        774        1288.815        16.518        0.000991;
    423        769        1436.251        15.815        0.001581;
    3          19         699.988         75.464        0.90236;
    3          28         134.544         129.544       0.110295;
    160        250        3427.912        56.613        0.024493;
    160        250        3751.772        54.451        0.029156;
    160        250        3918.78         54.736        0.024667;
    160        250        3379.58         58.034        0.016517;
    160        250        3345.296        55.981        0.026584;
    160        250        3138.754        61.52         0.00754;
    160        250        3453.05         58.635        0.01643;
    160        250        5119.3          44.647        0.045934;
    165        504        1898.415        71.584        0.000044;
    165        504        1898.415        71.584        0.000044;
    165        504        1898.415        71.584        0.000044;
    165        504        1898.415        71.584        0.000044;
    180        471        2473.39         85.12         0.002528;
    180        561        2781.705        87.682        0.000131;
    103        341        5515.508        69.532        0.010372;
    198        617        3478.3          78.339        0.007627;
    100        312        6240.909        58.172        0.012464;
    153        471        9960.11         46.636        0.039441;
    163        500        3671.997        76.947        0.007278;
    95         302        1837.383        80.761        0.000044;
    160        511        3108.395        70.136        0.000044;
    160        511        3108.395        70.136        0.000044;
    196        490        7095.484        49.84         0.018827;
    196        490        3392.732        65.404        0.010852;
    196        490        7095.484        49.84         0.018827;
    196        490        7095.484        49.84         0.018827;
    130        432        4288.32         66.465        0.03456;
    130        432        13813.001       22.941        0.08154;
    137        455        4435.493        64.314        0.023534;
    137        455        9750.75         45.017        0.035475;
    195        541        1042.366        70.644        0.000915;
    175        536        1159.895        70.959        0.000044;
    175        540        1159.895        70.959        0.000044;
    175        538        1303.99         70.302        0.001307;
    175        540        1156.193        70.662        0.000392;
    330        574        2118.968        71.101        0.000087;
    160        531        779.519         37.854        0.000521;
    160        531        829.888         37.768        0.000498;
    200        542        2333.69         67.983        0.001046;
    56         132        2028.954        77.838        0.13205;
    115        245        4412.017        63.671        0.096968;
    115        245        2982.219        79.458        0.054868;
    115        245        2982.219        79.458        0.054868;
    207        307        3174.939        93.966        0.014382;
    207        307        3218.359        94.723        0.013161;
    175        345        3723.822        66.919        0.016033;
    175        345        3551.405        68.185        0.013653;
    175        345        4332.615        60.821        0.028148;
    175        345        3493.739        68.551        0.01347;
    360        580        226.799         2.842         0.000064;
    415        645        382.932         2.946         0.000252;
    795        984        156.987         3.096         0.000022;
    795        978        154.484         3.04          0.000022;
    578        682        332.834         1.709         0.000203;
    615        720        326.599         1.668         0.000198;
    612        718        345.306         1.789         0.000215;
    612        720        350.372         1.815         0.000218;
    758        964        370.377         2.726         0.000193;
    755        958        367.067         2.732         0.000197;
    750        1007       124.875         2.651         0.000324;
    750        1006       130.785         2.798         0.000344;
    713        1013       878.746         1.595         0.00069;
    718        1020       827.959         1.503         0.00065;
    791        954        432.007         2.425         0.000233;
    786        952        445.606         2.499         0.000239;
    795        1006       467.223         2.674         0.000261;
    795        1013       475.94          2.692         0.000259;
    795        1021       899.462         1.633         0.000707;
    795        1015       1000.367        1.816         0.000786;
    94         203        1269.132        89.83         0.014355;
    94         203        1269.132        89.83         0.014355;
    94         203        1269.132        89.83         0.014355;
    244        379        4965.124        64.125        0.030266;
    244        379        4965.124        64.125        0.030266;
    244        379        4965.124        64.125        0.030266;
    95         190        2243.185        76.129        0.024027;
    95         189        2290.381        81.805        0.00158;
    116        194        1681.533        81.14         0.022095;
    175        321        6743.302        46.665        0.07681;
    2          19         394.398         78.412        0.953443;
    4          59         1243.165        112.088       0.000044;
    15         83         1454.74         90.871        0.072468;
    9          53         1011.051        97.116        0.000448;
    12         37         909.269         83.244        0.599112;
    10         34         689.378         95.665        0.244706;
    112        373        1443.792        91.202        0.000042;
    4          20         535.553         104.501       0.085145;
    5          38         617.734         83.015        0.524718;
    5          19         90.966          127.795       0.176515;
    50         98         974.447         77.929        0.063414;
    5          10         263.81          92.779        2.740485;
    42         74         1335.594        80.95         0.112438;
    42         74         1033.871        89.073        0.041529;
    41         105        1391.325        161.288       0.000911;
    17         51         4477.11         161.829       0.005245;
    7          19         57.794          84.972        0.234787;
    7          19         57.794          84.972        0.234787;
    26         40         1258.437        16.087        1.111878;];

%% [ Po          UR         DR ];
Data2= ...
[   98.4         30         120;
    134          30         120;
    141.5        60         60;
    183.3        60         60;
    125          150        150;
    91.3         150        150;
    401.1        180        300;
    329.5        180        300;
    386.1        300        510;
    427.3        300        510;
    412.2        300        510;
    370.1        300        510;
    301.8        600        600;
    368          600        600;
    301.9        600        600;
    476.4        600        600;
    283.1        600        600;
    414.1        600        600;
    328          600        600;
    389.4        600        600;
    354.7        600        600;
    262          600        600;
    461.5        600        600;
    371.6        600        600;
    462.6        300        300;
    379.2        300        300;
    530.8        360        360;
    391.9        360        360;
    480.1        180        180;
    319          180        180;
    329.5        600        600;
    333.8        600        600;
    390          600        600;
    432          600        600;
    402          660        660;
    428          900        900;
    178.4        180        180;
    194.1        180        180;
    474          600        600;
    609.8        600        600;
    17.8         210        210;
    6.9          366        366;
    224.3        702        702;
    210          702        702;
    212          702        702;
    200.8        702        702;
    220          702        702;
    232.9        702        702;
    168          702        702;
    208.4        702        702;
    443.9        1350       1350;
    426.0        1350       1350;
    434.1        1350       1350;
    402.5        1350       1350;
    357.4        1350       1350;
    423          720        720;
    220          720        720;
    369.4        2700       2700;
    273.5        1500       1500;
    336          1656       1656;
    432          2160       2160;
    220          900        900;
    410.6        1200       1200;
    422.7        1200       1200;
    351          1014       1014;
    296          1014       1014;
    411.1        1014       1014;
    263.2        1014       1014;
    370.3        1350       1350;
    418.7        1350       1350;
    409.6        1350       1350;
    412          1350       1350;
    423.2        780        780;
    428          1650       1650;
    436          1650       1650;
    428          1650       1650;
    425          1650       1650;
    497.2        1620       1620;
    510          1482       1482;
    470          1482       1482;
    464.1        1668       1668;
    118.1        120        120;
    141.3        180        180;
    132          120        180;
    135          120        180;
    252          120        180;
    221          120        180;
    245.9        318        318;
    247.9        318        318;
    183.6        318        318;
    288          318        318;
    557.4        18         18;
    529.5        18         18;
    800.8        36         36;
    801.5        36         36;
    582.7        138        204;
    680.7        144        216;
    670.7        144        216;
    651.7        144        216;
    921          48         48;
    916.8        48         48;
    911.9        36         54;
    898          36         54;
    905          30         30;
    846.5        30         30;
    850.9        30         30;
    843.7        30         30;
    841.4        36         36;
    835.7        36         36;
    828.8        36         36;
    846          36         36;
    179          120        120;
    120.8        120        120;
    121          120        120;
    317.4        480        480;
    318.4        480        480;
    335.8        480        480;
    151          240        240;
    129.5        240        240;
    130          120        120;
    218.9        180        180;
    5.4          90         90;
    45           90         90;
    20           300        300;
    16.3         162        162;
    20           114        114;
    22.1         120        120;
    125          1080       1080;
    10           60         60;
    13           66         66;
    7.5          12         6;
    53.2         300        300;
    6.4          6          6;
    69.1         60         60;
    49.9         60         60;
    91           528        528;
    41           300        300;
    13.7         18         30;
    7.4          18         30;
    28.6         72         120;];

%% Loss Co-efficients
        B1=zeros(140,140);
        B2=zeros(1,140);
        B3=0; 

%% [ Unit No.     e           f]
Data3 = ...
    [   5        700        0.080;
        10       600        0.055;
        15       800        0.060;
        22       600        0.050;
        33       600        0.043;
        40       600        0.043;
        52       1100       0.043;
        70       1200       0.030;
        72       1000       0.050;
        84       1000       0.050;
        119      600        0.070;
        121      1200       0.043;];
%% INITIALIZATIONS
Data1(Data3(:,1),6:7)= Data3(:,2:3);
Pmin = Data1(:,1)'; 
Pmax = Data1(:,2)';
a = Data1(:,5)';
b = Data1(:,4)';
c = Data1(:,3)';
e = Data1(:,6)';
f = Data1(:,7)';
%% [ Unit No.  Zone1min    Zone1max     Zone2min   Zone2max    Zone3min    Zone3max]
Data4 = ...
    [    8        250         280         305        335         420         450;
         32       220         250         320        350         390         420;
         74       230         255         365        395         430         455;
         136      50          75          85         95          0           0;];    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Initial_Generations = Data2(:,1)';
Up_Ramp = Data2(:,2)';
Down_Ramp = Data2(:,3)';
Up_Ramp_Limit = min(Pmax,Initial_Generations+Up_Ramp);
Down_Ramp_Limit = max(Pmin,Initial_Generations-Down_Ramp);

Data2(Data4(:,1),4:9)= Data4(:,2:end);
Prohibited_Operating_Zones_POZ = Data2(:,4:end)';
No_of_POZ_Limits = size(Prohibited_Operating_Zones_POZ,1);
POZ_Lower_Limits = Prohibited_Operating_Zones_POZ(1:2:No_of_POZ_Limits,:);
POZ_Upper_Limits = Prohibited_Operating_Zones_POZ(2:2:No_of_POZ_Limits,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATIONS
for i = 1:Pop_Size
    x = Input_Population(i,:);
    
    Power_Loss  = (x*B1*x') + (B2*x') + B3;
    Power_Loss  = round(Power_Loss *10000)/10000;
%%% Power Balance Penalty Calculation
    Power_Balance_Penalty = abs(Power_Demand + Power_Loss - sum(x));    
%%% Capacity Limits Penalty Calculation
    Capacity_Limits_Penalty = (abs(x-Pmin)-(x-Pmin)) + sum(abs(Pmax-x)-(Pmax-x));
%%% Ramp Rate Limits Penalty Calculation
    Ramp_Limits_Penalty = (abs(x-Down_Ramp_Limit)-(x-Down_Ramp_Limit)) + (abs(Up_Ramp_Limit-x)-(Up_Ramp_Limit-x));    
%%% Prohibited Operating Zones Penalty Calculation
    temp_x = repmat(x,No_of_POZ_Limits/2,1);
    POZ_Penalty = sum(sum((POZ_Lower_Limits<temp_x & temp_x<POZ_Upper_Limits).*min(temp_x-POZ_Lower_Limits,POZ_Upper_Limits-temp_x)));    
%%% Total Penalty Calculation
    Total_Penalty(i,1) = 1e7*Power_Balance_Penalty + 1e5*sum(Capacity_Limits_Penalty) + 1e7*sum(Ramp_Limits_Penalty) + 1e5*POZ_Penalty;
%%% Cost Calculation
    Cost(i,1) = sum( a.*(x.^2) + b.*x + c + abs(e.*sin(f.*(Pmin-x))) );
    Total_Cost(i,1) = Cost(i,1) + Total_Penalty(i,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (nargin>1)
        disp('----------------------------------------------------------------------------');
        disp(sprintf('140 UNIT SYSTEM'));
        disp(sprintf('Total_Power_Generation   : %17.8f ',sum(x))); 
        disp(sprintf('Power_Balance_Penalty    : %17.8f ',Power_Balance_Penalty));
        disp(sprintf('Capacity_Limits_Penalty  : %17.8f ',sum(Capacity_Limits_Penalty) ));
        disp(sprintf('Ramp_Limits_Penalty      : %17.8f ',sum(Ramp_Limits_Penalty)));
        disp(sprintf('POZ_Penalty              : %17.8f ',POZ_Penalty));
        disp(sprintf('Cost                     : %17.8f ',Cost(i,1)));
        disp(sprintf('Total_Penalty            : %17.8f ',Total_Penalty(i,1)));
        disp(sprintf('Total_Objective_Value    : %17.8f ',Total_Cost(i,1))); 
        disp('----------------------------------------------------------------------------');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end
end