% Downloaded from:
% http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2fepnsugan%2fPublicSite%2fShared%20Documents%2fCEC%202011%2d%20RWP&FolderCTID=&View=%7bDAF31868%2d97D8%2d4779%2dAE49%2d9CEC4DC3F310%7d


% INPUT FILE FOR GRAVER 6 BUS SYSTEM %%%
n=6;
    %      No  from  to   X      num_Line   Pijmax       cost  Overloads    %Pij
Linedata=[ 1    1    2   0.4        1        1           inf       1      ;
           2    1    4   0.6        1        0.8         inf       0    ;
           3    1    5   0.2        1        1           inf       0    ;
           4    2    3   0.2        1        1           inf       0    ;
           5    2    4   0.4        1        1           inf       1    ;
           6    3    5   0.2        1        1           20       1    ;
           7    6    2   0.3        1        1           30       1    ;];

Candidate=[ 1    1    2   0.4        1        1           40       1      ;
            2    1    3   0.38       1        1           38        1;
            3    1    4   0.6        1        0.8         60       0    ;   
            4    1    5   0.2        1        1           20       0    ;
            5    1    6   0.68      1        0.7          68        0;
            6    2    3   0.2        1        1           20       0    ;
            7    2    4   0.4       1        1            40       1    ;
            8    2    5   0.31      1        1            31       1    ;
            9    6    2   0.3        1        1           30        1 ;
            10   3    4   0.69      1          0.82       59        1;
            11    3    5   0.2        1        1           20        1 ;
            12    6    3   0.48        0       1           48       0;
            13    4    5   0.63        0       0.75        63       0;
            14    4    6   0.30        0        1          30       0;
            15    5     6  0.61        0        0.78       61       0;
                 
           ];
Pgen=[0.5 0 1.65 0 0 5.45];
Pload=[0.8 2.4 0.4 1.6 2.4 0 ];
