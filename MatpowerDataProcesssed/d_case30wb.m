% 01/10/12 File data originated from Matpower data file
% 

Bus.con = [ ...
      1      135        1        0    1    1;
      2      135        1        0    1    1;
      3      135        1        0    1    1;
      4      135        1        0    1    1;
      5      135        1        0    1    1;
      6      135        1        0    1    1;
      7      135        1        0    1    1;
      8      135        1        0    1    1;
      9      135        1        0    1    1;
     10      135        1        0    3    1;
     11      135        1        0    1    1;
     12      135        1        0    2    1;
     13      135        1        0    2    1;
     14      135        1        0    2    1;
     15      135        1        0    2    1;
     16      135        1        0    2    1;
     17      135        1        0    2    1;
     18      135        1        0    2    1;
     19      135        1        0    2    1;
     20      135        1        0    2    1;
     21      135        1        0    3    1;
     22      135        1        0    3    1;
     23      135        1        0    2    1;
     24      135        1        0    3    1;
     25      135        1        0    3    1;
     26      135        1        0    3    1;
     27      135        1        0    3    1;
     28      135        1        0    1    1;
     29      135        1        0    3    1;
     30      135        1        0    3    1;
   ];

SW.con = [ ...
      1      100      135        1        0       15       -2     1.05     0.95   0.2354 1 1  1;
   ];

PV.con = [ ...
   2      100      135   0.6097        1        6       -2     1.05     0.95 1  1;
  22      100      135   0.2159        1     6.25     -1.5     1.05     0.95 1  1;
  27      100      135   0.2691        1     4.87     -1.5     1.05     0.95 1  1;
  23      100      135    0.192        1        4       -1     1.05     0.95 1  1;
  13      100      135     0.37        1     4.47     -1.5     1.05     0.95 1  1;
   ];

PQ.con = [ ...
   2 100      135    0.217    0.127     1.05     0.95 0 1;
   3 100      135    0.024    0.012     1.05        0 0 1;
   4 100      135    0.076    0.016     1.05        0 0 1;
   7 100      135    0.228    0.109     1.05        0 0 1;
   8 100      135      0.3      0.3     1.05        0 0 1;
  10 100      135    0.058     0.02     1.05        0 0 1;
  12 100      135    0.112    0.075     1.05        0 0 1;
  14 100      135    0.062    0.016     1.05        0 0 1;
  15 100      135    0.082    0.025     1.05        0 0 1;
  16 100      135    0.035    0.018     1.05        0 0 1;
  17 100      135     0.09    0.058     1.05        0 0 1;
  18 100      135    0.032    0.009     1.05        0 0 1;
  19 100      135    0.095    0.034     1.05        0 0 1;
  20 100      135    0.022    0.007     1.05        0 0 1;
  21 100      135    0.175    0.112     1.05        0 0 1;
  23 100      135    0.032    0.016     1.05     0.95 0 1;
  24 100      135    0.087    0.067     1.05        0 0 1;
  26 100      135    0.035    0.023     1.05        0 0 1;
  29 100      135    0.024    0.009     1.05        0 0 1;
  30 100      135    0.106    0.019     1.05        0 0 1;
   5 100      135        0        0     1.05        0 0 1;
   6 100      135        0        0     1.05        0 0 1;
   9 100      135        0        0     1.05        0 0 1;
  11 100      135        0        0     1.05        0 0 1;
  25 100      135        0        0     1.05        0 0 1;
  28 100      135        0        0     1.05        0 0 1;
   ];

Shunt.con = [ ...
   5 100      135 60        0   0.0019 1;
  24 100      135 60        0   0.0004 1;
   ];

Line.con = [ ...
   1    2 100      135 60 0         0     0.02     0.06     0.03        0        0       13      1.3      1.3  1;
   1    3 100      135 60 0         0     0.05     0.19     0.02        0        0       13      1.3      1.3  1;
   2    4 100      135 60 0         0     0.06     0.17     0.02        0        0      6.5     0.65     0.65  1;
   3    4 100      135 60 0         0     0.01     0.04        0        0        0       13      1.3      1.3  1;
   2    5 100      135 60 0         0     0.05      0.2     0.02        0        0       13      1.3      1.3  1;
   2    6 100      135 60 0         0     0.06     0.18     0.02        0        0      6.5     0.65     0.65  1;
   4    6 100      135 60 0         0     0.01     0.04        0        0        0        9      0.9      0.9  1;
   5    7 100      135 60 0         0     0.05     0.12     0.01        0        0        7      0.7      0.7  1;
   6    7 100      135 60 0         0     0.03     0.08     0.01        0        0       13      1.3      1.3  1;
   6    8 100      135 60 0         0     0.01     0.04        0        0        0      3.2     0.32     0.32  1;
   6    9 100      135 60 0         0    1e-05     0.21        0        0        0      6.5     0.65     0.65  1;
   6   10 100      135 60 0         0    1e-05     0.56        0        0        0      3.2     0.32     0.32  1;
   9   11 100      135 60 0         0    1e-05     0.21        0        0        0      6.5     0.65     0.65  1;
   9   10 100      135 60 0         0    1e-05     0.11        0        0        0      6.5     0.65     0.65  1;
   4   12 100      135 60 0         0    1e-05     0.26        0        0        0      6.5     0.65     0.65  1;
  12   13 100      135 60 0         0    1e-05     0.14        0        0        0      6.5     0.65     0.65  1;
  12   14 100      135 60 0         0     0.12     0.26        0        0        0      3.2     0.32     0.32  1;
  12   15 100      135 60 0         0     0.07     0.13        0        0        0      3.2     0.32     0.32  1;
  12   16 100      135 60 0         0     0.09      0.2        0        0        0      3.2     0.32     0.32  1;
  14   15 100      135 60 0         0     0.22      0.2        0        0        0     0.16      1.6     0.16  1;
  16   17 100      135 60 0         0     0.08     0.19        0        0        0      1.6     0.16     0.16  1;
  15   18 100      135 60 0         0     0.11     0.22        0        0        0      1.6     0.16     0.16  1;
  18   19 100      135 60 0         0     0.06     0.13        0        0        0      1.6     0.16     0.16  1;
  19   20 100      135 60 0         0     0.03     0.07        0        0        0      3.2     0.32     0.32  1;
  10   20 100      135 60 0         0     0.09     0.21        0        0        0      3.2     0.32     0.32  1;
  10   17 100      135 60 0         0     0.03     0.08        0        0        0      3.2     0.32     0.32  1;
  10   21 100      135 60 0         0     0.03     0.07        0        0        0      3.2     0.32     0.32  1;
  10   22 100      135 60 0         0     0.07     0.15        0        0        0      3.2     0.32     0.32  1;
  21   22 100      135 60 0         0     0.01     0.02        0        0        0      3.2     0.32     0.32  1;
  15   23 100      135 60 0         0      0.1      0.2        0        0        0      1.6     0.16     0.16  1;
  22   24 100      135 60 0         0     0.12     0.18        0        0        0      1.6     0.16     0.16  1;
  23   24 100      135 60 0         0     0.13     0.27        0        0        0      1.6     0.16     0.16  1;
  24   25 100      135 60 0         0     0.19     0.33        0        0        0      1.6     0.16     0.16  1;
  25   26 100      135 60 0         0     0.25     0.38        0        0        0      1.6     0.16     0.16  1;
  25   27 100      135 60 0         0     0.11     0.21        0        0        0      1.6     0.16     0.16  1;
  28   27 100      135 60 0         0    1e-05      0.4        0        0        0      6.5     0.65     0.65  1;
  27   29 100      135 60 0         0     0.22     0.42        0        0        0      1.6     0.16     0.16  1;
  27   30 100      135 60 0         0     0.32      0.6        0        0        0      1.6     0.16     0.16  1;
  29   30 100      135 60 0         0     0.24     0.45        0        0        0      1.6     0.16     0.16  1;
   8   28 100      135 60 0         0     0.06      0.2     0.02        0        0      3.2     0.32     0.32  1;
   6   28 100      135 60 0         0     0.02     0.06     0.01        0        0      3.2     0.32     0.32  1;
   ];

Supply.con = [ ...
   1       100        0        8        0        0        0        2        0        0        0        0  1 0 1;
   2       100        0        8        0        0        0     1.75        0        0        0        0  1 0 1;
  22       100        0        5        0        0        0        1        0        0        0        0  1 0 1;
  27       100        0      5.5        0        0        0     3.25        0        0        0        0  1 0 1;
  23       100        0        3        0        0        0        3        0        0        0        0  1 0 1;
  13       100        0        4        0        0        0        3        0        0        0        0  1 0 1;
   ];

Bus.names = {...
      'Bus 1'; 'Bus 2'; 'Bus 3'; 'Bus 4'; 'Bus 5'; 
      'Bus 6'; 'Bus 7'; 'Bus 8'; 'Bus 9'; 'Bus 10'; 
      'Bus 11'; 'Bus 12'; 'Bus 13'; 'Bus 14'; 'Bus 15'; 
      'Bus 16'; 'Bus 17'; 'Bus 18'; 'Bus 19'; 'Bus 20'; 
      'Bus 21'; 'Bus 22'; 'Bus 23'; 'Bus 24'; 'Bus 25'; 
      'Bus 26'; 'Bus 27'; 'Bus 28'; 'Bus 29'; 'Bus 30'};

