    A = [1.0000         0        0
    1.0000    0.2000         0
    1.0000    0.4000         0
    1.0000         0    0.2000
    1.0000    0.2000    0.2000
    1.0000    0.4000    0.2000
    1.0000         0    0.4000
    1.0000    0.2000    0.4000
    1.0000    0.4000    0.4000];
B = [  0  
    1.0396
    0.0206
         0
    1.0396
    0.0206
         0
    1.0396
    0.0206];
C = A\B

 AA = [   1.0000         0    0.2000
    1.0000    0.2000    0.2000
    1.0000    0.4000    0.2000];
BB = [      0
    1.0396
    0.0206];
CC = AA\BB

 AAA = [   1.0000         0
    1.0000    0.2000
    1.0000    0.4000];
BBB = [0
    1.0396
    0.0206];
CCC = AAA\BBB

uc =     [0.6930
    0.6898
    0.6834
    0.6801
    0.6930
    0.6898
    0.6834
    0.6801];
A = [0.0200    0.0200    0.0200    0.0200    0.0200    0.0200    0.0200    0.0200];
c = [0.6914; 0.6817];
h = [0.2,0.2];

