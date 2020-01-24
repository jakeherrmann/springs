


clc
C = (1:100)' ;
R = (2/sqrt(3))*C + 0.5*(sqrt(3)-1) ;
sortrows( [ round(R) , C , mod(R,1) ] , [3,1,2] )