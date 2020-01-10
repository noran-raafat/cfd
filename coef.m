function [c11,c12,c22]=coef(ip,jp)
global x y imax jmax jair il it cord yal yau ps psp dx dy r d1 d2 omega Vinf cosa sina 

% '____________________________________________________________________
% '____________________________________________________________________________
% coef;  ' Subprogram determines the geometrical coefficients
%        ' -------------------------------------------------
%        ' given;
%        '    1]  x(ii,jj) and y(ii,jj) for the whole domain ,and
%        '    2]  The specific point  ip , jp
%        '        where            ip = 2*i-1   &  jp = 2*j-1
%        ' this subprogram calculates the metric terms and the jacobian.
      
         d1x = (x(ip + 1, jp) - x(ip - 1, jp)) / d1;
         d1y = (y(ip + 1, jp) - y(ip - 1, jp)) / d1;
         d2x = (x(ip, jp + 1) - x(ip, jp - 1)) / d2;
         d2y = (y(ip, jp + 1) - y(ip, jp - 1)) / d2;
       
         jaco = d1x * d2y - d1y * d2x;
         c11 = (d2x * d2x + d2y * d2y) / jaco;
         c12 = -(d1x * d2x + d1y * d2y) / jaco;
         c22 = (d1x * d1x + d1y * d1y) / jaco;
         