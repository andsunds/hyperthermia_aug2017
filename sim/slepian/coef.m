function [a] = coef(v1,v2,rot_order)
%The coefficient function that calculates the coefficients of the symmetric
%matrix
if v1==v2
    a=.5*((besselj(rot_order,v1)).^2-...
          besselj(rot_order-1,v1).*besselj(rot_order+1,v1));
else
    a=(v2*besselj(rot_order,v1)*besselj(rot_order-1,v2)-...
       v1*besselj(rot_order-1,v1)*besselj(rot_order,v2))/...
       (v1^2-v2^2);
end
a=a*sqrt(v1*v2);
end