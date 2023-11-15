function [y,X_0] = Chirp_parameter(a,C)

X_0 = 2*sqrt((1+C.^2)/12)*cos((1/3)*(acos(4.5*(a)*sqrt((1+C.^2)/12))));
y = 32*X_0.^2*C -8*X_0*(1+C.^2) - (a)*(1+C.^2)^2;

end


