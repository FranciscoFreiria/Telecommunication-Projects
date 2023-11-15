

beta2 = 0;
beta3 = 2;
beta_2 = beta2*1e-24;
beta_3 = beta3*1e-36;

L = 100;
a = ((beta_3)^2/(abs(beta_2)^3*L));

const = (beta_3)^(2/3)*L^(2/3);

C = -10:0.0001:10;
C_critico = zeros(length(a),length(C));
for i=1:length(C)
    
   [C_critico(i)]=Chirp_parameter(a,C(i));
    
end

[val, indx] = min(abs(C_critico-0));  
%C_cri = C(indx);
C_cri = 0;
p = 0.25*(1+C_cri^2);

Y_0 = p^(2/3);

fun = @(eps)p^(2/3)*(1+(1/2)*eps^2);
Y_1 = fun(1);

var = sqrt(Y_1)*(abs(beta_3))^(1/3)*(L)^(1/3);
Bit_rate = 1/(4*var);

Broad = sqrt((Y_1*const)/(Y_0*const));

function [y,Y_0] = Chirp_parameter(a,C)

p = 0.25*(1+C^2);
Y_0 = p^(2/3);
y = 32*Y_0.^2*C -8*Y_0*(1+C.^2) - (a)*(1+C.^2)^2;

end