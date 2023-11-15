%  Antenas
% Computation Work
%
% Francisco Freiria - 97236
% Joao Morais – 83916


function [M, X, I, A] = Fun(L, N, a,theta,Zl,Voc_calc)

format short
%format longg

lambda=1;
delta = L/N;
k0 = (2*pi)/lambda; 
E0=1;

%Sistemas de equações Ax=B
A = zeros(N, N);
B = zeros(N,1);
M = zeros (N+1,1);
X = zeros (N+1,1);

%Funções K(zm,z') e z(m) z´=x
K = @(x,zm) (1/(4*pi)*((exp(-1i*k0*sqrt(a.^2+(zm-x).^2)))/(sqrt(a.^2+(zm-x).^2))));
z = @(y) -L/2+(y+1/2)*delta;

%Loop da Matriz A e do vector B
for m =1:N
    for n = 1:N
        if(n ~= N-1 && n ~= N && n ~= (N-1)/2)
        A(m,n) = integral(@(x)K(x,z(m-1)), z(n)-delta/2, z(n)+delta/2, 'ArrayValued', true);
        elseif n == N-1
            A(m,n) = cos(k0*z(m-1));  
        elseif n == N
            A(m,n) = sin(k0*z(m-1));
        elseif n == (N-1)/2
            if Voc_calc == 0
                t= 1j*(Zl/(240*pi))*(sin(k0 * abs(z(m-1))));
                A(m,n) = integral(@(x)K(x,z(m-1) - t), z(n)-delta/2, z(n)+delta/2, 'ArrayValued', true); 
            elseif Voc_calc == 1
                Zl=1;
                t= 1j*(Zl/(240*pi))*(sin(k0 * abs(z(m-1))));
                A(m,n) = -t;
            end
        end  
    end
    B(m) = -((1i*E0)/(120*pi*k0*sin(theta)))*(exp(1i*k0*cos(theta)*z(m-1)));
end

I = pinv(A)*B;
disp(A);
disp(B);
disp(I);

for f =1:N
    if (f ~=(N-1) && f ~= N)
        M(f+1) = I(f); 
    elseif (f==(N-1))
        M(N)=0;
    elseif (f==(N))
        M(N)=0;   
    end
end


for g =1:N+1
    X(g)= z(g-1.5);
    
end

end

