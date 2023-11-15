
n = 960;
tr = 9.07;
di = 57447.9;
d0 = 13283.92;
N=525960;
infl = 0.02;
x = zeros(25,1);
Custo = zeros(25,1);

for i=1:25
    x(i)=i; 
    Custo(i) = (((d0/tr)+di)/(n*(0.2+0.02*i)*N/3))*(1+infl)^(i);
    
end   

figure(1)
hold on
plot(x,Custo)
grid on
xlabel('Anos')
ylabel('Custo (€)')
title('Evolução do custo médio de uma chamada de 3 minutos em 25 anos')
hold off