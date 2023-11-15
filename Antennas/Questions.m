%  Antenas
% Computation Work
%
% Francisco Freiria - 97236
% Joao Morais – 83916
%
Questao=1; %################### SELECIONAR A QUESTÃO [1,2,4,5,6,7] ###################%

%-------------------------Questão 1----------------------------------------
if (Questao == 1)
    
    clc;
    clear all;
    
    lambda = 1;
    L = 0.47*lambda;
    a = 0.005*lambda;
    theta= pi/2;
    Zl=0;
    
    N = 21;
    [Z1,X1]= Fun(L, N, a,theta,Zl,0);
    N=51;
    [Z2,X2]= Fun(L, N, a,theta,Zl,0);
    N=71;
    [Z3,X3]= Fun(L, N, a,theta,Zl,0);

    figure (1)
    stairs(X1, imag(Z1), 'red');
    hold on
    stairs(X2, imag(Z2), 'blue');
    stairs(X3, imag(Z3), 'green');
    grid on;
    xlabel('2Z/L');
    ylabel('imag.I(z)');
    legend({'N=21','N=51','N=71'});
    title(['Q1: Parte imaginária Corrente de destribuição de uma antena dipolo de  L/lamda = ',num2str(L/lambda)]);
    hold off

    figure (2)
    stairs(X1, real(Z1), 'red');
    hold on
    stairs(X2, real(Z2), 'blue');
    stairs(X3, real(Z3), 'green');
    grid on;
    xlabel('2Z/L');
    ylabel('|I(z)|');
    legend({'N=21','N=51','N=71'});
    title(['Q1: Parte real Corrente de destribuição de uma antena dipolo de  L/lamda = ',num2str(L/lambda)]);
    hold off
    
    figure (3)
    stairs(X1, abs(Z1), 'red');
    hold on
    stairs(X2, abs(Z2), 'blue');
    stairs(X3, abs(Z3), 'green');
    grid on;
    xlabel('2Z/L');
    ylabel('|I(z)|');
    legend({'N=21','N=51','N=71'});
    title(['Q1- Valor absoluto Corrente de destribuição de uma antena dipolo de  L/lamda = ',num2str(L/lambda)]);
    hold off

%-------------------------Questão 2----------------------------------------

%Fun(N, L, lambda, a)
elseif (Questao == 2)
    
    clc;
    clear all;
    
    lambda = 1;
    L = 1;
    a = 0.005;
    Zl=0;
    theta= pi/2;

    N = 21;
    [Z1,X1]= Fun(L, N, a,theta,Zl,0);
    N=51;
    [Z2,X2]= Fun(L, N, a,theta,Zl,0);
    N=71;
    [Z3,X3]= Fun(L, N, a,theta,Zl,0);

    figure (1)
    stairs(X1, imag(Z1), 'red');
    hold on
    stairs(X2, imag(Z2), 'blue');
    stairs(X3, imag(Z3), 'green');
    grid on;
    xlabel('2Z/L');
    ylabel('|I(z)|');
    legend({'N=21','N=51','N=71'});
    title(['Q2- Parte Imaginária Corrente de destribuição de uma antena dipolo de L/lamda = ',num2str(L/lambda)]);
    hold off

    figure (2)
    stairs(X1, real(Z1), 'red');
    hold on
    stairs(X2, real(Z2), 'blue');
    stairs(X3, real(Z3), 'green');
    grid on;
    xlabel('2Z/L');
    ylabel('|I(z)|');
    legend({'N=21','N=51','N=71'});
    title(['Q2- Parte Real Corrente de destribuição de uma antena dipolo de L/lamda = ',num2str(L/lambda)]);
    hold off
    
    figure (3)
    stairs(X1, abs(Z1), 'red');
    hold on
    stairs(X2, abs(Z2), 'blue');
    stairs(X3, abs(Z3), 'green');
    grid on;
    xlabel('2Z/L');
    ylabel('|I(z)|');
    legend({'N=21','N=51','N=71'});
    title(['Q2- Valor absoluto Corrente de destribuição de uma antena dipolo de L/lamda = ',num2str(L/lambda)]);
    hold off
    
    %------------------------------
    
    L = 1.5;
    N = 21;
    [Z1,X1,p,y]= Fun(L, N, a,theta,Zl,0);
    N=51;
    [Z2,X2]= Fun(L, N, a,theta,Zl,0);
    N=71;
    [Z3,X3]= Fun(L, N, a,theta,Zl,0);
    
    figure (4)
    stairs(X1, imag(Z1), 'red');
    hold on
    stairs(X2, imag(Z2), 'blue');
    stairs(X3, imag(Z3), 'green');
    grid on;
    xlabel('2Z/L');
    ylabel('|I(z)|');
    legend({'N=21','N=51','N=71'});
    title(['Q2- Parte Imaginária Corrente de destribuição de uma antena dipolo de L/lamda = ',num2str(L/lambda)]);
    hold off

    figure (5)
    stairs(X1, real(Z1), 'red');
    hold on
    stairs(X2, real(Z2), 'blue');
    stairs(X3, real(Z3), 'green');
    grid on;
    xlabel('2Z/L');
    ylabel('|I(z)|');
    legend({'N=21','N=51','N=71'});
    title(['Q2- Parte Real Corrente de destribuição de uma antena dipolo de L/lamda = ',num2str(L/lambda)]);
    hold off
    
    figure (6)
    stairs(X1, abs(Z1), 'red');
    hold on
    stairs(X2, abs(Z2), 'blue');
    stairs(X3, abs(Z3), 'green');
    grid on;
    xlabel('2Z/L');
    ylabel('|I(z)|');
    legend({'N=21','N=51','N=71'});
    title(['Q2- Valor absoluto Corrente de destribuição de uma antena dipolo de L/lamda = ',num2str(L/lambda)]);
    hold off

%-------------------------Questão 4----------------------------------------

elseif(Questao == 4)  

    clc;
    clear all;
    
    lambda = 1;
    L = [0.47,0.75,1.5]*lambda;
    N=51;
    a = 0.005*lambda;
    theta= 0.001:0.02:2*pi;
    Zl=1;
    E0=1;
    k0 = (2*pi)/lambda;

    h = zeros(N, 1);
    s= zeros(N,1);
    Voc_Teorico = zeros (length(theta),1);
    Voc = zeros (length(theta),1);

    for c =  1:length(L)
        for b = 1:length(theta)
            [M, X, I]= Fun(L(c), N, a, theta(b),Zl,1);
            Voc(b) = I((N-1)/2);
        
            s= @(x)(sin(k0*((L(c)/2)-abs(x)))*exp(2i*k0*x*cos(theta(b))));
            h(b)=  E0 * (((sin(theta(b)).^2)/ sin((k0*L(c))/2)));
            Voc_Teorico(b)=  h(b) * integral(@(x) s(x), -L(c)/2, L(c)/2,'ArrayValued', true);
        end
        
        figure(c)
        polarplot(theta, abs(Voc), 'blue');
        hold on;
        polarplot(theta, abs(Voc_Teorico), 'red');
        legend({'|Voc|','|Voc|Aproximação Sinusoidal'});
        title(['Q4 - Voc para L/lambda = ', num2str(L(c)/lambda)]);
        hold off;
        
    end
    
%-------------------------Questão 5----------------------------------------    

elseif(Questao == 5)
    
    clc;
    clear all;
    
    N=[21,51,71];
    lambda = 1;
    L = 0.47*lambda;
    a = 0.005*lambda;
    Zl=0;
    theta= pi/2;

    for c =  1:length(N)
        %I_zero = zeros(N(c),1)
        [M, X, I] = Fun(L, N(c), a,theta,Zl,0);
        I_zero=I((N(c)-1)/2);

        [M, X, I]= Fun(L, N(c), a, theta,Zl,1);
        Voc=I((N(c)-1)/2);

        Zin = Voc/I_zero;
        disp(Zin);
    end
    
%-------------------------Questão 6----------------------------------------

elseif(Questao == 6) 
    
    clc;
    clear all;
    
    Steps = 70;
    lambda = 1;
    a = 0.005*lambda;
    N = 51;
    Z = zeros(1,Steps);
    Ll= linspace(0.01,1.3,Steps);
    theta= pi/2;
    Zl=0;

    I_zero = zeros (Steps,1);
    Voc = zeros (Steps,1);
    Zin = zeros (Steps,1);

    for f = 1:Steps
        [M, X, I] = Fun(Ll(f), N, a,theta,Zl,0);
        I_zero(f)=I((N-1)/2);

        [M, X, I]= Fun(Ll(f), N, a, theta,Zl,1);
        Voc(f)=I((N-1)/2);

        Zin(f)=Voc(f)/I_zero(f);      
    end

    figure(1);
    plot(Ll, real(Zin),'blue');
    hold on;
    plot(Ll, imag(Zin), 'red');
    grid on
    xlabel('L/lambda');
    ylabel('Z [Ohm]');
    legend({'Resistencias','Reatancia'});
    title('Q6 - Impedância de entrada da antena dipolo para a = 0.005*lambda');


%-------------------------Questão 7----------------------------------------
elseif(Questao == 7)
    
    clc;
    clear all;
    
    Zl=73;
    E0=1;
    lambda = 1;
    k0 = (2*pi)/lambda;
    L=[0.47, 0.75];
    a = 0.005*lambda;
    N = 51;
    theta= 0.001:0.05:2*pi;
    delta = L/N;
    Z = zeros(1,length(theta));

    I_zero = zeros (length(theta),1);
    Voc = zeros (length(theta),1);
    Za = zeros (length(theta),1);
    V = zeros (length(theta),1);
    V_teorico = zeros (length(theta),1);
    h = zeros(length(theta), 1);
    s= zeros(length(theta),1);
    Voc_teorico = zeros (length(theta),1);
    Za_teorico = zeros (length(theta),1);

    for c = 1:length(L)
        for b = 1:length(theta)

        	[M, X, I] = Fun(L(c), N, a,theta(b),Zl,0);
            I_zero(b)=I((N-1)/2);

            [M, X, I]= Fun(L(c), N, a, theta(b),Zl,1);
            Voc(b)=I((N-1)/2);

            Za(b)=Voc(b)/I_zero(b);
    
            V(b)=(Zl/(Za(b)+Zl)*Voc(b));

            %------------------------

            s= @(x)(sin(k0*((L(c)/2)-abs(x)))*exp(2i*k0*x*cos(theta(b))));
            h(b)=  E0 * (((sin(theta(b)).^2)/ sin((k0*L(c))/2)));
            Voc_teorico(b)=  h(b) * integral(@(x) s(x), -L(c)/2, L(c)/2,'ArrayValued', true);
        
            Za_teorico(b)=Voc_teorico(b)/I_zero(b);
            V_teorico(b)= (Za_teorico(b)*I_zero(b))+Voc_teorico(b);
        end

    figure(c);
    polarplot(theta, abs(V),'blue');
    hold on;
    polarplot(theta, abs(V_teorico), 'red');
    legend({'|V|','|V| teórico'});
    title(['Q7 - Tensão induzida ao terminais da antena da antena dipolo para L = ' num2str(L(c)/lambda)]);
    hold  off
    end    
end
