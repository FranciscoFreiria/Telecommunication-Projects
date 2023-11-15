% Parte 2 Projeto

% Ponto 1.1

zeta=2;

C=6;
eta6=((1-zeta*C).^2+zeta.^2).^0.5;

C=-6;
etamenos6=((1-zeta*C).^2+zeta.^2).^0.5;

tau=linspace(-30,30,100);
quociente1=(exp(-(tau/eta6).^2))/eta6;

figure()
plot(tau, quociente1)

grid on
hold all

quociente2=(exp(-(tau/etamenos6).^2))/etamenos6;

plot(tau, quociente2)
title({'Profile of the chirped Gaussian pulse in the single-mode fiber'; 'for different values of the chirp parameter at the input and at zeta=2'});
xlabel('\tau');
ylabel('|A|^2 / A0');
legend('C=6', 'C=-6');

% Ponto 1.2

zeta = linspace(0, 2, 100);
C=6;
quociente3 = 1./((1-zeta*C).^2+zeta.^2).^0.5;
C=-6;
quociente4 = 1./((1-zeta*C).^2+zeta.^2).^0.5;

figure()
plot(zeta, quociente3)


grid on
hold all

plot(zeta, quociente4)
title({'Evolution along the fiber of the amplitude of the chirped Gaussian pulse'; 'in the single-mode fiber for different vales of the initial chirp parameter'});
xlabel('\zeta');
ylabel('|A1 / A0');
legend('C=6', 'C=-6');


% Ponto 1.3

C=6;
eta6novo=((1-zeta*C).^2+zeta.^2).^0.5;
C16=C-(1+C.^2)*zeta;

C=-6;
etamenos6novo=((1-zeta*C).^2+zeta.^2).^0.5;
C1menos6=C-(1+C.^2)*zeta;

f = figure;
p = uipanel('Parent',f,'BorderType','none'); 
subplot(1,2,1,'Parent',p) 
plot(zeta,eta6novo)
grid on
hold all
plot(zeta,etamenos6novo)
title({'Evolution along the fiber, for \beta2 < 0, of the broadening coefficient'});
xlabel('\zeta');
ylabel('\eta');
legend('C=6', 'C=-6');
subplot(1,2,2,'Parent',p) 
plot(zeta, C16)
hold all
grid on
plot(zeta, C1menos6)
title({'Evolution along the fiber, for \beta2 < 0, of the the total chirp parameter'});
xlabel('\zeta');
ylabel('C1');
legend('C=6', 'C=-6');

% Ponto 2

C=-6;
L=40;
beta2=-1;
beta3=0.1;

p=0.25*(1+C.^2);
aquadrado=beta3.^2./(L*abs(beta2).^3);

xpequeno = linspace(0, 15, 100);
xgrande =  xpequeno - C + p*(1+aquadrado*p./(2*xpequeno))./xpequeno;

figure()
plot(xpequeno, xgrande)
ylim([10,30]);
grid on
title({'Variation of the normalized output pulse width X'; 'with the normalized input pulse width x'});
xlabel('x');
ylabel('X');
legend('a^2 = 2.5*10^-4');


% Ponto 3

C=linspace(-10,10,100);
p=0.25*(1+C.^2);

xzero=2.*((p/3).^0.5).*cos(acos(4.5.*aquadrado.*((p/3).^0.5))./3);
xgrande =  xzero - C + p.*(1+aquadrado.*p./(2.*xzero))./xzero;

mu = (xgrande./xzero).^0.5;

figure()
hold all
plot(C, mu)

grid on
title({'Variation of the optimum broadening coefficient \mu'; 'with the chirp parameter for \beta2 = -1 and a^2=2.5*10^-4'});
xlabel('C');
ylabel('\mu');
legend('a^2=2.5*10^-4');


% Ponto 4

T=50;
C=-6;
L=40;
beta2=-1;
beta3=0.1;

t=linspace(-250,250,1000);
A1 = sech(t/T).*exp(-1i.*C.*((t/T).^2)./2);
A2 = exp(-(1+i.*C).*((t/T).^8)./2);

f = figure;
p = uipanel('Parent',f,'BorderType','none'); 
subplot(1,2,1,'Parent',p) 
plot(t,A1)
grid on
hold all
subplot(1,2,2,'Parent',p) 
plot(t, A2)
%xlim(-30,30);
hold all
grid on