function [x2ponto] = escravo(parametros,sinal,xescravo)

%x(1) eixo x
%x(2) eixo y 
%x(3) eixo z

%Parametros(1) a
%Parametros(2) b
%Parametros(3) c
%Parametros(4) eta

x2ponto(1) = -(xescravo(2)+xescravo(3))+parametros(4)*(sinal-xescravo(1));
x2ponto(2) = xescravo(1)+parametros(1)*xescravo(2);
x2ponto(3) = parametros(2)+xescravo(3)*(xescravo(1)-parametros(3));
