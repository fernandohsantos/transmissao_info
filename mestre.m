function [x1ponto] = mestre(parametros,xmestre)

%x(1) eixo x
%x(2) eixo y 
%x(3) eixo z

%Parametros(1) a
%Parametros(2) b
%Parametros(3) c

x1ponto(1) = -(xmestre(2)+xmestre(3));
x1ponto(2) = xmestre(1)+parametros(1)*xmestre(2);
x1ponto(3) = parametros(2)+xmestre(3)*(xmestre(1)-parametros(3));
