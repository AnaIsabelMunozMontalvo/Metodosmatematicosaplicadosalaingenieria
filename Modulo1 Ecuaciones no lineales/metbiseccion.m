%____copyright___="Copyright (C) 2022 A. Nolla, A.I. Mu�oz, E. Schiavi."
%____license____="GPL-3.0-only"

%Detalles sobre el empleo de los c�digos y ejercicios de aplicaci�n pueden
%encontrarse en las direcciones de los siguientes documentos elaborados por
%A. Nolla, A.I. Mu�oz, E. Schiavi:
%https://burjcdigital.urjc.es/handle/10115/20132
%https://burjcdigital.urjc.es/handle/10115/20134
%As� como en el fichero Readme

% La mayor parte de los c�digos de la colecci�n presentada en esta librer�a
%son adaptaciones de los publicados en el libro "C�lculo cient�fico con 
% Matlab y Octave" de A. Quarteroni y F. Saliery, 
% que se pueden obtener en https://mox.polimi.it/qs/.



function [sol,itera]=metbiseccion(fecu,a,b,errorper,maxitera)

%   este c�digo encuentra una aproximaci�n de una ra�z en el intervalo
%   [a,b], de la ecuaci�n fecu=0. 
%   errorper es el error permitido
%   maxitera es el n�mero m�ximo de iteraciones que se permiten realizar
%   sol es la soluci�n num�rica aproximada obtenida
%   itera es el n�mero de iteraciones que ha realizado para obtener sol

x = [a, (a+b)*0.5, b];
fx = fecu(x);
%comprobamos que el intervalo elegido satisface el teorema de Bolzano
%fecu(a)*fecu(b)<0, si no lo verifica, aparece un mensaje de error.
%Tambi�n puede ocurrir que a o b sean ya raices
if fx(1)*fx(3) > 0
  error([' El signo de la funcion en los extremos',...
   ' del intervalo (a,b) tiene que ser distinto']);
elseif fx(1) == 0
    zero = a; res = 0; niter = 0; return
elseif fx(3) == 0
    zero = b; res = 0; niter = 0; return
end
itera = 0;
I = (b - a)*0.5;
while I >= errorper & itera < maxitera
   itera = itera + 1;
   if fx(1)*fx(2) <  0
      x(3) = x(2);
      x(2) = x(1)+(x(3)-x(1))*0.5;
      fx = fecu(x);
      I = (x(3)-x(1))*0.5;
   elseif fx(2)*fx(3) < 0
      x(1) = x(2);
      x(2) = x(1)+(x(3)-x(1))*0.5;
      fx = fecu(x);
      I = (x(3)-x(1))*0.5;
   else
       x(2) = x(find(fx==0)); I = 0;
   end
end
if  (itera==maxitera & I > errorper)
 fprintf(['El m�todo no converge',...
   'se ha llegado al n�mero m�ximo de iteraciones ',...
   'sin estar por debajo del error m�ximo permitido']);
end
sol = x(2); x = x(2);

