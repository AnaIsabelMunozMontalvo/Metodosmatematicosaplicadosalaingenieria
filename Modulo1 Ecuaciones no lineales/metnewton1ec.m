%____copyright___="Copyright (C) 2022 A. Nolla, A.I. Mu�oz, E. Schiavi.
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



function [sol,itera]=metnewton1ec(fecu,dfecu,x0,errorper,maxitera)

%   Este c�digo obtiene una soluci�n aproximada de una ra�z
%   de la ecuaci�n fecu=0, con el m�todo de Newton-Raphson
%   tomando como semilla x0
%   fecu es la funci�n que definie la ecuaci�n de la cual queremos obtener
%   una ra�z, dfecu es su derivada primera
%   sol es la soluci�n nnum�rica obtenida
%   errorper es el m�ximo error permitido en la aproximaci�n
%   maxitera es el n�mero m�ximo de iteraciones permitidas
%   itera es el n�mero de iteracioes realizadas por el esquema para
%   alcanzar un error por debajo del error cometido. El error se mide con
%   el valor absoluto de fecu(x)/dfecu(x)

x = x0;
fx = fecu(x);
dfx = dfecu(x);
itera = 0; medidaerror = errorper+1;
while medidaerror >= errorper & itera < maxitera
   itera = itera + 1;      
   medidaerror = - fx/dfx;
   x = x + medidaerror;
   medidaerror= abs(medidaerror);
   fx = fecu(x);
   dfx = dfecu(x);
end
if (itera==maxitera & medidaerror > errorper)
  fprintf(['el m�todo no converge por alcanzar el n�mero m�ximo  ',...
   'de iteraciones permitidas ',...
   'sin llegar a un error por debajo del error m�ximo permitido']);
end
sol = x;
return
