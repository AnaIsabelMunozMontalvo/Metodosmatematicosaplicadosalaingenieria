%____copyright___="Copyright (C) 2022 A. Nolla, A.I. Muñoz, E. Schiavi.
%____license____="GPL-3.0-only"

%Detalles sobre el empleo de los códigos y ejercicios de aplicación pueden
%encontrarse en las direcciones de los siguientes documentos elaborados por
%A. Nolla, A.I. Muñoz, E. Schiavi:
%https://burjcdigital.urjc.es/handle/10115/20132
%https://burjcdigital.urjc.es/handle/10115/20134
%Así como en el fichero Readme

% La mayor parte de los códigos de la colección presentada en esta librería
%son adaptaciones de los publicados en el libro "Cálculo científico con 
% Matlab y Octave" de A. Quarteroni y F. Saliery, 
% que se pueden obtener en https://mox.polimi.it/qs/.



function [sol,itera]=metnewton1ec(fecu,dfecu,x0,errorper,maxitera)

%   Este código obtiene una solución aproximada de una raíz
%   de la ecuación fecu=0, con el método de Newton-Raphson
%   tomando como semilla x0
%   fecu es la función que definie la ecuación de la cual queremos obtener
%   una raíz, dfecu es su derivada primera
%   sol es la solución nnumérica obtenida
%   errorper es el máximo error permitido en la aproximación
%   maxitera es el número máximo de iteraciones permitidas
%   itera es el número de iteracioes realizadas por el esquema para
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
  fprintf(['el método no converge por alcanzar el número máximo  ',...
   'de iteraciones permitidas ',...
   'sin llegar a un error por debajo del error máximo permitido']);
end
sol = x;
return
