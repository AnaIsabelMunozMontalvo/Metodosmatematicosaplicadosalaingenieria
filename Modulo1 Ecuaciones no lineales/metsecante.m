%____copyright___="Copyright (C) 2022 A. Nolla, A.I. Muñoz, E. Schiavi."
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



function [sol,itera]=metsecante(fecu,a,b,errorper,maxitera)

%   este código encuentra una aproximación de una raíz de la ecuación 
% dada por fecu=0 utilizando el método de la secante
%sol es la solución numérica obtenida tomando como semillas a y b
% errorper es el error máximo permitido
% el error se mide sobre el valor absoluto de la diferencia de dos valores
% obtenidos en iteraciones consecutivas
% maxitera es el número máximo de iteraciones permitido
% itera es el número de iteraciones realizadas para obtener sol, 
% cumpliendo con el eror máximo permitido
fx0=feval(fecu,a);
fx1=feval(fecu,b);
x0=a;
x1=b;
r=fx1*(x1-x0)/(fx1-fx0);
x2=x1-r;
c=1;
while (abs(r)>errorper && c < maxitera)
    x0=x1;
    x1=x2;
    fx0=feval(fecu,x0);
    fx1=feval(fecu,x1);
    r=fx1*(x1-x0)/(fx1-fx0);
    x2=x1-r;   
    c=c+1;    
end
if c >= nmax
 fprintf(' No converge en el maximo',...
          'número de iteraciones\n ');
end
sol = x2; itera = c;
end

