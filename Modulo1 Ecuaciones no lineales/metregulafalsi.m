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



function [sol,itera]=metregulafalsi(fecu,a,b,errorper,maxitera)

%   este c�digo encuentra una aproximaci�n de una ra�z de la ecuaci�n 
% dada por fecu=0 utilizando el m�todo de regulafalsi
%sol es la soluci�n num�rica obtenida tomando como semillas a y b
% errorper es el error m�ximo permitido
% el error se mide sobre el valor absoluto de la diferencia de dos valores
% obtenidos en iteraciones consecutivas
% maxitera es el n�mero m�ximo de iteraciones permitido
% itera es el n�mero de iteraciones realizadas para obtener sol, 
% cumpliendo con el eror m�ximo permitido
fx0=feval(fecu,a);
fx1=feval(fecu,b);
x0=a;
x1=b;
r=fx1*(x1-x0)/(fx1-fx0);
x2=x1-r;
c=1;
fx2=feval(fecu,x2);
while (abs(r)>errorper && c < maxitera)
    if fx2*fx1<0
        x0=x2;
            else
        x1=x2;
    end
    fx0=feval(fecu,x0);
    fx1=feval(fecu,x1);
    r=fx1*(x1-x0)/(fx1-fx0);
    x2=x1-r; 
    fx2=feval(fecu,x2);
    c=c+1;    
end
if c >= maxitera
 fprintf(' No converge en el maximo',...
          'n�mero de iteraciones\n ');
end
sol = x2; itera = c;
end

