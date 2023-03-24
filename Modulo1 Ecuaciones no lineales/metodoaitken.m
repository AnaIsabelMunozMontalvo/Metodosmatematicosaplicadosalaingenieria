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



function [sol,itera]=metodoaitken(g,x0,errorper,maxitera)

%este c�digo utiliza el m�todo de Aitken para obtener una aproximaci�n
%num�rica al punto fjo de una funci�n g, utilizando como semilla x0.
%errorper es el m�ximo error permitido en la aproximaci�n
%el error se vs calculando como la diferencia entre los valores obtenidos
%en dos iteraciones consecutivas
%maxitera es el m�ximo n�mero de iteraciones que se permite realizar en el
%esqueme

x = x0;
difer = errorper + 1;
itera = 0;
while itera < maxiter & difer >= errorper
    gx = g(x);
    ggx = g(gx);
    xnew = (x*ggx-gx^2)/(ggx-2*gx+x);
    difer = abs(x-xnew);
    x = xnew;
    itera = itera  + 1;
end
if (itera==maxitera & difer>erroper)
    fprintf([' El m�todo no converge en el n�mero',...
            ' m�ximo de iteraciones fijado ']);
end
sol = x;
return
