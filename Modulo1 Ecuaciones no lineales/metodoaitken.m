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



function [sol,itera]=metodoaitken(g,x0,errorper,maxitera)

%este código utiliza el método de Aitken para obtener una aproximación
%numérica al punto fjo de una función g, utilizando como semilla x0.
%errorper es el máximo error permitido en la aproximación
%el error se vs calculando como la diferencia entre los valores obtenidos
%en dos iteraciones consecutivas
%maxitera es el máximo número de iteraciones que se permite realizar en el
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
    fprintf([' El método no converge en el número',...
            ' máximo de iteraciones fijado ']);
end
sol = x;
return
