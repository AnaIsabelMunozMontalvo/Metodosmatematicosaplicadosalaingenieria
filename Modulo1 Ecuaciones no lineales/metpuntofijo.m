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



function [sol,itera] = metpuntofijo(g,x0,errorper,maxitera)

% este método iterativo encuentra una aproximación del punto fijo 
% de la función g, sol
%sol es la solución numérica obtenida tomando como semilla x0
% nótese que f(x)=0 <=> g(x)=x
% errorper es el error máximo permitido
% el error se mide sobre el valor absoluto de la diferencia de dos valores
% obtenidos en iteraciones consecutivas
% maxitera es el número máximo de iteraciones permitido
% el valor absoluto de la diferencia entre los valores obtenidos en dos
% iteraciones consecutivas sea menor que tol
% itera es el número de iteraciones realizado para obtener sol, 
% cumpliendo con el eror máximo permitido

 x = x0; diferencia= errorper +1; itera =0;
 while itera <= maxitera & diferencia >= errorper
 gx = feval(g,x); %x(i+1)=g(x(i)), xnueva=g(xanterior)
 diferenciait=x-gx;%x(i+1)-x(i)
 diferencia=abs(diferenciait);
 xnueva=gx;
 x=xnueva;
 itera=itera+1; 
 end
 if itera >= maxitera
 fprintf('Se ha llegado al número máximo de iteraciones',...
          'sin estar por debajo del máximo error permitido ');
 end
 sol = x;
 return
endfunction