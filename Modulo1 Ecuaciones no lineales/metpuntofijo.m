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



function [sol,itera] = metpuntofijo(g,x0,errorper,maxitera)

% este m�todo iterativo encuentra una aproximaci�n del punto fijo 
% de la funci�n g, sol
%sol es la soluci�n num�rica obtenida tomando como semilla x0
% n�tese que f(x)=0 <=> g(x)=x
% errorper es el error m�ximo permitido
% el error se mide sobre el valor absoluto de la diferencia de dos valores
% obtenidos en iteraciones consecutivas
% maxitera es el n�mero m�ximo de iteraciones permitido
% el valor absoluto de la diferencia entre los valores obtenidos en dos
% iteraciones consecutivas sea menor que tol
% itera es el n�mero de iteraciones realizado para obtener sol, 
% cumpliendo con el eror m�ximo permitido

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
 fprintf('Se ha llegado al n�mero m�ximo de iteraciones',...
          'sin estar por debajo del m�ximo error permitido ');
 end
 sol = x;
 return
endfunction