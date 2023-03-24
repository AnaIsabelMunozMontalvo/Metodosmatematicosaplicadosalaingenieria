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


function [soluciont,soluciony]=eulerexplicito(f,intiempo,valorini,npasos)
%este código resuelve el problema de valor inicial dado por la ecuación
%diferencial y' = f(t,y) y el valor inicial, y(t0) denotado por valorini, 
%en el intervalo temporal [t0,T), definido en intiempo 
%con el método de Euler explícito
%npasos es el número de pasos que hay que dar, tomando un tamaño de 
%discretización constante h, para llegar a T, es decir, (T-t0)/h
%soluciont, es el vector que contiene los nodos temporales, t0,t1, ...tn=T
%soluciony, es el vector que contiene los valores numéricos aproximados de
%y(t), obtenidos para los distintos nodos temporales

h=(intiempo(2)-intiempo(1))/npasos;
%inicializamos el vector que albergará la solución
soluciont=ones(npasos+1,1);
soluciony=ones(npasos+1,1);
soluciont(1)=intiempo(1);%t0
soluciony(1)=valorini;%y(t0)
for i=2:npasos+1
 soluciont(i)=soluciont(i-1)+h;
 soluciony(i)=soluciony(i-1)+h*f(soluciont(i-1),soluciony(i-1));
 
end

return
