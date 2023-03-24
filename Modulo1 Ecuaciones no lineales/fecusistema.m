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



function F=fecusistema(x,y)

%en este código se muestra a modo de ejemplo, como definir la función
%vectorial asociada a un sistema de ecuaciones no lineales que se pretende 
% resolver con el código metnewtonsistema.m

F(1,1)=log(x.^2+y.^2)-sin(x.*y)-(log(2)+log(pi));
F(2,1)=exp(x-y)+cos(x.*y);
end 