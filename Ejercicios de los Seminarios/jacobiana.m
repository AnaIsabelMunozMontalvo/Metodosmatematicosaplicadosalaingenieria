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



function J=jacobiana(x,y)

%en este código se muestra a modo de ejemplo, como definir la jacobiana
%asociada a un sistema de ecuaciones no lineales que se pretende resolver con
%el código metnewtonsistema.m

J(1,1)=2.*x.*(x.^2+y^2).^(-1)-y.*cos(x.*y);
J(1,2)=2.*y.*(x.^2+y^2).^(-1)-x.*cos(x.*y);
J(2,1)=exp(x-y)-y.*sin(x.*y);
J(2,2)=-exp(x-y)-x.*sin(x.*y);
end