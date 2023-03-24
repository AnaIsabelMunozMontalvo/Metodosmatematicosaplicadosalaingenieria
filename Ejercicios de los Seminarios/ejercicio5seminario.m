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

%Este código resuelve el ejercicio 5 de los seminarios que aparecen en:
%https://burjcdigital.urjc.es/handle/10115/20132
% o en el fichero de la carpeta de documentación



f=@(t,y) -2.*t.*y.^2;
intiempo=[0 1];
npasos=100;
valorini=1;
[tt,uheun]=heun(f,intiempo,valorini,npasos);
[tt,urk3]=rungekuttao3(f,intiempo,valorini,npasos);
plot(tt,uheun,'ro',tt,urk3,'b*')
t1=min(find(tt>=0.1));% también t1=find(tt>=0.1)
t2=min(find(tt>=0.5));%t2=find(tt>=0.5)
t3=min(find(tt>=1));%t3=find(tt>=1)
uheun(t1),uheun(t2),uheun(t3)
urk3(t1),urk3(t2),urk3(t3)