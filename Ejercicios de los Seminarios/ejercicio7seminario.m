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

%Este código resuelve el ejercicio 7 de los seminarios que aparecen en:
%https://burjcdigital.urjc.es/handle/10115/20132
% o en el fichero de la carpeta de documentación



intespacio=[0 1];
intiempo=[0 1];
pasosespacio=20;pasostiempo=100;
theta=0;
c=1;
u0=@(x) 0.*x;
cc=@(t,x) (x.^2).*(x-1);
f=@(t,x) 2.*t-x.^2+10.*x.*t;
[x,u]=ecucalor(c,intespacio,intiempo,pasosespacio, ...
pasostiempo,theta,u0,cc,f);
figure;
plot(x,u)
x1=min(find(x>=0.85))
u(x1)
dx=0.05;
dt=(0.05)^2;
pasost=1/(0.05)^2;% no sale, con 800 pasos si 0.0013
pasosespacio=20;pasostiempo=800;
[x,u1]=ecucalor(c,intespacio,intiempo,pasosespacio, ...
pasostiempo,theta,u0,cc,f);
figure;
plot(x,u1)
x1=min(find(x>=0.85))
u1(x1)



