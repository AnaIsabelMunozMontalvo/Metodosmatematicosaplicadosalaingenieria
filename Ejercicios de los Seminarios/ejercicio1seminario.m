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

%Este código resuelve el ejercicio 1 de los seminarios que aparecen en:
%https://burjcdigital.urjc.es/handle/10115/20132
% o en el fichero de la carpeta de documentación




fecu=@(x) x.^3+4.*x.^2-10;
a=0;b=2;maxitera=100;
fecu(a)
fecu(b)
errorper=1e-3;
[sol1,itera1]=metbiseccion(fecu,a,b,errorper,maxitera)
dfecu=inline('3.*x.^2+8.*x','x');
x0=1;
[sol2,itera2]=metnewton1ec(fecu,dfecu,x0,errorper,maxitera)

