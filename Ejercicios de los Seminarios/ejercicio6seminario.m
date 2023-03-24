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

%Este código resuelve el ejercicio 6 de los seminarios que aparecen en:
%https://burjcdigital.urjc.es/handle/10115/20132
% o en el fichero de la carpeta de documentación



a=0;
b=2;
D=1;
V=-4;
q=0;
f=@(x) -16.*x.^3+34.*x-1;
ua=4;
ub=2;
numeronodos=17;
[xh,uh]=bvpdirichlet(a,b,numeronodos,D,V,q,f,ua,ub);
solexac=xh.^4-xh.^3-3.5.*xh.^2+2.*xh+4;
figure;
plot(xh,uh,'r',xh,solexac,'g')
nmax=max(uh)
nmin=min(uh)
emin=min(solexac)
emax=max(solexac)
figure;
err=uh-solexac;
plot(xh,err)
xm=find(uh>=nmax)
xmin=find(uh<=nmin)
xh(xm),xh(xmin)
