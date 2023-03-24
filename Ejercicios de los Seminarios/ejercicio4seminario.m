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

%Este c�digo resuelve el ejercicio 4 de los seminarios que aparecen en:
%https://burjcdigital.urjc.es/handle/10115/20132
% o en el fichero de la carpeta de documentaci�n




f=inline('y-sin(t)+cos(t)','t','y');
intiempo=[0 2];valorini=1;npasos=20;
[t,u1]=eulerexplicito(f,intiempo,valorini,npasos);
[t,u2]=eulerimplicito(f,intiempo,valorini,npasos);
[t,u3]=cranknicolson(f,intiempo,valorini,npasos);
solexac=exp(t)+sin(t);
figure;
plot(t,u1,'ro',t,u2,'b+',t,u3,'g*',t,solexac,'k^')
t1=min(find(t>=0.5));% tambi�n t1=find(t==0.5);
t2=min(find(t>=1.5));% tambi�n t2=find(t==1.5);
u1(t1),u2(t1),u3(t1)
u1(t2),u2(t2),u3(t2)

