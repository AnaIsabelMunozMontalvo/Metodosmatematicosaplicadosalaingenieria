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




function [soluciont,soluciony]=rungekuttao3(f,intiempo,valorini,npasos)
% este c�digo resuelve un problema de valor inicial y'=f(t,y), y(t0)=y0
%utlizando el m�todo Runge Kutta de orden 3
%f es la funci�n que define y', intiempo es el intervalo de tiempo donde se
%quiere resolver el problema, e valorini es el dato inicial
%npasos es el n�mero de pasos que hay que dar para llegar al tiempo final y
%depende del tama�o de discretizaci�n h
%soluciont y soluciony son los vectores soluci�n, e indican los tiempos
%intermedios considerados y los correspondientes valores aproximados de y

soluciont=linspace(intiempo(1),intiempo(2),npasos+1);
h=(intiempo(2)-intiempo(1))/npasos;
soluciony=ones(npasos+1:1);
soluciony(1)=valorini;

for i=2:npasos+1
    ti1=soluciont(i-1);
    ti2=soluciont(i-1)+0.5.*h;
    ti3=soluciont(i);
    
    yi1=soluciony(i-1);
    yi2=soluciony(i-1)+h.*f(ti1,yi1);
    yi3=soluciony(i-1)+h.*(-f(ti1,yi1)+2.*f(ti2,yi2));
    soluciony(i)=soluciony(i-1)+(1/6).*h.*(f(ti1,yi1)+4.*f(ti2,yi2)+...
        f(ti3,yi3));
end
