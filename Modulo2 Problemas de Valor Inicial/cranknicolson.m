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


function [soluciont,soluciony]=cranknicolson(f,intiempo,valorini,npasos)
% este c�digo resuelve un problema de valor inicial y'=f(t,y), y(t0)=y0
%utlizando el m�todo de crank nicolson
%f es la funci�n que define y', intiempo es el intervalo de tiempo donde se
%quiere resolver el problema, e valorini es el dato inicial
%npasos es el n�mero de pasos que hay que dar para llegar al tiempo final y
%depende del tama�o de discretizaci�n h
%soluciont y soluciony son los vectores soluci�n, e indican los tiempos
%intermedios considerados y los correspondientes valores aproximados de y
%utilizamos el comando fsolve para resolver la posible ecuaci�n no lineal
%en y que pueda aparecer al aplicar el esquema impl�cito

h=(intiempo(2)-intiempo(1))/npasos;
soluciont=(intiempo(1):h:intiempo(2));
soluciony=ones(npasos+1:1);
soluciony(1)=valorini; 
global glob_h glob_t glob_tt glob_y glob_f;
glob_h=h;
glob_y=valorini;
glob_f=f;

for i=2:npasos+1
    glob_t=soluciont(i);
    glob_tt=soluciont(i-1);
  w = fsolve(@(w) crank(w),glob_y);
  soluciony(i)=w;
  glob_y=w;
end   

clear glob_h glob_t glob_tt glob_y glob_odefun;
end
function [z]=crank(w)
  global glob_h glob_t glob_tt glob_y glob_f;
  z=w-glob_y-0.5.*glob_h*feval(glob_f,glob_t,w)- ...
      0.5.*glob_h*feval(glob_f,glob_tt,glob_y);
end


