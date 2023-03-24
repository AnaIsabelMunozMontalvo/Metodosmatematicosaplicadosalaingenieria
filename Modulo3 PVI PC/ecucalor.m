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



function [x,uf]=ecucalor(C,intespacio,intiempo,pasosespacio, ...
               pasostiempo,theta,u0,cc,f)
%este c�digo resuelve la ecuaci�n del calor en una dimesi�n espacial
%utilizando un theta-m�todo para la discretizaci�n temporal
%du/dt=C (d^2 u/dx^2) +f  (t0,T)x(a,b)
%intespacio=(a,b) intiempo=(t0,T)
%la condici�n inicial se deota por u0
%la condici�n de contorno tipo Dirichlet, la denotamos por cc
%theta es el valor de theta empleado en el theta m�todo
%theta=0 es para Euler expl�cito
%theta=1 es para Euler impl�cito
%theta=0.5 para Crank-Nicolson
%pasosespacio es el n�mero de pasos en espacio
%pasostiempo es el n�mero de pasos en tiempo
%x es el vector que contiene los nodos en espacio
%uf contiene los valores num�ricos obtenidos en t=T

h  = (intespacio(2)-intespacio(1))/pasosespacio;
dt = (intiempo(2)-intiempo(1))/pasostiempo;
N = pasosespacio+1;%n�mero de nodos en espacio
%perparamos la matriz de coeficientes, que variar� dependiendo del
%theta-m�todo empleado
e =ones(N,1); D = spdiags([-e 2*e -e],[-1,0,1],N,N);
I = speye(N); M = I+C*dt*theta*D/h^2;
Mn = I-C*dt*(1-theta)*D/h^2;
 M(1,:) = 0;
 M(1,1) = 1; M(N,:) = 0; M(N,N) = 1;
 x =linspace(intespacio(1),intespacio(2),N);
 x = x'; fn = feval(f,intiempo(1),x); 
 un =feval(u0,x);
 [L,U]=lu(M); 
 for t = intiempo(1)+dt:dt:intiempo(2)
    fn1 = feval(f,t,x);
    %lado derecho
    rhs = Mn*un+dt*(theta*fn1+(1-theta)*fn);
    temp = feval(cc,t,[intespacio(1),intespacio(2)]);
    rhs([1,N]) = temp;
    %el sistema se resuleve utilizando una factorizaci�n LU
    u = L\rhs;
    u = U\u;
    fn = fn1;
    un = u;
 
 end
uf=u;
return
