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



function [xh,uh]=bvpdirichlet(a,b,numeronodos,D,V,Q,f,ua,ub)
%este c�digo resuelve un problema de contorno unidimensional con 
%condiciones de contorno Dirichlet :%
%     -D u''+V*u'+Qu=f en (a,b)
%con condiciones de contorno u(a)=ua, u(b)=ub
%utilizando diferencias centradas
%numerosnodos es el n�mero de nodos inclu�dos los extremos del intervalo
%xh es el vector que contiene los nodos
%uh es la soluci�n num�rica

h = (b-a)/(numeronodos-1);%n�mero de intervalos=n�mero de nodos-1
xh = (linspace(a,b,numeronodos))';
nodosinteriores=numeronodos-2;
%construimos la matriz de coeficientes que resulta de aplicar el esquema
%centrado
hD = D/h^2;
hV = V/(2*h);
e =ones(nodosinteriores,1);
%matriz de coeficientes A
A = spdiags([-hD*e-hV (2*hD+Q)*e -hD*e+hV],...
    -1:1,nodosinteriores,nodosinteriores);
xi = xh(2:end-1);%valores de nodos interiores
%lado derecho del sistema ff
ff =feval(f,xi);
ff(1) =  ff(1)+ua*(hD+hV);
ff(end) = ff(end)+ub*(hD-hV);
uh = A\ff;
uh=[ua; uh; ub];
return
