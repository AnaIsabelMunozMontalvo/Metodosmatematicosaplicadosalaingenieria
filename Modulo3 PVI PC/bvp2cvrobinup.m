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




function [xh,uh]=bvp2cvrobinup(a,b,N,D,V,q,bvpfun,c11,c12,c21,c22,ua,ub,esquema,varargin)

%  BVP2CVROBINUP resuelve problemas de contorno 1-dim.
%  de la forma 
%     -D(x)*(d^2U/dX^2)+V(x)*dU/dX+q(x)*U=BVPFUN
%  en el intervalo (A,B) con las condiciones de contorno
%     c11*U'(A)+c12*U(A)=UA
%     c21*U'(B)+c22*U(B)=UB
%  mediante difirencias finitas en N nodos internos 
%  equiespaciado en (A,B). 
%  Parametros:
%     - BVPFUN, D, V, q  son funciones inline o anonimas
%     - esquema = 'C' se usan formulas centradas
%     - esquema = 'U' se usan formulas descentradas a 
%                     contracorriente ("upwind") para el 
%                     termino convectivo
%     - Varargin permite pasar paremetros adicionales a la
%                     funcion BVPFUN
%  Salida:
%     - XH contiene los nodos de la discretizacion
%     - UH contiene la solucion numerica


h = (b-a)/(N+1);
xh = (linspace(a,b,N+2))';
xi = xh(2:N+1);

if (c11 == 0 && c12 == 0)
  disp('Condicion de contorno en la frontera izquierda incorrecta')
  return
end 
if (c21 == 0 && c22 == 0)
  disp('Condicion de contorno en la frontera derecha incorrecta')
  return
end

%% Ecuaciones en la frontera
if c11 ~= 0
  TrIzq = [-D(a)/(h^2) 2*D(a)/(h^2)-V(a)*c12/c11+q(a) -D(a)/(h^2)];
  FlIzq = [-c11/(2*h) c12 c11/(2*h)];
end
if c21 ~= 0
  TrDcha = [-D(b)/(h^2) 2*D(b)/(h^2)-V(b)*c22/c21+q(b) -D(b)/(h^2)];
  FlDcha = [-c21/(2*h) c22 c21/(2*h)];    
end



if esquema == 'C'    % Esquema centrado

  nW = -D(xi)/(h^2)-V(xi)/(2*h);
  nC = 2*D(xi)/(h^2)+q(xi);
  nE = -D(xi)/(h^2)+V(xi)/(2*h);

  K = [nW nC nE];

elseif esquema == 'U'    % Esquema descentrado
%  disp('Esquema upwind para la parte convectiva')
  
  for i=2:N+1  % parametro rho
    if V(a+(i-1)*h) ~= 0
      rho(i-1) = 1/2+abs(V(a+(i-1)*h))./(2*V(a+(i-1)*h));
    else 
      rho(i-1) = 1/2;
    end
  end
  
  nW = -D(xi)/(h^2)-V(xi).*rho'/h;
  nC = 2*D(xi)/(h^2)+V(xi).*(2*rho'-1)/h+q(xi);
  nE = -D(xi)/(h^2)+V(xi).*(1-rho')/h;

  K = [nW nC nE];
  
else 
  disp('Error en el tipo de esquema. Las opciones son:')
  disp('C = centrado')
  disp('U = upwind')
  return
    end



% Construccion de las matrices para cada caso
if c11 ~= 0
%  disp('Cond. tipo Robin en Frontera izquierda')
  if c21 ~= 0
    %%% Caso Robin izq, Robin dcha
%    disp('Cond. tipo Robin en Frontera derecha')
    A = full([spdiags(FlIzq,[0 1 2],1,N+4);...
              spdiags(TrIzq,[0 1 2],1,N+4);...
              spdiags(K,[1 2 3],N,N+4);...
              spdiags(TrDcha,[N+1 N+2 N+3],1,N+4);...
              spdiags(FlDcha,[N+1 N+2 N+3],1,N+4)]);
    %B = [ua; bvpfun(a)-V(a)*ua/c11; feval(bvpfun,xi); bvpfun(b)-V(b)*ub/c21; ub];  %sin varargin
    B = [ua; bvpfun(a)-V(a)*ua/c11; feval(bvpfun,xi,varargin{:}); bvpfun(b)-V(b)*ub/c21; ub];  %con varargin
    uh = A\B;
    uh = [uh(2:N+3)];
  
  else 
    %%% Caso Robin izq, Dirichlet dcha
%    disp('Cond. tipo Dirichlet en Frontera derecha')
    cD = zeros(1,N+3); cD(1,end) = 1;
    A = full([spdiags(FlIzq,[0 1 2],1,N+3);...
              spdiags(TrIzq,[0 1 2],1,N+3);... 
              spdiags(K,[1 2 3],N,N+3); cD]);
    %B = [ua; bvpfun(a)-V(a)*ua/c11; feval(bvpfun,xi); ub/c22];  %sin varargin
    B = [ua; bvpfun(a)-V(a)*ua/c11; feval(bvpfun,xi,varargin{:}); ub/c22];  %con varargin
    uh = A\B;
    uh = [uh(2:N+2); ub];
  end
else 
%  disp('Cond. tipo Dirichlet en Frontera izquierda')
  if c21 ~= 0    
    %%% Caso Diriclet izq, Robin dcha
%    disp('Cond. tipo Robin en Frontera derecha')
    cI = zeros(1,N+3); cI(1,1) = 1;
    A = full([cI; spdiags(K,[0 1 2],N,N+3);...
               spdiags(TrDcha,[N N+1 N+2],1,N+3);...
               spdiags(FlDcha,[N N+1 N+2],1,N+3)]);
    %B = [ua/c12; feval(bvpfun,xi); bvpfun(b)-V(b)*ub/c21; ub];  %sin varargin
    B = [ua/c12; feval(bvpfun,xi,varargin{:}); bvpfun(b)-V(b)*ub/c21; ub];  %con varargin
    uh = A\B;
    uh = [ua; uh(2:N+2)];
  else 
    %%% Caso Diriclet izq, Dirichlet dcha
%    disp('Cond. tipo Dirichlet en Frontera derecha')
    cI = zeros(1,N+2); cI(1,1) = 1;
    cD = zeros(1,N+2); cD(1,end) = 1;
    A = full([cI; spdiags(K,[0 1 2],N,N+2); cD]);
    %B = [ua/c12; feval(bvpfun,xi); ub/c22];  %sin varargin
    B = [ua/c12; feval(bvpfun,xi,varargin{:}); ub/c22];  %con varargin
    uh = A\B;
    uh = [ua; uh(2:N+1); ub];
  end
end


return

