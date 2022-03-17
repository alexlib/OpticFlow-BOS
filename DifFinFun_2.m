function [Dm]=DifFinFun_2(Matx,Maty,m,a)

Matxt=Matx'; %Transpone las matrices del desplazamiento en x
Matyt=Maty'; %Transpone las matrices del desplazamiento en y

%Calcula las derivadas de Matx y Maty utilizando Diferencias Finitas

%% Realiza la operación (-/+)(Ex(i,j) - Ex(i+/-1,j)/dx), el signo depende 
%%de si está en el último cuadro(abajo) o no(arriba)
p=1;%(m/a);
% od=2;


for i=1+3*p:size(Matx,2)-3*p;
    for j=1+3*p:size(Matx,1)-3*p;
%         if i == size(Matx,2);
%              DvMatxt(i,j)=(Matxt(i,j)-Matxt(i-1,j))/(m/a);
%         else 
%              DvMatxt(i,j)=(Matxt(i+1,j)-Matxt(i,j))/(m/a);
%              DvMatxt(i,j)=(-1.5*Matxt(i,j)+2.0*Matxt(i+p,j)-0.5*Matxt(i+2*p,j))/(m/a);
               DvMatxt(i,j)=(-(1/60)*Matxt(i-3*p,j)+(3/20)*Matxt(i-2*p,j)-0.75*Matxt(i-p,j)+0.75*Matxt(i+p,j)-(3/20)*Matxt(i+2*p,j)+(1/60)*Matxt(i+3*p,j))/(m/a);
%         end
    end
end


%% Realiza la operación (-/+)(Ey(i,j) - Ey(i+/-1,j)/dy), el signo depende 
%%de si está en el último cuadro(abajo) o no(arriba)
for j=1+3*p:size(Maty,1)-3*p;
    for i=1+3*p:size(Maty,2)-3*p;
%         if j == size(Maty,1);
%              DvMatyt(i,j)=(Matyt(i,j)-Matyt(i,j-1))/(m/a);
%         else
%              DvMatyt(i,j)=(Matyt(i,j+1)-Matyt(i,j))/(m/a);
%              DvMatyt(i,j)=(-1.5*Matyt(i,j)+2.0*Matyt(i+1,j)-0.5*Matyt(i+2*p,j))/(m/a);
               DvMatyt(i,j)=(-(1/60)*Matyt(i,j-3*p)+(3/20)*Matyt(i,j-2*p)-0.75*Matyt(i,j-p)+0.75*Matyt(i,j+p)-(3/20)*Matyt(i,j+2*p)+(1/60)*Matyt(i,j+3*p))/(m/a);
%         end
    end
end

Dmx=DvMatxt'; %Matriz que tiene la diferencia finita de la derivada de la componente x del desplazamiento, respecto de x
Dmy=DvMatyt'; %Matriz que tiene la diferencia finita de la derivada de la componente y del desplazamiento, respecto de y

Dm=Dmx+Dmy; %Matriz de la divergencia de desplazamientos
end