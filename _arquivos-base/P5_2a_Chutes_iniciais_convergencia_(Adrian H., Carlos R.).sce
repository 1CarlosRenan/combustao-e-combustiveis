clc
clear

 //Criado em 12 de novembro de 2018
 //Engenheiro responsável: Adrian Beppu Hirata
 //Engenheiro verificador: Carlos Renan Cândido da Silva
 
//PARTE 2a - P5 - Equilíbrio Químico
//Programa para encontrar o vetor de convergência i_2 para um dado intervalo de temperatura

//ARGUMENTOS DE ENTRADA
//Para o metano com CH4 (dados retirados da P1):
a = 1;                      //nº de átomos de carbono
b = 4;                     //nº de átomos de hidrogênio
c = 0;                      //nº de átomos de oxigênio
d = 0;                      //nº de átomos de nitrogênio
as = a + b/4 - c/2;         //coeficiente estequiométrico

//Condições estabelecidas para o metano:
phi = 0.8;                          //Razão de equivalência ACs/AC
T = linspace(1000,3500,35);         //Temperatura em K
P = 50*0.986923;                    //Pressão em atm

//Matriz com as constantes para K
//         A                 B               C               D               E
M = [0.432168*10^00  -0.112464*10^05  0.267269*10^01 -0.745744*10^-04  0.242484*10^-08
     0.310805*10^00  -0.129540*10^05  0.321779*10^01 -0.738336*10^-04  0.344645*10^-08
    -0.141784*10^00  -0.213308*10^04  0.853461*10^00  0.355015*10^-04 -0.310227*10^-08
     0.150879*10^-01 -0.470959*10^04  0.646096*10^00  0.272805*10^-05 -0.154444*10^-08
    -0.752364*10^00   0.124210*10^05 -0.260286*10^01  0.259556*10^-03 -0.162687*10^-07
    -0.415302*10^-02  0.148627*10^05 -0.475746*10^01  0.124699*10^-03 -0.900227*10^-08];

for k = 1:length(T);

//Equação para o cálculo da constante Ki de equilíbrio químico
    K = 10^(M(:,1).*log(T(k)/1000)+M(:,2)./T(k)+M(:,3)+M(:,4).*T(k)+M(:,5).*(T(k)^2));
    K_2 = [K.^2];
    
//Função para o sistema de equações não-lineares

for i_2=0:0.00001:1;  //intervalo possível para i

function x = g(y)

    x(1) = b/a*(y(1)+y(5))-(2*y(2)+2*y(6)+y(7)+y(9));
    x(2) = (c+2*as/phi)/a*(y(1)+y(5))-(2*y(1)+y(2)+2*y(4)+y(5)+y(8)+y(9)+y(10));
    x(3) = (d+7.52*as/phi)/a*(y(1)+y(5))-(2*y(3)+y(10));
    x(4) = y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)-1;
    x(5) = K_2(1)*y(6)-y(7)*y(7)*P;
    x(6) = K_2(2)*y(4)-y(8)*y(8)*P;
    x(7) = K_2(3)*y(4)*y(6)-y(9)*y(9);
    x(8) = K_2(4)*y(4)*y(3)-y(10)*y(10);
    x(9) = K_2(5)*y(4)*y(6)*y(6)*P-y(2)*y(2);
    x(10) = K_2(6)*y(4)*y(5)*y(5)*P-y(1)*y(1);

endfunction

//Comando para a solução do sistema de equações

[y fv in] = fsolve([i_2,i_2,i_2,i_2,i_2,i_2,i_2,i_2,i_2,i_2],g,10^-15)

y=y'

//CONDIÇÕES LÓGICAS DE EQUILÍBRIO QUÍMICO PARA A SELEÇÃO DE FRAÇÃO MOLAR
//Condições caso y seja negativo (o que não pode)
if  y(1)<0
    y(:,1)=[]
    elseif y(2)<0
    y(:,1)=[]
    elseif y(3)<0
    y(:,1)=[]
    elseif y(4)<0
    y(:,1)=[]
    elseif y(5)<0
    y(:,1)=[]
    elseif y(6)<0
    y(:,1)=[]
    elseif y(7)<0
    y(:,1)=[]
    elseif y(8)<0
    y(:,1)=[]
    elseif y(9)<0
    y(:,1)=[]
    elseif y(10)<0
    y(:,1)=[]
end

//Condições na convergência de fsolve e soma das frações próxima de 1
if in == 1 & 0.99<=sum(y) & sum(y)<=1.1;

T_i_2 = [T(k) i_2] //Relação entre cada temperatura com o valor de convergência

disp(T_i_2) //Mostra o vetor a ser copiado no Console do Scilab

break

end

end

end
