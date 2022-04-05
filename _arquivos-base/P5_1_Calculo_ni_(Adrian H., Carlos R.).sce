clc
clear

 //Criado em 10 de novembro de 2018
 //Engenheiro responsável: Adrian Beppu Hirata
 //Engenheiro verificador: Carlos Renan Cândido da Silva
 
//PARTE 1 - P5 - Equilíbrio Químico
//Programa para calcular as composições ni de cada um dos seguintes produtos de combustão:
//CO2, H2O, N2, O2, CO, H2, H, O, OH e NO, respectivamente, a uma dada temperatura, pressão e
//razão de equivalência phi.

//ARGUMENTOS DE ENTRADA
//Para o exemplo da gasolina com C7H17 (dados retirados da P1):
a = 7;                              //nº de átomos de carbono
b = 17;                             //nº de átomos de hidrogênio
c = 0;                              //nº de átomos de oxigênio
d = 0;                              //nº de átomos de nitrogênio
as = a + b/4 - c/2;                 //coeficiente estequiométrico

//Condições estabelecidas para a gasolina (dados retirados da P4):
phi = 0.8;                    //Razão de equivalência ACs/AC
T = 3000;                     //Temperatura em K
P = 50*0.986923;              //Pressão em atm

//Matriz com as constantes para K
//         A                 B               C               D               E
M = [0.432168*10^00  -0.112464*10^05  0.267269*10^01 -0.745744*10^-04  0.242484*10^-08
     0.310805*10^00  -0.129540*10^05  0.321779*10^01 -0.738336*10^-04  0.344645*10^-08
    -0.141784*10^00  -0.213308*10^04  0.853461*10^00  0.355015*10^-04 -0.310227*10^-08
     0.150879*10^-01 -0.470959*10^04  0.646096*10^00  0.272805*10^-05 -0.154444*10^-08
    -0.752364*10^00   0.124210*10^05 -0.260286*10^01  0.259556*10^-03 -0.162687*10^-07
    -0.415302*10^-02  0.148627*10^05 -0.475746*10^01  0.124699*10^-03 -0.900227*10^-08];
  
//Equação para o cálculo da constante Ki de equilíbrio químico
    K = 10^(M(:,1).*log(T/1000)+M(:,2)./T+M(:,3)+M(:,4).*T+M(:,5).*(T^2));
    K_2 = [K.^2];

//SEÇÃO PARA ENCONTRAR VALOR DE CHUTES INICIAIS (i) QUE FAZEM O PROBLEMA CONVERGIR

for i=0:0.001:1;       //intervalo para os chutes iniciais

//Função para o sistema de equações não-lineares
function x = g(y);

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

//COMANDO PARA A SOLUÇÃO DO SISTEMA DE EQUAÇÕES
//y: valores de fração molar; fp: resultados ao substituir y nas equações;
// in: é o indicador de convergência, converge quando in = 1;

[y fp in] = fsolve([i,i,i,i,i,i,i,i,i,i],g,10^-15)

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

//Condições na convergência de fsolve e soma das frações = 1
if in == 1 & sum(y)==1;

//Cálculo do nº total de mols da mistura N
    N = a/(y(1)+y(5));  //[mols] 

//Cálculo da composição de cada um dos produtos de combustão
//CO2, H2O, N2, O2, CO, H2, H, O, OH e NO, respectivamente
    ni = N.*y;          //[mols] 

disp([ 'Fração molar' 'Composição molar [mol]'])
disp([       y                    ni])

break

end

end
