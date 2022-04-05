clc
clear

 //Criado em 12 de novembro de 2018
 //Engenheiro responsável: Adrian Beppu Hirata
 //Engenheiro verificador: Carlos Renan Cândido da Silva
 
//PARTE 2b - P5 - Equilíbrio Químico 
//Programa para gerar o gráfico de fração molar y x Temperatura T [K] para as curvas de:
//CO2, H2O, N2, O2, CO, H2, H, O, OH e NO, em um intervalo de temperatura (1000 a 3500 K),
//a pressão de 50 bar e a uma dada razão de equivalência phi para o metano, utilizando
//o vetor de convergência i_2 encontrado em 2a

//ARGUMENTOS DE ENTRADA
//Para o metano com CH4 (dados retirados da P1):
a = 1;                //nº de átomos de carbono
b = 4;               //nº de átomos de hidrogênio
c = 0;                //nº de átomos de oxigênio
d = 0;                //nº de átomos de nitrogênio

//Vetor de convergência i_2 encontrado na parte 2a
    i_2 =[0.00805 0.0149 0.00448 0.00444 0.00455 0.00452 0.00455 0.00931 0.00105 0.00135 0.00138 0.00182 0.00258 1.9E-4 3.0E-4 0.0448 0.00239 2.0E-4 4.2E-4 0.00184 2.0E-4 8.0E-5 0.1245 2.0E-4 2.3E-4 4.7E-4 3.1E-4 2.0E-5 2.6E-4 9.0E-5 4.7E-4 1.1E-4 0.0 2.0E-5 6.0E-5];
    
as = a+b/4-c/2;                             //Coeficiente estequiométrico
phi = 0.8;                                  //Razão de equivalência ACs/AC
T = linspace(1000,3500,35);                 //Intervalo de temperatura em K
P = 50*0.986923;                            //Pressão em atm

//Matriz M com as constantes para Ki
//         A                 B               C               D               E
M = [0.432168*10^00  -0.112464*10^05  0.267269*10^01 -0.745744*10^-04  0.242484*10^-08
     0.310805*10^00  -0.129540*10^05  0.321779*10^01 -0.738336*10^-04  0.344645*10^-08
    -0.141784*10^00  -0.213308*10^04  0.853461*10^00  0.355015*10^-04 -0.310227*10^-08
     0.150879*10^-01 -0.470959*10^04  0.646096*10^00  0.272805*10^-05 -0.154444*10^-08
    -0.752364*10^00   0.124210*10^05 -0.260286*10^01  0.259556*10^-03 -0.162687*10^-07
    -0.415302*10^-02  0.148627*10^05 -0.475746*10^01  0.124699*10^-03 -0.900227*10^-08];

//Função para o sistema de equações não-lineares
    
function[x] = g(y)
   
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

for p=1:length(T);

//Equação para o cálculo da constante Ki de equilíbrio químico 
    K = 10^(M(:,1).*log(T(p)/1000)+M(:,2)./T(p)+M(:,3)+M(:,4).*T(p)+M(:,5).*T(p).^2);
    K_2 = [K.^2];
//Comando para a solução do sistema de equações
   [y(p,:)] = [fsolve([i_2(p),i_2(p),i_2(p),i_2(p),i_2(p),i_2(p),i_2(p),i_2(p),i_2(p),i_2(p)],g)];
end

//Gráfíco de fração molar y por Temperatura T [K] em escala logarítimica no eixo de y

figure()
T=T'                                  //Intervalo da temperatura no eixo de x
plot2d1("onl",T,y)                    //Plot do gráfico
zoom_rect([1000,10^-3,3500,10^0])     //Zoom no gráfico

//Legendas para as curvas e eixos
h1 = legend(['CO2' 'H2O' 'N2' 'O2' 'CO' 'H2' 'H' 'O' 'OH' 'NO'],-4); 
xlabel("Temperatura T [K]","fontsize",4)
ylabel("Fração molar y","fontsize",4)

//Título do gráfico
titulo = ["Curvas para phi = 0.8"];
title(titulo,"fontsize",5) 

