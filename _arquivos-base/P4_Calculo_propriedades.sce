clc
clear

 //Criado em 25 de setembro de 2018
 //Engenheiro responsável: Adrian Beppu Hirata
 //Engenheiro verificador: Carlos Renan Cândido da Silva

//Programa para estimar as proprioedades termodinâmicas cp [kJ/kg/K], s [kJ/kg/K], 
//h [kJ/kg], v [m3/kg], u [kJ/kg] em função da temperatura T [K], pressão e composição (ni) 
//da mistura dos produtos de combustão

T = 3000;                       //Temperatura [K]
P = 5000;                       //Pressão [kPa]

ni = [0.0775;0.1064;0.7203;0.0359;0.0191;0.0036;0.0013;0.003;0.0133;0.0196]; 
N = sum(ni);                        //Número de mols da mistura
yi = [ni./N];                       //Fração molar
m = [44;18;28;32;28;2;1;16;17;30];  //Massa molar de cada produto [kg/kmol]
M = sum(yi.*m);                     //Massa molar da mistura [kg/kmol]
Ru = 8.314;                    //Constante dos gases ideais [kJ/kmol/K]
R = Ru/M;                      //Constante dos gases ideais para a mistura [kJ/kg/K]
pi = yi.*P;                    //Pressão parcial dos produtos da mistura [kPa]
po = 101.325;                  //Pressão atmosférica [kPa]

//Matriz com as constantes a1,a2,3,a4,a5,a6 e a7 para estimar as propriedades termodinâmcias
a = [0.446080*10^1 0.309817*10^-2 -0.123925*10^-5 0.227413*10^-9  -0.155259*10^-13 -0.489614*10^5 -0.986359
0.271676*10^1 0.294513*10^-2 -0.802243*10^-6 0.102266*10^-9  -0.484721*10^-14 -0.299058*10^5 0.663056*10^1
0.289631*10^1 0.151548*10^-2 -0.572352*10^-6 0.998073*10^-10 -0.652235*10^-14 -0.905861*10^3 0.616151*10^1
0.362195*10^1 0.736182*10^-3 -0.196522*10^-6 0.362015*10^-10 -0.289456*10^-14 -0.120198*10^4 0.361509*10^1
0.298406*10^1 0.148913*10^-2 -0.578996*10^-6 0.103645*10^-9  -0.693535*10^-14 -0.142452*10^5 0.634791*10^1
0.310019*10^1 0.511194*10^-3  0.526442*10^-7 -0.349099*10^-10 0.369453*10^-14 -0.877380*10^3 -0.196294*10^1
    0.25*10^1        0                  0           0                0          0.254716*10^5 -0.460117
0.254205*10^1 -0.275506*10^-4 -0.310280*10^-8 0.455106*10^-11 -0.436805*10^-15 0.292308*10^5 0.492030*10^1
0.291064*10^1 0.959316*10^-3 -0.194417*10^-6 0.137566*10^-10 0.142245*10^-15 0.393538*10^4 0.544234*10^1
  0.3189*10^1 0.133822*10^-2 -0.528993*10^-6 0.959193*10^-10 -0.648479*10^-14 0.982832*10^4 0.674581*10^1];

//Cálculo do Calor Específico à pressão constante
cp_bar = Ru*(a(:,1) + a(:,2).*T + a(:,3).*T^2 + a(:,4).*T^3 + a(:,5).*T^4);
cp = sum(cp_bar.*yi)/M;

//Cálculo da Entalpia
h_bar = Ru*T*(a(:,1) + a(:,2).*T/2 + a(:,3).*T^2/3 + a(:,4).*T^3/4 + a(:,5).*T^4/5 + a(:,6)/T);
h = sum(h_bar.*yi)/M;

//Cálculo da Entropia
s_bar = Ru*(a(:,1)*log(T) + a(:,2).*T + a(:,3).*T^2/2 + a(:,4).*T^3/3 + a(:,5).*T^4/4 + a(:,7));
si_bar = s_bar - Ru*log(pi./po);
s = sum(si_bar.*yi)/M;

//Cálculo do Volume Específico
v = Ru*T/P/M;

//Cálculo da Energia Interna
u = h - P*v;
