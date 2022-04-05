clc
clear

 //Criado em 10 de outubro de 2018
 //Engenheiro responsável: Adrian
 //Engenheiro verificador: Carlos Renan Cândido da Silva

//Programa para calcular a Temperatura Adiabática de chama

R = 8.315/(44+18+28);   //Constante da mistura
To = 298.15;            //Temperatura ambiente [K]
PCI = 191.759*4.184;    //Cálculo para conversão do PCI para kJ

a = [0.446080*10^1  0.309817*10^-2 -0.123925*10^-5  0.227413*10^-9  -0.155259*10^-13 -0.489614*10^5 -0.986359
     0.271676*10^1  0.294513*10^-2 -0.802243*10^-6  0.102266*10^-9  -0.484721*10^-14 -0.299058*10^5  0.663056*10^1
     0.289631*10^1  0.151548*10^-2 -0.572352*10^-6  0.998073*10^-10 -0.652235*10^-14 -0.905861*10^3  0.616151*10^1
     0.362195*10^1  0.736182*10^-3 -0.196522*10^-6  0.362015*10^-10 -0.289456*10^-14 -0.120198*10^4  0.361509*10^1
     0.298406*10^1  0.148913*10^-2 -0.578996*10^-6  0.103645*10^-9  -0.693535*10^-14 -0.142452*10^5  0.634791*10^1
     0.310019*10^1  0.511194*10^-3  0.526442*10^-7 -0.349099*10^-10  0.369453*10^-14 -0.877380*10^3 -0.196294*10^1
         0.25*10^1         0               0               0                0         0.254716*10^5 -0.460117
     0.254205*10^1 -0.275506*10^-4 -0.310280*10^-8  0.455106*10^-11 -0.436805*10^-15  0.292308*10^5  0.492030*10^1
     0.291064*10^1  0.959316*10^-3 -0.194417*10^-6  0.137566*10^-10  0.142245*10^-15  0.393538*10^4  0.544234*10^1
       0.3189*10^1  0.133822*10^-2 -0.528993*10^-6  0.959193*10^-10 -0.648479*10^-14  0.982832*10^4  0.674581*10^1];
       
//fuel(combustivel) = [a b c d], onde a é o número de atomos de Carbono, 
//b é o número de átomos de Hidrogênio, c é o número de átomos de Oxigênio e
//d é o número de átomos de Nitrogênio

fuel = [1 4 0 0];
phi = [1 inv(1.2) inv(1.5)]; //Razão de equivalência ACs/AC, onde AC é a relação ar-combustível
A = [1 0 0 0; 0 1 0 0; 2 1 0 2; 0 0 2 0]; //Matriz para a solução do sistema linear do cálculo estequiométrico
s = zeros(3,4);
Ta = zeros(5,3);

for p = 1:3
    //Coeficiente referente ao número de mols de ar
    a_s = fuel(1)+fuel(2)/4-fuel(3)/2;
    //Matriz para a solução do sistema linear do cálculo estequiométrico
    b = [fuel(1); fuel(2)/2; fuel(3)+(2*a_s/phi(p)); fuel(4)+(3.76*2*a_s/phi(p))]
    s(p,:) = A\b; //Solução do sistema linear /(p,:)
    ni = [s(p,1) s(p,2) s(p,3) s(p,4) 0 0 0 0 0 0];
    N = sum(ni);
    yi = [ni./N]';
    ay_co2 = yi(1)*a(1,:);
    ay_h2o = yi(2)*a(2,:);
    ay_n2 = yi(3)*a(3,:);
    ay_o2 = yi(4)*a(4,:);
    ay = [ay_co2;ay_h2o;ay_n2;ay_o2];
    AY = [sum(ay,1)];
    cte = -(R*AY(1)*To + R*AY(2)*To^2/2 + R*AY(3)*To^3/3 + R*AY(4)*To^4/4 + R*AY(5)*To^5/5) - PCI;
    f1 = [R*AY(5)/5 R*AY(4)/4 R*AY(3)/3 R*AY(2)/2 R*AY(1) cte];
    Ta(:,p) = roots(f1)
end

disp(Ta(5,:))
