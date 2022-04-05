 //Criado em 31 de agosto de 2018
 //Modificado em 16 de setembro de 2018
 //Engenheiro responsável: Adrian
 //Engenheiro verificador: Carlos Renan Cândido da Silva
 
 //Programa para: calcular os valores de n1, n2, n3 e n4, que são os coeficientes 
 //referentes à aplicação da conservação da massa para cada reação química de combustão;
 //calcular os valores de massas de cada componente; plotar gráficos em barras referentes 
 //às proporções de cada componente; aplicar o conceito de limites de inflamabilidade; 
 //avaliar as concentrações de produtos das reações na base seca.

//fuel(combustivel) = [a b c d], onde a é o número de atomos de Carbono, 
//b é o número de átomos de Hidrogênio, c é o número de átomos de Oxigênio e
//d é o número de átomos de Nitrogênio
fuel = [1 4 0 0 ;
        3 8 0 0 ;
        7 17 0 0 ;
        8 18 0 0 ;
        14.4 24.9 0 0 ;
        15 32 0 0 ;
        1 4 1 0 ;
        2 6 1 0 ;
        1 3 2 1 ;
        0 2 0 0 ;
        2 2 0 0 ;
        2 0 0 2 ;
        0 3 0 1 ;
        6 6 0 0 ;
        10 8 0 0 ;
        1 0 0 0 ;
        176 144 8 3 ];

//Razão de equivalência phi=ACs/AC, onde AC é a relação ar-combustível e o subscrito "s"
//é referente à relação ar-combustível estequiométrica
phi = 1; 
//Matriz para a solução do sistema linear do cálculo estequiométrico
A = [1 0 0 0; 0 1 0 0; 2 1 2 0; 0 0 0 2]; 
//Matriz Solução
x = zeros(4,17);

//Matriz nome, com o nome de cada combustível na posição correspondente
titulo = ['Metano';'Propano';'Gasolina';'Octano';'Diesel';'Pentadecano';
         'Metanol';'Etanol';'Nitrometano';'Hidrogênio';'Acetileno';'Cianogênio';
         'Amônia';'Benzeno';'Naftaleno';'Grafite';'Carvão']

//Matriz de Limites Inferiores de Inflamibilidade
LI = [.05;.021;.014;.01;.006;.0045;.067;.033;.073;.04;.025;.06;.15;.0135;.009;%nan;%nan];
//Matriz de Limites Superiores de Inflamibilidade
LS = [.15;.101;.076;.07;.075;.065;.36;.19;.222;.75;.81;.426;.28;.0665;.059;%nan;%nan];

for p = 1:17
    a_s(p) = fuel(p,1)+fuel(p,2)/4-fuel(p,3)/2; //Coeficiente referente ao número de mols de ar
    b = [fuel(p,1); fuel(p,2)/2; fuel(p,3)+(2*a_s(p)/phi); fuel(p,4)+...
    (3.76*2*a_s(p)/phi)] //Matriz para a solução do sistema linear do cálculo estequiométrico
    x(:,p) = A\b; //Solução do sistema linear
end

//Dados referentes a reação estequiométrica
N = x';//Onde N(1)=n1, N(2)=n2, N(3)=n3 e N(4)=n4
M = [12; 1; 16; 14]; //Massas molares do carbono, hidrogênio, oxigênio e nitrogênio [g/mol]
m_comb = fuel*M; //Massa de combustivel [g]
m_ar = a_s/phi*(M(3)*2+3.76*M(4)*2); //Massa de ar [g]
m_co2 = N(:,1)*(M(1)+M(3)*2); //Massa de CO2 [g]
m_h2o = N(:,2)*(M(2)*2+M(3)); //Massa de H2O [g]
m_o2 = N(:,3)*M(3)*2; //Massa de O2 [g]
m_n2 = N(:,4)*M(4)*2; //Massa de N2 [g]
m_r = m_comb + m_ar; //Massa dos reagentes [g]
m_p = m_co2 + m_h2o + m_o2 + m_n2; //Massa dos produtos [g]

//Dados referentes às reações nos Limites Inferiores e Superiores de Inflamibilidade

y_LI = (1 ./LI-1)/4.76; //número de mols para o LI
y_LS = (1 ./LS-1)/4.76; //número de mols para o LS
m_ar_LI = y_LI*(M(3)*2+3.76*M(4)*2); //Massa de ar para o LI [g]
m_ar_LS = y_LS*(M(3)*2+3.76*M(4)*2); //Massa de ar para o LS [g]
m_r_LI = m_comb + m_ar_LI; //Massa dos reagentes para o LI  [g]
m_r_LS = m_comb + m_ar_LS; //Massa dos reagentes para o LS [g]
l_LI = m_comb./m_r_LI; //Relação de massa de combustivel e massa dos reagentes para o LI
l_LS = m_comb./m_r_LS; //Relação de massa de combustivel e massa dos reagentes para o LS
m_comb_LI = m_r(:).*l_LI(:); //Massa de combustível para o LI ajustada a massa de reagentes do cálculo estequiométrico [g]
m_comb_LS = m_r(:).*l_LS(:); //Massa de combustível para o LS ajustada a massa de reagentes do cálculo estequiométrico [g]
m_ar_LI_2 = m_r(:) - m_comb_LI(:); //Massa de ar para o LI ajustada a massa de reagentes do cálculo estequiométrico [g]
m_ar_LS_2 = m_r(:) - m_comb_LS(:);  //Massa de ar para o LS ajustada a massa de reagentes do cálculo estequiométrico [g]

//Dados referentes às concentrações dos produtos das reações na base seca

prop_co2 = N(:,1)./(N(:,1)+N(:,3)+N(:,4))*100;   //Proporção de co2 (%)
prop_o2 = N(:,3)./(N(:,1)+N(:,3)+N(:,4))*100;    //Proporção de o2 (%)
prop_n2 = N(:,4)./(N(:,1)+N(:,3)+N(:,4))*100;    //Proporção de n2 (%)

P = [prop_co2 prop_o2 prop_n2];     //Matriz para as proporções dos produtos de todas as reações

//Plotagem dos gráficos
for p = 1:17
    scf(p)
    bar([1 2 3 4],[m_comb_LI(p,:),m_ar_LI_2(p,:),0,0,0,0;m_comb_LS(p,:),m_ar_LS_2(p,:),0,0,0,0;m_comb(p,:),m_ar(p,:),0,0,0,0;0,0,m_co2(p,:),m_h2o(p,:),m_n2(p,:),m_o2(p,:)],'stacked') 
    a=gca()
    a.x_ticks.labels=["LI";"LS";"Reagentes";"Produtos"] //Legenda para as barras empilhadas
    
    //Legenda para as cores das barras
    hl=legend(['Combustível';'Ar';'Dióxido de carbono';'Água';'Nitrogênio';'Oxigênio'],-1);
    title(titulo(p),"fontsize",5) //Título do gráfico
    ylabel("Massa [g]","fontsize",4) //Título para o eixo y
end
