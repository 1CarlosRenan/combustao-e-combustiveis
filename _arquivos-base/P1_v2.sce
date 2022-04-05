 //Criado em 31 de agosto de 2018
 //Engenheiro responsável: Adrian
 //Engenheiro verificador: Carlos Renan Cândido da Silva
 
 //Programa para calcular os valores de n1, n2, n3 e n4, que são os coeficientes 
 //referentes à aplicação da conservação da massa para cada reação química de combustão,
 //calculando valores de massas e plotando gráficos em barras

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
        
phi = 1; //Razão de equivalência ACs/AC, onde AC é a relação ar-combustível e o subscrito "s"
        //é referente à relação estequiométrica 
A = [1 0 0 0; 0 1 0 0; 2 1 2 0; 0 0 0 2];
x = zeros(4,17);

for p = 1:17
    a_s(p) = fuel(p,1)+fuel(p,2)/4-fuel(p,3)/2;
    b = [fuel(p,1); fuel(p,2)/2; fuel(p,3)+(2*a_s(p)/phi); fuel(p,4)+...
    (3.76*2*a_s(p)/phi)]
    x(:,p) = A\b;
end

N = x';//Onde N(1)=n1, N(2)=n2, N(3)=n3 e N(4)=n4
M = [12; 1; 16; 14]; //Massas molares do carbono, hidrogênio, oxigênio e nitrogênio [g/mol]
m_comb = fuel*M; //Massa de combustivel [g]
m_ar = a_s/phi*(M(3)*2+3.76*M(4)*2); //Massa de ar [g]]
m_co2 = N(:,1)*(M(1)+M(3)*2); //Massa de CO2 [g]
m_h2o = N(:,2)*(M(2)*2+M(3)); //Massa de H2O [g]
m_o2 = N(:,3)*M(3)*2; //Massa de O2 [g]
m_n2 = N(:,4)*M(4)*2; //Massa de N2 [g]
m_r = m_comb + m_ar; //Massa dos reagentes [g]
m_p = m_co2 + m_h2o + m_o2 + m_n2; //Massa dos produtos [g]

titulo = ['Metano';'Propano';'Gasolina';'Octano';'Diesel';'Pentadecano';
         'Metanol';'Etanol';'Nitrometano';'Hidrogênio';'Acetileno';'Cianogênio';
         'Amônia';'Benzeno';'Naftaleno';'Grafite';'Carvão']; //Matriz nome, com o nome de cada combustível na posição correspondente

for p = 1:17
    figure(p)
    bar([1 2],[m_comb(p,:),m_ar(p,:),0,0,0,0;0,0,m_co2(p,:),m_h2o(p,:),m_n2(p,:),m_o2(p,:)],'stacked') //Plotagem dos gráficos
    
    hl=legend(['Combustível';'Ar';'Dióxido de carbono';'Água';'Nitrogênio';'Oxigênio'],-1); //Legenda para as cores das barras
    
     title(titulo(p),"fontsize",5) //Título do gráfico
    ylabel("Massa [g]","fontsize",4) //Título para o eixo y

end
