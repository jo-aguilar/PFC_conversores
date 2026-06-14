 #include <iostream>
 #include <string>
 #include <cmath>
 #include <fstream>
 #define pi 3.1415
 using namespace std;

 int main() {
   cout << "Teste 1" << endl;
   float passo = 0.01;
   float percentual = 0.20;
   float ress_sup = 1.5e3;
   float ress_inf = 15e3;
   float tolerancia = 9e-3;
   float Co, Cf, Lo, Lf = 0;
   ofstream arquivo("Elementos de energia.txt");
  
   cout << "[!!!] Procurando..." << endl;
     for(float j = 0.01; j <= percentual; j = j + passo){
       for(float k = 0.01; k <= percentual; k = k + passo){
         for(float l = 0.01; l <= percentual; l = l + passo){
           for(float m = 0.01; m <= percentual; m = m + passo){
           // [!!!] VALORES PRESCRITOS
           // Valores mantêm a relação de potência
           float D = 0.5;
           float Ro = 10;
           float V_in = 50;
           float VCo = V_in*D; //valor ideal;
           float ILo = VCo/Ro;
           float ILf = ILo*D; //[???]
           float VCf = 50;
           float Ts = 1/50e3;
           float f_com = 1/Ts;
           float f_ress_i = ress_sup;
           float f_ress_o = ress_inf;
          
           //valores ajustáveis
           float Delta_iLo = ILo*j;
           float Delta_VCf = V_in*k;
           float Delta_vCo = VCo*l;
           float Delta_iLf = ILf*m;
          
           //valores iniciais dos elementos
           Cf = (ILf*Ts*(1- D))/(Delta_VCf);
           Lo = VCo*Ts*(1- D)/(Delta_iLo);
           Lf = (1/Cf)*pow((1/(2*pi*f_ress_i)),2);
           Co = (1/Lo)*pow((1/(2*pi*f_ress_o)),2);
          
           //ondulação da corrente de entrada
           float Delta_iLf_ver = Delta_VCf/(8*Lf*f_com);
           float Delta_vCo_ver = Delta_iLo/(8*Lo*f_com);
           //percentual de ondulação
          
           //frequências prescritas
           float freq_ress_ou = (1/(2*pi*sqrt(Lo*Co)));
           float freq_ress_in = (1/(2*pi*sqrt(Lf*Cf)));
          
           // Teste
           //Co = Delta_iLo*Ts/(8*Delta_vCo)
           //Lf = Delta_VCf*Ts/(8*Delta_iLf)
          
           freq_ress_ou = (1/(2*pi*sqrt(Lo*Co)));
           freq_ress_in = (1/(2*pi*sqrt(Lf*Cf)));
          
           //cout << "\n[!!!] PERCENTUAL DE ONDULAÇÃO" << endl;
           float P_iLf = Delta_iLf_ver*100/(ILf);
           float P_vCo = Delta_vCo_ver*100/(VCo);
          
           if (Cf < 5e-6 and
             Lo < tolerancia and
             Lf < tolerancia and
             (Co < 5e-6 and Co > 2e-7)){
             arquivo << "Cf: " << Cf << endl;
             arquivo << "Co: " << Co << endl;
             arquivo << "Lf: " << Lf << endl;
             arquivo << "Lo: " << Lo << endl;
             arquivo << "Delta_iLo[%]: " << j*100 << "%" << endl;
             arquivo << "Delta_VCf[%]: " << k*100 << "%" << endl;
             arquivo << "Delta_vCo[%]: " << l*100 << "%" << endl;
             arquivo << "Delta_iLf[%]: " << m*100 << "%" << endl;
             cout << endl;
             }
           }
         }
       }
     }
   arquivo.close();
 }
