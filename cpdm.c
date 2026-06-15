
#include <Stdlib.h>
#include <String.h>
#include <string.h>
#include <math.h>
#include <Psim.h>
#include <stdio.h>
#define BUFF 100
#define LC_MAX 10
//#define MATRIZ_NULA() matriz((double[BUFF]){0}, 0, 0)



//PROTÓTIPOS DE OPERAÇÃO DE MATRIZES
//typedef struct Matriz Matriz;

typedef struct Matriz {
	int colunas;
	int linhas;
	double* matriz;
	double (*ret)(Matriz*, int, int);
	double (*transposta)(Matriz*);
	void (*print)(Matriz*);
} Matriz;



void cop_mat(Matriz* m1, Matriz* m2);
Matriz matriz(double* mat, int linhas, int colunas);
Matriz matriz_nula();
void mat_pot(Matriz* mat, Matriz* res, int pot);
double ret(Matriz* self, int l, int c);
void transposta(Matriz* self, Matriz* res);
void print(Matriz *self);
void pm(Matriz* m1, Matriz* m2, Matriz* res, char* op);
void sm (Matriz* m1, Matriz* m2, Matriz* res);
void sub(Matriz* m1, Matriz* m2, Matriz* res);
void ident(Matriz* m, int t);
void ret_lin(Matriz* m, Matriz* dest, int l);
void mult(Matriz* m, Matriz* destino, double v);
void subst(Matriz* m_linha, Matriz* dest , int indice);
void extend(Matriz* m, Matriz* n, Matriz* dest);
void mat_mult(Matriz* a, Matriz* b, Matriz* res);
double int_prod(Matriz* l1, Matriz* l2);
void inv(Matriz *origin, Matriz* dest);
Matriz zeros(int l, int c); 

Matriz MATRIZ_NULA( ) {
	double dados[BUFF];
	int i;
	for(i = 0; i < BUFF; i++) dados[i] = 0;
	Matriz m = matriz(dados, 6, 6);
	return m;
}

//PROTÓTIPOS DE DMPC
void qp_hild(Matriz* H, Matriz* f, Matriz* A_cons, Matriz* b, Matriz* x_opt);
double malha_fechada(int Nc, int Np, Matriz* F, Matriz* Phi, Matriz* R_barra, 
		   Matriz* Rs_barra, double r_ki, double rw, Matriz* x_ki);
void func_deltau(Matriz* Phi, Matriz* R_barra, Matriz* Rs_barra, Matriz* F, double r_ki, Matriz* x_ki, Matriz* ret); //pronto
void func_Phi(Matriz* A, Matriz* B, Matriz* C, int Np, int Nc, Matriz* dest); //pronto
void func_F(Matriz* A, Matriz* C, Matriz* dest, int Np); //pronto
void ret_A(Matriz* Am, Matriz* Cm, Matriz* A); //pronto
void ret_B(Matriz* Bm, Matriz* Cm, Matriz* B); //pronto
void ret_C(Matriz* Cm, Matriz* C); //pronto 
void discret(double, Matriz*, Matriz*, Matriz*, Matriz*); //pronto
void zoh(Matriz* A, Matriz* B, double tempo, Matriz* Ad, Matriz* Bd);
void func_Kmpc_Ky(Matriz* Phi, Matriz* R_barra, Matriz* Rs_barra, Matriz* F, Matriz* K_mpc, double* K_y);

// PLACE GLOBAL VARIABLES OR USER FUNCTIONS HERE...


//===============================================================================//
//                                                   VARIÁVEIS GLOBAIS                                                  //
//===============================================================================//

static double   ilf = 0.0;
static double vcf = 0.0;
static double ilo = 0.0;
static double vco = 0.0;
static double vd  = 0.0;

static double u_anterior = NULL;

double Cf, Co, Lf, Lo, D, Ro;
Matriz A_planta, B_planta, C_planta, D_planta, Phi, F;
double Delta_T, Ts;
Matriz Ad, Bd, Cd;
Matriz A, B, C;
int Nc, Np, r_ki;
double rw;


static Matriz x_anterior;
static Matriz x_ki;

double rs_barra[BUFF]; Matriz Rs_barra;
double r_barra1[BUFF];  Matriz R_barra1;
double r_barra[BUFF];    Matriz R_barra; 
double deltau[BUFF]; Matriz Delta_U;

static double k_mpc_array[BUFF];
static Matriz K_mpc;
static double K_y = 0.0;




//===============================================================================//
//                                                            MAIN                                                                //
//===============================================================================//

/////////////////////////////////////////////////////////////////////

void SimulationStep(
		double t, double delt, double *in, double *out,
		 int *pnError, char * szErrorMsg,
		 void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr)
{

//essas duas têm de ser jogadas para o começo do código 


	x_ki.matriz[0] = in[0];
	x_ki.matriz[1] = in[1];
	x_ki.matriz[2] = in[2];
	x_ki.matriz[3] = in[3];
	x_ki.matriz[4] = in[4];


// Proteção do estado anterior
	if(u_anterior == NULL) { 
		u_anterior = 0; 
		cop_mat(&x_ki, &x_anterior); 
	}
	
	// 1. Multiplica o ganho constante pelo estado: M_kx = K_mpc * x_ki
	double buff_kx[BUFF]; Matriz M_kx = matriz(buff_kx, 1, 1);
	mat_mult(&K_mpc, &x_ki, &M_kx);
	
	// 2. Lei de controle puramente algébrica: Delta_u = Ky*r - Kmpc*x
	double delta_u_escalar = (K_y * r_ki) - M_kx.ret(&M_kx, 0, 0);
	
	// 3. Integra o controle (u(k) = u(k-1) + Delta_u)
	double du_calc = u_anterior + delta_u_escalar;
	
	// 4. Aplica as restrições físicas (saturação do Duty Cycle entre 0 e 0.5)
	if(du_calc > 0.5) { du_calc = 0.5; } 
	else if(du_calc < 0) { du_calc = 0; } 
	
	// 5. Atualiza variáveis para o próximo step e envia o sinal
	u_anterior = du_calc;
	out[0] = du_calc;


}



void SimulationBegin(
		const char *szId, int nInputCount, int nOutputCount,
		 int nParameterCount, const char ** pszParameters,
		 int *pnError, char * szErrorMsg,
		 void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr)
{

	double x_ant[BUFF]; x_anterior = matriz(x_ant, LC_MAX, LC_MAX);
	static double xki[5]; x_ki = matriz(xki, 5, 1);
	Cf = 4.62963e-06;
	Co = 2.00543e-07;
	Lf = 0.00379977;
	Lo = 0.000877193;
	D  = 0.5;
	Ro = 10;
	Nc = 5;
	Np = 10;
	r_ki = 25;
	rw = 0.5;
	Delta_T=1/300e3;
	Ts = Delta_T;

	double a_planta[] = {     0,  -1/Lf,         0,                 0, 
	                              1/Cf,        0,   -1/Cf,                 0,
								       0,   D/Lo,         0,           -1/Lo,
									   0,        0,     1/Co,   -1/(Ro*Co)};
	double b_planta[] = {0, -1/Cf, 1/Lo, 0};
	double c_planta[] = {0, 0, 0, 1};
	     
	A_planta = matriz(a_planta, 4, 4);
	B_planta = matriz(b_planta, 4, 1);
	C_planta = matriz(c_planta, 1, 4);
	
	//passa-se do contínuo para o discreto
	//void discret(double T, Matriz* Am, Matriz* Bm, Matriz* Ad, Matriz* Bd)
	double ad[BUFF];   Ad = matriz(ad, LC_MAX, LC_MAX);
	double bd[BUFF];   Bd = matriz(bd, LC_MAX, LC_MAX);
	double cd[] = {0, 0, 0, 1};   Cd = matriz(cd, 1, 4); 
	zoh(&A_planta, &B_planta, Delta_T, &Ad, &Bd);

	
	//estende-se o modelo
	double a[BUFF]; A = matriz(a, LC_MAX, LC_MAX);
	double b[BUFF]; B = matriz(b, LC_MAX, LC_MAX);
	double c[BUFF]; C = matriz(c, LC_MAX, LC_MAX);

	ret_A(&Ad, &Cd, &A);
	ret_B(&Bd, &Cd, &B);
	ret_C(&Cd, &C);

	Rs_barra  = matriz(rs_barra, Np, 1);
	int i;
	for (i = 0; i < Np; i++)
		Rs_barra.matriz[i] = 1;
	R_barra1 = matriz(r_barra1, LC_MAX, LC_MAX); ident(&R_barra1, Nc);
	R_barra = matriz(r_barra, LC_MAX, LC_MAX);
	mult(&R_barra1, &R_barra, rw); //VERIFICAR

	Delta_U = matriz(deltau, LC_MAX, LC_MAX);

	
	//cria=se as matrizes de base
	double phi[BUFF]; Phi = matriz(phi, LC_MAX, LC_MAX);
	double f[BUFF];   F = matriz(f, LC_MAX, LC_MAX);
	func_Phi(&A, &B, &C, Np, Nc, &Phi);
	func_F(&A, &C, &F, Np);

	
	
	// Instanciação da Matriz K_mpc global como matriz linha (1 x n_estados)
	K_mpc = matriz(k_mpc_array, 1, A.linhas);
	
	// Calcula os ganhos de forma offline antes do simulador rodar
	func_Kmpc_Ky(&Phi, &R_barra, &Rs_barra, &F, &K_mpc, &K_y);

}


void SimulationEnd(const char *szId, void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr)
{
}

//===============================================================================//
//                                                    FIM DE MAIN                                                            //
//===============================================================================//


//###############################################################################//
//                                                             DMPC                                                              //
//###############################################################################//

void zoh(Matriz* A, Matriz* B, double tempo, Matriz* Ad, Matriz* Bd) {
	//Cálculo de Ad
	double m0[BUFF]; Matriz M0 = matriz(m0, LC_MAX, LC_MAX); ident(&M0, A->linhas); //I
	double m1[BUFF]; Matriz M1 = matriz(m1, LC_MAX, LC_MAX); mult(A, &M1, tempo);   //A*T
	double m2[BUFF]; Matriz M2 = matriz(m2, LC_MAX, LC_MAX); mat_mult(&M1, &M1, &M2);
	double m3[BUFF]; Matriz M3 = matriz(m3, LC_MAX, LC_MAX); mult(&M2, &M3, 0.5);		//(A*T)^2/2!
	double m4[BUFF]; Matriz M4 = matriz(m4, LC_MAX, LC_MAX); sm(&M0, &M1, &M4);
	double m5[BUFF]; Matriz M5 = matriz(m5, LC_MAX, LC_MAX); sm(&M4, &M3, &M5);
	cop_mat(&M5, Ad);

	//Cálculo de Bd
	double m_bd1[BUFF]; Matriz M_bd1 = matriz(m_bd1, LC_MAX, LC_MAX); 
	mult(&M0, &M_bd1, tempo); 

	// 2. M_bd2 = A * (T^2)/2. Note que M1 já é (A*T). 
	//    Logo, (A*T) * (T/2) resulta em A*(T^2)/2
	double m_bd2[BUFF]; Matriz M_bd2 = matriz(m_bd2, LC_MAX, LC_MAX); 
	mult(&M1, &M_bd2, tempo / 2.0); 

	// 3. M_bd3 = (I*T + A*(T^2)/2)
	double m_bd3[BUFF]; Matriz M_bd3 = matriz(m_bd3, LC_MAX, LC_MAX); 
	sm(&M_bd1, &M_bd2, &M_bd3); 

	// 4. Bd = M_bd3 * B
	double m_bd4[BUFF]; Matriz M_bd4 = matriz(m_bd4, LC_MAX, LC_MAX); 
	mat_mult(&M_bd3, B, &M_bd4); 
	
	cop_mat(&M_bd4, Bd);
}

void func_Kmpc_Ky(Matriz* Phi, Matriz* R_barra, Matriz* Rs_barra, Matriz* F, Matriz* K_mpc, double* K_y) {
	// Phi_T = transposta(Phi)
	double phi_t[BUFF]; Matriz Phi_T = matriz(phi_t, LC_MAX, LC_MAX);
	transposta(Phi, &Phi_T);
	
	// Phi_T_Phi = Phi_T * Phi
	double phi_t_phi[BUFF]; Matriz Phi_T_Phi = matriz(phi_t_phi, LC_MAX, LC_MAX);
	mat_mult(&Phi_T, Phi, &Phi_T_Phi);
	
	// M0 = Phi_T_Phi + R_barra (Matriz Hessiana)
	double m0[BUFF]; Matriz M0 = matriz(m0, LC_MAX, LC_MAX);
	sm(&Phi_T_Phi, R_barra, &M0);
	
	// M1 = inv(M0)
	double m1[BUFF]; Matriz M1 = matriz(m1, LC_MAX, LC_MAX);
	inv(&M0, &M1);
	
	// M2 = M1 * Phi_T
	double m2[BUFF]; Matriz M2 = matriz(m2, LC_MAX, LC_MAX);
	mat_mult(&M1, &Phi_T, &M2);
	
	// -----------------------------------------------------------------
	// Cálculo do K_mpc (M2 * F)
	// -----------------------------------------------------------------
	double kmpc_full[BUFF]; Matriz K_mpc_full = matriz(kmpc_full, LC_MAX, LC_MAX);
	mat_mult(&M2, F, &K_mpc_full);
	
	// K_mpc recebe apenas a PRIMEIRA LINHA de K_mpc_full (índice 0)
	ret_lin(&K_mpc_full, K_mpc, 0);
	
	// -----------------------------------------------------------------
	// Cálculo do K_y (M2 * Rs_barra)
	// -----------------------------------------------------------------
	double ky_full[BUFF]; Matriz K_y_full = matriz(ky_full, LC_MAX, LC_MAX);
	mat_mult(&M2, Rs_barra, &K_y_full);
	
	// K_y recebe apenas o PRIMEIRO ELEMENTO de K_y_full (linha 0, col 0)
	*K_y = K_y_full.ret(&K_y_full, 0, 0);
}

void func_deltau(Matriz* Phi, Matriz* R_barra, Matriz* Rs_barra, Matriz* F, double r_ki, Matriz* x_ki, Matriz* ret){
//Obtenção do incremento do sinal de controle Delta_U
	double phi_t[BUFF]; Matriz Phi_T     = matriz(phi_t, LC_MAX, LC_MAX); transposta(Phi, &Phi_T);
	double phi_t_phi[BUFF]; Matriz Phi_T_Phi = matriz(phi_t_phi, LC_MAX, LC_MAX); mat_mult(&Phi_T, Phi, &Phi_T_Phi);
	double m0[BUFF]; Matriz M0        = matriz(m0, LC_MAX, LC_MAX); sm(&Phi_T_Phi, R_barra, &M0);
	double m1[BUFF]; Matriz M1        = matriz(m1, LC_MAX, LC_MAX); inv(&M0, &M1);
	double m2[BUFF]; Matriz M2        = matriz(m2, LC_MAX, LC_MAX); mat_mult(&M1, &Phi_T, &M2);
	double rs_barra[BUFF]; Matriz Rs_barra_r_ki = matriz(rs_barra, LC_MAX, LC_MAX); mult(Rs_barra, &Rs_barra_r_ki, r_ki);
	double fx_ki[BUFF]; Matriz Fx_ki     = matriz(fx_ki, LC_MAX, LC_MAX); mat_mult(F, x_ki, &Fx_ki);
	double m3[BUFF]; Matriz M3        = matriz(m3, LC_MAX, LC_MAX); sub(&Rs_barra_r_ki, &Fx_ki, &M3);
	mat_mult(&M2, &M3, ret);
	//inv(Phi_T*Phi + R_barra)*Phi^T*(Rs_barra*r_ki - F*x_ki)
}


void ret_A(Matriz* Am, Matriz* Cm, Matriz* A){
//Entrega da versão estendida de A para CPDM
	Matriz M0 = zeros(Am->linhas, 1);
	Matriz M1   = zeros(Am->linhas, 1);
	double m2[BUFF]; Matriz M2   = matriz(m2, LC_MAX, LC_MAX); extend(Am, &M1, &M2);              
	double cmam[BUFF]; Matriz CmAm = matriz(cmam, LC_MAX, LC_MAX); mat_mult(Cm, Am, &CmAm);        
	double intg[BUFF];     Matriz Intg = matriz(intg, LC_MAX, LC_MAX); ident(&Intg, 1);
	double m3[BUFF];      Matriz M3   = matriz(m3, LC_MAX, LC_MAX); extend(&CmAm, &Intg, &M3);
	double m4[BUFF];      Matriz M4   = matriz(m4, LC_MAX, LC_MAX); transposta(&M2, &M4);  
	double m5[BUFF];      Matriz M5   = matriz(m5, LC_MAX, LC_MAX); transposta(&M3, &M5); 
	double m6[BUFF];      Matriz M6   = matriz(m6, LC_MAX, LC_MAX); extend(&M4, &M5, &M6); 
	transposta(&M6, A);
}


void ret_B(Matriz* Bm, Matriz* Cm, Matriz* B){
//Entrega da versão estendida de B para CPDM
	double cmbm[BUFF]; Matriz CmBm = matriz(cmbm, LC_MAX, LC_MAX); mat_mult(Cm, Bm, &CmBm);
	double m0[BUFF];     Matriz M0   = matriz(m0, LC_MAX, LC_MAX); transposta(Bm, &M0);
	double m1[BUFF];     Matriz M1   = matriz(m1, LC_MAX, LC_MAX); extend(&M0, &CmBm, &M1);
	transposta(&M1, B);
}


void ret_C(Matriz* Cm, Matriz* C){
//Entrega da versão estendida de C para CPDM
	double i[BUFF]; Matriz I = matriz(i, LC_MAX, LC_MAX);; ident(&I, 1);
	Matriz M0 = zeros(1, Cm->colunas);
	extend(&M0, &I, C);
}


void func_Phi(Matriz* A, Matriz* B, Matriz* C, int Np, int Nc, Matriz* dest){
//Matriz deslizada de Toeplitz formada a partir de um vetor de fatores de
//potência a serem deslizados


	double phi[BUFF]; Matriz Phi = matriz(phi, 1, Np);
	double phi_prim[BUFF]; Matriz Phi_primario = matriz(phi_prim, 1, Np);
	double temp1[BUFF]; Matriz Temp1 = matriz(temp1, LC_MAX, LC_MAX);
	double temp2[BUFF]; Matriz Temp2 = matriz(temp2, LC_MAX, LC_MAX);
	double temp3[BUFF]; Matriz Temp3 = matriz(temp3, LC_MAX, LC_MAX);
	int indice = Np - 1;
	
	
	//Crição do vetor (1xNp) a ser deslizado para a criação de Toeplitz
	int i;
	for(i = 1; i <= Np; i++){
		if((Np-i) != 0) {
			int pot = Np - i;
			mat_pot(A, &Temp1, pot);     //Temp1.print(&Temp1);
			mat_mult(C, &Temp1, &Temp2); //Temp2.print(&Temp2);
			mat_mult(&Temp2, B, &Temp3); //Temp3.print(&Temp3);
			  
		} else {
			mat_mult(C, B, &Temp3); 
			}
		Phi_primario.matriz[indice--] = Temp3.ret(&Temp3, 0, 0);
	}
	
	dest->linhas = Np;
	dest->colunas = Nc;
	int indice2 = 0;
	int j;
	for(j = 0; j < Np; j++){
		int indice3 = j;
		int k;
		for(k = 0; k < Nc; k++){
			if(indice3>=0){
				dest->matriz[indice2++] = 
					Phi_primario.ret(&Phi_primario, 0, indice3--); 
			} else {
				dest->matriz[indice2++] = 0;
			}
		}
	}

}


void func_F(Matriz* A, Matriz* C, Matriz* dest, int Np){
//Entrega a matriz compacta F (representativa das operações inerentes
//ao sistema
	double m0[BUFF]; Matriz M0 = matriz(m0, LC_MAX, LC_MAX);
	double m1[BUFF]; Matriz M1 = matriz(m1, LC_MAX, LC_MAX);
	double m2[BUFF]; Matriz M2 = matriz(m2, LC_MAX, LC_MAX);
	double f[BUFF];    Matriz F   = matriz(f,    LC_MAX, LC_MAX); 
	dest->linhas = Np;
	dest->colunas = 2;

	int indice = 0;
	int i;
	for(i = 1; i <= Np; i++){
		mat_pot(A, &M0, i);    
		mat_mult(C, &M0, &M1); 
		transposta(&M1, &M2);  
		dest->matriz[indice++] = M2.ret(&M2, 0, 0);
		dest->matriz[indice++] = M2.ret(&M2, 0, 1); 
	}
}



/*
void qp_hild(Matriz* H, Matriz* f, Matriz* A_cons, Matriz* b, Matriz* x_opt){
	const int n1 = A_cons->linhas;
	const int m1 = A_cons->colunas;
	Matriz eta = MATRIZ_NULA();
	Matriz inv_H = MATRIZ_NULA(); inv(H, &inv_H);
	Matriz neg_inv_H = MATRIZ_NULA(); mult(&inv_H, &neg_inv_H, -1);
	mat_mult(&neg_inv_H, f, &eta);
	cop_mat(&eta, x_opt);
	
	double max(double a, double b){
	//função de retorno de máximo de 2 números
		if(a > b) return a;
		else { return b; }
	}

	int kk = 0;
	for(int i = 0; i < n1; i++){
		Matriz A_cons_linha = MATRIZ_NULA();
		Matriz A_cons_linha_eta = MATRIZ_NULA();
		ret_lin(A_cons, &A_cons_linha, i);
		mat_mult(&A_cons_linha, &eta, &A_cons_linha_eta);
		if(A_cons_linha_eta.ret(&A_cons_linha_eta, 0, 0) 
		   > b->ret(b, i, 1))
		   	kk++;
		else { kk = kk + 0; }
	}
	if(kk==0) { 
		printf("Solução ótima encontrada.\nRetornando...\n");
		return;
	}//solução ótima encontrada inicialmente
	
	Matriz P = MATRIZ_NULA(),
	       A_cons_tran = MATRIZ_NULA(),
	       H_A_cons = MATRIZ_NULA(),
	       H_f = MATRIZ_NULA(),
	       A_cons_H_f = MATRIZ_NULA(),
	       d = MATRIZ_NULA();
	transposta(A_cons, &A_cons_tran); //A_cons'
	mat_mult(&inv_H, &A_cons_tran, &H_A_cons); // H^(-1)*A_cons'
	mat_mult(A_cons, &H_A_cons, &P);  //A_cons*(H^(-1)*A_cons')

	mat_mult(&inv_H, f, &H_f); //H^(-1)*f
	mat_mult(A_cons, &H_f, &A_cons_H_f); //A_cons*H^(-1)*f
	sm(&A_cons_H_f, b, &d);
	
	const int n = d.linhas;
	const int m = d.colunas;
	Matriz x_ini = zeros(n, m);
	Matriz lambda   = MATRIZ_NULA(); cop_mat(&x_ini, &lambda);
	double al = 10; //número qualquer acima de 0 para evitar problemas
	                //de convergência
	for(int i = 0; i < 38; i++){ //1 SOMENTE PARA TESTES
	//38 iterações como procura máxima para se evitar recursão infinita
		Matriz lambda_p = MATRIZ_NULA();
		cop_mat(&lambda, &lambda_p);
		for(int j = 0; j < n; j++){
			Matriz P_ret = MATRIZ_NULA(); ret_lin(&P, &P_ret, j);
			double w = int_prod(&P_ret, &lambda) 
			           - P.ret(&P, j, j)*lambda.ret(&lambda, j, 0);
			//primeiro resultado é igual a 0.000 por lambda = |0
			w = w + d.ret(&d, j, 0);
			double la = -1*w/P.ret(&P, j, j);
			lambda.matriz[j] = max(0, la);
		}

		Matriz M1 = MATRIZ_NULA(); sub(&lambda, &lambda_p, &M1);
		Matriz M0 = MATRIZ_NULA(); transposta(&M1, &M0);
		double al = int_prod(&M0, &M1);
		
		if(al < 10e-8) {
			printf("al menor que limite.\nParando...\n");
			break;
		}
		
	
	}
	Matriz M2 = MATRIZ_NULA(), 
	       M3 = MATRIZ_NULA(),
	       M4 = MATRIZ_NULA(),
	       H_2 = MATRIZ_NULA();
	mat_mult(&neg_inv_H, f, &M2); //-H^(-1)*f 
	mat_mult(&neg_inv_H, &A_cons_tran, &M3); //-H^(-1)*A_cons'
	mat_mult(&M3, &lambda, &M4); //-H^(-1)*A_cons'*lambda
	sm(&M2, &M4, x_opt);
}

double malha_fechada (int Nc, int Np, Matriz* F, Matriz* Phi, Matriz* R_barra, 
		      Matriz* Rs_barra, double r_ki, double rw, Matriz* x_ki){
//Malha fechada para execução real, não manipulando entrada por meio de ganhos,
//enviando para fora da função somente o sinal de controle
//[!!!] DEVE SER TESTADO POSTERIORMENTE	
	static double u_anterior = NAN;
	static Matriz* x_anterior;
	if(isnan(u_anterior)){
		x_anterior = x_ki;
		u_anterior = 0;
	}
	Matriz Delta_U = MATRIZ_NULA(); 
	func_deltau(Phi, R_barra, Rs_barra, F, r_ki, x_ki, &Delta_U);
	double du_calc = u_anterior + Delta_U.ret(&Delta_U, 0, 0);
	
	if(du_calc >0.5)
		du_calc = 0.5;
	else if(du_calc < 0)
		du_calc = 0;
	
	u_anterior = du_calc;
	//não se está fazendo update de x_ki, aparentemente
	return du_calc;
}

*/



void discret(double T, Matriz* Am, Matriz* Bm, Matriz* Ad, Matriz* Bd){
//Modifica os valores internos das matrizes Bd, e Ad a partir da discretização
//de Tustin feita sobre as Matrizes Am, Bm e sobre o valor de tempo amostrado T

	double mi[BUFF];         Matriz I    = matriz(mi, LC_MAX, LC_MAX);  ident(&I, Am->colunas);
	double a_modif[BUFF]; Matriz A_modif = matriz(a_modif, LC_MAX, LC_MAX); mult(Am, &A_modif, T/2);
	double m0[BUFF];        Matriz M0 = matriz(m0, LC_MAX, LC_MAX); sub(&I, &A_modif, &M0);
	double m1[BUFF];        Matriz M1 = matriz(m1, LC_MAX, LC_MAX); inv(&M0, &M1);          
	double m2[BUFF];        Matriz M2 = matriz(m2, LC_MAX, LC_MAX); sm(&I, &A_modif, &M2);  
	mat_mult(&M1, &M2, Ad);

	double m3[BUFF];        Matriz M3 = matriz(m3, LC_MAX, LC_MAX);; mult(Bm, &M3, T);
	mat_mult(&M1, &M3, Bd);
}



//###############################################################################//
//                                            OPEAÇÕES MATRICIAIS                                                    //
//###############################################################################//
void cop_mat(Matriz* m1, Matriz* m2){
	m2->linhas = m1->linhas;
	m2->colunas = m1->colunas;
	int contador = 0;
	int linhas;
	int colunas;
	for(linhas = 0; linhas < m2->linhas; linhas++){
		for(colunas = 0; colunas < m2->colunas; colunas++){
			m2->matriz[contador++] = m1->ret(m1, linhas, colunas);
		}
	}
}

Matriz matriz(double* mat, int linhas, int  colunas) {
//Construtor do pseudo-objeto Matriz
	Matriz m;
	m.matriz = mat;
	m.colunas = colunas;
	m.linhas = linhas;
	m.ret = ret;
	m.print = print;
	return m;
}


void mat_pot(Matriz* mat, Matriz* res, int pot) {
//Entrega uma potência qualquer de uma matriz fornecida pelo usuário
	if(pot == 0){
		ident(res, mat->linhas);
	}

	 else if  (pot == 1){
		cop_mat(mat, res);
	}
	
	else {
		double buff[BUFF];
		Matriz r = matriz(buff, 6, 6);
		double temp_buff[BUFF];
		Matriz temp = matriz(temp_buff, mat->linhas, mat->colunas);
		cop_mat(mat, &r);
		int i;
		for(i = 1; i < pot; i++){
			mat_mult(mat, &r, &temp);
			cop_mat(&temp, &r);
		}
		cop_mat(&r, res);
	}

}

Matriz zeros(int l, int c){
//Retorna uma matriz de zeros do tamanho especificado pelo
//usuário.
	double* mat = malloc(l*c*sizeof(double));
	memset(mat, 0, l*c*sizeof(double));
	return matriz(mat, l, c);
} 

void inv(Matriz *origin, Matriz* dest){
//Utiliza Gauss-Jordan para a obtenção da matriz inversa de uma 
//matriz quadrada fornecida pelo usuário
	int linhas  = origin->linhas;
	int colunas = origin->colunas; 
	double buff_id[BUFF];    Matriz id  = matriz(buff_id, 6, 6); ; ident(&id, linhas);
	double buff_est[BUFF];  Matriz est = matriz(buff_est, 6, 6); extend(origin, &id, &est); 
	double buff_l1[BUFF]; Matriz l1 = matriz(buff_l1, 6, 6);
	double buff_l2[BUFF]; Matriz l2 = matriz(buff_l2, 6, 6);
	double buff_novo_l1[BUFF]; Matriz novo_l1 = matriz(buff_novo_l1, 6, 6);
	double buff_novo_l2[BUFF]; Matriz novo_l2 = matriz(buff_novo_l2, 6, 6);
	double buff_temp[BUFF];    Matriz temp = matriz(buff_temp, 6, 6);
	int i;
	for(i = 0; i < linhas; i++){
	//Opera Gauss-Jordan sobre uma matriz estendida [A|B], transformando
	//A em uma matriz identidade e B na inversa
		int n = i;
		double fator = 1/ret(&est, n, n);
		ret_lin(&est, &l1, n);
		mult(&l1, &temp, fator);
		subst(&temp, &est, i);
		
		int j;
		for(j = 0; j < linhas; j++){
			if(i!=j){
				ret_lin(&est, &l2, j);
				double novo_fator = ret(&est, j, i);
				ret_lin(&est, &l1, i);
				mult(&l1, &novo_l1, novo_fator);
				sub(&l2, &novo_l1, &novo_l2);
				subst(&novo_l2, &est, j);
			}
		}
	}
	int contador = 0;
	double buff_lin_temp[BUFF]; Matriz Lin_temp = matriz(buff_lin_temp, 6, 6);
	int i2;
	for(i2 = 0; i2 < linhas; i2++){
	//Extrai a matriz B como inversa e passa seu resultado para a
	//matriz de destino
		ret_lin(&est, &Lin_temp, i2);
		int j2;
		for(j2 = colunas; j2 < 2*colunas; j2++)
			dest->matriz[contador++] = Lin_temp.ret(&Lin_temp, 0, j2);
	}
	dest->linhas = origin->linhas;
	dest->colunas = origin->colunas;
}

double int_prod(Matriz* l1, Matriz *l2){
//A partir de duas matrizes-linha, retorna o produto
//interno de ambas
	double prod = 0;
	int i;
	for(i=0; i < l1->colunas; i++)
		prod = prod + (l1->ret(l1, 0, i))*(l2->ret(l2, 0,i));
	return prod;
}


void mat_mult(Matriz* a, Matriz* b, Matriz* res){
//Multiplica duas matrizes fornecidas pelo usuário, contanto que 
//ambas possuam as mesmas dimensões linhas/colunas, colunas/linhas
	
	int linhas = a->linhas;
	int colunas = b->colunas;
	double dados[BUFF];
	Matriz b_transp = matriz(dados, 6, 6);
	transposta(b, &b_transp);
	int indice = 0;
	res->linhas = linhas;
	res->colunas = colunas;
	
	
	int i;
	for(i = 0; i < linhas; i++){
		double dados2[BUFF];
		Matriz l1 = matriz(dados2, 6, 6); 
		ret_lin(a, &l1, i);
		int j;
		for(j = 0; j < colunas; j++){
			double dados3[BUFF];
			Matriz l2 = matriz(dados3, 6, 6);
			ret_lin(&b_transp, &l2, j);
			res->matriz[indice++] = int_prod(&l1, &l2);
		}
	}
}


void extend(Matriz* m, Matriz*n, Matriz* dest){
//Entrega uma matriz estendida [A|B] a partir de duas matrizes 
//A e B fornecidas pelo usuário (com a mesma quantidade de linhas)
	if(m->linhas!=n->linhas){
		printf("ERRO!!! Dimensões de matrizes não conformantes.\nRetornando...\n");
		return;
	}

	int colunas = m->colunas;
	int linhas  = m->linhas;
	int contador = 0;
	int i;
	for(i = 0; i < linhas; i++){
		int j1;
		for(j1 = 0; j1 < m->colunas; j1++)
			dest->matriz[contador++] = m->ret(m, i, j1);
		int j2;
		for( j2 = 0; j2 < n->colunas; j2++)
			dest->matriz[contador++] = n->ret(n, i, j2);
	}
	dest->linhas = linhas;
	dest->colunas = m->colunas + n->colunas;
}




void mult(Matriz* m, Matriz* destino, double v){
//Devolve uma linha qualquer completa de uma matriz multiplicada 
//por um fator v fornecido pelo usuário
	destino->linhas  = m->linhas;
	destino->colunas = m->colunas;
	int i;
	for(i = 0; i < m->colunas*m->linhas; i++)
		destino->matriz[i] = v*m->matriz[i];
}


void subst(Matriz* m_linha, Matriz* dest, int indice) {
//Substitui uma linha arbitrária de uma matriz por uma nova linha
	if(m_linha->colunas != dest->colunas) {
		printf("ERRO!!! Colunas de linha e de destino de tamanhos diferentes"
			"\nRetornando em subst(...)...\n");
	} 
	else {
		int colunas = m_linha->colunas;
		//printf("%f\n", retorno);
		int i;
		for(i = 0; i < colunas; i++) 
			{ dest->matriz[(indice)*colunas + i] = m_linha->ret(m_linha, 0, i); }
	}
}


void ret_lin(Matriz* m, Matriz* dest, int l){
//Retorna uma linha qualquer completa de dentro da matriz fornecida
	int colunas = m->colunas;
	dest->linhas = 1;
	dest->colunas = colunas;
	int i;
	for(i = 0; i < colunas; i++)
		dest->matriz[i] = m->ret(m, l, i);
}


void ident(Matriz* m,  int t){
//Retorna uma matriz identidade de dimensão txt
	m->linhas = t;
	m->colunas = t;
	//memset(m->matriz, 0, sizeof(m->matriz));
	int i;
	for(i = 0; i < t; i++) { m->matriz[(t+1)*i] = 1.0;}
}


double ret(Matriz* self, int l, int c){
//Retorna um elemento definido pelo índice (a,b) de uma matriz
//de qualquer dimensão
	return self->matriz[l*self->colunas + c];
}

void print(Matriz *self){
//Apresenta uma matriz formatada para a visualização do usuário
	const int l = self->linhas;
	const int c = self->colunas;
	double* mat = self->matriz;
	printf("[");
	
	int i;
	for(i = 0; i < l; i++){
		int j;
		for(j = 0; j < c; j++){
			if(j == c-1) printf("%3lg", self->ret(self, i, j));
			else printf("%8.3lg \t", self->ret(self, i, j));
		}
		if(i == l-1) printf("]\n");
		else printf("\n");
	}	
	printf("\n");
}


void transposta(Matriz* self, Matriz* res){
//Retorna uma matriz transposta de uma matriz retangular
//de qualquer dimensão
	const int l = self->linhas;
	const int c = self->colunas;
	res->linhas  = c;
	res->colunas = l;
	int i;
	for(i = 0; i < c; i++){
		int j;
		for(j = 0; j < l; j++){
			res->matriz[i*l + j] = self->ret(self, j, i);
		}
	}
}

void pm (Matriz* m1, Matriz* m2, Matriz* res, char* op){
//Função superior para soma e subtração de matrizes
	res->linhas = m2->linhas;
	res->colunas = m2->colunas;
	if(m1->linhas == m2->linhas && m1->colunas==m2->colunas){
		const int l = m1->linhas;
		const int c = m1->colunas;
		int i;
		for(i = 0; i < l; i++){
			int j;
			for(j = 0; j < c; j++){
				if(strcmp(op, "mais")==0){
					res->matriz[i*c + j] =
						m1->ret(m1, i, j) + m2->ret(m2, i, j);}
				else if (strcmp(op, "menos")==0){
					res->matriz[i*c + j] =
						m1->ret(m1, i, j) - m2->ret(m2, i, j);}
				}	
			}
	}
	else {
		//Erro: Matrizes de tamanhos diferentes
		//Retorna uma matriz |0
		return;
	}				       
}

void sm (Matriz* m1, Matriz* m2, Matriz* res){
//Retorna a soma entre duas matrizes
	pm(m1, m2, res, "mais");
}

void sub(Matriz* m1, Matriz* m2, Matriz* res){
//Retorna a subtração entre duas matrizes
	pm(m1, m2, res, "menos");
}


