
#include <Stdlib.h>
#include <String.h>
#include <string.h>
#include <math.h>
#include <Psim.h>
#include <stdio.h>
#define BUFF 100
#define LC_MAX 10


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
void zoh(Matriz* A, Matriz* B, double tempo, Matriz* Ad, Matriz* Bd);

//PROTÓTIPOS DE DMPC COM LAGUERRE
void ret_A(Matriz* Am, Matriz* Cm, Matriz* A); //pronto
void ret_B(Matriz* Bm, Matriz* Cm, Matriz* B); //pronto
void ret_C(Matriz* Cm, Matriz* C); //pronto
void zoh(Matriz* A, Matriz* B, double tempo, Matriz* Ad, Matriz* Bd);
void Al_func (double a, int Nc, Matriz* Al);
void L0_func(double a, int Nc, Matriz* L0);
void Sc1_func(Matriz* B, Matriz* L0, Matriz* Res);
void phi_m_t(Matriz* A, Matriz* Sc_anterior,  Matriz* Sc1, Matriz* Al, Matriz* Res, int m);
void Omega_func(Matriz* A, Matriz* B, Matriz* C, Matriz* Om, int Np,  int Nc, double rw);
void Psi_func(Matriz* A, Matriz* B, Matriz* C, Matriz* Psi, int Np,  int Nc);
void Kmpc_func (Matriz* A, Matriz* B, Matriz* C, Matriz* Kmpc, int Np,  int Nc, int n_estados, double rw, double a);


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
Matriz Kmpc;
Matriz A, B, C;
int Nc, Np, r_ki;
double rw;
int n_estados;
double a;

static Matriz x_anterior;
static Matriz x_ki;
static Matriz Du;
Matriz Neg_kmpc;


//===============================================================================//
//                                                            MAIN                                                                //
//===============================================================================//

void SimulationStep(double t, double delt, double *in, double *out,  int *pnError, char * szErrorMsg,
		 void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr) {


	x_ki.matriz[0] = in[0];
	x_ki.matriz[1] = in[1];
	x_ki.matriz[2] = in[2];
	x_ki.matriz[3] = in[3];
	x_ki.matriz[4] = in[4];

	double kx[1]; Matriz Kmpc_xki = matriz(kx, 1, 1);


	mat_mult(&Neg_kmpc, &x_ki, &Kmpc_xki);
	double du = Kmpc_xki.matriz[0];

	if(du > 0.5) {du = 0.5;}
	else if(du < 0) {du = 0;}
	else {du = du;}

	out[0] = du;

 }

void SimulationBegin( const char *szId, int nInputCount, int nOutputCount, int nParameterCount, const char ** pszParameters,
		 int *pnError, char * szErrorMsg, void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr)
{
	double x_ant[BUFF]; x_anterior = matriz(x_ant, LC_MAX, LC_MAX);
	static double xki[5]; x_ki = matriz(xki, 5, 1);
	n_estados = x_ki.linhas;
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
	Delta_T=1.0/300e3;
	Ts = Delta_T;
	a = 0.5;

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
	double ad[BUFF];   Ad = matriz(ad, LC_MAX, LC_MAX);
	double bd[BUFF];   Bd = matriz(bd, LC_MAX, LC_MAX);
	double cd[] = {0, 0, 0, 1};   Cd = matriz(cd, 1, 4); 
	zoh(&A_planta, &B_planta, Ts, &Ad, &Bd); 
	Bd.print(&Bd);
	
	//estende-se o modelo
	double aa[BUFF]; A = matriz(aa, LC_MAX, LC_MAX);
	double bb[BUFF]; B = matriz(bb, LC_MAX, LC_MAX);
	double cc[BUFF]; C = matriz(cc, LC_MAX, LC_MAX);

	ret_A(&Ad, &Cd, &A);
	ret_B(&Bd, &Cd, &B);
	ret_C(&Cd, &C);

	Bd.print(&Bd);

	double kmpc[BUFF]; Kmpc = matriz(kmpc, 1, n_estados);
	double n_k[BUFF]; Neg_kmpc = zeros(Kmpc.linhas, Kmpc.colunas);
	Kmpc_func(&A, &B, &C, &Kmpc, Np, Nc, n_estados, rw, a);
	
	mult(&Kmpc, &Neg_kmpc, -1);
}


void SimulationEnd(const char *szId, void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr)
{
}

//===============================================================================//
//                                                    FIM DE MAIN                                                            //
//===============================================================================//


//###############################################################################//
//                                             DMPC COM LAGUERRE                                                    //
//###############################################################################//

void Psi_func(Matriz* A, Matriz* B, Matriz* C, Matriz* Psi, int Np,  int Nc) {
	//n_estados deve ser passado como argumento, ou ele já existe na forma de outro argumento?
	//Matrizes padrão
	int n_estados = A->linhas;

	B->print(B);
	double buff1[BUFF]; Matriz Al = matriz(buff1, Nc, Nc);			   	   	     Al_func(0.5, Nc, &Al);      //mudar primeiro argumento
	double buff2[BUFF]; Matriz L0 = matriz(buff2, Nc, 1);                        L0_func(0.5, 5, &L0);     //mudar primeiro e segundo argumento
	double buff3[BUFF]; Matriz Sc1 = matriz(buff3, LC_MAX, LC_MAX);  Sc1_func(B, &L0, &Sc1);
	double sc_ant[BUFF]; Matriz Sc_anterior = matriz(sc_ant, n_estados, Nc); cop_mat(&Sc1, &Sc_anterior);
	double c_t[BUFF]; Matriz C_T = matriz(c_t, LC_MAX, LC_MAX); transposta(C, &C_T);                                                                           
	double q[BUFF]; Matriz Q = matriz(q, n_estados, n_estados); mat_mult(&C_T, C, &Q);      //Q 
	double som[BUFF]; Matriz Somatorio = matriz(som, Nc, n_estados);

	int mm;
	for (mm = 1; mm<=Np; mm++) { 
		double phi_mt[BUFF]; Matriz Phi_m_t = matriz(phi_mt, n_estados, Nc);

		if (mm == 1) cop_mat(&Sc1, &Phi_m_t); 
		else  phi_m_t(A, &Sc_anterior, &Sc1, &Al, &Phi_m_t, mm); //é?
		
		double phi_m[BUFF];  Matriz Phi_m = matriz(phi_m, Nc, n_estados);
		transposta(&Phi_m_t, &Phi_m);                                                              //phi(m)
		double m0[BUFF]; Matriz M0        = matriz(m0, Nc, n_estados); mat_mult(&Phi_m, &Q, &M0);
		double a_pot[BUFF]; Matriz A_pot = matriz(a_pot, A->linhas, A->colunas); mat_pot(A, &A_pot, mm);
		double m1[BUFF]; Matriz M1        = matriz(m1, Nc, n_estados); mat_mult(&M0, &A_pot, &M1);
		double m2[BUFF]; Matriz M2        = matriz(m2, Nc, n_estados); sm(&M1, &Somatorio, &M2);
		cop_mat(&M2, &Somatorio);
		cop_mat(&Phi_m_t, &Sc_anterior);                                                          //Sc(m-1)
	}
	cop_mat(&Somatorio, Psi);
	
}

void Kmpc_func (Matriz* A, Matriz* B, Matriz* C, Matriz* Kmpc, int Np,  int Nc, int n_estados, double rw, double a) {
	double buff1[BUFF]; Matriz L0 = matriz(buff1, Nc, 1);    L0_func(a, Nc, &L0); //mudar os argumentos
	double psi[BUFF]; Matriz Psi = matriz(psi, Nc, n_estados);
	double omega[BUFF]; Matriz Omega = matriz(omega, Nc, Nc);
	Omega_func(A, B, C, &Omega, Np,  Nc, rw);
	
	
	Psi_func(A, B, C, &Psi, Np, Nc);
	double omega_inv[BUFF]; Matriz Omega_inversa = matriz(omega_inv, Nc, Nc);
	inv(&Omega, &Omega_inversa);	
	double lo_t[BUFF];  Matriz L0_T = matriz(lo_t, 1, Nc); transposta(&L0, &L0_T);
	double m0[BUFF];  Matriz M0 = matriz(m0, 1, n_estados); mat_mult(&Omega_inversa, &Psi, &M0);
	mat_mult(&L0_T, &M0, Kmpc);
	//printf("Valor de Psi\n");
	//Psi.print(&Psi);
}

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

void Al_func (double a, int Nc, Matriz* Al) {
	//printf("Teste dentro de Al\n");
	Matriz Z = zeros(Nc, Nc);
	Matriz Z2 = zeros(1, Nc);
	double B = 1 - a*a;
	int i;
	for (i=0; i<Nc; i++) {
		if (i == 0)  Z2.matriz[i] = a;
		else Z2.matriz[i] = pow(-1, i-1)*pow(a, i-1)*B;   
	}

	int j;
	int k;
	int mult = 0;
	int contador = 0;
	for(j = 0; j <Nc; j++) {
		for(k = j; k >=0; k--) {
			Z.matriz[Nc*j + k] = Z2.matriz[contador++];
		}
		contador=0;
	}
	cop_mat(&Z, Al);
}

void L0_func(double a, int Nc, Matriz* L0) {
	Matriz Z = zeros(Nc, 1);
	int i;
	double B = sqrt(1 - a*a);
	for(i = 0; i < Nc; i++) {
		Z.matriz[i] = pow(-a, i)*B;
	}
	cop_mat(&Z, L0);
}

void Sc1_func(Matriz* B, Matriz* L0, Matriz* Res) {
		double m0[BUFF]; Matriz M0 = matriz(m0, LC_MAX, LC_MAX);
		transposta(L0, &M0);
		double m1[BUFF]; Matriz M1 = matriz(m1, LC_MAX, LC_MAX);
		mat_mult(B, &M0, &M1);
		cop_mat(&M1, Res);
}

void phi_m_t(Matriz* A, Matriz* Sc_anterior,  Matriz* Sc1, Matriz* Al, Matriz* Res, int m) {
	double buff0[BUFF]; Matriz M0 = matriz(buff0, LC_MAX, LC_MAX);
	double buff1[BUFF]; Matriz M1 = matriz(buff1, LC_MAX, LC_MAX);
	double buff2[BUFF]; Matriz M2 = matriz(buff2, LC_MAX, LC_MAX);
	double buff3[BUFF]; Matriz M3 = matriz(buff3, LC_MAX, LC_MAX);
	double buff4[BUFF]; Matriz M4 = matriz(buff4, LC_MAX, LC_MAX);
	
	mat_mult(A, Sc_anterior, &M0); //A*Sc(m-1)
	mat_pot(Al, &M1, m-1);                //Al^{m-1}
	transposta(&M1, &M2);             //(Al*{m-1})^T
	mat_mult(Sc1, &M2,&M3);        //Sc1*(Al*{m-1})^T
	sm(&M0, &M3, &M4);               //A*Sc(m-1) + Sc1*(Al*{m-1})^T
	cop_mat(&M4, Res);
}

void Omega_func(Matriz* A, Matriz* B, Matriz* C, Matriz* Om, int Np,  int Nc, double rw) {
	//n_estados deve ser passado como argumento, ou ele já existe na forma de outro argumento?
	//Matrizes padrão
	int n_estados = A->linhas;
	double buff1[BUFF]; Matriz Al = matriz(buff1, Nc, Nc);			   	   	     Al_func(0.5, Nc, &Al);      //mudar primeiro argumento
	double buff2[BUFF]; Matriz L0 = matriz(buff2, Nc, 1);                        L0_func(0.5, 5, &L0);     //mudar primeiro e segundo argumento
	double buff3[BUFF]; Matriz Sc1 = matriz(buff3, n_estados, Nc);  Sc1_func(&B, &L0, &Sc1);
	double sc_ant[BUFF]; Matriz Sc_anterior = matriz(sc_ant, n_estados, Nc); cop_mat(&Sc1, &Sc_anterior);
	double c_t[BUFF]; Matriz C_T = matriz(c_t, LC_MAX, LC_MAX); transposta(C, &C_T);                                                                           
	double q[BUFF]; Matriz Q = matriz(q, n_estados, n_estados); mat_mult(&C_T, C, &Q);      //Q 
	double som[BUFF]; Matriz Somatorio = matriz(som, Nc, Nc);

	int mm;
	for (mm = 1; mm<=Np; mm++) { 
	//Somatório. Talvez o limite superior esteja errado
		//phi_m_t(Matriz* A, Matriz* Sc_anterior,  Matriz* Sc1, Matriz* Al, Matriz* Res, int m) 
		double phi_mt[BUFF]; Matriz Phi_m_t = matriz(phi_mt, n_estados, Nc);

		if (mm == 1) cop_mat(&Sc1, &Phi_m_t); 
		else  phi_m_t(A, &Sc_anterior, &Sc1, &Al, &Phi_m_t, mm); //é?
		
		double phi_m[BUFF];  Matriz Phi_m = matriz(phi_m, Nc, n_estados);
		transposta(&Phi_m_t, &Phi_m);                                                              //phi(m)
		double m0[BUFF]; Matriz M0 = matriz(m0, Nc, n_estados); mat_mult(&Phi_m, &Q, &M0);
		double m1[BUFF]; Matriz M1 = matriz(m1, Nc, Nc); mat_mult(&M0, &Phi_m_t, &M1);
		double m2[BUFF]; Matriz M2 = matriz(m2, Nc, Nc); sm(&M1, &Somatorio, &M2);
		cop_mat(&M2, &Somatorio);
		cop_mat(&Phi_m_t, &Sc_anterior);                                                          //Sc(m-1)
	}

	double i[BUFF]; Matriz I = matriz(i, Nc, Nc); ident(&I, Nc);
	double rl[BUFF]; Matriz RL = matriz(rl, Nc, Nc); mult(&I, &RL, rw);  //RL
	double soma[BUFF]; Matriz Sm = matriz(soma, Nc, Nc); sm(&Somatorio, &RL, &Sm);
	cop_mat(&Sm, Om);
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



