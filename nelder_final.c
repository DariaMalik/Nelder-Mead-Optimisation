
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define alpha 1
#define gamma 2
#define rho 0.5
#define sigma 0.5
#define pi 3.1415


//-----------------------------------------------------------------------------//
//---------------------------------TYPEDEF-------------------------------------//
//-----------------------------------------------------------------------------//


typedef float point[3];
typedef point simplex[4];
typedef int valeurs_f[4];


//-----------------------------------------------------------------------------//
//--------------------------PROTOTYPES DES FONCTIONS---------------------------//
//-----------------------------------------------------------------------------//


// remplir un tableau u[] de 10^a à 10^b avec n points échantillonnés logarithmiquement
void logspace(float a, float b, int n, double *u);

// calculer le module d'un nombre complexe
float module (float rel, float freq, float omega_x);

// calculer le gain d'une fonction de transfert aux valeurs des composants données A, w0, w1, à une fréquence donné f
float gain(float f, float A,float w0,float w1);

// calculer le nombre de points hors gabarit pour la fonction de transfert aux valeurs des composants données
int hors_gabarit (point X, double *frequence, int N);

// remplit un tableau V[4] avec des nombres de points hors gabarit pour chaque sommets du simplex
// sachant que chaque sommet du simplex contient A, w0, w1 auxquels est associée une unique fonction de transfert
void fonctiona_du_simplex (simplex S, valeurs_f V, double *frequence, int N);

// trier les nombres de points hors gabarit dans l'ordre croissant (le meilleur point a l'indice 0, le pire - 3 (4 points)
// et trier en même temps les coordonnées des sommets correspondants à ces valeurs
void sorting_hat (simplex S, valeurs_f V);

// calcul le centre de gravite du simplex sans prendre en compte le pire point qui va être modifié
void gravity (simplex S, point Xg);

// genere un nombre aleatoire entre 0 et 20 pour les valeurs de A du simplex initial (dans les limites physique possibles de ces valeurs)
int randa(int min,int max);

//genere un nombre aleatoire entre 0 et 1000 pour les valeurs de w0 et w1 du simplex initial (dans les limites physique possibles de ces valeurs)
int randw(int min,int max);

// on ne peut pas avoir A, w0, w1 négatifs, ainsi que w0 et w1 nulles car se sont les pulsations
// on va donc mettre une contrainte pour Nelder-Mead et recommencer la recherche d'une solution tant qu'une des valeurs dans les sommets du simplex soit négative ou nulle
// (A n'est pas nul non plus car le gabarit passe par le gain positif qu'on ne peut obtenir qu'avec A positif
// la fonction retourne 1 si tous les A, w0, w1 du simplex sont non-nuls et positifs
// 0 sinon
int no_negatif (simplex S);

// creation du fichier avec des frequences et valeurs du gain associées pour 0 points hors gabarit afin de tracer Bode sur Octave
void fichier_bode (double *frequence, simplex S, int n);

// création des sommets initiaux du simplex de départ pour commencer Nelder-Mead
void creation_sommets (simplex S);

// calcul Nelder_Mad
void Nelder_Mead (simplex S, double *frequences, int n);


//-----------------------------------------------------------------------------//
//-----------------------------CORPS DES FONCTIONS-----------------------------//
//-----------------------------------------------------------------------------//


void logspace(float a, float b, int n, double *u){
    
	double c;
    int i;

    // le pas "logarithmique"
    c = (b - a)/(n - 1);

    for(i = 0; i < n -1; i++){
        u[i] = pow(10., a + i*c);
        //printf("%lf\n",u[i]);

    }
    
    // le dernier point
    u[n - 1] = pow(10., b);
    //printf("%lf\n",u[n-1]);

}


float module (float rel, float freq, float omega_x){
	if (rel!=0){
	//printf ("w (%f Hz) = %f\n", freq, 2*pi*freq);
		float b = (2*pi*freq)/omega_x;
		return sqrt((rel*rel) + (b*b));
	}
	else{
		return ((2*pi*freq)/omega_x);
	}
}


float gain(float f, float A,float w0,float w1){

	double T;
	double G;

	T = (A * module(0,f,w1)) / (module(1,f,w1) * module(1,f,w0));
	G = 20*log10(T);
	//printf ("|T| = %lf\n", T);	
	//printf ("G = %f\n", 20*log10(T));
	return G;
	
}


int hors_gabarit (point X, double *frequence, int N){
	
	int i;
    int nbpts=0; // nbpts dans gabarit
    
    double *G;
    G = (double*)malloc(N*sizeof(double));
	
    for(i=0;i<N;i++)
	{
		G[i]=gain(frequence[i],X[0],X[1],X[2]);
		//printf("%.2f Hz	=>	%.2f Db\n", frequence[i], G[i]);
	}

    for(i=0; i<N; i++)
    {
        if((frequence[i]<100 && G[i]<-15)||(frequence[i]<7000 && frequence[i] > 2000 && G[i]>5)||(frequence[i]>200000 && G[i]<-20)){
            nbpts+=1;
        }
		if( (frequence[i]>=100 && frequence[i]<=2000) || (frequence[i]>=7000 && frequence[i]<=200000) ){
			nbpts+=1;
		}
    }
    
    free (G);
    
    return (N-nbpts);
    
}


void fonctiona_du_simplex (simplex S, valeurs_f V, double *frequence, int N){
	int i;
	for (i=0; i<4; i++){
		V[i] = hors_gabarit( S[i], frequence, N );
	}
}


void sorting_hat (simplex S, valeurs_f V) {
	int i, j, k, pass, decal;
	valeurs_f poubelle;
	simplex coord;
	
	for (i=0;i<4; i++){
		pass = 0;
		for (j=0; j<4; j++){
			if (j!=i){
				if (V[i]>V[j]) pass++;
			}
		}
		if (poubelle[pass]!=V[i]){
			poubelle[pass] = V[i];
			for (k=0; k<3; k++) coord[pass][k] = S[i][k];
		}
		else{
			decal=pass;
			while (poubelle[decal]==V[i]){
				decal++;
			}
			poubelle[decal] = V[i];
			for (k=0; k<3; k++) coord[decal][k] = S[i][k];
		}		
	}
    for (i=0; i<4; i++){
    	V[i] = poubelle[i];
    	for (j=0; j<3; j++) S[i][j] = coord[i][j] ;
	}
}


void gravity (simplex S, point Xg){
	int i, j;

	for (j=0; j<3; j++){
		for (i=0;i<3;i++){
		//printf ("Xg = %f\n", Xg[j]);
			Xg[j] += S[i][j];
		}
		//printf ("%f\n", Xg[j]);
		Xg[j] = Xg[j]/3;
		//printf (" coord = %f\n ", Xg[j]);
		//puts ("");
	}
}


int randa(int min,int max){
    int nbgen=rand()%(max-min+1)+min;
    return nbgen;
}


int randw(int min,int max){
	int nbgen=rand()%(max-min+1)+min;
    return nbgen;
}


int no_negatif ( simplex S){
	
	int res=1, i, j;
	
	for (i=0; i<4; i++){
		for (j=0; i<3; i++){
			if (S[i][j]<=0) res=0;
		}
	}
	
	return res;
}


void fichier_bode (double *frequence, simplex S, int n){

	int i;
	float *gain_super;	
	gain_super = (float*)malloc(n*sizeof(float));
	
	FILE *trace;
	trace = fopen("gain.txt", "w");
	
	for (i=0;i<n;i++){
		gain_super[i] = gain (frequence[i], S[0][0], S[0][1], S[0][2]);
		fprintf (trace, "%f\t%f\n", frequence[i], gain_super[i]);
	}
	
	fclose(trace);
	free (gain_super);

}


void creation_sommets (simplex S){
	int i;
	for (i=0; i<4; i++){
		S[i][0]=randa(0,100);													//A
		S[i][1]=randw(0,10000);												//w0
		do{
			S[i][2]=randw(0,10000);											//w1
		}while(S[i][2]<S[i][1]);
	
		//printf ("A=%f w0=%f w1=%f\n", S[i][0], S[i][1], S[i][2]);
	}
}


void Nelder_Mead (simplex S, double *frequences, int n){

	point G={0}, R={0}, E={0}, C={0};																					// R - reflechi, E - expansion, C - contraction
	int iter, i, j;
	float Fr;
	valeurs_f F;

	do{

		creation_sommets (S);

		iter = 0;
		
		while (iter<20) {

			//printf ("\n");

			fonctiona_du_simplex (S, F, frequences, n);																	// etape 2 : valeurs de la fonction aux sommets du simplex
			sorting_hat (S, F);																							// etape 2 : trie des valeurs et des coordonnées correspondantes
			
			//for (i=0; i<4; i++) printf ("%d\n", F[i]);																// possible d'afficher les nombres de points hors gabarits pour chaque sommets du simplex
			
			gravity (S, G);																								// etape 3 : calcul du centre de gravité du simplex sans le 4ième point
			//printf ("le centre de gravite est %f %f %f\n", G[0], G[1], G[2]);
		
			for (i=0; i<3; i++) R[i] = G[i] + alpha*(G[i] - S[3][i]);													// etape 4 : calcul des coordonées du point de réflexion
			//printf ("le point reflexion est %f %f %f\n", R[0], R[1], R[2]);
		
			Fr = hors_gabarit (R, frequences, n);																		// nombre de points hors gabarit pour le point réflechi
			//printf ("f (R) = %f\n", Fr);

			if (Fr>=F[0] && Fr<F[2]){
				for (i=0; i<3; i++) S[3][i] = R[i];
				//printf ("\nReflection\n");
			} // RETOUR A L'ETAPE 2 SI EFFECTUÉ

			if (Fr<F[0]){
				for (i=0; i<3; i++) E[i] = G[i] + gamma*(R[i] - G[i]);
				if (hors_gabarit(E, frequences, n)<=Fr){
					for (i=0; i<3; i++) S[3][i] = E[i];
					//printf ("\nExpansion\n");
				}
				else{
					for (i=0; i<3; i++) S[3][i] = R[i];
					//printf ("\nReflection\n");
				}
			} // RETOUR A L'ETAPE 2 SI EFFECTUÉ

			if (Fr>=F[2]){
				for (i=0; i<3; i++) C[i] = G[i] + rho*(S[3][i] - G[i]);
				if (hors_gabarit(C, frequences, n)<F[3]){
					for (i=0; i<3; i++) S[3][i] = C[i];
					//printf ("\nContraction\n");
				} // RETOUR A L'ETAPE 2 SI EFFECTUÉ

				else{
					for (i=1; i<4; i++){
						for (j=0; j<3; j++){
							//printf ("%.1f ---> ", S[i][j]);
							S[i][j] = S[0][j] + sigma*(S[i][j]-S[0][j]);
							//printf ("%.1f\n", S[i][j]);
						}
						//puts("");
					}
				}
			} // RETOUR A L'ETAPE 2 SI EFFECTUÉ
			
			iter++;
			
		}

	}while (no_negatif(S)==0 || hors_gabarit(S[0], frequences, n)!=0 || S[0][2]<S[0][1]);

}


//-----------------------------------------------------------------------------//
//-----------------------------------MAIN--------------------------------------//
//-----------------------------------------------------------------------------//


int main (void){
	
	//*************************************************************************************************************************************************//
	puts("\n*************************************************************************");
	simplex sommets;																										
	int i, j, n;
	float a, b;																												
	double *freq;														

	//*************************************************************************************************************************************************//

	srand(time(NULL)); 																									// initialisation de rand
	
	//*************************************************************************************************************************************************//
	
	a = 1; 																												// puissance de la première decade
	b = 6;																												// puissance de la dernière decade
	printf ("Combien de points voulez vous?\n");
	scanf ("%d",&n);

	freq = (double*)malloc(n*sizeof(double));

    logspace(a,b,n,freq);
    
	//*************************************************************************************************************************************************//
	
	/*for (i=0; i<4; i++){
		printf ("Coordonnes du sommet separes par espace (A w0 w1): ");
		scanf ("%f %f %f", &sommets[i][0], &sommets[i][1], &sommets[i][2]);
		printf ("%f %f %f\n", sommets[i][0], sommets[i][1], sommets[i][2]);
	}*/
	
	//*************************************************************************************************************************************************//

	Nelder_Mead (sommets, freq, n);

	//*************************************************************************************************************************************************//
	
	puts("*************************************************************************");
	printf ("Le nombre min de points hors gabarit = %d / %d\n", hors_gabarit(sommets[0], freq, n), n);
	puts("*************************************************************************");
	printf ("A = %f\nw0 = %f\nw1 = %f\n", sommets[0][0],sommets[0][1],sommets[0][2]);
	
	//*************************************************************************************************************************************************//

	/*int ilow, imed1, imed2, ihigh;
	
	ilow = 0;
	while (freq[ilow]<100){
		ilow++;
	}
	imed1 = ilow;
	while (freq[imed1]<2000){
		imed1++;
	}
	imed2 = imed1;
	while (freq[imed2]<7000){
		imed2++;
	}
	ihigh = imed2;
	while (freq[ihigh]<200000){
		ihigh++;
	}*/
	
	puts("*************************************************************************");
	//printf("100 Hz -> %f\n2/7 kHz -> %f %f\n200 kHz -> %f\n", gain (freq[ilow], sommets[0][0], sommets[0][1], sommets[0][2]), gain (freq[imed1], sommets[0][0], sommets[0][1], sommets[0][2]), gain (freq[imed2], sommets[0][0], sommets[0][1], sommets[0][2]), gain (freq[ihigh], sommets[0][0], sommets[0][1], sommets[0][2]));
	
	
	//*************************************************************************************************************************************************//
	
	fichier_bode (freq, sommets, n);	
	
	//*************************************************************************************************************************************************//

	free(freq);
	
	//*************************************************************************************************************************************************//
	
	return 0;
}
