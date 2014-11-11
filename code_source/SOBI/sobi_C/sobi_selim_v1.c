#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

main(int argc, char** argv)
{
	//variables d'increment
	int k,i,j = 0;
	
	//ETAPE 1 - Récuperation des données
        
	/** Valeurs de test / debuggage
	double Melange1[6] = {1,2,3,1,1,1};
	double Melange2[6] = {4,5,6,2,1,5};
	double Signal1[6] = {4,5,6,9,0,8};
	double Signal2[6] = {4,5,6,2,3,3};

	const int tailleSignal = 6; //longueur des vecteurs Melange
	//taille des vecteurs enfants
	const int sizeChild = 2;
	//nombre de vecteurs enfants
	const int nbX = tailleSignal - sizeChild + 1;
	printf("nbX : %d\n", nbX);
	**/

	//calcul du temps d'execution du programme
	clock_t start_t, end_t;
	double total_t;
	start_t = clock();
/// debut de chargement des donnees
        printf("\n*** importation des fichiers signaux ***\n");
	
	FILE *fr1 = fopen("../../data/signalReference1.txt","r");
	FILE *fr2 = fopen("../../data/signalReference2.txt","r");
	FILE *fm1 = fopen("../../data/signalMelange1.txt","r");
	FILE *fm2 = fopen("../../data/signalMelange2.txt","r");

        if(fr1 == NULL || fr2 == NULL || fm1 == NULL || fm2 == NULL)
        {
		printf("erreur lors du chargement des fichiers...\n");
		return 1;
	}

	//on verifie que tous les fichiers contiennent le meme nombre d'echantillons

        fseek(fr1, 0, SEEK_END);
        fseek(fr2, 0, SEEK_END);
        fseek(fm1, 0, SEEK_END);
        fseek(fm2, 0, SEEK_END);
        const int NbEchr1 = ftell(fr1)/sizeof(double);
        const int NbEchr2 = ftell(fr2)/sizeof(double);
        const int NbEchm1 = ftell(fm1)/sizeof(double);
        const int NbEchm2 = ftell(fm2)/sizeof(double);

	if(NbEchr1 != NbEchr2 || NbEchr2 != NbEchm1 || NbEchm1 != NbEchm2)
	{
        	printf("Erreur : les fichiers ne contiennent pas tous le meme nombre d'echantillons...\n");
		return 1;
	}
	
	//longueur des vecteurs Melange
	const int tailleSignal = NbEchr1; 
	//taille des vecteurs enfants = CHOIX de la taille de fenetre
	const int sizeChild = 40;
	printf("\ntaille de la fenetre glissante : %d\n\n", sizeChild);
	//nombre de vecteurs enfants
	const int nbX = tailleSignal - sizeChild + 1;
	
	double Melange1[tailleSignal];
	double Melange2[tailleSignal];
	double Signal1[tailleSignal];
	double Signal2[tailleSignal];

        rewind(fr1);
        rewind(fr2);
        rewind(fm1);
        rewind(fm2);
        fread(Signal1, sizeof(double), tailleSignal, fr1);
        fread(Signal2, sizeof(double), tailleSignal, fr2);
        fread(Melange1, sizeof(double), tailleSignal, fm1);
        fread(Melange2, sizeof(double), tailleSignal, fm2);

        close(fr1);
        close(fr2);
        close(fm1);
        close(fm2);
        
/// fin de chargement des donnees

	

	//Matrices contenant des vecteurs enfants
	double** X1 = (double**)malloc(sizeChild*sizeof(double*));
	double** X2 = (double**)malloc(sizeChild*sizeof(double*));

	for(i = 0; i < sizeChild; i++)
	{
		X1[i] = (double*)malloc(nbX*sizeof(double));
		X2[i] = (double*)malloc(nbX*sizeof(double));
		memset(X1[i] , 0 , nbX*sizeof(double));
		memset(X2[i] , 0 , nbX*sizeof(double));
	}
	
	//remplissage des matrices contenant les enfants de Melange1 et Melange2
	for(k = 0; k < nbX ; k++)
	{
		for(i = 0; i < sizeChild ; i++)
		{
			X1[i][k] = Melange1[k+i];
			X2[i][k] = Melange2[k+i];
		}
	}
	

	//matrices d'intercovariance
	double** R11 = (double**)malloc(sizeChild*sizeof(double*));
	double** R22 = (double**)malloc(sizeChild*sizeof(double*));
	double** R12 = (double**)malloc(sizeChild*sizeof(double*));

	for(i = 0; i < sizeChild; i++)
	{
		R11[i] = (double*)malloc(sizeChild*sizeof(double));
		R22[i] = (double*)malloc(sizeChild*sizeof(double));
		R12[i] = (double*)malloc(sizeChild*sizeof(double));
		memset(R11[i] , 0 , sizeChild*sizeof(double)) ;
		memset(R22[i] , 0 , sizeChild*sizeof(double)) ;
		memset(R12[i] , 0 , sizeChild*sizeof(double)) ;
	}
	
	for(k = 0; k < nbX ; k++)
	{
		for(j = 0; j < sizeChild ; j++)
		{
			for(i = 0; i < sizeChild ; i++)
			{
				R11[i][j] = R11[i][j] + X1[j][k]*X1[i][k];
				R22[i][j] = R22[i][j] + X2[j][k]*X2[i][k];
				R12[i][j] = R12[i][j] + X1[j][k]*X2[i][k];
			}
		}
	}

	// calcul de l'esperance des matrice d'intercovariance 
	for(k = 0; k < sizeChild; k++)
	{
		for(i = 0; i < sizeChild; i++)
		{
			R11[i][k] /= nbX;
			R22[i][k] /= nbX;
			R12[i][k] /= nbX;
		}
	}

	//calcul de F et T
	double T1, T2, T12, F1, F2, F12 = 0;
	for(i = 0; i < sizeChild; i++)
	{
		T1 += R11[i][i]; 
		T2 += R22[i][i]; 
		T12 += R12[i][i]; 
	}

	for(i = 0; i < sizeChild ; i++)
		for(j = 0; j < sizeChild ; j++)
		{
			F1 += R11[i][j];
			F2 += R22[i][j];
			F12 += R12[i][j];
		}

	F1 = 1/((double)sizeChild*((double)sizeChild - 1))*(F1 - T1);
	F2 = 1/((double)sizeChild*((double)sizeChild - 1))*(F2 - T2);
	F12 = 1/((double)sizeChild*((double)sizeChild - 1))*(F12 - T12);

	//calcul des autres constantes
	double alpha = 2*F12*T12 - (F1*T2 +F2*T1);
	double beta = 2*(pow((T12),2) - T1*T2);
	double gamma2 = pow((F1*T2 - F2*T1),2) + 4*(F12*T2-T12*F2)*(F12*T1-T12*F1);
	double gamma = sqrt(gamma2);
	double d1 = alpha - gamma;
	double d2 = alpha + gamma;

	//finalement la matrice de melange estime [a b;c d]
	double a, b, c, d;
	a =  beta*F1-T1*d1;	
	b = beta*F12-T12*d2;	
	c = beta*F12 - T12*d1;	
	d = beta*F2 - T2*d2;	

	//inversion de la matrice de melange
	double detA = a*d - b*c;
	double invA[2][2];
	invA[0][0] = (1/detA)*d ;
	invA[0][1] = (1/detA)*(-b) ;
	invA[1][0] = (1/detA)*(-c) ;
	invA[1][1] = (1/detA)*a ;		

	//vecteurs estimes :
	double s1Estime[tailleSignal];
	double s2Estime[tailleSignal];
	
	for(k = 0; k < tailleSignal; k++)
	{
		s1Estime[k] = invA[0][0]*Melange1[k] + invA[0][1]*Melange2[k];
		s2Estime[k] = invA[1][0]*Melange1[k] + invA[1][1]*Melange2[k];
	}


	//calcul du temps d'execution du programme
	end_t = clock();
	total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
	printf("temps d'execution du programme : %f\n", total_t);
	



	//calcul des produits scalaires et normes
	double prod1, prod2, nEstime1, nEstime2, nSignal1, nSignal2;

	for(k = 0; k < tailleSignal; k++)
	{
		prod1 += s1Estime[k]*Signal1[k];
		prod2 += s2Estime[k]*Signal2[k];
		nEstime1 += pow(s1Estime[k],2);
		nEstime2 += pow(s2Estime[k],2);
		nSignal1 += pow(Signal1[k],2);
		nSignal2 += pow(Signal2[k],2);
	}

	nEstime1 = sqrt(nEstime1);
	nEstime2 = sqrt(nEstime2);
	nSignal1 = sqrt(nSignal1);
	nSignal2 = sqrt(nSignal2);

	//calcul des erreurs :
	double err1 = 10*log10(1 - pow(prod1/(nEstime1*nSignal1),2));
	double err2 = 10*log10(1 - pow(prod2/(nEstime2*nSignal2),2));

	printf("\n*** CALCUL DES ERREURS ***\n");
	printf("erreur sur le signal 1 = %.2lf dB\n", err1);
	printf("erreur sur le signal 2 = %.2lf dB\n", err2);
	printf("erreur sur le signal 1+2 = %.2lf dB\n\n", err1+err2);


	//enregistrement des signaux estimes sur fichier texte, facilement exportables
	//vers matlab pour ecoute (fichier csv)
	//potential problem : csvread cutting precision on matlab
	FILE* fp1 = fopen("../../data/ESTIME1","w");
	FILE* fp2 = fopen("../../data/ESTIME2","w");

	if(fp1 == NULL || fp2 == NULL)
	{
		printf("\nproblem creating files...\n");
		return 1;
	}
	else
	{
        printf("*** Enregistrement des signaux estimes dans les fichiers ESTIME1 et ESTIME2 (dossier data)\n\n");
		for(k = 0; k < tailleSignal; k++)
		{
			if(k < tailleSignal - 1)
			{
				fprintf(fp1, "%lf,", s1Estime[k]);	
				fprintf(fp2, "%lf,", s2Estime[k]);	
			}
			else
			{
				fprintf(fp1, "%lf", s1Estime[k]);	
				fprintf(fp2, "%lf", s2Estime[k]);	
			}
			
		}
		fclose(fp1);
		fclose(fp2);
	}

	//liberation de l'espace alloué pour les matrices
	//X1 et X2 :
	for(i = 0; i < sizeChild; i++)
	{
		free((void*)X1[i]);
		free((void*)X2[i]);
	}
	free((void*)X1);
	free((void*)X2);
	
	//R11 , R22 et R12 :
	for(i = 0; i < sizeChild; i++)
	{
		free((void*)R11[i]);
		free((void*)R22[i]);
		free((void*)R12[i]);
	}
	free((void*)R11);
	free((void*)R22);
	free((void*)R12);

	return 0;
}

