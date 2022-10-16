#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define A 1.0
#define B 1.0
#define C 1.0
#define D 1.0
#define E 1.0
#define F 1.0

#define MAX_ITER 100

double g(double x, double y);
double p(double x, double y);
void assignValues(double x[], double y[], double** u, double** old_u, double h, int n);
void calculateSolution(double** s, double x[], double y[], int n);
void printMatrix(double** u, int n);
double matrixInfNorm(double** m1, double** m2, int n);
int jacobi(double** u, double** old_u, double x[], double y[], double a[9], int n, double error);
double calculateErrorBound(double delta, double old_delta);

int main(){
    /*
    Valors a canviar segons el que es vulgui aconseguir
    */
    const double error = 0.01;           //Valor de l'error que accepta
    const bool showMatrices = false;     //Valor que determina si es mostren les matrius a cada iteració
    const int method_iter = 10;           //Nombre d'iteracions amb n diferents que farà el programa

    //valor de n, es canvia més endavant
    int n;
    //nombre d'iteracions que ha fet el mètode iteratiu en cada iteració
    int total_iter;

    printf("\n---------------------------------------------------------\n");
    printf("| n i h \t\t|");
    printf("| Norma infinit|");
    printf("| Iteracions \t|");

    for(int n_iter=0; n_iter<method_iter; n_iter++){
        n = (n_iter + 1) * 10;    //valor de n segons la iteració (n = 10,20,30..., 500)

        /*
        Crea els vectors i matrius a usar
        */
        //vectors que representen els punts x, y on avaluem la solució
        double x[n+2];
        double y[n+2];
        //vector solució en forma de matriu i old_u co a iteracio per jacobi; i s la solució exacta
        double** s = (double**)calloc((n+2),sizeof(double*));
        double** u = (double**)calloc((n+2),sizeof(double*));   //representa la solució trobada en la iteració actual
        double** old_u = (double**)calloc((n+2),sizeof(double*));  //representa la solució en la iteració anterior
        for(int i=0; i<(n+2); i++){
            u[i] = (double*)calloc((n+2),sizeof(double));
            old_u[i] = (double*)calloc((n+2),sizeof(double));
            s[i] = (double*)calloc((n+2),sizeof(double));
        } 

        double h = 1/(double)(n+1); //pas de discretització

        //assigna valors als vectors
        assignValues(x, y, u, old_u, h, n);

        //Per no haver de refer els càlculs a cada iteració de Jacobi, guarda els coeficients:
        double a[9] = {B/(4*h*h), C/(h*h)-E/(2*h), -B/(4*h*h), A/(h*h)-D/(2*h), (-2*A-2*C)/(h*h)+F, A/(h*h)+D/(2*h), -B/(4*h*h), C/(h*h)+E/(2*h), B/(4*h*h)};

        //resol per Jacobi
        total_iter = jacobi(u, old_u, x, y, a, n, error); 

        //Calcula la solució esperada
        calculateSolution(s, x, y, n);
        if(showMatrices){
            printf("\nSolució esperada:\n");
            //Imprimeix la solució esperada
            printMatrix(s, n+2);
            printf("\nSolució calculada:\n");
            //Imprimeix la solució calculada
            printMatrix(u, n+2);
        }
        
        //Imprimeix els resultats d'aquesta iteració
        printf("\n---------------------------------------------------------\n");
        printf("| n = %d, h = %lf\t|", n, h);
        printf("| %lf \t|", matrixInfNorm(s, u, n+2));
        printf("| %d  \t|", total_iter);

        /*Abans de la següent iteració, alliberem els vectors
        Els haurem de tornar a crear, ja que tindran mides diferents*/
        free(u);
        free(old_u);
        free(s);
    } 
    printf("\n---------------------------------------------------------\n");
    return 0;
}

/*
---------------------------------------------------------------------------
FUNCIONS AUXILIARS
---------------------------------------------------------------------------
*/

/*
Funció p: valors coneguts de la frontera
*/
double p(double x, double y){
    return x*y;
}

/*
Funció g: valors coneguts
*/
double g(double x, double y){
    return x*y + y + x + 1;
}

/*
Funció que implementa el mètode de Jacobi. 
Modifica u fins trobar la solució aproximada i retorna el nombre d'iteracions que ha fet
La n aquí representa el nombre fixat per aquesta iteració de Jacobi
*/
int jacobi(double** u, double** old_u, double x[], double y[], double a[9], int n, double error){
    int iter = 0;   //Controla en quina iteració es troba
    double delta, old_delta = 0;
    double error_aprox;
    do{
        //Càlcul de la nova solució
        for(int i=1; i<(n+1); i++){
            for(int j=1; j<(n+1); j++){
                u[i][j] = (g(x[i], y[j]) - a[0]*old_u[i-1][j-1] - a[1]*old_u[i][j-1] - a[2]*old_u[i+1][j-1] - a[3]*old_u[i-1][j] - a[5]*old_u[i+1][j] - a[6]*old_u[i-1][j+1] - a[7]*old_u[i][j+1] - a[8]*old_u[i+1][j+1])/a[4];
            }
        }
        //Càlcul de la fita de l'error i actualització de deltes
        delta = matrixInfNorm(u, old_u, n+2);    
        error_aprox = calculateErrorBound(delta, old_delta);
        old_delta = delta;

        //actualització de la solució de l'a "iteració anterior"
        for(int i=1; i<(n+1); i++){
            for(int j=1; j<(n+1); j++){
                old_u[i][j] = u[i][j];
            }
        }
        iter++;
    }while((iter < MAX_ITER) && (error_aprox>error));
    return iter;
}

/*
Funció que assigna valors als vectors x, y, u i old_u
n és el nombre de punts, fixat 
*/
void assignValues(double x[], double y[], double** u, double** old_u, double h, int n){
    for(int i=0; i<(n+2); i++){
        x[i] = i*h;
        y[i] = i*h;
        //Asignem els valors de la frontera a mesura que calculem x_i, y_i amb la funcio p
        if(i==0){
            u[0][0] = p(x[0], y[0]);
            old_u[0][0] = p(x[0], y[0]);
        }
        else if(i==(n+1)){
            for(int j=0; j<n+2; j++){
                u[j][i] = p(x[j], y[i]);
                u[i][j] = p(x[i], y[j]);
                old_u[j][i] = p(x[j], y[i]);
                old_u[i][j] = p(x[i], y[j]);
            }
        }
        else{
            u[0][i] = p(x[0], y[i]);
            u[i][0] = p(x[i], y[0]);
            old_u[0][i] = p(x[0], y[i]);
            old_u[i][0] = p(x[i], y[0]);
        }
    }
}

/*
Guarda a s la solució esperada per segons la p escollida
n és el nombre de punts de l'interior (el que varia en cada iteració del programa)
*/
void calculateSolution(double** s, double x[], double y[], int n){
    for(int i=0; i<(n+2); i++){
        for(int j=0; j<(n+2); j++){
            s[i][j] =  p(x[i], y[j]);
        }
    }
}

/*
Funció que calcula una cota de l'error segons vist a teoria donades delta^(k+1) i delta(k)
*/
double calculateErrorBound(double delta, double old_delta){
    return fabs(delta * delta / (old_delta - delta));
}

/*
Imprimeix una matriu nxn
n és la mida de la matriu
*/
void printMatrix(double** a, int n){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
        printf("%lf ", a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

/*
Funció que calcula la norma infinit de 2 vectors de dimensió n
*/
double vectorInfNorm(double* v1, double* v2, int n){
    double max_abs_val = 0;
    double current;
    for(int i=0; i<n; i++){
        current = fabs(v1[i]-v2[i]);
        if(max_abs_val<current){
            max_abs_val = current;
        }
    }
    return max_abs_val;
}

/*
Funció que calcula la norma infinit de dos vectors escrits com a matriu nxn
n és la mida de la matriu
*/
double matrixInfNorm(double** m1, double** m2, int n){ 
    //TODO es pot començar un despres i acabar un abans perque nomes l'interior afecta
    double max_abs_val = 0;
    double current;
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            current = fabs(m1[i][j]-m2[i][j]);
            if(current>max_abs_val){
                max_abs_val = current;
            }
        }
    }
    return max_abs_val;
}

/*
---------------------------------------------------------------------------
RESPOSTES A LES PREGUNTES DE LA PRÀCTICA
---------------------------------------------------------------------------
*/

/*
Resultat fixant 100 iteracions per controlar l'error, amb error tolerat de 0.01
---------------------------------------------------------
| n i h                 || Norma infinit|| Iteracions   |
---------------------------------------------------------
| n = 10, h = 0.090909  || 0.010207     || 100          |
---------------------------------------------------------
| n = 20, h = 0.047619  || 0.161963     || 100          |
---------------------------------------------------------
| n = 30, h = 0.032258  || 0.297176     || 100          |
---------------------------------------------------------
| n = 40, h = 0.024390  || 0.397221     || 100          |
---------------------------------------------------------
| n = 50, h = 0.019608  || 0.472266     || 100          |
---------------------------------------------------------
| n = 60, h = 0.016393  || 0.531572     || 100          |
---------------------------------------------------------
| n = 70, h = 0.014085  || 0.577665     || 100          |
---------------------------------------------------------
| n = 80, h = 0.012346  || 0.616011     || 100          |
---------------------------------------------------------
| n = 90, h = 0.010989  || 0.647087     || 100          |
---------------------------------------------------------
| n = 100, h = 0.009901 || 0.673599     || 100          |
---------------------------------------------------------
*/

/*Canvia els valors de la frontera als valors que té u
    for(int i=0; i<(N+2); i++){
            old_u[i][0] = u[i][0];
            old_u[i][N+1] = u[i][N+1];
            old_u[0][i] = u[0][i];
            old_u[N+1][i] = u[N+1][i];
        }
    */