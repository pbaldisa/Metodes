#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define A 1.0
#define B 0.0
#define C 1.0
#define D 0.0
#define E 0.0
#define F 0.0

double g(double x, double y);
double p(double x, double y);
void assignValues(double x[], double a[9], double** u, double h, int n);
void assignValuesZero(double x[], double** u, double h, int n);
void calculateSolution(double** s, double x[], int n);
void printMatrix(double** u, int n);
double matrixInfNorm(double** m1, double** m2, int n);
int sor(double** u, double x[], double a[9], int n, double w, double error, int maxIter);
double calculateErrorBound(double delta, double old_delta);
void writeSolutionFile(double** u, double x[], int n);
double checkRealError(double** u, double x[], int n);

int main(){

    //Llegeix la n que vol l'usuari
    int n;              //Sistema nxn (punts de l'interior discretitzats)
    printf("n: ");
    scanf("%d", &n);

    //Llegeix l'error que vol l'usuari
    double error;    
    printf("Precisió (error): ");
    scanf("%lf", &error);

    //llegeix el nombre d'iteracions màxim permès
    int maxIter;
    printf("Màxim d'iteracions: ");
    scanf("%d", &maxIter);

    //legeix el factor de relaxació
    double w;
    printf("Factor de relaxació: ");
    scanf("%lf", &w);

    /*
    Crea els vectors i matrius a usar
    */
    //vectors que representen els punts x, y on avaluem la solució
    double x[n+2];
    //vector solució en forma de matriu co a iteracio per sor
    double** u = (double**)calloc((n+2),sizeof(double*));   //representa la solució trobada en la iteració actual
    for(int i=0; i<(n+2); i++){
        u[i] = (double*)calloc((n+2),sizeof(double));
    } 

    double h = 1/(double)(n+1); //pas de discretització

    //Per no haver de refer els càlculs a cada iteració de sor, guarda els coeficients:
    double a[9] = {B/(4*h*h), C/(h*h)-E/(2*h), -B/(4*h*h), A/(h*h)-D/(2*h), (-2*A-2*C)/(h*h)+F, A/(h*h)+D/(2*h), -B/(4*h*h), C/(h*h)+E/(2*h), B/(4*h*h)};

    //assigna valors als vectors
    //assignValues(x, a, u, h, n);
    assignValuesZero(x, u, h, n);

    //resol per sor
    printf("El nombre d'iteracions necessari ha estat: %d\n", sor(u, x, a, n, w, error, maxIter));

    //Escriu la solució en un fitxer
    writeSolutionFile(u, x, n);

    //Calcula l'error real, comparant amb la funció coneguda
    printf("L'error real ha estat de %E\n", checkRealError(u, x, n+2));

    //allibera els vectors
    for(int i=0; i<(n+2); i++){
        free(u[i]);
    } 
    free(u);
    
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
    return pow(x,5)*pow(y,4);
}

/*
Funció g: valors coneguts
*/
double g(double x, double y){
    return 12*pow(x,5)*pow(y,2) + 20*pow(x,3)*pow(y,4);
}

/*
Funció que implementa el mètode de sor. 
Modifica u fins trobar la solució aproximada i retorna el nombre d'iteracions que ha fet
La n aquí representa el nombre fixat per aquesta iteració de sor
*/
int sor(double** u, double x[], double a[9], int n, double w,  double error, int maxIter){
    int iter = 0;   //Controla en quina iteració es troba
    double error_aprox=0, old_error_aprox=1;
    double current;
    double error_q;
    do{
        error_aprox = 0;
        //Càlcul de la nova solució
        for(int i=1; i<(n+1); i++){
            for(int j=1; j<(n+1); j++){
                current = u[i][j];
                u[i][j] = (1-w)*u[i][j] + w*(g(x[i], x[j]) - a[0]*u[i-1][j-1] - a[1]*u[i][j-1] - a[2]*u[i+1][j-1] - a[3]*u[i-1][j] - a[5]*u[i+1][j] - a[6]*u[i-1][j+1] - a[7]*u[i][j+1] - a[8]*u[i+1][j+1])/a[4];
                //Si l'error en aquesta component és major, actualitzem l'error de la iteració
                if(error_aprox<fabs(current-u[i][j])){
                    error_aprox = fabs(current - u[i][j]);
                }
            }
        }
        //Càlcul de la fita de l'error i actualització de deltes   
        error_q = error_aprox / old_error_aprox;
        old_error_aprox = error_aprox;

        iter++;
        //Per no carregar la pantalla, imprimeix el quocient d'error cada 1000 iteracions
        if(iter%1000==0){
            printf("\nQuocient d'error a iteració %d: %lf\n", iter, error_q);
        }
    }while((iter < maxIter) && (error_aprox>error));

    printf("\nQuocient d'error a iteració %d: %lf\n", iter, error_q);

    return iter;
}

/*
Funció que assigna valors als vectors x, u 
Inicia el vector solucio amb x_i = b_i / a_{ii}
n és el nombre de punts, fixat 
*/
void assignValues(double x[], double a[9], double** u, double h, int n){
    //valors del vector x i la frontera de la "matriu" solució
    //posant valors a fora manualment ens estalviem 4 assignacions repetides
    x[0] = 0;
    x[n+1] = 1;
    u[0][0] = p(0, 0);
    u[0][n+1] = p(0, 1);
    u[n+1][0] = p(1, 0);
    u[n+1][n+1] = p(1, 1);
    for(int i=1; i<(n+1); i++){
        x[i] = i*h;
        u[i][0] = p(x[i], 0);
        u[i][n+1] = p(x[i], 1);
        u[0][i] = p(0, x[i]);
        u[n+1][i] = p(1, x[i]);
    }
    
    //assigna els valors de l'interior segons x_i = b_i / a_{ii}
    for(int i=1; i<(n+1); i++){
        for(int j=1; j<(n+1); j++){
            if(i == 1){
                if(j==1){
                    //cantonada esquerra inferior
                    u[i][j] = (g(x[i], x[j]) - a[0]*u[0][0] - a[3]*u[0][1] - a[6]*u[0][2] - a[1]*u[1][0] - a[2]*u[2][0])/a[4];
                }
                else if(j==n){
                    //cantonada esquerra superior
                    u[i][j] = (g(x[i], x[j]) - a[0]*u[0][n-1] - a[3]*u[0][n] - a[6]*u[0][n+1] - a[7]*u[1][n+1] - a[8]*u[n+1][n+1])/a[4];
                }
                else{
                    //costat esquerra
                    u[i][j] = (g(x[i], x[j]) - a[0]*u[0][j-1] - a[3]*u[0][j] - a[6]*u[0][j+1])/a[4];
                }
            }
            else if(i==n){
                if(j==1){
                    //cantonada dreta inferior
                    u[i][j] = (g(x[i], x[j]) - a[0]*u[n-1][0] - a[1]*u[n][0] - a[2]*u[n+1][0] - a[5]*u[n+1][1] - a[8]*u[n+1][2])/a[4];
                }
                else if(j==n){
                    //cantonada dreta superior
                    u[i][j] = (g(x[i], x[j]) -a[2]*u[n+1][n-1] - a[5]*u[n+1][n] - a[8]*u[n+1][n+1] - a[7]*u[n][n+1] - a[6]*u[n-1][n+1])/a[4];
                }
                else{
                    //costat dret
                    u[i][j] = (g(x[i], x[j]) -a[2]*u[n+1][j-1] - a[5]*u[n+1][j] - a[8]*u[n+1][j+1])/a[4];
                }
            }
            else{
                if(j==1){
                    //costat inferior
                    u[i][j] = (g(x[i], x[j]) - a[0]*u[i-1][0] - a[1]*u[i][0] - a[2]*u[i+1][0])/a[4];
                }
                else if(j==n){
                    //costat superior
                    u[i][j] = (g(x[i], x[j]) - a[8]*u[i+1][n+1] - a[7]*u[i][n+1] - a[6]*u[i-1][n+1])/a[4];
                }
                else{
                    //punts de l'interior
                    u[i][j] = g(x[i], x[j])/a[4]; 
                }
            }
        }
    }
    
}

/*
Funció que assigna valors als vectors x, u 
El vector solució queda iniciat a zeror a tot arreu menys la frontera
n és el nombre de punts, fixat 
*/
void assignValuesZero(double x[], double** u, double h, int n){
    for(int i=0; i<(n+2); i++){
        x[i] = i*h;
        //Asignem els valors de la frontera a mesura que calculem x_i, y_i amb la funcio p
        if(i==0){
            u[0][0] = p(x[0], x[0]);
        }
        else if(i==(n+1)){
            for(int j=0; j<n+2; j++){
                u[j][i] = p(x[j], x[i]);
                u[i][j] = p(x[i], x[j]);
            }
        }
        else{
            u[0][i] = p(x[0], x[i]);
            u[i][0] = p(x[i], x[0]);
        }
    }
}

/*
Funció que escriu la solució trobada en un fitxer
*/
void writeSolutionFile(double** u, double x[], int n){
    FILE *fptr;
    fptr = fopen("solution_sor.txt", "w");
    fprintf(fptr,"x, y, u\n");
    for(int i=1; i<n+1; i++){
        for(int j=1; j<n+1; j++){
            fprintf(fptr, "%lf, %lf, %+.6le\n", x[i], x[j], u[i][j]);
        }
    }
    fclose(fptr);
}

/*
Funció que retorna l'error real comparant-ho amb la funció coneguda
*/
double checkRealError(double** u, double x[], int n){
    double max_abs_val = 0;
    double current;
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            current = fabs(u[i][j]-p(x[i], x[j]));
            if(current>max_abs_val){
                max_abs_val = current;
            }
        }
    }
    return max_abs_val;
}

/*
Guarda a s la solució esperada per segons la p escollida
n és el nombre de punts de l'interior (el que varia en cada iteració del programa)
*/
void calculateSolution(double** s, double x[], int n){
    for(int i=0; i<(n+2); i++){
        for(int j=0; j<(n+2); j++){
            s[i][j] =  p(x[i], x[j]);
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