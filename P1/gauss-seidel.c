#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define A 1.0
#define B 1.0
#define C 0.0
#define D 0.0
#define E 0.0
#define F 0.0

double g(double x, double y);
double p(double x, double y);
void assignValues(double x[], double y[], double** u, double h, int n);
void calculateSolution(double** s, double x[], double y[], int n);
void printMatrix(double** u, int n);
double matrixInfNorm(double** m1, double** m2, int n);
int gs(double** u, double x[], double y[], double a[9], int n, double error, int maxIter);
double calculateErrorBound(double delta, double old_delta);
void writeSolutionFile(double** u, double x[], double y[], int n);

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

    //nombre d'iteracions que ha fet el mètode iteratiu 
    int totalIter;

    /*
    Crea els vectors i matrius a usar
    */
    //vectors que representen els punts x, y on avaluem la solució
    double x[n+2];
    double y[n+2];
    //vector solució en forma de matriu co a iteracio per gs; i s la solució exacta
    double** s = (double**)calloc((n+2),sizeof(double*));
    double** u = (double**)calloc((n+2),sizeof(double*));   //representa la solució trobada en la iteració actual
    for(int i=0; i<(n+2); i++){
        u[i] = (double*)calloc((n+2),sizeof(double));
        s[i] = (double*)calloc((n+2),sizeof(double));
    } 

    double h = 1/(double)(n+1); //pas de discretització

    //assigna valors als vectors
    assignValues(x, y, u, h, n);

    //Per no haver de refer els càlculs a cada iteració de gs, guarda els coeficients:
    double a[9] = {B/(4*h*h), C/(h*h)-E/(2*h), -B/(4*h*h), A/(h*h)-D/(2*h), (-2*A-2*C)/(h*h)+F, A/(h*h)+D/(2*h), -B/(4*h*h), C/(h*h)+E/(2*h), B/(4*h*h)};

    //resol per gs
    totalIter = gs(u, x, y, a, n, error, maxIter); 

    //Escriu la solució en un fitxer
    writeSolutionFile(u, x, y, n);

    printf("El nombre d'iteracions necessari ha estat: %d\n", totalIter);

    //Calcula la solució esperada
    calculateSolution(s, x, y, n);

    printf("L'error real ha estat de %E", matrixInfNorm(s, u, n+2));

    //allibera els vectors
    for(int i=0; i<(n+2); i++){
        free(u[i]);
        free(s[i]);
    } 
    free(u);
    free(s);
    
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
Funció que implementa el mètode de gs. 
Modifica u fins trobar la solució aproximada i retorna el nombre d'iteracions que ha fet
La n aquí representa el nombre fixat per aquesta iteració de gs
*/
int gs(double** u, double x[], double y[], double a[9], int n, double error, int maxIter){
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
                u[i][j] = (g(x[i], y[j]) - a[0]*u[i-1][j-1] - a[1]*u[i][j-1] - a[2]*u[i+1][j-1] - a[3]*u[i-1][j] - a[5]*u[i+1][j] - a[6]*u[i-1][j+1] - a[7]*u[i][j+1] - a[8]*u[i+1][j+1])/a[4];
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
    }while((iter < maxIter) && (error_aprox>error));

    printf("\nEl quocient d'error ha estat %lf\n", error_q);

    return iter;
}

/*
Funció que assigna valors als vectors x, y, u 
n és el nombre de punts, fixat 
*/
void assignValues(double x[], double y[], double** u, double h, int n){
    for(int i=0; i<(n+2); i++){
        x[i] = i*h;
        y[i] = i*h;
        //Asignem els valors de la frontera a mesura que calculem x_i, y_i amb la funcio p
        if(i==0){
            u[0][0] = p(x[0], y[0]);
        }
        else if(i==(n+1)){
            for(int j=0; j<n+2; j++){
                u[j][i] = p(x[j], y[i]);
                u[i][j] = p(x[i], y[j]);
            }
        }
        else{
            u[0][i] = p(x[0], y[i]);
            u[i][0] = p(x[i], y[0]);
        }
    }
}

/*
Funció que escriu la solució trobada en un fitxer
*/
void writeSolutionFile(double** u, double x[], double y[], int n){
    FILE *fptr;
    fptr = fopen("solution_gs.txt", "w");
    fprintf(fptr,"x, y, u\n");
    for(int i=1; i<n+1; i++){
        for(int j=1; j<n+1; j++){
            fprintf(fptr, "%lf, %lf, %lf\n", x[i], y[j], u[i][j]);
        }
    }
    fclose(fptr);
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
