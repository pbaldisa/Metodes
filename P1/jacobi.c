#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A 1.0
#define B 1.0
#define C 1.0
#define D 1.0
#define E 1.0
#define F 1.0

#define N 4

double g(double x, double y);
double p(double x, double y);
void assignValues(double x[], double y[], double** u, double h);
void printSolution(double** u);

int main(){
    const double error;
    //vector solució en forma de matriu
    double** u = (double**)calloc((N+2),sizeof(double*));
    for(int i=0; i<(N+2); i++){
        u[i] = (double*)calloc((N+2),sizeof(double));
    } 

    //vectors que representen els punts x, y on avaluem la solució
    double x[N+2];
    double y[N+2];

    double h = 1/(double)(N+1); //pas de discretització

    assignValues(x, y, u, h);

    //printSolution(u);

    int iter = 0;

    double** uk = (double**)calloc((N+2),sizeof(double*));
    for(int i=0; i<(N+2); i++){
        uk[i] = (double*)calloc((N+2),sizeof(double));
    } 

    //Canvia els valors de la frontera als valors que té u
    for(int i=0; i<(N+2); i++){
            uk[i][0] = u[i][0];
            uk[i][N+1] = u[i][N+1];
            uk[0][i] = u[0][i];
            uk[N+1][i] = u[N+1][i];
        }

    //Per no haver de refer els càlculs cada cop, guardem els coeficients:
    double a[9] = {B/(4*h*h), C/(h*h)-E/(2*h), -B/(4*h*h), A/(h*h)-D/(2*h), (-2*A-2*C)/(h*h)+F, A/(h*h)+D/(2*h), -B/(4*h*h), C/(h*h)+E/(2*h), B/(4*h*h)};

    //TODO posar condicions d'error
    while(iter < 2000){
        for(int i=1; i<(N+1); i++){
            for(int j=1; j<(N+1); j++){
                uk[i][j] = 1/a[4]*(g(x[i], y[j]) - a[0]*u[i-1][j-1] - a[1]*u[i][j-1] - a[2]*u[i+1][j-1] - a[3]*u[i-1][j] - a[5]*u[i+1][j] - a[6]*u[i-1][j+1] - a[7]*u[i][j+1] - a[8]*u[i+1][j+1]);
            }
        }
        for(int i=1; i<(N+1); i++){
            for(int j=1; j<(N+1); j++){
                u[i][j] = uk[i][j];
            }
        }
        iter++;
    }  
    
    printf("Solució esperada:\n");
    for(int i=0; i<(N+2); i++){
        for(int j=0; j<(N+2); j++){
            printf("%lf ", g(x[i], y[j]));
        }
        printf("\n");
    }
    
    printSolution(u);

    free(uk);
    free(u);
}

/*
---------------------------------------------------------------------------
FUNCIONS AUXILIARS
---------------------------------------------------------------------------
*/

double p(double x, double y){
    return x + y;
}

double g(double x, double y){
    return x + y;
}

void assignValues(double x[], double y[], double** u, double h){
    for(int i=0; i<(N+2); i++){
        x[i] = i*h;
        y[i] = i*h;
        //Asignem els valors de la frontera a mesura que calculem x_i, y_i
        if(i==0){
            u[0][0] = g(x[0], y[0]);
        }
        else if(i==(N+1)){
            for(int j=0; j<N+2; j++){
                u[j][i] = g(x[j], y[i]);
                u[i][j] = g(x[i], y[j]);
            }
        }
        else{
            u[0][i] = g(x[0], y[i]);
            u[i][0] = g(x[i], y[0]);
        }
    }
}

void printSolution(double** u){
    printf("u:\n");
    for(int i=0; i<N+2; i++){
        for(int j=0; j<N+2; j++){
        printf("%lf ", u[i][j]);
        }
        printf("\n");
    }
}
