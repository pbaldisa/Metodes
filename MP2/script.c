#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define TOL 1e-5
#define PRECISION 1e-8
#define POINTS_IN_CURVE 30000
#define NEWTON_ITER 5
#define DELTA 5e-3

//define type for Point
typedef struct Point{
    double x;
    double y;
} Point;

void readPoints(Point points[4]);
Point baricentre(Point points[4]);
double dist(Point p1, Point p2);
Point graddist(Point p1, Point p2);
double f(Point X, Point points[4]);
Point gradf(Point X, Point points[4]);
double h(Point X, Point X0, Point points[4]);
Point gradh(Point X, Point points[4]);
Point perpendicular(Point X);
void computeInitialPoints(Point B, Point points[4], Point initialPoints[16]);
Point sumPoints(Point points[], int n);
Point sumTwoPoints(Point p1, Point p2);
Point mulPoint(Point p, double k);
Point divPoint(Point p, double k);
double norm(Point p);
Point Fun(Point X, Point X0, Point points[4], Point centre);

void predictorCorrector(Point initial, Point points[4], FILE **f);
Point predictor(Point initial, Point points[4]);
Point corrector(Point initial, Point centre, Point X0, Point points[4]);


int main(){
    //LLegeix els punts del pla que introdueix l'usuari
    printf("Introdueix els 4 punts del pla:\n");
    //Conjunt de 4 punts a partir del NIUB de l'usuari
    Point points[4];
    readPoints(points);

    //Calcula el baricentre dels punts
    Point B = baricentre(points);
    printf("El baricentre dels punts es: (%lf, %lf)\n", B.x, B.y);

    //Calcula els 16 punts inicials segons la definició de l'enunciat
    Point initialPoints[16];
    computeInitialPoints(B, points, initialPoints);
    
    // Per cada punt inicial, troba la corba associada
    for(int i=0; i<16; i++){
        //printf("Gradient pel punt inicial (%lf, %lf):\n", initialPoints[i].x, initialPoints[i].y);
        //Point aux = gradh(initialPoints[i], points);
        FILE *f = fopen("points.txt", "w");
        fprintf(f, "ind x y");
        //Punt inicial per trobar la corba
        Point initial = initialPoints[i];

        //Calcula predicció-correcció per aquest punt inicial (escriu els punts en un fitxer i va imprimint per pantalla)
        predictorCorrector(initial, points, &f);
        fclose(f);
    }
}

/*
    Funció que donat un punt inicial i els 4 punts del pla escriu un fitxer amb els punts de la corba que
    va trobant amb el mètode preidctor-corrector.
*/
void predictorCorrector(Point initial, Point points[4], FILE **f){
    Point X0 = initial;
    
    for(int index=0; index<POINTS_IN_CURVE; index++){
        // Calcula la predicció del punt següent
        Point X1 = predictor(X0, points);
        // Calcula la correcció del punt següent
        Point X2 = corrector(X1, X0, initial, points);

        // Comprova que no ha fet la volta
        if (dist(initial, X1) < DELTA){
            printf("Ha fet la volta\n");
            break;
        }
        // Escriu el punt trobat al fitxer       
        fprintf(f, "%d %lf %lf", index, X2.x, X2.y);
        // Actualitza el punt inicial per a la següent iteració
        X0 = X2;
    }
}

/*
    Funció predictora.
    Donat el punt x_i, calcula el gradient de la funció h (amb els punts corresponents) i el perpendicular.
    Amb això calcula la predicció del punt x_i+1, a una distància delta del punt x_i.
*/
Point predictor(Point initial, Point points[4]){
    Point grad = gradh(initial, points);
    Point perp = perpendicular(grad);
    if (abs(perp.x) < TOL && abs(perp.y) < TOL){
        printf("El gradient és zero\n");
        exit (EXIT_FAILURE);
    }
    // delta divided by the norm of the perpendicular vector
    double norm_perp = DELTA / norm(perp);
    Point X = sumTwoPoints(initial, mulPoint(perp, norm_perp));
    //imprimeix la predicció per pantalla
    printf("#pred: x=(%lf, %lf)\n", X.x, X.y);
    return X;
}

// Funció que volem fer zero al corrector amb el mètode de Newton
Point Fun(Point X, Point X0, Point points[4], Point centre){
    Point F = {h(X, X0, points), dist(X, centre) - pow(DELTA,2)};
    return F;
}

//Corrector function of the predictor-corrector algorithm
//It implements newton's method to finda point that makes h zero and is in the circumference of center initial point and radius delta
/*
    Funció correctora.
    Donat el punt predit per la funció predictora, utilitza el mètode de Newton per trobar un punt que faci que h sigui zero i estigui a una circumferència de radi delta i centre el punt centre.
*/
Point corrector(Point pred, Point centre, Point X0, Point points[4]){
    Point X_prev = pred;
    int i = 0;
    for (int i=0; i<NEWTON_ITER; i++){
        Point grad = gradh(X_prev, points);
        Point F = Fun(X_prev, X0, points, centre);
        // Determinant de DF (per la inversa de la matriu)
        double det = 2 * (grad.x * (X_prev.y - pred.y) - grad.y * (X_prev.x - pred.x));
        // Producte de la inversa de DF per F
        Point aux = {2 * (X_prev.y - pred.y) * F.x - grad.y * F.y, -2 * (X_prev.x - pred.x) * F.x - grad.x * F.y};
        aux = divPoint(aux, det);
        // X_i+1 = X_i - DF^-1 * F
        Point X = sumTwoPoints(X_prev, mulPoint(aux, -1));
        //imprimeix la correcció per pantalla
        printf("#corr: k=%d x=(%lf, %lf) dlf=%lf\n", i, X.x, X.y, dist(X, X_prev));
        // Si la distància entre X_i i X_i+1 és menor que la precisió desitjada, retornem X_i+1
        if (dist(X, X_prev) < PRECISION){
            return X;
        }
        X_prev = X;
    }
    //S'han acabat les iteracions i no ha trobat un punt
    exit (EXIT_FAILURE);
    
    /*while (h1 != 0 && i < NETWONITER){
        X = sumTwoPoints(X0, mulPoint(gradh(X0, points), -h0/dist(X0, initial)));
        h1 = h(X, initial, points);
        i++;
    }
    return X;*/
}

//Read 4 pairs of numbers from stdin delimited by a space in 4 different lines
//and assign an array of points
void readPoints(Point points[4]){
    for(int i=0; i<4; i++){
        scanf("%lf %lf", &points[i].x, &points[i].y);
    }
}

//compute distance between two points
double dist(Point p1, Point p2){
    return sqrt(pow(p2.x-p1.x, 2) + pow(p2.y-p1.y, 2));
}

double norm(Point p){
    return sqrt(pow(p.x, 2) + pow(p.y, 2));
}

//compute gradient of distance between two points
Point graddist(Point p1, Point p2){
    Point grad;
    grad.x = (p2.x-p1.x)/dist(p1, p2);
    grad.y = (p2.y-p1.y)/dist(p1, p2);
    return grad;
}

//Compute function f
double f(Point X, Point points[4]){
    return dist(X, points[0]) * dist(X, points[1]) * dist(X, points[2]) * dist(X, points[3]);
}

//compute gradient of f
Point gradf(Point X, Point points[4]){
    Point grad = {0,0};
    for(int i=0; i<4; i++){
        grad = sumTwoPoints(grad, mulPoint(graddist(X, points[i]), f(X,points)/dist(X, points[i])));
    }
    //maybe this has better error??
    /*printf("gradf: (%lf, %lf)\n", grad.x, grad.y);

    grad.x = dist(X, points[0]) * dist(X, points[1]) * dist(X, points[2]) * graddist(X, points[3]).x
            + dist(X, points[0]) * dist(X, points[1]) * dist(X, points[3]) * graddist(X, points[2]).x
            + dist(X, points[0]) * dist(X, points[2]) * dist(X, points[3]) * graddist(X, points[1]).x
            + dist(X, points[1]) * dist(X, points[2]) * dist(X, points[3]) * graddist(X, points[0]).x;
    grad.y = dist(X, points[0]) * dist(X, points[1]) * dist(X, points[2]) * graddist(X, points[3]).y
            + dist(X, points[0]) * dist(X, points[1]) * dist(X, points[3]) * graddist(X, points[2]).y
            + dist(X, points[0]) * dist(X, points[2]) * dist(X, points[3]) * graddist(X, points[1]).y
            + dist(X, points[1]) * dist(X, points[2]) * dist(X, points[3]) * graddist(X, points[0]).y;

    printf("gradf: (%lf, %lf)\n", grad.x, grad.y);*/
    return grad;
}

//Compute function h
double h(Point X, Point X0, Point points[4]){
    return f(X, points) - f(X0, points);
}

//Compute gradient of h
Point gradh(Point X, Point points[4]){
    return gradf(X, points);
}

//Compute perpendicular to a vector
Point perpendicular(Point X){
    Point perp;
    perp.x = -X.y;
    perp.y = X.x;
    return perp;
}

//Sum points
Point sumPoints(Point points[], int n){
    Point sum;
    sum.x = 0;
    sum.y = 0;
    for(int i=0; i<n; i++){
        sum.x += points[i].x;
        sum.y += points[i].y;
    }
    return sum;
}

//sum two points
Point sumTwoPoints(Point p1, Point p2){
    Point sum;
    sum.x = p1.x + p2.x;
    sum.y = p1.y + p2.y;
    return sum;
}

//Multiply point by a scalar
Point mulPoint(Point p, double k){
    Point mul;
    mul.x = p.x * k;
    mul.y = p.y * k;
    return mul;
}

//Divide point by a scalar
Point divPoint(Point p, double k){
    Point div;
    div.x = p.x / k;
    div.y = p.y / k;
    return div;
}

//Compute baricentre of 4 points
Point baricentre(Point points[4]){
    Point B;
    B = sumPoints(points, 4);
    B = divPoint(B, 4);
    return B;
}

//Compute initial points
void computeInitialPoints(Point B, Point points[4], Point initialPoints[16]){
    for(int i=0; i<4; i++){
        initialPoints[i] = sumTwoPoints(mulPoint(points[i], 1.2), mulPoint(B, -0.2));
        initialPoints[i+4] = sumTwoPoints(mulPoint(points[i], 1.1), mulPoint(B, -0.1));
        initialPoints[i+8] = sumTwoPoints(mulPoint(points[i], 0.95), mulPoint(B, 0.05));
        initialPoints[i+12] = sumTwoPoints(mulPoint(points[i], 0.85), mulPoint(B, 0.15));
    }
}