#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

// Definició de constants com a l'enunciat
#define TOL 1e-5                // Tolerància pel zero
#define PRECISION 1e-8          // Precisió desitjada en la correccio
#define POINTS_IN_CURVE 30000   // Nombre màxim de punts que volem trobar per cada corba
#define NEWTON_ITER 5           // Nombre màxim d'iteracions de Newton
#define DELTA 5e-3              // Distància entre els punts de la corba

// Els punts com a tipus de dades
typedef struct Point{
    double x;
    double y;
} Point;

// Funcions que utilitzarem
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
Point substractTwoPoints(Point p1, Point p2);
Point mulPoint(Point p, double k);
Point divPoint(Point p, double k);
double norm(Point p);
Point Fun(Point X, Point X0, Point points[4], Point centre);

// Funcions principals
void predictorCorrector(Point initial, Point points[4], FILE *f);
Point predictor(Point initial, Point points[4]);
Point corrector(Point initial, Point centre, Point X0, Point points[4]);


int main(){
    //LLegeix els punts del pla que introdueix l'usuari
    printf("Introdueix els 4 punts del pla:\n");
    //Conjunt de 4 punts a partir del NIUB de l'usuari

    Point points[4];
    readPoints(points);
    //Point points[4] = {(Point){2.0, 0.0}, (Point){1.0, 5.0}, (Point){0.0, 2.0}, (Point){7.0, 0.0}};

    //Calcula el baricentre dels punts
    Point B = baricentre(points);
    printf("El baricentre dels punts es: (%lf, %lf)\n", B.x, B.y);

    //Calcula els 16 punts inicials segons la definició de l'enunciat
    Point initialPoints[16];
    computeInitialPoints(B, points, initialPoints);
    
    // Obre el fitxer on escriure els punts de la corba
    FILE *f = fopen("points.res", "w");
    fprintf(f, "#ind x y\n");

    // Per cada punt inicial, troba la corba associada
    for(int i=0; i<16; i++){
        // Escriu un espai abans de començar una corba al fitxer
        fprintf(f, "\n");
        
        //Punt inicial per trobar la corba
        Point initial = initialPoints[i];

        // Escriu per pantalla el nombre de corba que està calculant i el seu punt inicial
        printf("#Corba %d\n", i+1);
        printf("#Punt inicial: (%lf, %lf)\n", initial.x, initial.y);
        //Calcula predicció-correcció per aquest punt inicial (escriu els punts en un fitxer i va imprimint per pantalla)
        predictorCorrector(initial, points, f);
        
    }
    // Tanca el fitxer i acaba el programa
    fclose(f);
    return 0;
}

/*
------------------------------------------------------------------------------------------------------------------
    Funcions principals
------------------------------------------------------------------------------------------------------------------
*/

/*
    Funció que donat un punt inicial i els 4 punts del pla escriu un fitxer amb els punts de la corba que
    va trobant amb el mètode preidctor-corrector.
*/
void predictorCorrector(Point initial, Point points[4], FILE *f){
    Point X0 = initial; //fa 0 la corba h
    Point X1, X2;
    // Escriu el punt inicial al fitxer
    fprintf(f, "%d %lf %lf\n", 0, initial.x, initial.y);
    for(int index=0; index<POINTS_IN_CURVE; index++){
        // Calcula la predicció del punt següent
        X1 = predictor(X0, points);
        // Calcula la correcció del punt següent
        X2 = corrector(X1, X0, initial, points);
        
        // Escriu el punt trobat al fitxer       
        fprintf(f, "%d %lf %lf\n", index + 1, X2.x, X2.y);

        // Comprova que no ha fet la volta
        if ((dist(initial, X2) < DELTA) && (index != 0)){
            printf("Ha fet la volta!\n\n");
            // Per mostrar el gràfic tancat, ens assegurem que el punt inicial i el final són iguals
            fprintf(f, "%d %lf %lf\n", index, initial.x, initial.y);
            break;
        }

        // Actualitza el punt inicial per a la següent iteració
        X0 = X2;
    }
    printf("\n");
}

/*
    Funció predictora.
    Donat el punt x_i, calcula el gradient de la funció h (amb els punts corresponents) i el perpendicular.
    Amb això calcula la predicció del punt x_i+1, a una distància delta del punt x_i.
*/
Point predictor(Point initial, Point points[4]){
    // Calcula el gradient de la funció h
    Point grad = gradh(initial, points);
    Point perp = perpendicular(grad);
    // Comprova que el gradient no sigui zero
    if (fabs(perp.x) < TOL && fabs(perp.y) < TOL){
        printf("El gradient és zero\n");
        // Si el gradient és 0 sortim del programa
        exit (EXIT_FAILURE);
    }
    // Per tal de tenir un vector de norma delta, el normalitzem
    double norm_perp = DELTA / norm(perp);
    // Calcula el punt predit: z_k+1 = z_k + delta * normalized(perp)
    Point X = sumTwoPoints(initial, mulPoint(perp, norm_perp));
    //  Imprimeix la predicció per pantalla
    printf("#pred: x=(%.10lf, %.10lf)\n", X.x, X.y);

    // Retorna el punt predit
    return X;
}

/*
    Funció correctora.
    Donat el punt predit per la funció predictora, utilitza el mètode de Newton per trobar un punt que faci que h sigui zero i estigui a una circumferència de radi delta i centre el punt centre.
*/
Point corrector(Point pred, Point centre, Point X0, Point points[4]){
    // Punt anterior
    Point X_prev = pred;
    // Punts gradient, avaluat a la funció i actual
    Point grad, F, X;
    // Diferència entre punts i determinant de la diferencial
    double dif, det;
    for (int i=0; i<NEWTON_ITER; i++){
        // Gradient de la funció h (primera fila de la matriu DF)
        grad = gradh(X_prev, points);
        // Avaluació de la funció Fun
        F = Fun(X_prev, X0, points, centre);

        // Determinant de DF (per la inversa de la matriu)
        det = 2.0 * (grad.x * (X_prev.y - centre.y) - grad.y * (X_prev.x - centre.x));
        // Producte de la inversa de DF per F
        Point aux = {2.0 * (X_prev.y - centre.y) * F.x - grad.y * F.y, -2.0 * (X_prev.x - centre.x) * F.x + grad.x * F.y};
        aux = divPoint(aux, det);
        // X_i+1 = X_i - DF^-1 * F
        X = substractTwoPoints(X_prev, aux);

        // Calcula la distància entre X_i i X_i+1
        dif = dist(X, X_prev);
        // Imprimeix la correcció per pantalla
        printf("#corr: k=%d x=(%.10lf, %.10lf) dif=%e\n", i, X.x, X.y, dif);

        // Si la distància entre X_i i X_i+1 és menor que la precisió desitjada, retornem X_i+1
        if (dif < PRECISION){
            return X;
        }

        // Actualitza el punt anterior
        X_prev = X;
    }
    // S'han acabat les iteracions i no ha trobat un punt
    // Es deixa d'executar el programa (no podem confiar que continui correctament)
    exit (EXIT_FAILURE);
}

/*
------------------------------------------------------------------------------------------------------------------
    Funcions auxiliars
------------------------------------------------------------------------------------------------------------------
*/

/*
    Funció que volem fer zero al corrector amb el mètode de Newton
    Trobar zeros de Fun correspon a trobar solucions del sistema 
    h(x) = 0
    (x-x0)^2 + (y-y0)^2 - delta^2 = 0
    Corresponent a la correcció: punt de la corba a distància delta del punt anterior
*/
Point Fun(Point X, Point X0, Point points[4], Point centre){
    Point F = {h(X, X0, points), pow(X.x-centre.x,2) + pow(X.y-centre.y,2) - pow(DELTA,2)};
    return F;
}

/*
    Funció que llegeix 4 punts per pantalla (4 línies amb dos nombres separats per un espai)
*/
void readPoints(Point points[4]){
    for(int i=0; i<4; i++){
        scanf("%lf %lf", &points[i].x, &points[i].y);
    }
}

/*
    Funció que calcula la distància euclidiana entre dos punts
*/
double dist(Point p1, Point p2){
    return sqrt(pow(p1.x-p2.x, 2) + pow(p1.y-p2.y, 2));
}

/*
    Funció que calcula la norma d'un punt
*/
double norm(Point p){
    return sqrt(pow(p.x, 2) + pow(p.y, 2));
}

/*
    Funció que calcula el gradient de la funció distància
*/
Point graddist(Point p1, Point p2){
    Point grad;
    grad.x = (p1.x-p2.x)/dist(p1, p2);
    grad.y = (p1.y-p2.y)/dist(p1, p2);
    return grad;
}

/*
    Funció que calcula la funció f de l'enunciat (multilicació de distàncies)
*/
double f(Point X, Point points[4]){
    return dist(X, points[0]) * dist(X, points[1]) * dist(X, points[2]) * dist(X, points[3]);
}

/*
    Funció que calcula el gradient de f en un punt X
*/
Point gradf(Point X, Point points[4]){
    Point grad = {0,0};
    for(int i=0; i<4; i++){
        // La derivada és la suma de 4 termes
        grad = sumTwoPoints(grad, mulPoint(graddist(X, points[i]), f(X,points)/dist(X, points[i])));
        // El terme quan i=1 és:
        // grad_x = dist(x, p0) * dist(x, p2) * dist(x, p3) * dist_x(x,p1)
        // grad_y = dist(x, p0) * dist(x, p2) * dist(x, p3) * dist_y(x,p1)
    }
    return grad;
}

/*
    Funció que calcula la funció h de la pràctica. f(x) - f(x0)
*/
double h(Point X, Point X0, Point points[4]){
    return f(X, points) - f(X0, points);
}

/*
    Funció que calcula el gradient de la funció h en un punt X
*/
Point gradh(Point X, Point points[4]){
    // Coincideix amb el de f perquè f(x0) és una constant
    return gradf(X, points);
}

/*
    Donat un vector (punt), retorna un vector perpendicular
*/
Point perpendicular(Point X){
    Point perp;
    perp.x = -X.y;
    perp.y = X.x;
    return perp;
}

/*
    Funció que suma els punts d'un vector de n punts
*/
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

/*
    Funció que suma 2 punts
*/
Point sumTwoPoints(Point p1, Point p2){
    Point sum;
    sum.x = p1.x + p2.x;
    sum.y = p1.y + p2.y;
    return sum;
}

/*
    Funció que resta un punt de l'altre (p1 - p2)
*/
Point substractTwoPoints(Point p1, Point p2){
    Point sub;
    sub.x = p1.x - p2.x;
    sub.y = p1.y - p2.y;
    return sub;
}

/*
    Funció que multiplica un punt per un escalar
*/
Point mulPoint(Point p, double k){
    Point mul;
    mul.x = p.x * k;
    mul.y = p.y * k;
    return mul;
}

/*
    Funció que divideix un punt per un escalar
*/
Point divPoint(Point p, double k){
    Point div;
    div.x = p.x / k;
    div.y = p.y / k;
    return div;
}

/*
    Funció que calcula el baricentre de 4 punts passats en un vector
*/
Point baricentre(Point points[4]){
    Point B;
    B = sumPoints(points, 4);
    B = divPoint(B, 4);
    return B;
}

/*
    Funció que calcula els 16 punts inicials segons la definició de l'enunciat
*/
void computeInitialPoints(Point B, Point points[4], Point initialPoints[16]){
    for(int i=0; i<4; i++){
        initialPoints[i] = sumTwoPoints(mulPoint(points[i], 1.2), mulPoint(B, -0.2));
        initialPoints[i+4] = sumTwoPoints(mulPoint(points[i], 1.1), mulPoint(B, -0.1));
        initialPoints[i+8] = sumTwoPoints(mulPoint(points[i], 0.95), mulPoint(B, 0.05));
        initialPoints[i+12] = sumTwoPoints(mulPoint(points[i], 0.85), mulPoint(B, 0.15));
    }
}