
 
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <string>
#define    MAX    160801

using namespace std;

double V_new[MAX],V_old[MAX];

int N,N1;                    // V array N1*N1 dim
double delta;
double deltaVmin;
double dV;                   // deltaV jacobi

double V_in,V_out;           // potenziale faccie prisma

int P1,P2;                   // quadrato interno

int i,j;
int iarray;                  // iarray = i + N1*J

unsigned int n_iterations=0;
unsigned int n_iterations_rid=0;
void iterate();
void init();
void initprismacomp();
int update(double V1[MAX], double V0[MAX]);
void print_V(double V[MAX]);
void print_Efield(double V[MAX]);
void simmetry(double V[]);
void iterate_rid();
void init_rid();
void init_ridprismacomp();
int update_rid(double V1[MAX], double V0[MAX]);

int main()
{
    
    cout <<"Sistema ridotto:\n";
    init_rid();
    iterate_rid();
    simmetry(V_old);
    
    print_V(V_old);
    print_Efield(V_old);
    
    cout <<"Sistema completo:\n";
    
    init();
    iterate();
    print_V(V_old);
    print_Efield(V_old);
    
    double a = update_rid(V_new,V_old);
    double b = update(V_new,V_old);
    
    double rapporto = a/b;
    
    //cout << "iterazionirid " << update_rid(V_new,V_old) << endl;
    //cout << "interazionicomp" << update(V_new,V_old) << endl;
    
    cout << "Il rapporto fra il numero di iterazioni del sistema ridotto e quello completo vale: " << rapporto << endl;
    
    return 0;
}


void iterate()
{
    int count=0;
    while(dV>deltaVmin){
        
        update(V_new,V_old);
        cout << dV << "\t" << deltaVmin << "\t" << count << "\t" << n_iterations << endl;
        
        count++;
        
        update(V_old,V_new);
        cout << dV << "\t" << deltaVmin << "\t" << count << "\t" << n_iterations <<endl;
        count++;
    }
    
    
    return;
}


void init()
{
    
    cout << "Inserisci N della griglia:" << endl;
    cin >> N;
    N1 = N + 1;
    delta = 2.0/N;
    P1 = N/4;    // -0.5
    P2 = 3*P1;   // +0.5
    
    cout << "Inserisci il valore di soglia deltaVmin:" << endl;
    cin>> deltaVmin;
    dV = deltaVmin + 1.0; // inizializzo dV e faccio partire l'iterazione
    
    if (N1*N1>MAX) cout << "Array troppo grande!" << endl;
    
    
    initprismacomp();
    
    
    return;
}

void initprismacomp(){   /* V initialization for prism capacitor */
    
    V_in = 1.0;    // condizioni al bordo
    V_out = -1.0;
    
    // inizializzo V_new e V_old e uno aggiorna l'altro
    
    for(i=0;i<N1;i++){
        for(j=0;j<N1;j++){
            iarray = i + N1*j;
            if (i==0||i==N||j==0||j==N){
                V_old[iarray] = V_out;                  //valore di bordo sulla faccia esterna per i due array
                V_new[iarray] = V_out;
            }
            else if (i>=P1&&i<=P2&&j>=P1&&j<=P2){
                V_old[iarray] = V_in;                   //valore di bordo sulla faccia interna per i due array
                V_new[iarray] = V_in;
                //  V_in e' assegnato a tutti i punti dentro il quadrato interno
                
            }
            else V_old[iarray] = 0.0;                 // condizione iniziale per i punti di bulk
        }
    }
    
    return;
    
}

int update(double V1[MAX], double V0[MAX]){
    
    int l_center,l_left,l_right,l_up,l_down,n_grid;
    
    dV = 0.0;
    n_grid = 0;
    
    for(i=0;i<N1;i++){ // loop su tutte le colonne
        if(i>0&&i<N){    // filtro prima e ultima colonna, sempre di bordo
            
            
            for(j=0;j<N1;j++){  // loop su tutte le righe
                if(j>0&&j<N){     // filtro prima e ultima riga, di bordo
                    
                    
                    if (i<P1||i>P2||j<P1||j>P2)  // condizione del bulk
                        
                        
                    {
                        
                        //Jacobi
                        
                        l_center = i + N1*j;
                        l_left = l_center - 1;
                        l_right = l_center + 1;
                        l_up = l_center - N1;
                        l_down = l_center + N1;
                        V1[l_center] = (V0[l_left]+V0[l_right]+V0[l_up]+V0[l_down])/4.0;          // Laplaciano
                        dV += fabs(V1[l_center]-V0[l_center]);          //calcolo dV per la convergenza poi divido per n_grid
                        n_iterations++;
                        n_grid ++;
                        
                    }
                }
            }
        }
    }
    
    dV /= n_grid;
    
    return n_iterations;
}

void print_V(double V[MAX]){
    
    double x,y;
    
    string s;
    
    cout << "Scegli il nome del file per V:\n";
    cin >>s;
    
    ofstream out (s,std::ios_base::out);
    
    for(j=0;j<N1;j++){
        y = 1.0 - j*delta;
        for(i=0;i<N1;i++){
            x = -1.0 + i*delta;
            iarray = i + j*N1;
            out << x <<"\t" <<y << "\t" << V[iarray] << endl;
        }
    }
    
    return;
}

void print_Efield(double V[MAX]){
    
    double x,y,E_x,E_y,dvec;
    int l_center,l_left,l_right,l_up,l_down,n_grid;
    dvec = 0.5;
    
    string g;
    
    cout << "Scegli il nome del file per E:\n";
    cin >>g;
    
    ofstream out (g,std::ios_base::out);
    
    for(i=0;i<N1;i++){
        x = -1.0 + i*delta;
        for(j=0;j<N1;j++){
            y = 1.0 - j*delta;
            l_center = i + N1*j;
            l_left = l_center - 1;
            l_right = l_center + 1;
            l_up = l_center - N1;
            l_down = l_center + N1;
            if (j==0) E_y = -(V[l_center] - V[l_down])/delta;      /* downward y-derivative */
            else if (j==N) E_y = -(V[l_up] - V[l_center])/delta;   /* upward y-derivative   */
            else E_y = -0.5*(V[l_up] - V[l_down])/delta;           /* bulk y-derivative     */
            if (i==0) E_x = -(V[l_right] - V[l_center])/delta;     /* forward x-derivative  */
            else if (i==N) E_x = -(V[l_center] - V[l_left])/delta; /* backward x-derivative */
            else E_x = -0.5*(V[l_right] - V[l_left])/delta;        /* bulk x-derivative     */
            out << x-E_x*delta*dvec/2 <<"\t"<< y-E_y*delta*dvec/2 <<"\t"<< E_x*delta*dvec <<"\t"<< E_y*delta*dvec << endl;
        }
    }
    
    return;
}



void iterate_rid()
{
    int count=0;
    while(dV>deltaVmin){
        
        update_rid(V_new,V_old);
        cout << dV << "\t" << deltaVmin << "\t" << count << "\t" << n_iterations_rid << endl;
        
        count++;
        
        update_rid(V_old,V_new);
        cout << dV << "\t" << deltaVmin << "\t" << count << "\t" << n_iterations_rid <<endl;
        count++;
    }
    
    
    return;
}


void init_rid()
{
    
    cout << "Inserisci N della griglia:" << endl;
    cin >> N;
    N1 = N + 1;
    delta = 2.0/N;
    P1 = N/4;    // -0.5
    P2 = 3*P1;   // +0.5
    
    cout << "Inserisci il valore di soglia deltaVmin:" << endl;
    cin>> deltaVmin;
    dV = deltaVmin + 1.0; // inizializzo dV e faccio partire l'iterazione
    
    if (N1*N1>MAX) cout << "Array troppo grande!" << endl;
    
    
    init_ridprismacomp();
    
    
    return;
}

void init_ridprismacomp(){   /* V init_ridialization for prism capacitor */
    
    V_in = 1.0;    // condizioni al bordo
    V_out = -1.0;
    
    // inizializzo V_new e V_old e uno aggiorna l'altro
    
    for(i=0;i<N1;i++){
        for(j=0;j<N1;j++){
            iarray = i + N1*j;
            if (i==0||i==N||j==0||j==N){
                V_old[iarray] = V_out;                  //valore di bordo sulla faccia esterna per i due array
                V_new[iarray] = V_out;
            }
            else if (i>=P1&&i<=P2&&j>=P1&&j<=P2){
                V_old[iarray] = V_in;                   //valore di bordo sulla faccia interna per i due array
                V_new[iarray] = V_in;
                //  V_in e' assegnato a tutti i punti dentro il quadrato interno
                
            }
            else V_old[iarray] = 0.0;                 // condizione iniziale per i punti di bulk
        }
    }
    
    return;
    
}

int update_rid(double V1[MAX], double V0[MAX]){
    
    int l_center,l_left,l_right,l_up,l_down,n_grid, k1,k2,k3,k4,k5,k6,k7,k8;
    
    dV = 0.0;
    n_grid = 0;
    
    for(i=0;i<N1;i++){ // loop su tutte le colonne
        if(i>0&&i<N){    // filtro prima e ultima colonna, sempre di bordo
            
            
            for(j=0;j<N1;j++){  // loop su tutte le righe
                if(j>0&&j<N){     // filtro prima e ultima riga, di bordo
                    
                    
                    if (i>P2&&j<=i&&j>=(N/2))  // condizione del bulk
                        
                        
                    {
                        
                        //Jacobi
                        
                        l_center = i + N1*j;
                        l_left = l_center - 1;
                        l_right = l_center + 1;
                        l_up = l_center + N1;
                        l_down = l_center - N1;
                        if (j==N/2) V1[l_center] = (V0[l_left]+V0[l_right]+2*V0[l_up])/4.0;
                        else
                        {
                            if(j==i) V1[l_center] = (2*V0[l_right]+2*V0[l_down])/4.0;
                            else V1[l_center] = (V0[l_left]+V0[l_right]+V0[l_up]+V0[l_down])/4.0;
                        }
                        
                        
                        // Laplaciano
                        dV += fabs(V1[l_center]-V0[l_center]);          //calcolo dV per la convergenza poi divido per n_grid
                        n_iterations_rid++;
                        n_grid ++;
                        
                        
                    }
                }
            }
        }
    }
    
    dV /= n_grid;
    
    return n_iterations_rid;
}

void simmetry (double V1[MAX])
{
    int l_center, k1,k2,k3,k4,k5,k6,k7,k8;
    for(i=0;i<N1;i++){ // loop su tutte le colonne
        if(i>0&&i<N){    // filtro prima e ultima colonna, sempre di bordo
            
            
            for(j=0;j<N1;j++){  // loop su tutte le righe
                if(j>0&&j<N){     // filtro prima e ultima riga, di bordo
                    
                    
                    if (i>P2&&j<=i&&j>=(N/2))  // condizione del bulk
                        
                        
                    {
                        
                        //Jacobi
                        
                        l_center = i + N1*j;
                        
                        k1=N1*i+j;
                        V1[k1]=V1[l_center];
                        
                        k2=N-i + N1*j;
                        V1[k2]=V1[l_center];
                        
                        k3=i + N1*(N-j);
                        V1[k3]=V1[l_center];
                        
                        k4=N-i+N1*(N-j);
                        V1[k4]=V1[l_center];
                        
                        k5=N1*(N-i)+j;
                        V1[k5]=V1[l_center];
                        
                        k6=N1*i +(N -j);
                        V1[k6]=V1[l_center];
                        
                        k7=N1*(N-i)+(N-j);
                        V1[k7]=V1[l_center];
                        
                        
                        
                    }
                }
            }
        }
    }
    
    
}


