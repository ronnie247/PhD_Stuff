//Kinetic Monte Carlo Project Part 2 Version 1
//C++ Version

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <stdlib.h>
// #include <time>
#include <chrono>
// #include "masses.h"

using namespace std;

//make a class of a given state (Lattice)
//variables : no of atoms, charges, geometry, energy
//functions : constructor/destructor function
class State 
{
    int NAtoms;
    int Zs[NAtoms]; // get Atomic numbers
    double masses[NAtoms]; //get masses from masses.h from crawford
    double charges[NAtoms]; //stores charges for each atom
    double geom[NAtoms][3]; //stores the xyz coordinates
//     double rates[3][3]; //rates in different direction
//     double wait; //rate of waiting = rate[1][1]
//     double elecfield[3]; //net electric field at point (x,y,z)
//     double force[3]; //net force acting on a particle at (x,y,z)
//     double energy; //net energy at (x,y,z)
    
    void State(int N, G[][], Q[], A[]){
        cout << "Initializing State";
        NAtoms = N;
        for(int i=0; i<N; i++){
            charges[i]=Q[i];
            geom[i][0]=G[i][0];
            geom[i][1]=G[i][1];
            geom[i][2]=G[i][2];
            Zs[i]=A[i];
//             masses[i]=masses(Zs[i]);
        }
        //do not calculate anything in the constructor function
        cout << "Initialized State S";
    }  
    
    public void State(){
        cout << "Initialized NULL State ";
    } 
    
    ~State();  
};

//make a class of a given atom that moves 
//variables : position, charge
//functions : constructor/destructor function
class Particle
{
    int Z; //get Atomic Number
    double M; //get mass from masses.h
    double Q; //get charge
    double pos[3]; //get cartesian coords
    
    void Particle(P[], Q1, A){
        cout << "Initializing Particle";
        Z = A;
//         M = masses(Z);
        pos[0] = P[0];
        pos[1] = P[1];
        pos[2] = P[2];
        Q = Q1;
        cout << "Initialized Particle";
    }
    
    public void Particle(){
        cout << "Initialized NULL Particle ";
    }
    
    double * distance(double p, double G[]);
    //calculates particle distance from each charge in the State
    double * elecfield(double r[][3] double C[], double k = 1.0); 
    //calculate electric field component on particle
    double elecpot(double r[][3], double field[]); 
    //calculate potential of the particle
    double elecenergy(double q, double p);
    //calculate electric potential energy of the particle
    double * rate(double q, double f[], double T , double step = 1e-12, double D = 1); 
    //calculate rate component in each direction
    
    ~Particle();
};

//Create a function to calculate compound index, if it might be needed
//Usage:
// input two indices - get the compound index
double CompIndex (int a, int b){
    double ab = a*(a+1)/(2+b);
    return ab;
}

//Usage:
//double **pdistx;
//double **pdisty;
//double **pdistz;
// double x = pos[0];
// double y = pos[1];
// double z = pos[2];
// double Gx[N];
// double Gy[N];
// double Gz[N];
// for(int i = 0; i<N; i++){
//     Gx[i] = object.geom[i][0];
//     Gy[i] = object.geom[i][1];
//     Gz[i] = object.geom[i][2];
// }
// pdistx = object.distance(x, Gx);
// pdisty = object.distance(y, Gy);
// pdistz = object.distance(z, Gz);
// double NetR[N][3];
// for(int i = 0; i<N; i++){
//         NetR[i][0] = *(pdistx+i);
//         NetR[i][1] = *(pdisty+i);
//         NetR[i][2] = *(pdistz+i);
//     } 
// }  
double * Particle::distance(double p, double G[]){
    len = sizeof(G) / sizeof(G[0]);
    double R[len];
    for(int i = 0; i<len; i++){
            R[i] = G[i] - p;
    }
    cout << "Calculated Distances";
    return R;
}

//Usage:
//double *pfield;
//pfield = object.elecfield(particle r[][3], geometry[N][3], charges[N], dielectric const);
//double NetE[3];
//NetE[0] = *pfield;
//NetE[1] = *(pfield+1);
//NetE[2] = *(pfield+2);  
double * Particle::elecfield(double r[][3], double C[], double k = 1.0){
    double E[3]; //stores the three components
    double pi = 3.14159265358979323846 ;//value of pi
    double e0 = 8.8541878128e-12 ;//permittivity of free space
    int len = sizeof(G) / sizeof(G[0]);
    for(int i = 0; i<len; i++){
        for(int j = 0; j<3; j++){
            E[0] += ((C[i])/(r[i][j]*r[i][j]))/(4*pi*e0*k);//C1/xi2
            }
    }
    cout << "Calculated Fields";
    return E;
}

//Usage:
//Potential = object.elecpot(particle r[][3], elecfield[3])
double Particle::elecpot(double r[][3], double field[]){
    double E = 0.0;
    int len = sizeof(r) / sizeof(r[0]);
    for(int j = 0; j<len; j++){
        for(int i = 0; i<3; i++){
            E += r[j][i]*field[i];
        }
    }
    cout << "Calculated Potentials";
    return E;
}

//Usage:
//Energy = object.elecenergy(particle q, elecpot)
double Particle::elecenergy(double q, double p){
    double E = q*p;
    cout << "Calculated Energies";
    return E;
}

//Usage:
//double *prate;
//prate = object.rate(particle q, net field[3], temp = 273.16 K, step size = 1ps, diffusion const = 1);
//double NetR[3];
//NetR[0] = *prate;
//NetR[1] = *(prate+1);
//NetR[2] = *(prate+2); 
double * Particle::rate(double q, double f[], double T = 273.16, double step = 1e-12, double D = 1.0){
    double r[3];
    double k = 1.380649e-23;
    for(int i = 0; i<3; i++){
            r[i] = (q*D*f[i])/(k*T*step);
        }
    cout << "Calculated Rates";
    return r;
}

//Usage:
//Disorder = object.disorder(rate[3])
double Particle::disorder(double r[]){
    double r1 = rand();
    double T;
    for(int i = 0; i<3; i++){
            T += r[i];
        }
    D = (-1)*(log(r1))/T; //log in cmath is ln
    cout << "Calculated disorder";
    return D;
}

//MAIN FUNCTION STARTS HERE
int main(){
    //read lattice site/protein state/anything from xyz/pdb
    //read atom position from pdb
    
    //initialize system to 0 - run constructor function on a new object
    State Si = State(N, GS, QS, AS);
    Particle Pi = Particle(PPi, QP, AP);
    //calculate initial distances
    double **pdistx;
    double **pdisty;
    double **pdistz;
    double x = PPi[0];
    double y = PPi[1];
    double z = PPi[2];
    double Gx[N];
    double Gy[N];
    double Gz[N];
    for(int i = 0; i<N; i++){
        Gx[i] = Si.GSi[i][0];
        Gy[i] = Si.GSi[i][1];
        Gz[i] = Si.GSi[i][2];
    }
    pdistx = Pi.distance(x, Gx);
    pdisty = Pi.distance(y, Gy);
    pdistz = Pi.distance(z, Gz);
    double NetRi[N][3]; // This will be needed later
    for(int i = 0; i<N; i++){
        NetR[i][0] = *(pdistx+i);
        NetR[i][1] = *(pdisty+i);
        NetR[i][2] = *(pdistz+i);
        } 
    }
    //calculate initial elec fields
    double *pfield;
    pfield = Pi.elecfield(NetRi, GSi, QS); //Assuming Dielectric Const = 1 now
    double NetEi[3];
    NetEi[0] = *pfield;
    NetEi[1] = *(pfield+1);
    NetEi[2] = *(pfield+2);
    //calculate initial potential
    NetPoti = Pi.elecpot(NetRi, NetEi);
    //calculate initial energy
    NetEni = Pi.elecenergy(QP, NetPoti);
    //calculate inital rate
    double *PRatei;
    PRatei = Pi.rate(QP, NetEi);
    double NetRatei[3];
    NetRatei[0] = *PRatei;
    NetRatei[1] = *(PRatei+1);
    NetRatei[2] = *(PRatei+2);
    double RateSumI = NetRatei[0]+NetRatei[1]+NetRatei[2];
    double Disorder = 0.0;//no initial disorder

    //Iteration State and Particle
//     double GSf[N][3];
//     double QSf[N];
//     double ASf[N];
    double PPf[3];
    double NetRf[N][3];
    
    for(int i = 0; i<N, i++){
        for(int j = 0; j<3; j++){
//             GSf[i][j] = GS[i][j];
            NetRf[i][j] = NetRi[i][j];
        }
//         QSf[i] = QSi[i];
//         ASf[i] = ASi[i];
    }
    for(int j = 0; j<3; j++){
        PPf[j] = PPi[j];  
    }
    
    int NSteps = 1000; //for now
    State Sf = State(N, GS, QS, AS);
    Particle Pf = Particle(PPf, QP, AP);
    auto t_start = std::chrono::high_resolution_clock::now(); //system start time
    double time_start = 0.0; //this is the simulation time
    double time_end = 0.0;
    double del_t = 0.0;
    double *pfieldi;
    double NetEf[3];
    double *PRatef;
    double NetRatef[3];
    double RateSumF;
    double disp[NSteps];
    double Traj[NSteps][3];
    for(int i = 0; i<3; i++){
        disp[0] += (PPf[i]-PPi[i])*(PPf[i]-PPi[i]);
        Traj[0][i] = PPf[i]; 
    }
    double msd;
    double Diff_const = 0.0;
    
    //DEFINE TERMINATION CONDITIONS
    for(int w = 0; w<NSteps; w++){           //WRITE THIS AS A WHILE LOOP
        //iteration particle initialize
        Pf = Particle(PPf, QP, AP);
        //calculate the fields
        pfieldi = Pf.elecfield(NetRf, GS, QS);
        NetEf[0] = *pfieldi;
        NetEf[1] = *(pfieldi+1);
        NetEf[2] = *(pfieldi+2);
        //calculate the potentials
        NetPotf = Pf.elecpot(NetRf, NetEf);
        //calculate the energy
        NetEnf = Pf.elecenergy(QS, NetPotf);
        //make list of all possible rates = eqn 5 from the proposal
        PRatef = Pf.rate(QP, NetEf);
        NetRatef[0] = *PRatef;
        NetRatef[1] = *(PRatef+1);
        NetRatef[2] = *(PRatef+2);
        //calculate net rate as sum of all rates acting on the system
        RateSumF = NetRatef[0]+NetRatef[1]+NetRatef[2];
        //also add the rates for waiting
        //generate a random number
        double rn1 = rand();
        
        //find event to carry out - waiting
        //determine what direction and magnitude to hop
        //find event to carry out - hopping
        
        //carry out the event
        //change the PPf values - and store them in the trajectory - add to disp
        
        //recalculate all rates - this is done at start of every iteration
        //generate another random number
        double rn2 = rand();
        //update the simulation time t+=del t where del t = -ln(random no)/net R
        del_t = (-1)*(log(rn2))/RateSumF;
        time_end += del_t;
        
        //CHECK FOR VANISHING DERIVATIVE (but that's only at steady-state)?????
        //go back with new state
    } //end of iteration loop
    //once the iterations are done, calculate diffusion coeff
    //D = avg of x^2 / 6* time - time = total time
    for(int i = 0; i<NSteps; i++){
        msd += disp[i];
    }
    msd = msd/NSteps;
    Diff_const = msd/(6*time_end);
    double rmsd = sqrt(msd);
    
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
    cout << "Total System Time Elapsed" << elapsed_time_ms;
    
    double rTraj[NSteps];
    double xTraj[NSteps];
    double yTraj[NSteps];
    double zTraj[NSteps];
    
    for(int i = 0; i<NSteps; i++){
        xTraj[i] += Traj[i][0];
        yTraj[i] += Traj[i][1];
        xTraj[i] += Traj[i][2];
        rTraj[i] += sqrt(((xTraj[i])*(xTraj[i]))+((yTraj[i])*(yTraj[i]))+((zTraj[i])*(zTraj[i])));
    }
    //plot rTraj and disp
    return 0;
}





























