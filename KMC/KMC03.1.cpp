//Kinetic Monte Carlo Project Part 3 Version 1
//C++ Version

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <stdlib.h>
#include <chrono>

using namespace std;

//make a class of a given state (Lattice)
//variables : no of atoms, geometry, force term
//functions : constructor/destructor function
class State 
{
    int NLattice;
    int NAtoms;
    double forceterm[3]; //stores force at a point in each direction
    double occupancy[NLattice][NLattice][NLattice];
    // double lattice[NLattice][NLattice][NLattice];
    // double * field_init_G(int N); //initalizes a normal Gaussian field
    // double * field_init_R(int N); //initialzes a random(0,1) field

    void State(int NL, int NA){
        cout << "Initializing State";
        NLattice = NL;
        NAtoms = NA;
        temp_NA = NA;
        for(int i=0; i<NL; i++){
            for(int j=0; j<NL; j++){
                for(int k=0; k<NL; k++){
                    if ((i+j+k)!=0){ //not the origin
                        if (temp_NA > 0){
                            int occ = (int) round(rand()/RAND_MAX);
                            occupancy[i][j][k] = 1;
                            temp_NA -=1;
                            }
                        else occupancy[i][j][k] = 0;
                        }
                    }
                }
            }
        }
        //do not calculate anything in the constructor function
        cout << "Initialized State S";
    }  
    
    public void State(){
        cout << "Initialized NULL State ";
    } 
    ~State();  
};

//Create a function to calculate compound index, if it might be needed
//Usage:
// input two indices - get the compound index
double CompIndex (int a, int b){
    double ab = a*(a+1)/(2+b);
    return ab;
};

//Create a distance file
//Usage:
//double *r;
//r = dist(pos1[], pos2[]);
//double distance[4];
//for(int i = 0; i<4; i++){
//     distance[i] = *(r+i);
// }
double * dist (int A[], int B[]){
    double d[4] = {0.0,0.0,0.0,0.0};
    d[0] = A[0]-B[0];
    d[1] = A[1]-B[1];
    d[2] = A[2]-B[2];
    int x = pow(d[0],2);
    int y = pow(d[1],2);
    int z = pow(d[2],2);
    d[3] = sqrt(x+y+z);
    return d;
}

//Create a function to assign normally distributed value to the field term
//Usage
//double * pF;
//pF = object.field_init_G(NAtoms);
//double F[N];
//for(int i = 0; i<N; i++){
//     F[i] = *(pF+i);
// }
// double * State::field_init_G(int N){
//     double nG_field[N];
//     std::default_random_engine generator;
//     std::normal_distribution<double> distribution(0.0,1.0);
//     for(int i =0; i<N; i++){
//         nG_field[i] = distribution(generator);
//     }
//     return nG_field;
// };

//Create a function to assign a random(0,1) value to the field term
//Usage
//double * pF;
//pF = object.field_init_R(NAtoms);
//double F[N];
//for(int i = 0; i<N; i++){
//     F[i] = *(pF+i);
//  }
// double * State::field_init_R(int N){
//     double nR_field[N];
//     for(int i =0; i<N; i++){
//         nR_field[i] = (double) rand()/RAND_MAX;
//     }
//     return nR_field;
// };

//MAIN FUNCTION STARTS HERE
int main(){
    //read NL and NA from user
    int t_NL = 0;
    int t_NA = 0;
    cout << "\n Enter Box length ";
    cin >> t_NL;
    cout << "\n Enter No of Atoms < length cubed ";
    cin >> t_NA;
    
    //initialize system to 0 - run constructor function on a new object
    State S = State(t_NL,t_NA);
    int pos[3] = {0,0,0};



