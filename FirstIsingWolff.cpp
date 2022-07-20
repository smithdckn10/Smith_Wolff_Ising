#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <string>

#define L 100
#define SIZE L*L
#define UPDATES 10000
#define J 0.01  //If J < 0, then the model is a ferromagnet, if J > 0, then the model is an anti-ferromagnet, but I've got negative signs in the code right now so don't make it negative here
                //As of right now, this J isn't really doing much, since I modified the code to run the sim with multiple values of J
                //This next variable, numberOfTemps, tells the program how many different values of J it will use, 
                //and to change these, there is an array lower down where you can enter whatever J values you'd like, as long as you enter the correct amount (numberOfTemps)
#define numberOfTemps 3
using namespace std;

//create the array of "particles",


//This initializes the array with all the spins pointing up
void orderedinit(int spins[]);

//This initializes the array with random spins
void randominit(int spins[]);


//Get values for neighborhood of each site
void neighborhood(int size[], int neighbors[]);


//This outputs the spins to terminal for testing
void output(int spins[]);

//Checks if a lattice site should be flipped, and returns true or false
bool checkflip(double Jval, int spins[], int pos);

double getTotEnergy(int spins[], int neighbors[]);

double getChangeLocEnergy(double Jval, int spins[], int pos);

double getMagnetization(int spins[]);

//The simpleWriteEnergy function just outputs the energies and magnetizations at the end of each sweep, and leaves the data analysis to python
void simpleWriteEnergy(double Jval, double energies[], double magnets[]);

void getNeighborPositions(int spins[], int neighborPos[]);

//void wolffy(int initPos, int spins[], int clusterPositions[], int checkedSites[]);

void wolf(int initPos, int spins[], int neighborPos[], int checked[], double Jval);

void update(int spins[], int checked[]);



int main(){
    cout << "WOLFF TIME" << endl;
    //seeds the random number generator
    srand((unsigned) time(NULL));

    //Stores spins of each lattice site, with 1 meaning up, 0 meaning down
    int spins[SIZE];

    //Stores whether or not the sites have been checked by wolff; -1 = unchecked, 0 = checked not included, 1 = included
    int checked[SIZE];

    //Stores spins of each site's neighbors, so you can calc probability of whether or not the site will flip
    int neighbors[SIZE * 4];

    //Stores the positions of each sites neighbors, since periodic boundries are silly. 
    int neighborPos[SIZE*4];

    getNeighborPositions(spins, neighborPos);

    //Boolean that decides whether or not the spin is flipped
    bool YoN = 0;

    //Array that stores the total energy of the system after each sweep
    double energies[UPDATES];

    //Array that stores the total magnetization of the system after each sweep
    double magnets[UPDATES];
    

    double Jvalues[numberOfTemps] = {0.4, 0.45, 0.50};
    double Jval = 0.01;

    int initialSite = 0;

    cout << "SIZE: " << SIZE << "  UPDATES: " << UPDATES << "  TEMPS: " << Jvalues[0];
    for(int i = 0; i < numberOfTemps; i++){
        cout << ", " << Jvalues[i];
    }
    cout << endl;

    //This is for writing data to the file, not necessary rn
    ofstream EnergyText;
    string fileName = "Smith_Wolff_Energies.txt";
    EnergyText.open(fileName);
    EnergyText << SIZE << "," << UPDATES;
    for(int i = 0; i < numberOfTemps; i++){
        EnergyText << "," << Jvalues[i];
    }
    EnergyText << endl;
    EnergyText.close();
    

    //I think this should be the beginning of the loop for multiple temps
    for(int t = 0; t < numberOfTemps; t++){                                                 //CHANGE THE ONE TO numberOfTemps!!!!!!!!!!!!!!!!!!!!!!!!!!
        cout << "Current Temp (J) : " << Jvalues[t] << endl;

        for(int i = 0; i < UPDATES; i++){ //initialize the energies and magnetizations arrays with zeroes
            energies[i] = 0;
            magnets[i] = 0;
        }

        // This sets the spin values to 1
        //orderedinit(spins);

        //This initializes the spins randomly
        randominit(spins);
        

        for(int upd = 0; upd < UPDATES; upd++){ //This is where Wolff occurs

            for(int i = 0; i < SIZE; i++){ //Reset the checked array so everyone is unchecked
                checked[i] = -1;
            }

            initialSite = rand() % int(SIZE); //Get a new initial position

            wolf(initialSite, spins, neighborPos, checked, Jvalues[t]); //Run the recursive wolff

            update(spins, checked);

            /*
            output(spins);
            cout << endl;
            */


            neighborhood(spins, neighbors);

            energies[upd] = getTotEnergy(spins, neighbors); //Add the energy of this sweep to the array
            
            magnets[upd] = getMagnetization(spins); //Add the magnetization of this sweep to the array
            
        }




        simpleWriteEnergy(Jvalues[t] ,energies, magnets); //Output the energy after every update for each Jval to a text file, using simpleWriteEnergy

    }

    EnergyText.close();
    cout << "Done!" << endl;
    return 0;
}

void orderedinit(int spins[]){
    for(int i = 0; i < SIZE; i++){
        spins[i] = 1;
    }
}

void randominit(int spins[]){
    int rando = 0;
    for(int i = 0; i < SIZE; i++){
        rando = rand() % 2;
        if(rando == 1){
            spins[i] = 1;
        }
        else{
            spins[i] = -1;
        }
    }
}

bool checkflip(double Jval, int spins[], int pos){
    double prob = 0;
    double dE = getChangeLocEnergy(Jval, spins, pos);
    double Bweight = exp(-dE);
    double check = rand();
    check = check / RAND_MAX;
    //cout << "dE: " << dE << " , check: " << check << endl;
    prob = min(1.0, Bweight);
    if(check <= prob){
        return true;
    }
    else{
        return false;
    }

}

double getTotEnergy(int spins[], int neighbors[]){
    double Energy = 0; //Stores Total Energy/Hamiltonian I think?
    double extMagF = 0; //Contribution from external magnetic field
    double neighborF = 0; //Contribution from neighbor sites
    int B = 0; //External Magnetic Field
    for(int i = 0; i < SIZE; i++){
        neighborF = neighborF + spins[i] * neighbors[4*i]; //Get the top neigbor term
            //cout << "i: " << i << "  neighborF: " << neighborF << endl;
        neighborF = neighborF + spins[i] * neighbors[4*i + 3]; //Get the right neighbor term
            //cout << "i2: " << i << "  neighborF: " << neighborF << endl;

        //extMagF = extMagF - spins[i]; //Get the external magnet term 
    }    
    //Energy = (float(J) * neighborF) + (B * extMagF);
    Energy = (neighborF) + (B * extMagF);

    //LEts try this ig? I think this is what PrimerLi is doing
    Energy = Energy / float(SIZE);
    return Energy;
}

double getChangeLocEnergy(double Jval, int spins[], int pos){
    //We need to find all the neighbors of the site to calc the energy, so uses code from the neighborhood subroutine
    int site = spins[pos];
    int antisite = site * -1;
    int top = 0;
    int bottom = 0;
    int left = 0;
    int right = 0;
    double change = 0;
    double before = 0;
    double after = 0;
    double energyI = 0;

    if(pos % L == 0){             //Check if its on the top side
            
        before = site * spins[pos + L - 1];
        after = antisite * spins[pos + L - 1];
    }
    else{
        before = site * spins[pos- 1];
        after = antisite * spins[pos - 1];
    }
    if(pos % L == L - 1){      //Check if its on the bottom side
        before = before + site * spins[pos - L + 1];
        after = after + antisite * spins[pos - L + 1];
    }
    else{
        before = before + site * spins[pos + 1];
        after = after + antisite * spins[pos + 1];
    }
    if(pos < L){              //Check if its on the left
        before = before + site * spins[SIZE - L + pos];
        after = after + antisite * spins[SIZE - L + pos];
    }
    else{
        before = before + site * spins[pos - L];
        after = after + antisite * spins[pos - L];
    }
    if(pos >= SIZE - L){              //Check if its on the right
        before = before + site * spins[pos % L];
        after = after + antisite * spins[pos % L];
    }
    else{
        before = before + site * spins[pos + L];
        after = after + antisite * spins[pos + L];
    }
    //ENERGY/HAMILTONIAN EQUATION    
    change = -1*Jval * ((after) - (before));
    //change = -J / 2 * ((after) - (before));

    return change;

}

double getMagnetization(int spins[]){
    double magnetization = 0;
    for(int i = 0; i < SIZE; i++){
    magnetization = magnetization + spins[i];
    }
    return sqrt(magnetization * magnetization) / float(SIZE);
}

void neighborhood(int spins[], int neighbors[]){ //This creates an array that stores each site's neighbor's spins
    int top = 0;
    int bottom = 0;
    int left = 0;
    int right = 0;
    int counter;

    for(int i = 0; i < SIZE; i++){
        
                 
        //Do the stuff to get periodic boundries, might make this a switch if it really matters

        if(i % L == 0){             //Check if its on the top side
            top = i + L - 1;
        }
        else{
            top = i - 1;
        }

        if(i % L == L - 1){      //Check if its on the bottom side
            bottom = i - L + 1;
        }
        else{
            bottom = i + 1;
        }

        if(i < L){              //Check if its on the left
            left = SIZE - L + i;
        }
        else{
            left = i - L;
        }

        if(i >= SIZE - L){              //Check if its on the right
            right = i % L;
        }
        else{
            right = i + L;
        }

        neighbors[4*i] = spins[top]; //Get the top neighbor
        neighbors[4*i + 1] = spins[bottom]; //Get the bottom neighbor
        neighbors[4*i + 2] = spins[left]; //Get the left neighbor
        neighbors[4*i + 3] = spins[right]; //Get the right neighbor


    }



}

void output(int spins[]){
    string space = "   ";
    for(int i = 0; i < L; i++){
        for(int j = 0; j < L; j++){

                //Enable this if you want to show the positions
            space = "   ";
            if(((L*j + i) > 9) && ((L*j + i) < 100)){
                space = "  "; 
            }
            if(((L*j + i) >= 100) && ((L*J + i) < 1000)){
                space = " ";
            }
            
            cout << "[" << (int(L) *j + i) << "]";
            

            cout << spins[L*j + i] << space; //This just marks the position before the spins, not super necessary
            if(spins[i + L * j] > 0){
                cout << " ";
            }
        }
        cout << endl;
    }
}

void simpleWriteEnergy(double Jval, double energies[], double magnets[]){
    ofstream EnergyText;
    string filename = "Smith_Wolff_Energies.txt";
    EnergyText.open(filename, ios::app);
    for(int i = 0; i < UPDATES; i++){
        EnergyText << Jval << ",[" << i << "]," << energies[i] << "," << magnets[i] << endl;
    }
}

void getNeighborPositions(int spins[], int neighborPos[]){ //This is like neighborhood, except for each site, it just stores the position of the neighbors, so for site 0 in a 4x4 lattice , neighbors[0 + 0] = 3 (the top neighbor of site 0 is site 3)
    int top = 0;    //You should only need to call this function once, at the very beginning 
    int bottom = 0;
    int left = 0;
    int right = 0;
    int counter;

    for(int i = 0; i < SIZE; i++){
        
                 
        //Do the stuff to get periodic boundries, might make this a switch if it really matters

        if(i % L == 0){             //Check if its on the top side
            top = i + L - 1;
        }
        else{
            top = i - 1;
        }

        if(i % L == L - 1){      //Check if its on the bottom side
            bottom = i - L + 1;
        }
        else{
            bottom = i + 1;
        }

        if(i < L){              //Check if its on the left
            left = SIZE - L + i;
        }
        else{
            left = i - L;
        }

        if(i >= SIZE - L){              //Check if its on the right
            right = i % L;
        }
        else{
            right = i + L;
        }

        neighborPos[4*i] = top; //Get the top neighbor
        neighborPos[4*i + 1] = bottom; //Get the bottom neighbor
        neighborPos[4*i + 2] = left; //Get the left neighbor
        neighborPos[4*i + 3] = right; //Get the right neighbor


    }
}

/* wolffy function, but now I've moved on to wolf now
void wolffy(int initPos, int spins[], int clusterPositions[], int checkedSites[]){ 
    //Every element of clusterPositions should be one before being passed to this function
    //initPos should be a randomly selected site, unless your calling the function within the function (recursion?), where it should be the neighbor of a site in the cluster, and it should have been added to the cluster already
    //Every element of checkedSites should be -1 before being passed to this function, unless your calling from within the function
    //This is probably not my best work, but this function is hard to optimize

    // int checkedSites[SIZE] This array keeps track of which sites have already been checked to see if they join the cluster, 1=true

    for(int i = 0; i < SIZE; i++){ //this was implemented before checkedSites was defined outside this function, and omitted for recursion purposes
        checkedSites[i] = -1;
    }
    
    
    int rando = rand();
    rando = rando % SIZE; //randomly select a site at position rando
    

    clusterPositions[initPos] = -1; //The first site added to the cluster is the site at position: initPos
    checkedSites[initPos] = 0;

    //Here we get the positions of the neighbors
    int t = 0;
    int b = 0;
    int l = 0;
    int r = 0;
    double chance = rand() / RAND_MAX;


    getNeighborPositions(initPos, t, b, l, r);
    
    //Now we have the positions of all the neighboring sites, and we can check to see if they get added to the cluster

    //I think in this model I've absorbed the beta term into the J term, but I'm not certain
    chance = rand() / RAND_MAX;
    if(checkedSites[t] == -1){ //If the top neighbor hasn't been checked yet
        double pAdd = 0;
        double ePower = min(0.0 , (4 * J * spins[initPos] * spins[t])); //In the textbook there should be a beta term in there
        pAdd = 1 - exp(ePower);
        if(pAdd >= chance){
            clusterPositions[t] = -1; //This means that the top neighbor is in the cluster, but we haven't checked all of that sites neighbors yet
            checkedSites[t] = 0;    //This means that the top neighbor has been checked and made it into the cluster, but none of its neighbors have been checked yet
        }
        else{
            checkedSites[t] = -2; //This means that the top neighbor has been checked, and it didn't make it into the cluster
        }
    }

    chance = rand() / RAND_MAX;
    if(checkedSites[l] == -1){ //If the left neighbor hasn't been checked yet
        double pAdd = 0;
        double ePower = min(0.0 , (4 * J * spins[initPos] * spins[l])); //In the textbook there should be a beta term in there
        pAdd = 1 - exp(ePower);
        if(pAdd >= chance){
            clusterPositions[l] = -1; //This means that the left neighbor is in the cluster, but we haven't checked all of that sites neighbors yet
            checkedSites[l] = 0;    //This means that the left neighbor has been checked and made it into the cluster, but none of its neighbors have been checked yet
        }
        else{
            checkedSites[l] = -2; //This means that the left neighbor has been checked, and it didn't make it into the cluster
        }
    }

    chance = rand() / RAND_MAX;
    if(checkedSites[b] == -1){ //If the bottom neighbor hasn't been checked yet
        double pAdd = 0;
        double ePower = min(0.0 , (4 * J * spins[initPos] * spins[b])); //In the textbook there should be a beta term in there
        pAdd = 1 - exp(ePower);
        if(pAdd >= chance){
            clusterPositions[b] = -1; //This means that the bottom neighbor is in the cluster, but we haven't checked all of that sites neighbors yet
            checkedSites[b] = 0;    //This means that the bottom neighbor has been checked and made it into the cluster, but none of its neighbors have been checked yet
        }
        else{
            checkedSites[b] = -2; //This means that the bottom neighbor has been checked, and it didn't make it into the cluster
        }
    }

    chance = rand() / RAND_MAX;
    if(checkedSites[r] == -1){ //If the right neighbor hasn't been checked yet
        double pAdd = 0;
        double ePower = min(0.0 , (4 * J * spins[initPos] * spins[r])); //In the textbook there should be a beta term in there
        pAdd = 1 - exp(ePower);
        if(pAdd >= chance){
            clusterPositions[r] = -1; //This means that the right neighbor is in the cluster, but we haven't checked all of that sites neighbors yet
            checkedSites[r] = 0;    //This means that the right neighbor has been checked and made it into the cluster, but none of its neighbors have been checked yet
        }
        else{
            checkedSites[r] = -2; //This means that the right neighbor has been checked, and it didn't make it into the cluster
        }
    }

    checkedSites[initPos] = 1;

    if(checkedSites[t] == 0){
        wolffy(t, spins, clusterPositions, checkedSites); //run this function again, except for the top neighbor
    }
    if(checkedSites[l] == 0){
        wolffy(l, spins, clusterPositions, checkedSites); //run this function again, except for the left neighbor
    }
    if(checkedSites[b] == 0){
        wolffy(b, spins, clusterPositions, checkedSites); //run this function again, except for the bottom neighbor
    }
    if(checkedSites[r] == 0){
        wolffy(r, spins, clusterPositions, checkedSites); //run this function again, except for the right neighbor
    }


}
*/
    

void wolf(int initPos, int spins[], int neighborPos[], int checked[], double Jval){
    //each site has an attribute checked, so -1 = unchecked, 0 = checked and not in, 1 = in the cluster
    double chance;
    double prob;
    checked[initPos] = 1;

    //Get the neighbor positions
    int t = neighborPos[4*initPos];
    int b = neighborPos[4*initPos + 1];
    int l = neighborPos[4*initPos + 2];
    int r = neighborPos[4*initPos + 3];
    //check if the top neighbor gets in
    if(spins[t] == spins[initPos]){ //Check if top neighbor has the same spin
        if(checked[t] == -1){   //Test if the top neighbor is unchecked
            prob = 1 - exp(min(0.0 , (-4*Jval))); //NORMALLY THIS WOULD BE A MORE COMPLEX FUNCTION, but for the Ising system, both sites need to be the same spin for the new site to join the cluster. 
                                                //Since the spins can only be -1 or 1, then sites will only join when prob = 1 - 4*J*(1)(1) or prob = 1-4*J*(-1)(-1)
                                                //These are the same values, but when switching to XY model, fix this up 

            chance = rand() / float(RAND_MAX);
            //cout << "site[ " << initPos << "] = " << spins[initPos] << ", t: " << t << endl;
            //cout << "prob = " << prob << ", chance = " << chance << endl;

            if(prob > chance){  //Test if the top neighbor gets in
                checked[t] = 1;
                wolf(t, spins, neighborPos, checked, Jval); //Run wolff again if the top neighbor gets in
            }
            else{
                checked[t] = 0;
            }
        }  
    }

    if(spins[b] == spins[initPos]){ //Check if bottom neighbor has the same spin
        if(checked[b] == -1){   //Test if the bottom neighbor is unchecked
            //prob = 1 - exp(min(0.0 , (-4*Jval))); //NORMALLY THIS WOULD BE A MORE COMPLEX FUNCTION, but for the Ising system, both sites need to be the same spin for the new site to join the cluster.
            prob = 1 - exp(min(0.0, (-2*Jval)));
            chance = rand() / float(RAND_MAX);
            if(prob > chance){  //Test if the bottom neighbor gets in
                checked[b] = 1;
                wolf(b, spins, neighborPos, checked, Jval); //Run wolff again if the bottom neighbor gets in
            }
            else{
                checked[b] = 0;
            }
        }  
    }

    if(spins[l] == spins[initPos]){ //Check if left neighbor has the same spin
        if(checked[l] == -1){   //Test if the left neighbor is unchecked
            prob = 1 - exp(min(0.0 , (-4*Jval))); //NORMALLY THIS WOULD BE A MORE COMPLEX FUNCTION, but for the Ising system, both sites need to be the same spin for the new site to join the cluster.
            chance = rand() / float(RAND_MAX);
            if(prob > chance){  //Test if the left neighbor gets in
                checked[l] = 1;
                wolf(l, spins, neighborPos, checked, Jval); //Run wolff again if the left neighbor gets in
            }
            else{
                checked[l] = 0;
            }
        }  
    }

    if(spins[r] == spins[initPos]){ //Check if right neighbor has the same spin
        if(checked[r] == -1){   //Test if the right neighbor is unchecked
            prob = 1 - exp(min(0.0 , (-4*Jval))); //NORMALLY THIS WOULD BE A MORE COMPLEX FUNCTION, but for the Ising system, both sites need to be the same spin for the new site to join the cluster.
            chance = rand() / float(RAND_MAX);
            if(prob > chance){  //Test if the right neighbor gets in
                checked[r] = 1;
                wolf(r, spins, neighborPos, checked, Jval); //Run wolff again if the right neighbor gets in
            }
            else{
                checked[r] = 0;
            }
        }  
    }

}

void update(int spins[], int checked[]){
    for(int i = 0; i < SIZE; i++){
        if(checked[i] == 1){
            spins[i] = spins[i]*-1;
        }
    }
}