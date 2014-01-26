//
//  simulateRib.cpp
//
//  Created by Jasmine on 8/28/13.
//
//

#include "simulateRib.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>

using namespace std;

struct singleStepOutput {
	vector<int> lattice;
	double arrRate;
	vector<int> pos;
};

struct burnInOutput {
	vector<int> lattice;
	vector<int> pos;
};

int initiateQ( vector<int> lattice, int w){
	int rightMaxPos = min( ((int) lattice.size()-1), w );
	int occ = 0;
	for ( int j = 0; j <= rightMaxPos; ++j ){
		occ += lattice.at(j);
	}
	if ( occ > 0 ){
		return 0;
	}
	else {
		return 1;
	}
}

int moveToRightQ( vector<int> lattice, int pos, int w){
	int rightMaxPos	= min( ((int) lattice.size()-1), pos + w );
	int occ = 0;
	for ( int j = pos+1; j <= rightMaxPos; ++j ){
		occ += lattice.at(j);
	}
	if ( occ > 0 ){
		return 0;
	}
	else {
		return 1;
	}
}

singleStepOutput simOneTimeStep( vector<int> lattice, vector<int> occupiedPos, double alpha, double beta, vector<int> gamma, int w){
	vector<int> newLattice(lattice);

	/* old version using arrays instead of vectors
	 
	 for ( int l = 0; l <= L; ++l ){
	 //	newLattice[l] = lattice[l]
	 //	if (lattice[l] > 0 && l != L){
	 }
	 }
	 
	 */
	
	// calculate total arrival rate
	int canEnter = initiateQ( lattice, w );
	double sumRates = 0;
	int sizeProp = 1;
	vector<int> flexRib(0,0);
	if ( canEnter == 1 ){
		sumRates += alpha;
	}
	
	for ( int j=0; j < ((int) occupiedPos.size()); ++j){
		int canMove = moveToRightQ( lattice, occupiedPos.at(j), w );
		if ( occupiedPos.at(j) != lattice.size() - 1){
			if ( canMove == 1) {
				sumRates += gamma.at(occupiedPos.at(j));
				sizeProp += 1;
				flexRib.insert(flexRib.end(), j);
			}
			else {
				sumRates += 0;
			}
		}
		else {
			flexRib.insert(flexRib.end(), j);
			sumRates += beta;
			sizeProp += 1;
		}
	}

	vector<double> prop(sizeProp,0); // propensity function
	
	// add on initiation
	if ( canEnter == 1 ) {
		prop.at(0) = alpha/sumRates;
	}
	else{
		prop.at(0) = 0;
	}
	// cumulative sum over elongation rates
	for ( int j = 0; j < flexRib.size(); ++j){
		if ( occupiedPos.at(flexRib.at(j)) < lattice.size() - 1 ){
			prop.at(j+1) = prop.at(j) + gamma.at(occupiedPos.at(flexRib.at(j)))/sumRates;
		}
		else {
			prop.at(j+1) = prop.at(j) + beta/sumRates;
		}
	}

	// which ribosome will jump?

	double u = ((double) rand() / (RAND_MAX));
	int whichMove = 0;
	
	while ( prop.at(whichMove) < u ){
		whichMove++;
	}
	
	/* figure out how to work with upper_bound
	vector<int>::iterator whichMove;

	whichMove = upper_bound ( prop.begin(), prop.end(), u );
	
	*/
	
	// now the jump! (or not)
	if ( whichMove == 0 ){
		newLattice.at(0) = 1;
		occupiedPos.insert(occupiedPos.begin(), 0);
	}
	else {
		if ( occupiedPos.at(flexRib.at(whichMove-1)) < lattice.size() - 1 ) {
			newLattice.at(occupiedPos.at(flexRib.at(whichMove-1))) = 0;
			newLattice.at(occupiedPos.at(flexRib.at(whichMove-1))+1) = 1;
			occupiedPos.at(flexRib.at(whichMove-1)) += 1;
		}
		else if ( occupiedPos.at(flexRib.at(whichMove-1)) == lattice.size() - 1 )  {
			flexRib.pop_back();
			occupiedPos.pop_back();
			newLattice.at( newLattice.size() - 1 ) = 0;
		}
	}
	
	singleStepOutput result;
	result.lattice = newLattice;
	result.arrRate = sumRates;
	result.pos = occupiedPos;
	
	return result;
	
}

burnInOutput burnIn( int prot4BurnIn, int L, double alpha, double beta, vector<int> gamma, int w){

	int nProt = 0; // how many proteins are made

	vector<int> lattice(L,0); // mRNA is empty to start
	vector<int> occupiedPos(0,0);
	
	while (nProt < prot4BurnIn){
		int prevLast = lattice.at(lattice.size() - 1);
		singleStepOutput result = simOneTimeStep(lattice, occupiedPos, alpha, beta, gamma, w);
		lattice = result.lattice;
		occupiedPos = result.pos;
		if ( prevLast == 1 && lattice.at(L-1) == 0 ){
			nProt ++;
		}
	}
	
	burnInOutput burnOut;
	burnOut.lattice = lattice;
	burnOut.pos = occupiedPos;
	
	return burnOut;
}

// pass in pointers

vector<int> simulation( ofstream &outfile, int maxTime, vector<int> occupiedPos, vector<int> lattice, double alpha, double beta, vector<int> gamma, int w){
	
	double timePassed = 0;
	
	int iter = 0;
	int nProt = 0;
	
	while (timePassed < maxTime){
		int prevLast = lattice.at(lattice.size()-1);
		singleStepOutput result = simOneTimeStep(lattice, occupiedPos, alpha, beta, gamma, w);
		lattice = result.lattice;
		occupiedPos = result.pos;
		double totalRate = result.arrRate;
		
		// draw number from exponential distribution to determine time step
		double timeStep = (log(1 - ((double) rand() / (RAND_MAX))))/(-totalRate);
		timePassed += timeStep;
		iter++;
		if ( prevLast == 1 && lattice.at(((int) lattice.size()-1 )) == 0) {
			nProt++;
		}
	}
	
	double protMade = ((double) nProt / timePassed);
	
	outfile << "Number of iterations: " << iter << endl;
	outfile << "Number of proteins produced: " << nProt << endl;
	outfile << "Protein production per unit time: " << protMade << endl;
	
	return lattice;
}


/////////////////////////////////////////////////////////////////////////////

int main(){
	// output to this file
	ofstream outfile;
	outfile.open("ribTransSimData.txt");
	
	// seed for random number generator
	srand(time(NULL));
	
	// variables regarding lattice size and rates
	// eventually these should all be input from another function
	
	int L = 300; // number of codons
	int w = 9; // width of ribosomal footprint in codons
	double alpha = 1; // initiation rate
	double beta = 1; // termination rate
	vector<int> gamma(L-1,1); // elongation rate vector (constant for now)
	
	// variables regarding time and output of simulation
	
	int prot4BurnIn = 1000; // number of proteins to end burn-in period
	int maxTime = 1000; // time for simulation after burn-in
	
	burnInOutput burnOut = burnIn(prot4BurnIn, L, alpha, beta, gamma, w);
	vector<int> lattice = burnOut.lattice;
	vector<int> occupiedPos = burnOut.pos;
	
	int rib0 = 0; // initial number of ribosomes on mRNA	
	for( int j = 0; j < L; ++j ){
		rib0 += lattice.at(j);
	}
	double initDens = ((double) rib0/L); // initial density

	lattice = simulation(outfile, maxTime, occupiedPos, lattice, alpha, beta, gamma, w);
	int ribEnd = 0; // number of ribosomes at end of simulation
	for( int j = 0; j < L; ++j ){
		ribEnd += lattice.at(j);
	}
	double endDens = ((double) ribEnd/L); // end density
	
	outfile << "Initial number of ribosomes on the mRNA: " << rib0;
	outfile << ", Density: " << initDens << endl;
	outfile << "Final number of ribosomes on the mRNA: " << ribEnd << ", Density: " <<  endDens << endl;
	
}