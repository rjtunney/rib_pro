//
//  simulateRib.cpp
//
//  Created by Jasmine on 8/28/13.
//
//

#include "simulateRib.h"
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>

using namespace std;

struct singleStepOutput {
	vector<int> lattice;
	double arrRate;
	vector<int> pos;
	vector<int> states;
};

struct burnInOutput {
	vector<int> lattice;
	vector<int> pos;
	vector<int> states;
};

struct simOutput {
	vector<int> lattice;
	double prot;
};

int initiateQ( vector<int> lattice, int w){
	int rightMaxPos = min( ((int) lattice.size()-1), w );
	int occ = 0;
	for ( int j = 0; j < rightMaxPos; ++j ){
		occ += lattice.at(j);
	}
	if ( occ > 0 ){
		return 0;
	}
	else {
		return 1;
	}
}

int moveToRightQ( vector<int> lattice, int trna, int pos, int w){
	int rightMaxPos	= min( ((int) lattice.size()-1), pos + w );
	int occ = 0;
	for ( int j = pos+1; j <= rightMaxPos; ++j ){
		occ += lattice.at(j);
	}
	if ( occ == 0 && trna == 1){
		return 1;
	}
	else {
		return 0;
	}
}

singleStepOutput simOneTimeStep( vector<int> lattice, vector<int> occupiedPos, vector<int> ribState, double alpha, double beta, vector<double> gamma, vector<int> delta, int w){
	vector<int> newLattice(lattice);
	
	// calculate total arrival rate
	int canEnter = initiateQ( lattice, w );
	double sumRates = 0;
	int sizeProp = occupiedPos.size() + 1;
	if ( canEnter == 1 ){
		sumRates += alpha;
	}
		
	for ( int j=0; j < occupiedPos.size(); ++j){
		if ( ribState.at(j) == 0 ) {
		//	cout << "i have no trna " << endl;
	//		cout << occupiedPos.at(j) << " " << ribState.at(j) << endl;
			sumRates += gamma.at(occupiedPos.at(j));
	//		sizeProp += 1;
		}
		else {
			int canMove = moveToRightQ( lattice, ribState.at(j), occupiedPos.at(j), w );
		//	cout << "can I move " << canMove << endl;
			if ( occupiedPos.at(j) != lattice.size() - 1){
				if ( canMove == 1 ) {
					sumRates += delta.at(occupiedPos.at(j));
	//				sizeProp += 1;
				}
				else {
					sumRates += 0;
				}
			}
			else {
				if ( ribState.at(j) == 1 ) {
					sumRates += beta;
	//				sizeProp += 1;
				}
				else {
					sumRates += gamma.at(occupiedPos.at(j));
	//				sizeProp += 1;
				}
			}
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
	for ( int j = 0; j < occupiedPos.size(); ++j ){
		if ( occupiedPos.at(j) < lattice.size() - 1 ){
			if ( ribState.at(j) == 0 ) {
				prop.at(j+1) = prop.at(j) + gamma.at(occupiedPos.at(j))/sumRates;
			}
			else if ( ribState.at(j) == 1 ) {
				int canMove = moveToRightQ( lattice, ribState.at(j), occupiedPos.at(j), w );
				if ( canMove == 1 ) {
					prop.at(j+1) = prop.at(j) + delta.at(occupiedPos.at(j))/sumRates;
				}
				else {
					prop.at(j+1) = prop.at(j);
				}
			}
		}
		else {
			if ( ribState.at(j) == 0 ) {
				prop.at(j+1) = prop.at(j) + gamma.at(occupiedPos.at(j))/sumRates;
			}
			else if ( ribState.at(j) == 1 ) {
				prop.at(j+1) = prop.at(j) + beta/sumRates;
			}
		}
	}

	// which ribosome will jump?

	double u = ((double) rand() / (RAND_MAX));
	int whichMove = 0;
	
	while ( prop.at(whichMove) < u ){
		whichMove++;
	}
	
//	cout << "ribosome " << whichMove << " moves" << endl;
	
	// now the jump! (or not)
	if ( whichMove == 0 ){
	//	cout << "enter " << endl;
		newLattice.at(0) = 1;
		occupiedPos.insert(occupiedPos.begin(), 0);
		ribState.insert(ribState.begin(), 1);
	}
	else {
	//	cout << "it is at " << occupiedPos.at(whichMove-1) << endl;
	//	cout << "it is in state " << ribState.at(whichMove-1) << endl;
		if ( occupiedPos.at(whichMove-1) < lattice.size() - 1 ) {
			if ( ribState.at(whichMove-1) == 1 ) {
				newLattice.at(occupiedPos.at(whichMove-1)) = 0;
				newLattice.at(occupiedPos.at(whichMove-1)+1) = 1;
				occupiedPos.at(whichMove-1) += 1;
				ribState.at(whichMove-1) = 0;
	//			cout << "now it is at " << occupiedPos.at(whichMove-1) << endl;
	//			cout << "see " << newLattice.at(occupiedPos.at(whichMove-1)) << endl;
	//			cout << "it is at state " << ribState.at(whichMove-1) << endl;
			}
			else {
				ribState.at(whichMove-1) = 1;
			}
		}
		else if ( occupiedPos.at(whichMove-1) == lattice.size() - 1 )  {
			if ( ribState.at(whichMove-1) == 1 ){
				occupiedPos.pop_back();
				ribState.pop_back();
				newLattice.at( newLattice.size() - 1 ) = 0;
			}
			else {
	//			cout << "i should come here " << endl;
				ribState.at(whichMove-1) = 1;
			}
		}
	}
			
	singleStepOutput result;
	result.lattice = newLattice;
	result.arrRate = sumRates;
	result.pos = occupiedPos;
	result.states = ribState;
	
	return result;
	
}

burnInOutput burnIn( int prot4BurnIn, int L, double alpha, double beta, vector<double> gamma, vector<int> delta, int w){

	int nProt = 0; // how many proteins are made

	vector<int> lattice(L,0); // mRNA is empty to start
	vector<int> occupiedPos(0,0);
	vector<int> ribState(0,0);
	
	while (nProt < prot4BurnIn){
		int prevLast = lattice.at(lattice.size() - 1);
		singleStepOutput result = simOneTimeStep(lattice, occupiedPos, ribState, alpha, beta, gamma, delta, w);
		lattice = result.lattice;
		occupiedPos = result.pos;
		ribState = result.states;
		if ( prevLast == 1 && lattice.at(L-1) == 0 ){
			nProt ++;
	//		cout << "protein made" << endl;
		}
	}
	
	burnInOutput burnOut;
	burnOut.lattice = lattice;
	burnOut.pos = occupiedPos;
	burnOut.states = ribState;
	
	return burnOut;
}


simOutput simulation( ofstream &outfile, int maxTime, vector<int> occupiedPos, vector<int> ribState, vector<int> lattice, double alpha, double beta, vector<double> gamma, vector<int> delta, int w){
	
	double timePassed = 0;
	
	int iter = 0;
	int nProt = 0;
	
	while (timePassed < maxTime){
		int prevLast = lattice.at(lattice.size()-1);
		singleStepOutput result = simOneTimeStep(lattice, occupiedPos, ribState, alpha, beta, gamma, delta, w);
		// draw number from exponential distribution to determine time step
		double totalRate = result.arrRate;
		double timeStep = (log(1 - ((double) rand() / (RAND_MAX))))/(-totalRate);
		timePassed += timeStep;
		if (timePassed < maxTime) {
			iter++;
			lattice = result.lattice;
			occupiedPos = result.pos;
			ribState = result.states;
			if ( prevLast == 1 && lattice.at(((int) lattice.size()-1 )) == 0) {
				nProt++;
			}
		}
	}
	
	double protMade = ((double) nProt / timePassed);
	
	outfile << "Number of iterations: " << iter << endl;
	outfile << "Number of proteins produced: " << nProt << endl;
	outfile << "Protein production per unit time: " << protMade << endl;
	
	simOutput simOut;
	simOut.lattice = lattice;
	simOut.prot = protMade;
	
	return simOut;
}


/////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]){
	// make sure we have info file
	if (argc < 4) {
		cerr << "Usage: " << argv[0] << " info_file.txt gene_wanted init_rate output_file (optional)" << endl;
		return 1;
	}
	
	// find info for the gene we want
	fstream fin(argv[1]);
	string gene, tai, seq;
	while (fin){
		fin >> gene >> tai;
		if (gene == argv[2]){
			break;
		}
	}
	//assign tRNA capture rate vector gamma
	vector<double> gamma; //  tRNA capture rate vector
	
	stringstream ss(tai);
	float curr_tai;
	while (ss >> curr_tai){
		if (ss.peek()==','){
			ss.ignore();
		}
		gamma.insert(gamma.begin(), 1000*curr_tai); //16.7735*curr_tai);
	}
	int L = gamma.size(); // number of codons
	
	ofstream outfile;
	if (argc == 5) {
		// output to this file
		outfile.open(argv[4]);
	}
	// seed for random number generator
	srand(time(NULL));
	
	// variables regarding lattice size and rates
	// eventually delta should be read from another function
	vector<int> delta(L-1,1); // translocation rate vector (constant for now)
	int w = 1; // width of ribosomal footprint in codons
	double beta = atof(argv[3]); // initiation rate
	double alpha = 1; //35; // termination rate
	
	// variables regarding time and output of simulation
	int totalLatticePop = 0;
	int prot4BurnIn = 10; // number of proteins to end burn-in period
	int maxTime = 2000; // time for simulation after burn-in
	int num_runs = 10; // number of runs
	vector<double> latticePop(L,0); // population of mRNAs
	vector<double> avgDensity(num_runs,0);
	vector<double> avgProt(num_runs,0);
	double overallAve = 0;
	double overallAvgProt = 0;
	for (int n=0; n < num_runs; ++n){
		if (argc == 5) {
			outfile << "ITERATION " << n+1 << endl;
		}
		
		burnInOutput burnOut = burnIn(prot4BurnIn, L, alpha, beta, gamma, delta, w);
		vector<int> lattice = burnOut.lattice;
		vector<int> occupiedPos = burnOut.pos;
		vector<int> ribState = burnOut.states;

		int rib0 = 0; // initial number of ribosomes on mRNA	
		for( int j = 0; j < L; ++j ){
			rib0 += lattice.at(j);
		}
		double initDens = ((double) rib0/L); // initial density

		simOutput simOut = simulation(outfile, maxTime, occupiedPos, ribState, lattice, alpha, beta, gamma, delta, w);
		lattice = simOut.lattice;
		overallAvgProt += simOut.prot;
		avgProt[n] = simOut.prot;
	
		int ribEnd = 0; // number of ribosomes at end of simulation
		for( int j = 0; j < L; ++j ){
			ribEnd += lattice.at(j);
		}
		avgDensity[n] = ((double) ribEnd/L); // end density
		
		if (argc == 5){
			outfile << "Initial number of ribosomes on the mRNA: " << rib0 << ",Density: " << initDens << endl;
			outfile << "Final number of ribosomes on the mRNA: " << ribEnd << ", Density: " <<  avgDensity[n] << endl;
			outfile << "Ribosome distribution: " ;
		}
		for ( int i=0; i < lattice.size(); ++i ){
			latticePop[i] += lattice[i];
			totalLatticePop += lattice[i];
			if (argc == 5) {
				outfile << lattice[i] << " " ;
			}
		}
		overallAve += avgDensity[n];
		
		if (argc == 5) {
			outfile << "\n" << endl;
		}
	}
	overallAve /= num_runs;
	overallAvgProt /= num_runs;
	
	double varProt = 0;
	double varDensity = 0;
	for (int n = 0; n < num_runs; ++n){
		varDensity += pow((avgDensity[n] - overallAve),2.0);
		varProt += pow((avgProt[n] - overallAvgProt),2.0);
	}
	varDensity /= num_runs;
	varProt /= num_runs;
	
	for (int i = 0; i < latticePop.size(); ++i) {
		latticePop[i] /= totalLatticePop;
	}

	if (argc == 5) {
		outfile << "\n" << "Population ribosome distribution: ";
		for ( int i=0; i < latticePop.size(); ++i){
			outfile << latticePop[i] << " ";
		}
	}
	else {
		for ( int i=0; i < latticePop.size(); ++i){
			cout << latticePop[i] << " ";
		}
	}
	if (argc == 5) {
		outfile << "\n" << "Mean ribosome density: " << overallAve << endl;
		outfile << "Variance in ribosome density: " << varDensity << endl;
		outfile << "Average protein production: " << overallAvgProt << endl;
		outfile << "Variance in protein production: " << varProt << endl;
	}
	else {
		cout << "\n" << "avg " << overallAve << endl;
		cout << "varD " << varDensity << endl;
		cout << "nprot " << overallAvgProt << endl;
		cout << "varP " << varProt << endl;
	}

}