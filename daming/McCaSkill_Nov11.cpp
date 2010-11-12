#include <cstdio>
#include <iostream>
#include <cmath>
#include <cstring>
#include <algorithm>

using namespace std;

#define NMAX 600

static char RNA[NMAX]; // RNA sequence input
static double Q_m_Mtx[NMAX][NMAX]; // Q_m Matrix for dynamic programming
static double Q_b_Mtx[NMAX][NMAX]; // Q_b Matrix for dynamic programming
static double Q_Mtx[NMAX][NMAX]; // Q_b Matrix for dynamic programming

double k,T,a,b,c; // Parameters

int n; // actual RNA sequence length, get from input

//------------------------------------------
bool auPair(char x, char y){
  return ( x=='A'&&y=='U' || x=='U'&&y=='A' );
}

bool cgPair(char x, char y){
  return ( x=='C'&&y=='G' || x=='G'&&y=='C' );
}

bool basePair(int i, int j){ // Need Modification, use 2-D table
	char x = RNA[i];
	char y = RNA[j];
  return( auPair(x,y) || cgPair(x,y) );
}
//------------------------------------------
// initialize n*n matrix to zeros
void Zeros(double Mtx[NMAX][NMAX], int n){
	int i,j;
	for (i=0;i<n;i++) 
		for (j=0;j<n;j++) 
		 Mtx[i][j]=0;
}
//------------------------------------------
// toy F1() and F2()
double F1(int i,int j){
	return 2.71828;
}
double F2(int i,int j,int h,int l){
	return 3.14159;
}
//------------------------------------------

double Q_b(int,int); // function declaration

double Q_m(int i, int j){

	if(Q_m_Mtx[i][j]!=0) return Q_m_Mtx[i][j];

	double sigma=0.0;
	// Eq.(8)
	int h,l;
	for(l=j; l>i; l--){  
		for(h=l-1; h>=i; h--){
			if(basePair(h,l)==true){
				
				sigma += (exp(-(c*(h-i-1))/(k*T))+Q_m(i,h-1))*Q_b(h,l)*exp(-((b+c*(j-l-1))/(k*T)));

			}
		}
	}
	//update matrix
	Q_m_Mtx[i][j] = sigma;
	return sigma;
}


double Q_b(int i, int j){ 

	if(Q_b_Mtx[i][j]!=0) return Q_b_Mtx[i][j];

	double result = exp(-(F1(i,j)/(k*T)));
	
	double sigma_1=0,sigma_2=0;
	
	int h,l;
	// Eq.(7): Calculate Sigma 1&2 
	for(l=j; l>i; l--){  
		for(h=l-1; h>=i; h--){
			if(basePair(h,l)==true){
				sigma_1 += exp(-(F2(i,j,h,l)/(k*T)));
				
				sigma_2 += Q_m(i+1,h-1) * Q_b(h,l) * exp( -((a+b+c*(j-l-1))/(k*T)) );

				result += (sigma_1 + sigma_2);

			}
		}
	} // EndOf: Eq.(7): Calculate Sigma 1&2

	// update matrix
	Q_b_Mtx[i][j] = result;
	return result;
}


double Q(int i, int j){ // Eq.(5)

	if(Q_Mtx[i][j]!=0) return Q_Mtx[i][j];
	if(i==j) return 1.0; // Qii=1.0

	// ???
	// if(j-1==1) return 1.0;

	double sigma=1.0;
	int h,l;

	for(l=j; l>i; l--){
		for(h=l-1; h>=i; h--){
			if(basePair(h,l)==true){
				sigma += Q_b(h,l) * Q(i,h-1);

			}
		}
	}
	//update matrix
	Q_Mtx[i][j]=sigma;
	return sigma;
}


void Find_Solution(double Mtx[NMAX][NMAX]){
	// more stuff
}

int main(int argc, char *argv[])
{
	FILE* inFile;
	inFile = fopen("rna.txt","r");	
	fscanf(inFile, "%s", RNA);
  n = strlen(RNA);

	Zeros(Q_Mtx,n);
	Zeros(Q_b_Mtx,n);
	Zeros(Q_m_Mtx,n);

	int i,j,step_length;
	
	for(step_length=1; step_length<n-1; step_length++){
	// starting from the shortest segments (step_length==1)
	// iteratively calculate Q_b(i,j) and Q(i,j) to get Q(1,n)
		for(i=0;i<n;i++){
			j = i + step_length;
			Q(i,j);
		}
	}

	Find_Solution(Q_Mtx); // this function is to back track matrix Q_Mtx to find the solution
	// The detail needs to be clarified
		
	return 0;
}
