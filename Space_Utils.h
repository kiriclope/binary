#ifndef __RINGUTILS__
#define __RINGUTILS__

#include "librairies.h"
#include "Net_Utils.h"

#define BL .25
#define m0 .25

///////////////////////////////////////////////////////////////////////

double Gaussian_1D(double mu, double sigma) {
  return exp(-mu*mu/2./sigma/sigma)/sqrt(2.*M_PI)/sigma ;
}

///////////////////////////////////////////////////////////////////////

double Wrapped_Gaussian(double mu, double sigma, int klim) {
  double sum = 0 ; 
  for(int k=-klim;k<=klim;k++)
    sum += Gaussian_1D(mu+M_PI*(double)k,sigma) ;
  return sum ;
}

///////////////////////////////////////////////////////////////////////

void CreateDir_Ring(int nbpop,string &path,int N,double Crec) {
  
  string mkdirp = "mkdir -p " ;
 
  char cCrec[10] ;
  sprintf(cCrec,"%0.2f",Crec) ;
  string sCrec = string(cCrec) ;
  
  path += "/Ring/Crec"+sCrec ; 
  mkdirp += path ; 

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) {
    cout << "error creating directories" << endl ;
  }
  
  cout << "Created directory : " ;
  cout << path << endl ;
}

///////////////////////////////////////////////////////////////////////

void CreateDir_Ring_Ca(int nbpop,string &path,int N,double *Crec) {
  
  string mkdirp = "mkdir -p " ;
  
  string popList[4] = {"E","I","S","V"} ;  
  if(nbpop==1) 
    popList[0] = "I" ;
  string strpop ;
  
  char cCrec[10] ;
  
  path += "/Ring/" ;
  for(int i=0;i<nbpop;i++) {
    strpop = popList[i] ;
    sprintf(cCrec,"%0.2f",Crec[i]) ; 
    path += "Crec"+ strpop + string(cCrec) ; 
  }
  mkdirp += path ; 

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) {
    cout << "error creating directories" << endl ;
  }
  
  cout << "Created directory : " ;
  cout << path << endl ;
}

void CreateDir_Ring_Cff(int nbpop,string &path,int N,double Cff) {
  
  string mkdirp = "mkdir -p " ;
 
  char cCff[10] ;
  sprintf(cCff,"%0.2f",Cff) ;
  string sCff = string(cCff) ;
  
  path += "Cff"+sCff ; 
  mkdirp += path ; 

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) {
    cout << "error creating directories" << endl ;
  }
  
  cout << "Created directory : " ;
  cout << path << endl ;
}

///////////////////////////////////////////////////////////////////////

void Ring_Connection_Probability(int nbpop,int N,int* Nk,vector<int> Cpt,double K,double Crec,vector<vector<double> > &c ) {
  
  for(int i=0;i<nbpop;i++)
    for(int k=0;k<Nk[i];k++)
      for(int j=0;j<nbpop;j++) 
  	for(int l=0;l<Nk[j];l++) 
  	  c[k+Cpt[i]][l+Cpt[j]] = K / ( (double) Nk[j] ) *(1.+2.*Crec*cos(2.*(phi(k,Nk[i])-phi(l,Nk[j]) ) ) ) ;
}

///////////////////////////////////////////////////////////////////////

void Ring_Connection_Probability_Ka(int nbpop,int N,int* Nk,vector<int> Cpt,double *K,double Crec,vector<vector<double> > &c ) {
  
  for(int i=0;i<nbpop;i++)
    for(int k=0;k<Nk[i];k++)
      for(int j=0;j<nbpop;j++) 
  	for(int l=0;l<Nk[j];l++) 
  	  c[k+Cpt[i]][l+Cpt[j]] = K[j] / ( (double) Nk[j] ) *(1.+2.*Crec*cos(2.*(phi(k,Nk[i])-phi(l,Nk[j]) ) ) ) ;
  
}

///////////////////////////////////////////////////////////////////////

void Ring_Connection_Probability_Ca(int nbpop,int N,int* Nk,vector<int> Cpt,double K,double *Crec,vector<vector<double> > &c ) {
  
  for(int i=0;i<nbpop;i++)
    for(int k=0;k<Nk[i];k++)
      for(int j=0;j<nbpop;j++) 
  	for(int l=0;l<Nk[j];l++) 
  	  c[k+Cpt[i]][l+Cpt[j]] = K / ( (double) Nk[j] ) *(1.+2.*Crec[j]*cos(2.*(phi(k,Nk[i])-phi(l,Nk[j]) ) ) ) ;
  
}

///////////////////////////////////////////////////////////////////////

void External_Input(int nbpop,int N,int* Nk,int *Cpt,double K,double Cff,double* Iext,double* IextBL,int ndI,vector<double> &IextFF,double phi0) {

  cout << "Feedforward Input" ;
  
  /* cout << "Standard Tuned FF to E and I" << endl ; */
  /* for(int i=0;i<nbpop;i++) */
  /*   for(int k=0;k<Nk[i];k++)  */
  /* IextFF[k+Cpt[i]] = Iext[i]*( 1. + Cff * cos( 2.*( phi(k,Nk[i]) - phi0 ) ) ) ; */
  /* } */

  /////////////////////////
  //Orthogonal FF Input //
  ///////////////////////
  /* cout << " Orthogonal" << endl ; */
  /* for(int i=0;i<nbpop;i++) */
  /*   for(int k=0;k<Nk[i];k++) */
  /*     if(i==0) */
  /* 	IextFF[k+Cpt[i]] = Iext[i]*(1.+Cff*cos(2.*( phi(k,Nk[i]) - phi0 - M_PI/2. ) ) ) ; */
  /*     else */
  /* 	IextFF[k+Cpt[i]] = Iext[i]*(1.+sqrt(K)*.01*cos(2.*( phi(k,Nk[i]) - phi0 ) ) ) ; */
  
  /////////////////////////
  // Perturbation of I ///
  ///////////////////////
  
  cout << " Perturbation of I" << endl ;
  
  double dI ;
  dI = Iext[ndI] - IextBL[ndI] ;

  cout << "Baseline " ;
  for(int i=0;i<nbpop;i++)
    cout << IextBL[i]/sqrt(K)/m0 << " " ;
  cout << endl ;

  cout << "Perturbation " ;
  for(int i=0;i<nbpop;i++)
    cout << Iext[i]/sqrt(K)/m0 - IextBL[i]/sqrt(K)/m0 << " " ;
  cout << endl ;

  for(int i=0;i<nbpop;i++)
    for(int k=0;k<Nk[i];k++)
      if(i==ndI)
	/* if(i==1 || (i==3 & nbpop==4) ) */
  	IextFF[k+Cpt[i]] = IextBL[i] + dI*( 1. + Cff*cos(2.*( phi(k,Nk[i]) - phi0 ) ) ) ;
      else
  	IextFF[k+Cpt[i]] = IextBL[i] ;  
}

#endif
