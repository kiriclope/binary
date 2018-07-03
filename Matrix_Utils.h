#ifndef __MATRIXUTILS__
#define __MATRIXUTILS__

#include "librairies.h"

using namespace::std ;

///////////////////////////////////////////////////////////////////    

void Create_Path(int nbpop,string &path,int N,double K) {
  
  string mkdirp = "mkdir -p " ;
  path += "Connectivity/"+to_string(nbpop)+"pop/N" + to_string(N) ;

  char cK[10] ;
  sprintf(cK,"%0.0f",K) ;
  string sK = string(cK) ;
  
  path += "/K"+sK; 
  mkdirp += path ; 

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) {
    cout << "error creating directories" << endl ;
  }
  
  cout << "Created directory : " ;
  cout << path << endl ;
}

///////////////////////////////////////////////////////////////////    

void Save_Param(int nbpop,string path,double K,int* &Nk,vector<vector<double> > &c) {

  string fparam = path + "/param.txt" ;
  ofstream param(fparam.c_str(), ios::out | ios::ate) ;
  
  cout << "N " ; 
  for(int i=0;i<nbpop;i++)
    cout << Nk[i] << " " ;
  cout << "K " << K << endl ;
  cout << "Connection Probability " << endl ;
  for(int i=0;i<nbpop;i++) {
    for(int j=0;j<nbpop;j++)
      cout << c[i][j] << " " ;
    cout << endl ;
  }
  
  param << "N " ;
  for(int i=0;i<nbpop;i++)
    param << Nk[i] << " " ;
  param << "K " << K << endl ;
  param << "Connection Probability " << endl ;
  for(int i=0;i<nbpop;i++) {
    for(int j=0;j<nbpop;j++)
      param << c[i][j] << " " ;
    param << endl ;
  }
  param.close() ;
}

///////////////////////////////////////////////////////////////////    

void NumberPres(int nbpop,string path,int* &Nk,vector<vector<vector<int> > > &nbPreSab,bool IF_Nk) {
  
  vector<string> List(4) ;
  List[0] = 'E' ;
  List[1] = 'I' ;
  List[2] = 'S' ;
  List[3] = 'V' ;

  string sCon = path ;

  for(int i=0;i<nbpop;i++) {	
    for(int j=0;j<nbpop;j++) {
      
      if(IF_Nk)
	sCon = path + "/Con"+List[i]+List[j]+"_Nk.txt" ; 
      else
	sCon = path + "/Con"+List[i]+List[j]+".txt" ; 

      // cout << sCon << endl ;
      ofstream fCon(sCon.c_str(), ios::out | ios::ate);

      fCon << nbPreSab[i][j][0] ; 
      for(int k=1;k<Nk[i];k++) 
	fCon << "\t" << nbPreSab[i][j][k] ; 
      fCon << endl ;

      sCon.clear();
      fCon.close();
    }
  }

  List.clear() ;
}

///////////////////////////////////////////////////////////////////    

void CheckPres(int nbpop,int* &Nk,vector<vector<vector<int> > > &nbPreSab) {

  double meanPreS=0 ;

  cout << "Average nbPreS : " ;
  for(int i=0;i<nbpop;i++) 
    for(int j=0;j<nbpop;j++) { 
      meanPreS=0 ;
      for(int k=0;k<Nk[i];k++)
	meanPreS += nbPreSab[i][j][k] ;
      cout << meanPreS/(double)Nk[i] << " " ;      
    }  
  cout << endl ;

}
///////////////////////////////////////////////////////////////////    

void WritetoFile(string path,int N,vector<int> &IdPost,vector<int> &nbPost,vector<unsigned long int> &idxPost,bool IF_Nk) {

  cout <<" Writing to Files :" << endl ;

  string nbpath = path ;
  string idxpath = path ;
  string  Idpath = path ;

  if(IF_Nk) {
    nbpath += "/nbPost_Nk.dat" ;
    idxpath += "/idxPost_Nk.dat" ;
    Idpath += "/IdPost_Nk.dat" ;
  }
  else {
    nbpath += "/nbPost.dat" ;
    idxpath += "/idxPost.dat" ;
    Idpath += "/IdPost.dat" ;
  }

  FILE *fIdPost, *fnbPost, *fidxPost ;
  
  fIdPost = fopen(Idpath.c_str(), "wb");
  fwrite(&IdPost[0], sizeof IdPost[0] , IdPost.size() , fIdPost);
  fclose(fIdPost);

  cout << Idpath << endl ;

  for(int i=1;i<N;i++)
    idxPost[i] = idxPost[i-1] + nbPost[i-1] ;

  fidxPost = fopen(idxpath.c_str(), "wb") ;
  fwrite(&idxPost[0], sizeof idxPost[0] , idxPost.size() , fidxPost); 
  fclose(fidxPost);

  cout << idxpath << endl ;

  fnbPost = fopen(nbpath.c_str(), "wb") ;
  fwrite(&nbPost[0], sizeof nbPost[0] , nbPost.size() , fnbPost) ;
  fclose(fnbPost);

  cout << nbpath << endl ;

}

void WriteMatrix(string path,int N,vector<int> &IdPost,vector<int> &nbPost,vector<unsigned long int> &idxPost) {

  cout << "Writing Cij Matrix to : " ;
  string fmatrix= path + "/Cij_Matrix.dat" ;
  cout << fmatrix << endl ;
  
  FILE *Out;
  Out = fopen(fmatrix.c_str(),"wb");
  
  int **M ;
  M = new int*[N] ;
  for(int i=0;i<N;i++) 
    M[i] = new int[N]() ;
  
  for(int i=0;i<N;i++) 
    for(int l=idxPost[i]; l<idxPost[i]+nbPost[i]; l++) 
      M[IdPost[l]][i] = 1 ;
  
    int dum ;
    for (int i=0; i<N; i++) 
      dum = fwrite(M[i], sizeof M[i][0], N, Out) ;
    
    fclose(Out) ;
    delete [] M ;
}

#endif
