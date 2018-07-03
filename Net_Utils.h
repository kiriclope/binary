#ifndef __BINUTILS__
#define __BINUTILS__

#include "librairies.h"

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

using namespace:: std ; 

void printProgress (double percentage) {
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}

int check_updates(vector<int> v) {
  int sum = 0 ;
  for(int i=0;i<v.size();i++)
    if(v[i]>0)
      sum ++ ;
  return sum ;
}

double sum(vector<double> v,  int a, int b) {
  long double sum = 0 ;
  for(int i=a;i<b;i++)
    sum += v[i] ;
  return sum ;
}

int sumf(vector<double> v,int a,int b) {
  double sum = 0 ;
  for(int i=a;i<b;i++)
    sum += v[i] ;
  return sum ;
}

///////////////////////////////////////////////////////////////////////

double phi(int k, int N) {
  return ( (double) k ) * M_PI / (double) N ;
}

double Heavyside(double x) {
  if(x>0)
    return 1 ;
  else
    return 0 ;
}

///////////////////////////////////////////////////////////////////////

void Set_Parameters(int argc , char** argv, string &dir, int &nbpop, int &N, double &K, double &duration, double &Tst, double &Tw, double &g, double* &Iext, double* &IextBL) {

  if(argv[1] != NULL) {
    dir = argv[1] ;
    nbpop = (int) atoi(argv[2]) ; // total number of neurons prefactor
    N = (int) atoi(argv[3]) ;
    K = (double) atof(argv[4]) ;
    duration = (double) atof(argv[5]) ;
    Tst = (double) atof(argv[6]) ;
    Tw = (double) atof(argv[7]) ;
    g = (double) atof(argv[8]) ;

    Iext = new double [nbpop]() ;
    IextBL = new double [nbpop]() ;

    string stream_path =  "../Parameters/"+to_string(nbpop)+"pop/"+ dir + "/Param.txt" ;
    cout << stream_path << endl ;
    ifstream stream(stream_path.c_str(), ios::in);
    string line=" ";
    
    getline(stream,line);
    /* cout << "line " << line << endl ; */
    line.erase(0,5) ;
    /* cout << "line " << line << endl ; */
    char *end;
    
    IextBL[0] = strtod(line.c_str(),&end);
    int i=1 ;
    while(i<nbpop) {
      /* cout << i << " " ; */
      IextBL[i] = strtod(end,&end);
        /* cout << Iext[i] << endl; */
      i++ ;
    }
    /* cout << endl; */
    if(argc<9) {
      for(int i=0;i<nbpop;i++)
	Iext[i] = IextBL[i] ;
    }
    else
      for(int i=0;i<nbpop;i++)
	Iext[i] = (double) atof(argv[i+9]) ;
  }
  else {
    cout << "Directory ? " ;
    cin >> dir ;
    cout << "nbpop ? " ;
    cin >> nbpop ;
    cout << "total number of neurons ? " ;
    cin >> N ;
    cout << "K ? " ;
    cin >> K ;
    cout << "Duration (s) ? " ;
    cin >> duration ;
    cout << "time before steady Tst ? " ;
    cin >> Tst ;
    cout << "time window Tw ? " ;
    cin >> Tw ;
    cout << "gain g ? " ;
    cin >> g ;
    cout << "External Inputs ? " ;
    Iext = new double [nbpop] ; 
    for(int i=0;i<nbpop;i++)
      cin >> Iext[i] ;
  }
}

///////////////////////////////////////////////////////////////////////

void Import_Connectivity_Parameters(int nbpop,double** &J,string dir) {

  J = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    J[i] = new double[nbpop] ;
      
  cout << "Reading Connectivity from : " ;
  string Jparam = "../" + to_string(nbpop)+"pop/Parameters/"+ dir +"/Jparam.txt" ;
  cout << Jparam << endl;

  struct stat buffer;   
  if (stat (Jparam.c_str(), &buffer) == 0) {
    FILE *Jfile ;
    Jfile = fopen(Jparam.c_str(),"r");
      
    int dum = 0 ;
    for(int i=0;i<nbpop;i++) 
      for(int j=0;j<nbpop;j++) 
	dum = fscanf(Jfile, "%lf", &J[i][j]) ; 
    
    fclose(Jfile) ;
  }
  else 
    cout << "ERROR ... Jparam.txt not found" << endl ;
}

///////////////////////////////////////////////////////////////////////

void Import_Connectivity_Matrix_Pres(int nbpop,int N,double K,double Crec,int* &nbPreS,unsigned long int* &idxPreS,int* &IdPreS, bool IF_RING ){

  string Jpath = "../" ;

  char cK[10] ;
  sprintf(cK,"%0.0f",K) ;
  string sK = string(cK) ;

  char cCrec[10] ;
  sprintf(cCrec,"%0.2f",Crec) ;
  string sCrec = string(cCrec) ;

  Jpath += "Connectivity/"+to_string(nbpop)+"pop/N"+to_string(N)+"/K"+ sK ;
  if(IF_RING) 
    Jpath += "Ring/Crec" + sCrec ;

  cout << "Generating Connectivity Matrix from :" << endl ; 
 
  string nbpath = Jpath ;
  string idxpath = Jpath ;
  string Idpath = Jpath ;

  if(IF_RING) {
    nbpath += "/nbPreS_ring.dat" ;
    idxpath += "/idxPreS_ring.dat" ;
    Idpath += "/IdPreS_ring.dat" ;
  }
  else {
    nbpath += "/nbPreS.dat" ;
    idxpath += "/idxPreS.dat" ;
    Idpath += "/IdPreS.dat" ;
  }

  N = N*1000 ;
  
  cout << nbpath << endl ;
  cout << idxpath << endl ;
  cout << Idpath << endl ;

  nbPreS = new int [N] ;
  idxPreS = new unsigned long int [N] ;

  // vector<int> nbPreS(N) ;
  // vector<unsigned long int> idxPreS(N) ;
  
  FILE *fnbPreS, *fidxPreS, *fIdPreS ;

  int dum ;
  fnbPreS = fopen(nbpath.c_str(), "rb") ;
  dum = fread(&nbPreS[0], sizeof nbPreS[0], N, fnbPreS);  
  fclose(fnbPreS);

  fidxPreS = fopen(idxpath.c_str(), "rb") ;
  dum = fread(&idxPreS[0], sizeof idxPreS[0], N, fidxPreS);
  fclose(fidxPreS);

  unsigned long int nbPreStot = 0 ;
  for(int j=0 ; j<N; j++)
    nbPreStot += nbPreS[j] ;

  IdPreS = new int [nbPreStot] ;

  fIdPreS = fopen(Idpath.c_str(), "rb");
  dum = fread(&IdPreS[0], sizeof IdPreS[0], nbPreStot, fIdPreS); 
  fclose(fIdPreS);
}

///////////////////////////////////////////////////////////////////////

void Import_Ring_Orientations(int nbpop, int N, double K, vector<double> &phi){

  string Jpath = "../" ;

  char cK[10] ;
  sprintf(cK,"%0.0f",K) ;
  string sK = string(cK) ;

  Jpath += +"Connectivity/"+to_string(nbpop)+"pop/N"+to_string(N/1000)+"/K"+ sK ;

  cout << "Importing Orientations from :" ;

  string str_phi = Jpath + "/ring.dat" ;
  cout << str_phi << endl ;

  struct stat buffer;   
  if (stat (str_phi.c_str(), &buffer) == 0) {
    FILE *fphi ;
    fphi = fopen(str_phi.c_str(), "rb");
    
    int dum = 0 ;
    dum = fread(&phi[0], sizeof phi[0], phi.size(), fphi);  
    
    fclose(fphi);
  }
  else
    cout << "Ring.dat not found ..." << endl ;
}

///////////////////////////////////////////////////////////////////////

void Create_Cij_Matrix(int nbpop,int N,double K,int** &M){

  string Jpath = "../" ;

  char cK[10] ;
  sprintf(cK,"%0.0f",K) ;
  string sK = string(cK) ;

  Jpath += "Connectivity/"+to_string(nbpop)+"pop/N"+to_string(N)+"/K"+ sK ;

  cout << "Generating Connectivity Matrix from :" ;
  cout << Jpath << endl ;
 
  string matrixpath = Jpath ;

  matrixpath += "/Cij_Matrix.dat" ;

  N = N*1000 ;
   
  M = new int*[N] ;
  for(int i=0;i<N;i++) 
    M[i] = new int[N]() ;

  FILE *fmatrix;
  fmatrix = fopen(matrixpath.c_str(),"rb");

  int dum ;
  for (int i=0; i<N; i++) 
    dum = fread(M[i], sizeof M[i][0], N, fmatrix) ;
  
  fclose(fmatrix);
}

///////////////////////////////////////////////////////////////////////

void Create_Cij_Matrix_Ring(int nbpop,int N,double K,int** &M){

  string Jpath = "../" ;

  char cK[10] ;
  sprintf(cK,"%0.0f",K) ;
  string sK = string(cK) ;

  Jpath += "Connectivity/"+to_string(nbpop)+"pop/N"+to_string(N)+"/K"+ sK ;

  cout << "Generating Connectivity Matrix from :" ;
  cout << Jpath << endl ;
 
  string matrixpath = Jpath ;

  matrixpath += "/Cij_Matrix_Ring.dat" ;

  N = N*1000 ;
   
  M = new int*[N] ;
  for(int i=0;i<N;i++) 
    M[i] = new int[N]() ;

  FILE *fmatrix;
  fmatrix = fopen(matrixpath.c_str(),"rb");

  int dum ;
  for (int i=0; i<N; i++) 
    dum = fread(M[i], sizeof M[i][0], N, fmatrix) ;
  
  fclose(fmatrix);
}

///////////////////////////////////////////////////////////////////////

void Import_Connectivity_Matrix_Post(int nbpop,int N,double K,int* &nbPost,unsigned long int* &idxPost,int* &IdPost, double Crec, bool IF_RING,bool IF_Nk) {

  char cK[10] ;
  sprintf(cK,"%0.0f",K) ;
  string sK = string(cK) ;

  char cCrec[10] ;
  sprintf(cCrec,"%0.2f",Crec) ;
  string sCrec = string(cCrec) ;

  string Jpath = "../" ;
  Jpath += "Connectivity/"+to_string(nbpop)+"pop/N"+to_string(N)+"/K"+ sK ;

  if(IF_RING) 
    Jpath += "/Ring/Crec" + sCrec ;

  N = N*10000 ;
  nbPost = new int [N] ;
  idxPost = new unsigned long int [N] ;
  
  cout << "Importing Connectivity from : " << endl ;
 
  string nbpath = Jpath ;
  string idxpath = Jpath ;
  string Idpath = Jpath ;

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

  cout << nbpath << endl ;
  cout << idxpath << endl ;
  cout << Idpath << endl ;

  struct stat buffer;   
  if (stat (nbpath.c_str(), &buffer) == 0) {
    FILE *fnbPost, *fidxPost, *fIdPost ;
    
    fnbPost = fopen(nbpath.c_str(), "rb") ;
    fidxPost = fopen(idxpath.c_str(), "rb") ;
    
    int dum1, dum2 ;
    dum1 = fread(&nbPost[0], sizeof nbPost[0], N , fnbPost);  
    dum2 = fread(&idxPost[0], sizeof idxPost[0], N , fidxPost);
    
    fclose(fnbPost);
    fclose(fidxPost);
    
    unsigned long int nbposttot = 0 ;
    for(int j=0 ; j<N; j++)
      nbposttot += nbPost[j] ;

    IdPost = new int [nbposttot] ;    
    fIdPost = fopen(Idpath.c_str(), "rb");
    int dum3 ;
    dum3 = fread(&IdPost[0], sizeof IdPost[0], nbposttot , fIdPost); 
    fclose(fIdPost);
  }
  else
    cout << "ERROR : nbPost.dat or idxPost.dat or IdPost.dat not found" << endl ;
}

///////////////////////////////////////////////////////////////////////

void Import_Connectivity_Matrix_Post_Ca(int nbpop,int N,double K,int* &nbPost,unsigned long int* &idxPost,int* &IdPost, double *Crec, bool IF_RING, bool IF_Nk) {

  char cK[10] ;
  sprintf(cK,"%0.0f",K) ;
  string sK = string(cK) ;

  char cCrec[10] ;

  string Jpath = "../" ;
  Jpath += "Connectivity/"+to_string(nbpop)+"pop/N"+to_string(N)+"/K"+ sK +"/";

  string popList[4] = {"E","I","S","V"} ;  
  if(nbpop==1) 
    popList[0] = "I" ;
  string strpop ;

  if(IF_RING) 
    Jpath += "Ring/" ;
    for(int i=0;i<nbpop;i++) {
      strpop = popList[i] ;
      sprintf(cCrec,"%0.2f",Crec[i]) ; 
      Jpath += "Crec"+ strpop + string(cCrec) ; 
    }
  
  N = N*10000 ;
  nbPost = new int [N] ;
  idxPost = new unsigned long int [N] ;
  
  cout << "Importing Connectivity from : " << endl ;
 
  string nbpath = Jpath ;
  string idxpath = Jpath ;
  string Idpath = Jpath ;

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

  cout << nbpath << endl ;
  cout << idxpath << endl ;
  cout << Idpath << endl ;

  struct stat buffer;   
  if (stat (nbpath.c_str(), &buffer) == 0) {
    FILE *fnbPost, *fidxPost, *fIdPost ;
    
    fnbPost = fopen(nbpath.c_str(), "rb") ;
    fidxPost = fopen(idxpath.c_str(), "rb") ;
    
    int dum1, dum2 ;
    dum1 = fread(&nbPost[0], sizeof nbPost[0], N , fnbPost);  
    dum2 = fread(&idxPost[0], sizeof idxPost[0], N , fidxPost);
    
    fclose(fnbPost);
    fclose(fidxPost);
    
    unsigned long int nbposttot = 0 ;
    for(int j=0 ; j<N; j++)
      nbposttot += nbPost[j] ;

    IdPost = new int [nbposttot] ;    
    fIdPost = fopen(Idpath.c_str(), "rb");
    int dum3 ;
    dum3 = fread(&IdPost[0], sizeof IdPost[0], nbposttot , fIdPost); 
    fclose(fIdPost);
  }
  else
    cout << "ERROR : nbPost.dat or idxPost.dat or IdPost.dat not found" << endl ;
}

///////////////////////////////////////////////////////////////////////

void Check_Connectivity_Matrix(int nbpop,int N,double K,double Crec,bool IF_RING, bool IF_Nk) {

  string Jpath = "../" ;

  char cK[10] ;
  sprintf(cK,"%0.0f",K) ;
  string sK = string(cK) ;

  char cCrec[10] ;
  sprintf(cCrec,"%0.2f",Crec) ;
  string sCrec = string(cCrec) ;

  Jpath += "Connectivity/"+to_string(nbpop)+"pop/N"+to_string(N)+"/K"+ sK ;

  vector<string> List(4) ;
  List[0] = 'E' ;
  List[1] = 'I' ;
  List[2] = 'S' ;
  List[3] = 'V' ;

  string sCon = Jpath ;
  double nbPres = 0.;
  int tot = 0 ;
  double val;
  struct stat buffer;   

  for(int i=0;i<nbpop;i++) {	
    for(int j=0;j<nbpop;j++) {
      
      if(IF_RING) {
	if(IF_Nk) {
	  sCon = Jpath + "/Crec"+sCrec + "/Con"+List[i]+List[j]+"_ring_Nk.txt" ; 
	}
	else {
	  sCon = Jpath + "/Crec"+sCrec + "/Con"+List[i]+List[j]+"_ring.txt" ; 
	}
      }
      else
	sCon = Jpath + "/Con"+List[i]+List[j]+".txt" ; 
      // cout << sCon << endl ;
      
      if (stat (sCon.c_str(), &buffer) == 0) {
	ifstream fCon(sCon.c_str(), ios::in);
	
	nbPres = 0. ;
	tot = 0 ;
	while(fCon>>val ) {
	  nbPres +=val ;
	  tot += 1 ;
	}
	
	string nbPresij = "nbPres"+List[i]+List[j] ;
	cout << nbPresij << " " << nbPres/(double)tot << " " ;
	
	sCon.clear() ;
	fCon.close();
      }
      else
	cout << " error " << sCon << " not found !" << endl ; 
    }
    cout << endl;
  }
}

///////////////////////////////////////////////////////////////////////

void nbNeurons(int nbpop, int &N, int* &Nk, bool IF_Nk) {
  N = N*10000 ;
  
  Nk = new int [nbpop]() ;
  for(int i=0;i<nbpop;i++)
    Nk[i] = N/nbpop ;

  if(IF_Nk) {
    if(nbpop==2) {
      Nk[0]= int(N*70./100.) ;
      Nk[1]= int(N*30./100.) ;
    }
    if(nbpop==3) {
      Nk[0]= int(N*70./100.) ;
      Nk[1]= int(N*20./100.) ;
      Nk[2]= int(N*10./100.) ;
    }
    if(nbpop==4) {
      Nk[0]= int(N*60./100.) ;
      Nk[1]= int(N*20./100.) ; 
      Nk[2]= int(N*10./100.) ; 
      Nk[3]= int(N*10./100.) ; 
    }
  }

}

///////////////////////////////////////////////////////////////////////

void Save_Parameters(string dir, int nbpop, int N, int *Nk, double K, double** J, double duration, int Tw, double *Iext, string path) {

  string fparam = path + "/Param.txt" ;
  ofstream param(fparam.c_str(), ios::out | ios::ate);

  vector<char> List(4) ;
  List[0] = 'E' ;
  List[1] = 'I' ;
  List[2] = 'S' ;
  List[3] = 'V' ;

  param << "Number of neurons" << endl ;
  param << N ;
  for(int i=0;i<nbpop;i++)
    param <<" "<< "N"<< List[i] <<"="<< Nk[i] ;
  param << endl;

  param << "Number of Inputs K " << endl ;
  param << K << endl ;

  param << "duration=" << duration << " " ;
  param << "Tw=" << Tw << endl ; 

  param << "External Inputs" << endl ;
  for(int i=0;i<nbpop;i++)
    param << List[i]<<"="<< Iext[i] << " " ;
  param << endl ;
 
  param << "Synaptic Strength" << endl ;
  for(int i=0;i<nbpop;i++) {
    param << "//" ;
    for(int j=0;j<nbpop;j++) 
      param << J[i][j] << " " ;
    param << "//" << endl;
  }  
  param << endl;
        
  param.close();
}

///////////////////////////////////////////////////////////////////////

void Create_Dir(string dir, int nbpop, int N, double K, double g, string &path) {

  string mkdirp = "mkdir -p " ;
  path += "Simulations/"+ to_string(nbpop) +"pop/"+ dir + "/N" + to_string(N) ;

  char cK[10] ;
  sprintf(cK,"%0.0f",K) ;
  string sK = string(cK) ;

  char cg[10] ;
  sprintf(cg,"%0.2f",g) ;
  string sg = string(cg) ;

  path += "/K"+ sK + "/g" + sg ;  
  mkdirp += path ;

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) 
    cout << "error creating directories" << endl ;

  cout << "Created directory : " ;
  cout << path << endl ;

}

void CreateDir_Iext(int nPrtr, double I, string &path) {

  char cI[10] ;
  sprintf(cI,"%0.4f",I) ;

  string sI = string(cI) ;

  string popList[4] = {"E","I","S","V"} ;
  
  string spop = popList[nPrtr] ;

  string mkdirp = "mkdir -p " ;
  path = path + "/Prtr_" + spop + "/Iext_" + spop + sI ;
  mkdirp += path ;

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) 
    cout << "error creating directories" << endl ;

  cout << "Created directories : " ; 
  cout << path << endl ;

}

void CreateDir_PHI(double phi0, string &path) {

  char cphi0[10] ;
  sprintf(cphi0,"%0.2f",phi0) ;

  string sphi0 = string(cphi0) ;

  string mkdirp = "mkdir -p " ;
  path = path + "/PHI_" + sphi0 ;
  mkdirp += path ;

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) 
    cout << "error creating directories" << endl ;

  cout << "Created directories : " ; 
  cout << path << endl ;

}

void Create_Dir_PE(int nPrtr, double I, double Imax, double dI, string &path) {

  char cI[10] ;
  sprintf(cI,"%0.4f",I) ;
  char cImax[10] ;
  sprintf(cImax,"%0.4f",Imax) ;
  char cdI[10] ;
  sprintf(cdI,"%0.4f",dI) ;

  string sI = string(cI) ;
  string sImax = string(cImax) ;
  string sdI = string(cdI) ;

  string popList[4] = {"E","I","S","V"} ;
  
  string spop = popList[nPrtr] ;

  string mkdirp = "mkdir -p " ;
  //  path = path + "/PE" + "/Imin" + sI + "Imax" + sImax + "dI" + sdI ;
  path = path + "/Prtr_" + spop + "/Imin" + sI + "Imax" + sImax + "dI" + sdI ;
  mkdirp += path ;

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) 
    cout << "error creating directories" << endl ;

  cout << "Created directories : " ; 
  cout << path << endl ;

}

#endif
