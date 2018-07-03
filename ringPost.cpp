// LAST GOOD VERSION OF THE RING MODEL
#include "Space_Utils.h"
#include "Net_Utils.h"

#define IF_Nk false
#define IF_Iext false
#define IF_PHI false
#define ndI 1

clock_t t1=clock();

int main(int argc , char** argv) {

  ///////////////////////////////////////////////////////////////////
  // parameters 
  ///////////////////////////////////////////////////////////////////
  
  string dir ;
  int nbpop, N ;
  double K,duration,**J,g,*Iext,*IextBL,Tst,Tw,Crec,Cff,phi0;
  
  Set_Parameters(argc,argv,dir,nbpop,N,K,duration,Tst,Tw,g,Iext,IextBL) ;

  Crec = (double) atof(argv[nbpop+9]) ;
  Cff = (double) atof(argv[nbpop+10]) ;

  if(IF_PHI)
    phi0 = (double) atof(argv[nbpop+11]) ;
  else
    phi0 = 0 ;

  ///////////////////////////////////////////////////////////////////    
  
  cout << "External Inputs : " ;
  for(int i=0;i<nbpop;i++) 
    cout << Iext[i] << " ";
  cout << endl ;

  cout << "Synaptic Strength : " << endl ;
  Import_Connectivity_Parameters(nbpop,J,dir) ;    
  for(int i=0;i<nbpop;i++) {
    for(int j=0;j<nbpop;j++) 
      cout << J[i][j] << " ";
    cout << endl ;
  }

  ///////////////////////////////////////////////////////////////////    

  string path = "../" ;
  Create_Dir(dir, nbpop, N, K, g, path) ;
  CreateDir_Ring( nbpop, path, N, Crec) ;
  CreateDir_Ring_Cff(nbpop, path,N, Cff) ;

  if(IF_Iext)
    CreateDir_Iext(1,Iext[1],path) ;
  if(IF_PHI)
    CreateDir_PHI(phi0,path) ;


  ///////////////////////////////////////////////////////////////////    
  

  int *nbPost ;
  unsigned long  int *idxPost ;
  int *IdPost ;

  Import_Connectivity_Matrix_Post(nbpop,N,K,nbPost,idxPost,IdPost,Crec,true,IF_Nk) ;
  // Check_Connectivity_Matrix(nbpop,N,K,Crec,true,IF_Nk) ;

  ///////////////////////////////////////////////////////////////////    

  int* Nk ;
  nbNeurons(nbpop,N,Nk,IF_Nk) ;
  cout << "Nk " ;
  for(int i=0;i<nbpop;i++) 
    cout << Nk[i] << " ";
  cout << endl ;

  ///////////////////////////////////////////////////////////////////    

  Save_Parameters(dir, nbpop, N, Nk, K, J, duration, Tw, Iext, path) ;

  ///////////////////////////////////////////////////////////////////    

  for(int i=0;i<nbpop;i++) {
    Iext[i] = sqrt(K)*Iext[i]*m0 ;
    IextBL[i] = sqrt(K)*IextBL[i]*m0 ;
    for(int j=0;j<nbpop;j++) 
      J[i][j] = g*J[i][j]/sqrt(K) ;
  }
  
  ///////////////////////////////////////////////////////////////////    
  // Offset between the different populations
  ///////////////////////////////////////////////////////////////////    

  int* Cpt ; // compteur Cpt0 = Ne, cpt1 = Ne+Ni ...
  Cpt = new int [nbpop+1]() ;
  
  for(int i=0;i<nbpop+1;i++) { 
    for(int j=0;j<i;j++) 
      Cpt[i] += Nk[j] ; 
    cout <<"Cpt="<< Cpt[i] << " ";
  }
  cout << endl;

  ///////////////////////////////////////////////////////////////////    
  // External Input Selectivity
  ///////////////////////////////////////////////////////////////////    
  
  vector<double> IextFF(N) ;
  External_Input(nbpop,N,Nk,Cpt,K,Cff,Iext,IextBL,ndI,IextFF,phi0) ;
  
  ///////////////////////////////////////////////////////////////////    
  // Dynamics of the network : Binary neurons
  ///////////////////////////////////////////////////////////////////    

  vector<double> S(N) ; //State variable of the neurons
  vector<double> nbspk(N), total_nbspk(N), inputs(N) ; //State variable of the neurons
  vector<double> Isyn(N) ; // Synaptic input to the neurons

  double Inet = 0 ;
  ///////////////////////////////////////////////////////////////////    

  string strPop_Activity = path + "/Mean.txt" ;
  ofstream Pop_Activity(strPop_Activity.c_str(), ios::out | ios::ate);
 
  string strIdv_Activity = path + "/IdvRates.txt" ;
  ofstream Idv_Activity(strIdv_Activity.c_str(), ios::out | ios::ate);

  string strIdv_Inputs = path + "/IdvInputs.txt" ;
  ofstream Idv_Inputs(strIdv_Inputs.c_str(), ios::out | ios::ate);

 ///////////////////////////////////////////////////////////////////    

  cout << "Initialization" << endl;

  //Using C++11
  random_device rd ;
  default_random_engine gen( rd() ) ;

  ///////////////////////////////////////////////////////////////////    

  uniform_real_distribution<double> unif_one(0,1) ;
  uniform_int_distribution<int> unif_pop(0,Nk[0]) ;

  int i=0,j=0,k=0, updpop=0, updneuron=0 ;
  unsigned long int l=0 ;

  vector<unsigned long int> sum_spk(nbpop) ;

  cout << "Check random seed " ;
  for(i=0;i<10;i++)
    cout << unif_one(gen) << " " ;
  cout << endl ;

  ///////////////////////////////////////////////////////////////////    
  //initial conditions
  ///////////////////////////////////////////////////////////////////    

  for(i=0;i<nbpop;i++) {
    
    for(j=0;j<Nk[i];j++)
      if(unif_one(gen)<.5) {
	S[j+Cpt[i]] = 1. ;
	
	for(l=idxPost[j+Cpt[i]]; l<idxPost[j+Cpt[i]]+nbPost[j+Cpt[i]];l++)
	  for(k=0;k<nbpop;k++) 
	    if(IdPost[l]>=Cpt[k] && IdPost[l]<Cpt[k+1]) 
	      Isyn[IdPost[l]] += J[k][i] ;
      }
  }

  for(i=0;i<nbpop;i++)
      cout << " " << sum(S,Cpt[i],Cpt[i+1]) / ( (double) Nk[i]) ;
  cout << endl ;

  ///////////////////////////////////////////////////////////////////    

  vector<double> tau(nbpop) ;
  vector<double> Pb(nbpop) ;
  tau[0] = 3. ;
  tau[1] = 2. ;
  if(nbpop>=3)
    tau[2] = 2. ;
  if(nbpop>=4)
    tau[3] = 2. ;

  double tw = 0 ;
  double dt = 0 ;

  dt = (tau[0] * tau[1]) / ( ( (double) Nk[0] ) * tau[1] + ( (double) Nk[1] ) * tau[0]) ;

  cout << "nbSteps " << duration/dt << endl ;

  for (int i=0;i<nbpop;i++)
    Pb[i] = dt * ((double) Nk[i]) / tau[i] ;

  cout << "Pb " ;
  for (int i=0;i<nbpop;i++)
    cout << Pb[i] << " " ;
  cout << endl ;
  
  ///////////////////////////////////////////////////////////////////    

  double percentage=0,rdn=0 ;
  long double nsteps = duration/dt ;
  cout << "Main loop :" ;
  cout << " duration " << duration << " | dt " << dt ; 
  cout << " | Tst " << Tst << " | Tw " << Tw << endl ;

  for(double t=0;t<duration;t+=dt) {
    
    rdn = unif_one(gen) ;
    if(nbpop=+2)
      if(rdn<Pb[0]) 
	updpop = 0 ;
      else
	updpop = 1 ;
    else {
      if( (unif_one(gen)<.5 && nbpop==3) || (unif_one(gen)<.33 && nbpop==4) )
	updpop = 1 ;
      else
	if(unif_one(gen)<.33 || nbpop==3)
	  updpop = 2 ;
	else
	  updpop = 3 ;
    }

    // if(rdn<Pb[0]) 
    //   // updpop = unif(gen) ;
    //   if(unif_one(gen)<.5)
    // 	updpop = 0 ; // FF ou E1
    //   else
    // 	if(nbpop<4)
    // 	  updpop = 1 ; // E ou I1
    // 	else
    // 	  updpop = 2 ; // I ou E2
    // else
    //   if(unif_one(gen)<.5 || nbpop<4)
    // 	updpop = 1 ; // I 
    //   else
    // 	updpop = 3 ; // S ou I2
    

    unif_pop.param(uniform_int_distribution<int>::param_type(0,Nk[updpop])) ;
   
    updneuron = unif_pop(gen) ;
    updneuron += Cpt[updpop] ;

    Inet = IextFF[updneuron] + Isyn[updneuron] ;
    // Inet = Iext[updpop] + Isyn[updneuron] ;
    
    if(Inet>=1) {
      if(S[updneuron] == 0) {// update only if flip 0 -> 1	
	for(l=idxPost[updneuron]; l<idxPost[updneuron]+nbPost[updneuron];l++) 
	  for(i=0;i<nbpop;i++) 
	    if(IdPost[l]>=Cpt[i] && IdPost[l]<Cpt[i+1]) 
	      Isyn[IdPost[l]] += J[i][updpop] ;
	S[updneuron] = 1 ;
      }
    }
    else {
      if(S[updneuron]==1) {// update only if flip 1 -> 0
	for(l=idxPost[updneuron]; l<idxPost[updneuron]+nbPost[updneuron]; l++) 
	  for(i=0;i<nbpop;i++) 
	    if(IdPost[l]>=Cpt[i] && IdPost[l]<Cpt[i+1]) 
	      Isyn[IdPost[l]] -= J[i][updpop] ;
	S[updneuron] = 0 ;
      }
    }

    percentage = t/duration ;

    if(tw>=Tw) {

      cout << " t " << t-Tst << " Rates";
      for(int i=0;i<nbpop;i++) 
      	cout << " " << sum(nbspk,Cpt[i],Cpt[i+1]) * (dt/Tw) / ( (double) Nk[i]) ;
      // cout << endl ;
      
      cout.flush();

      Pop_Activity << t-Tst ;
      for(int i=0;i<nbpop;i++) 
      	Pop_Activity << ";" << sum(nbspk,Cpt[i],Cpt[i+1]) * (dt/Tw) / ( (double) Nk[i]) ;
      Pop_Activity << endl ;
      
      Idv_Activity << t-Tst ;
      for(int k=0; k<N; k++)
      	Idv_Activity << ";" << ( (long double) nbspk[k] ) * (dt/Tw)  ;
      Idv_Activity << endl ;

      Idv_Inputs << t-Tst ;
      for(int k=0; k<N; k++)
      	Idv_Inputs << ";" << IextFF[k] + ( (long double) inputs[k] ) * (dt/Tw)  ;
      Idv_Inputs << endl ;
      
      tw = 0 ;
      
      fill(nbspk.begin(), nbspk.end(), 0) ;
      fill(inputs.begin(), inputs.end(), 0) ;
    }

    if(t>=Tst) {
      tw += dt ;
      transform (nbspk.begin(), nbspk.end(), S.begin(), nbspk.begin(), plus<double>());
      transform (inputs.begin(), inputs.end(), Isyn.begin(), inputs.begin(), plus<double>());
      transform (total_nbspk.begin(), total_nbspk.end(), S.begin(), total_nbspk.begin(), plus<double>());
    }

    printProgress(percentage) ;

  } // Endloop
  
  ///////////////////////////////////////////////////////////////////

  delete[] nbPost ;
  delete[] IdPost ;
  delete[] idxPost ;

  delete [] J ;
  delete [] Iext ;
  delete [] Nk ;
  delete[] Cpt ;
  
  S.clear() ;
  inputs.clear() ;
  nbspk.clear() ;
  total_nbspk.clear();
  
  Isyn.clear() ;
  
  tau.clear() ; 
  Pb.clear() ;
  
  Pop_Activity.close() ;
  Idv_Activity.close() ;
  Idv_Inputs.close() ;

  // Rates.close();

  ///////////////////////////////////////////////////////////////////    

  cout << "Simulation Done !" << endl ;

  clock_t t2=clock();  
  int HOURS=0,MIN=0,SEC=0;
  string str_TIME = path + "/CPU_TIME.txt" ; 
  ofstream TIME(str_TIME.c_str(), ios::out | ios::ate);

  SEC = (t2-t1)/CLOCKS_PER_SEC ;
  HOURS = SEC/3600 ;
  MIN = SEC/60 ;
  SEC = SEC % 60 ;
  cout << "Elapsed Time = " << HOURS << "h " << MIN << "m " << SEC << "s" << endl;
  TIME << "Elapsed Time = " << HOURS << "h " << MIN << "m " << SEC << "s" << endl;
  TIME.close() ;
  return 0;

}
