#include "highPrec.h"

using namespace std;

//mass and width of phi resonance 
const double Mphi=1.0; //GeV
const double W_phi0=0.005; //GeV
//mass and width of phi_prime resonanace
const double Mphi_p= 1.5; //GeV
const double W_phi_p0=0.005; //GeV
//mass of Ds meson
double Mds= 1.967; //GeV
//photon 3-momentum
const double ph_3mom= 0.4; //GeV
const double gamma_phi= 1.0/(sqrt( 1.0 - pow(ph_3mom,2)/(pow(Mphi,2) + pow(ph_3mom,2))));
const double gamma_phi_p= 1.0/(sqrt( 1.0 - pow(ph_3mom,2)/(pow(Mphi_p,2) + pow(ph_3mom,2))));
//width of phi nad phi_prime in motion
double W_phi= W_phi0/gamma_phi;
double W_phi_p= W_phi_p0/gamma_phi_p;
//temporal extent of the lattice [lattice units]
int Th= 200; 
//set lattice spacing [very fine lattice ]
double a=0.04/0.197327; // inverse GeV



//params to plot spec dens and FF
double step_size_erg=0.001; //1 MeV
double Emax=40; //40 GeV
int Nsteps= (int)(Emax/step_size_erg);
double xk_max=1.0;
double xk_step=0.01;
int Nsteps_xk=(int)(xk_max/xk_step);

int prec=256;



PrecFloat Get_exact_gauss(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0) {
  PrecFloat e = exp( -0.5*(E-m)*(E-m)/(s*s));
  PrecFloat norm= s*( 2 + 0.0*erf( (m-E0)/(s*sqrt(PrecFloat(2)))))*sqrt(precPi()/PrecFloat(2)) ;
  return e/norm;
}


PrecFloat BaseFunc(const PrecFloat& E, int t, int T) { return exp(-E*t) + exp( -E*T+ E*t);  }

PrecFloat aE0(const PrecFloat &E0,const PrecFloat &t) {
   return exp(-E0*t)/t; 
}

void Get_Atr(PrecMatr& Atr, const PrecFloat &E0, int T, int tmin, int tmax)  {

  Atr.resize(tmax-tmin+1, tmax-tmin+1);

  for(int t=tmin;t<= tmax; t++)
    for(int r=tmin; r<= tmax; r++) Atr(t-tmin,r-tmin) = aE0(E0, PrecFloat(t+r)) + aE0(E0, PrecFloat(T - t +r)) + aE0(E0,PrecFloat(T+t -r)) + aE0(E0,PrecFloat(2*T -t -r));
    
  
  return;
}


void Get_ft(PrecVect& ft, const PrecFloat &E0,  const PrecFloat &m, const PrecFloat &s,  int T, int tmin, int tmax, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&, const PrecFloat &)> &f) {

  ft.resize(tmax-tmin+1);

 
  for(int t=tmin;t<=tmax;t++) {

    cout<<"In Get_ft t: "<<t<<endl;
    const auto ftT=
      [&f, &m,&s,&E0, &t, &T](const PrecFloat& x) -> PrecFloat
      {
    
	return f(x,m,s, E0)*(exp(-x*t) + exp(-x*(T-t))) ;
      };

    ft(t-tmin) =   integrateUpToInfinite(ftT,E0.get());
  }
  

  return;
}

void Get_bt(PrecVect& bt,const PrecFloat &E,  int T, int tmin, int tmax) {

  bt.resize(tmax-tmin+1);

  for(int t=tmin;t<=tmax;t++) bt(t-tmin) = BaseFunc(E, t, T);

  return;


}








int main(int narg, char** argv) {

  if(narg != 1) exit(-1);

  PrecFloat::setDefaultPrecision(prec);

  //HLT parameters
  const PrecFloat E0= 0.00; //00MeV
  const PrecFloat sigma=0.01; //10MeV

  //test
  cout<<"Using precision: "<<prec<<endl<<flush;

  /*
  PrecVect ft_RE_test;
  auto f_RE_func= [](PrecFloat E, PrecFloat Eg, PrecFloat s, PrecFloat E0) { return (E-Eg)/( pow(E-Eg,2) + s*s);};
  PrecFloat Eg_prec_test= sqrt( pow(PrecFloat(ph_3mom),2) + pow(PrecFloat(Mds)*PrecFloat(0.1),2));
  cout<<"Testing with Eg= "<<Eg_prec_test.get()<<" sigma: "<<sigma.get()<<endl;
  Get_ft(ft_RE_test, E0, Eg_prec_test, sigma, 2*Th, 1, Th, f_RE_func);
  cout<<"test RE part passed! "<<endl;
  */

  




  
 
  auto rho = [](double E) {

	       double Erg_phi= sqrt( Mphi*Mphi + ph_3mom*ph_3mom);
	       double Erg_phi_p= sqrt( Mphi_p*Mphi_p + ph_3mom*ph_3mom);

	       double func=  (W_phi/M_PI)*(1.0/( pow((Erg_phi-E),2) + pow(W_phi,2))) + (W_phi_p/M_PI)*(1.0/( pow((Erg_phi_p-E),2) + pow(W_phi_p,2)));

	       return func;

	     };

 


  auto RE_FF = [&rho](double Eg) {

		 double Erg_phi= sqrt( Mphi*Mphi + ph_3mom*ph_3mom);
		 double Erg_phi_p= sqrt( Mphi_p*Mphi_p + ph_3mom*ph_3mom);

		 
		 double RE_FF_phi= (1.0/(2.0*M_PI))*(Erg_phi-Eg)*( M_PI + 2.0*atan2(Erg_phi,W_phi) + W_phi*log( (pow(W_phi,2) + pow(Erg_phi,2))/pow(Eg,2)))/(pow(W_phi,2) + pow( Erg_phi-Eg,2));

		 double RE_FF_phi_p= (1.0/(2.0*M_PI))*(Erg_phi_p-Eg)*( M_PI + 2.0*atan2(Erg_phi_p,W_phi_p) + W_phi_p*log( (pow(W_phi_p,2) + pow(Erg_phi_p,2))/pow(Eg,2)))/(pow(W_phi_p,2) + pow( Erg_phi_p-Eg,2));

		 return RE_FF_phi+ RE_FF_phi_p;

	       };



  auto IM_FF = [&rho](double Eg) {

		 return -M_PI*rho(Eg);

	       };



  cout<<"Computing correlator..."<<endl;

  vector<PrecFloat> V(Th-1,0.0);

  #pragma omp parallel for schedule(dynamic)
  for(int t=1;t<Th;t++) {

    cout<<"t: "<<t<<endl;

  
    auto Int_corr_prec = [&t](PrecFloat E) -> PrecFloat {

			 PrecFloat W_phi_prec= PrecFloat(W_phi);
			 PrecFloat W_phi_p_prec= PrecFloat(W_phi_p);
			 PrecFloat Erg_phi= sqrt( PrecFloat(Mphi)*PrecFloat(Mphi) + PrecFloat(ph_3mom)*PrecFloat(ph_3mom));
			 PrecFloat Erg_phi_p= sqrt( PrecFloat(Mphi_p)*PrecFloat(Mphi_p) + PrecFloat(ph_3mom)*PrecFloat(ph_3mom));

			 PrecFloat func= ( exp(-t*E*PrecFloat(a)) + exp( -(2*Th-t)*E*PrecFloat(a)))*( (W_phi_prec/precPi())*(1.0/( pow((Erg_phi-E),2) + pow(W_phi_prec,2))) + (W_phi_p_prec/precPi())*(1.0/( pow((Erg_phi_p-E),2) + pow(W_phi_p_prec,2))));

			 return func;
			  };

    auto Corr_func = [](int t) -> PrecFloat {

		       PrecFloat W_phi_prec= PrecFloat(W_phi);
		       PrecFloat W_phi_p_prec= PrecFloat(W_phi_p);
		       PrecFloat Erg_phi= sqrt( PrecFloat(Mphi)*PrecFloat(Mphi) + PrecFloat(ph_3mom)*PrecFloat(ph_3mom));
		       PrecFloat Erg_phi_p= sqrt( PrecFloat(Mphi_p)*PrecFloat(Mphi_p) + PrecFloat(ph_3mom)*PrecFloat(ph_3mom));

		       PrecFloat PH_1= W_phi_prec/Erg_phi;
		       PrecFloat PH_2 = W_phi_p_prec/Erg_phi_p;

		       PrecFloat MOD_1= sqrt( pow(W_phi_prec,2) + pow(Erg_phi,2))*PrecFloat(t)*PrecFloat(a);
		       PrecFloat MOD_2= sqrt( pow(W_phi_p_prec,2) + pow(Erg_phi_p,2))*PrecFloat(t)*PrecFloat(a);

		       PrecFloat Exp_1 = exp(-Erg_phi*PrecFloat(t)*PrecFloat(a))*( cos( W_phi_prec*PrecFloat(t)*PrecFloat(a)) + (1.0/precPi())*ExpEiComplexSum(MOD_1, PH_1, W_phi*PrecFloat(t)*PrecFloat(a), 0)  );
		       PrecFloat Exp_2 = exp(-Erg_phi_p*PrecFloat(t)*PrecFloat(a))*( cos( W_phi_p_prec*PrecFloat(t)*PrecFloat(a)) + (1.0/precPi())*ExpEiComplexSum(MOD_2, PH_2, W_phi_p*PrecFloat(t)*PrecFloat(a), 0)   );

		       PrecFloat MOD_3= sqrt( pow(W_phi_prec,2) + pow(Erg_phi,2))*PrecFloat((2*Th-t))*PrecFloat(a);
		       PrecFloat MOD_4= sqrt( pow(W_phi_p_prec,2) + pow(Erg_phi_p,2))*PrecFloat((2*Th-t))*PrecFloat(a);
		       

		       PrecFloat Exp_3 = exp(-Erg_phi*PrecFloat((2*Th-t))*PrecFloat(a))*( cos( W_phi_prec*(2*PrecFloat(Th-t))*PrecFloat(a)) + (1.0/precPi())*ExpEiComplexSum(MOD_3, PH_1, W_phi*(2*Th-t)*PrecFloat(a), 0)  );
		       PrecFloat Exp_4 = exp(-Erg_phi_p*PrecFloat((2*Th-t))*PrecFloat(a))*( cos( W_phi_p_prec*(2*PrecFloat(Th-t))*PrecFloat(a)) + (1.0/precPi())*ExpEiComplexSum(MOD_4, PH_2, W_phi_p*(2*Th-t)*PrecFloat(a), 0)   );

		       return Exp_1 + Exp_2 + Exp_3 + Exp_4;
		      

		     };

    V[t-1] = integrateUpToInfinite(Int_corr_prec, 0.0);

    cout.precision(  (int)(prec/3.32192809489));

    cout<<"t: "<<t<<" Corr[t]_int : "<<V[t-1]<<endl;
    cout<<"t: "<<t<<" Corr[t]_fun : "<<Corr_func(t)<<endl;
    
          
  }

  cout<<"Correlator computed!"<<endl;

  //Print vector correlator

  ofstream Print_corr("corr.dat");

  for(int t=1;t<Th;t++) Print_corr<<t<<" "<<V[t-1].get()<<endl;

  Print_corr.close();



  //Print spectral density

  ofstream Print_density("dens.dat");



  for(int e=0;e<Nsteps;e++) {

    double Erg= e*step_size_erg;

    Print_density<<Erg<<" "<<rho(Erg)<<endl;

  }

  Print_density.close();


  

  //Print real and immaginary part of the form factors and the spectral density

  ofstream Print_FF("FF.dat");

  PrecMatr Atr;
  Get_Atr(Atr, E0, 2*Th, 1, Th);
  const PrecMatr Atr_inv = Atr.inverse();
  
  #pragma omp parallel for schedule(dynamic)
  for(int ik=0;ik<Nsteps_xk;ik++) {
    
    double xk= xk_step*ik;
    double Eg= sqrt( ph_3mom*ph_3mom + pow(Mds*xk,2));

    cout<<"Inverting offshellness xk: "<<xk<<endl;

    PrecFloat Eg_prec= sqrt( pow(PrecFloat(ph_3mom),2) + pow(PrecFloat(Mds)*PrecFloat(xk),2));
    

    auto f_RE= [](PrecFloat E, PrecFloat Eg, PrecFloat s, PrecFloat E0) { return (E-Eg)/( pow(E-Eg,2) + s*s);};
    auto f_IM= [](PrecFloat E, PrecFloat Eg, PrecFloat s, PrecFloat E0) { return -precPi()*Get_exact_gauss(E, Eg, s, E0);};
    
    PrecVect ft_RE, ft_IM;
    Get_ft(ft_IM, E0, Eg_prec, sigma, 2*Th, 1, Th, f_IM);
    cout<<"fIM computed"<<endl;
    Get_ft(ft_RE, E0, Eg_prec, sigma, 2*Th, 1, Th, f_RE);
    cout<<"fRE computed"<<endl;
   

    PrecVect g_RE = Atr_inv*ft_RE;
    PrecVect g_IM = Atr_inv*ft_IM;
    PrecFloat FF_pred_RE=0.0;
    PrecFloat FF_pred_IM=0.0;
    for(int t=1;t<Th;t++) {
      FF_pred_RE += a*V[t-1]*g_RE(t-1);
      FF_pred_IM += a*V[t-1]*g_IM(t-1);
    }
    

    Print_FF<<xk<<" "<<Eg<<" "<<RE_FF(Eg)<<" "<<IM_FF(Eg)<<" "<<FF_pred_RE.get()<<" "<<FF_pred_IM.get()<<endl<<flush;
  

  }

  Print_FF.close();
  

  

  return 0;
}
