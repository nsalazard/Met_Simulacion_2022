#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
#include <algorithm>
#include <fstream>
#include <set>
using namespace std;

//---- declarar constantes ---
const double K=1.0e4;
const double Lx=60, Ly=120;
const int Nx=5, Ny=5, N=Nx*Ny;
const double m = 1.0, Ep = 1.0, r_0 = 10.0, R0=2.5;

const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//--- declarar clases -----
class Cuerpo;
class Colisionador;

//---- interface e implementacion de clases ----
//---- clase cuerpo ---
class Cuerpo{
private:
  vector3D r,V,F; double m,R; double theta,omega,tau,I;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,
	      double theta0,double omega0,double m0,double R0);
  void BorreFuerza(){F.load(0,0,0);};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
	double GetVx(void){return V.x();}; //inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,
		    double theta0,double omega0,double m0,double R0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0;
  theta=theta0; omega=omega0; I=2.0/5*m*R*R;
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt); theta+=omega*(Coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m); omega+=tau*(Coeficiente*dt/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)"; // ","
    //<<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";
}

//--- clase Colisionador ----
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Molecula);
  void CalculeFuerzaEntre(Cuerpo & Molecula1, Cuerpo & Molecula2);
};

void Colisionador::CalculeFuerzas(Cuerpo * Molecula){
   int i,j;
	vector3D Faux;
  //Borrar fuerzas
  for(i=0;i<N;i++){Molecula[i].BorreFuerza();}
  //Calcular las fuerzas entre todas las parejas de Moleculas
  for(i=0;i<N;i++){
		double h;
		// Fuerza entre paredes
		//cout << Molecula[i].Getx() << endl;
		if (Molecula[i].Getx() > Lx -R0 ){
		 h = Molecula[i].Getx()-(Lx -R0);
		Faux.load((-1)*K*pow(abs(h),1.5),0,0);
		Molecula[i].AdicioneFuerza(Faux); 
		
		}
		if (Molecula[i].Getx() < R0 ){
		h = R0- Molecula[i].Getx() ;
		Faux.load(K*pow(abs(h),1.5),0,0);
		Molecula[i].AdicioneFuerza(Faux); 
		//cout << "Fuerza x" <<Molecula[i].F.x() <<"\t"<< "Fuerza Aux " << pow(abs(h),1.5) << endl;
		}
		if (Molecula[i].Gety()  < R0 ){
		h = R0 - Molecula[i].Gety();
		Faux.load(0,K*pow(abs(h),1.5),0);
		Molecula[i].AdicioneFuerza(Faux); 
		}
		if (Molecula[i].Gety()  > (Ly - R0) ){
		h = Molecula[i].Gety() -(Ly - R0);
		Faux.load(0,(-1)*K*pow(abs(h),1.5),0);
		Molecula[i].AdicioneFuerza(Faux); 
		}
		
    for(j=i+1;j<N;j++){CalculeFuerzaEntre(Molecula[i],Molecula[j]);}
		}
	}
	
void Colisionador::CalculeFuerzaEntre(Cuerpo & Molecula1, Cuerpo & Molecula2){
   vector3D r21,n,F1; double d21,F;
  r21=Molecula2.r-Molecula1.r; d21=r21.norm(); n=r21/d21;
  F=(12*Ep/(d21*d21))*((pow((r_0/d21),12))-(pow(r_0/d21,6)));
	F1 = F*n;
  Molecula2.AdicioneFuerza(F1);   Molecula1.AdicioneFuerza(F1*(-1));
  }   


//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
   cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'e.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly+10<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
    cout<<endl;
}
// -------- 
void Grafica(void){
  cout<<"set term pdf"<<endl; 
  cout<<"set out 'HistogramaVelocidad.pdf'"<<endl;
	cout<<"set title 'Histograma de Vx'"<<endl;
  cout<<"set ylabel 'Numero de particulas'"<<endl;
	cout<<"set xlabel 'Vx'"<<endl;
	cout<<"set autoscale"<<endl;
	cout<<"set key"<<endl;
	cout<<"set font ',7'"<<endl;
	cout<<"plot 'HistoV.txt' u 1:2 with impulses"<<endl;
 
}

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Molecula[N]; //N+4
  Colisionador Hertz;
  Crandom ran64(1);
  double m0=1.0, R0=2.5, kT=10, V0=sqrt(2*kT/m0);
  int i,ix, iy;
  double t,tdibujo,dt=1e-3,tmax=500.0,tcuadro=tmax/200;
  double dx=Lx/(Nx+1), dy=Lx/(Ny+1);
  double Theta, OmegaMax=1.0;
	int ti = 0;
	double Vel = 0.0;
	set<double> svel;
  //Crear doc para los datos
	ofstream FVel;
  FVel.open ("FVel.txt");
  
  for(ix=0;ix<Nx;ix++)
    for(iy=0;iy<Ny;iy++){
      Theta=2*M_PI*ran64.r();
      //-----------------------(x0,y0,Vx0,Vy0, theta0,omega0  ,m0,R0)
      Molecula[Nx*iy+ix].Inicie((ix+1)*dx,(iy+1)*dy, V0,0,Theta,0,m0,R0);//OJO
    }
		for(t=0,tdibujo=0 ; t<tmax ; t+=dt,tdibujo+=dt){
    //Dibujar
    if(tdibujo>tcuadro){
			if(ti >=80){
      	for(i=0;i<N;i++){ 		//Guardar Velocidades
					svel.insert(Molecula[i].GetVx());
				}
			}
      tdibujo=0; ti+=1;
    }
    //--- Muevase por PEFRL ---
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,epsilon);
    Hertz.CalculeFuerzas(Molecula);
    for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Molecula);
    for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,chiepsilon);
    Hertz.CalculeFuerzas(Molecula);
    for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Molecula);
    for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,epsilon);  

  }    
	// Calcular la desviaci??n estandar
  double vprom,v2prom,sigma_v;
  //Calculo xprom y x2prom
	for(auto it=svel.begin(); it != svel.end(); it++){
		vprom+=*it;
    v2prom+=pow(*it,2);
  }
	vprom/=svel.size();
	v2prom/=svel.size();
  //Calculo sigma_v
  sigma_v=sqrt(v2prom-vprom*vprom);
  //Imprimir sigma_v
  FVel<<"sigma_v="<<sigma_v<<endl;
 
	for(auto it=svel.begin(); it != svel.end(); it++){FVel << *it << "\n";}
	
	FVel.close();
	//-------- Histograma-----------
	ofstream HistoV;
  HistoV.open ("HistoV.txt");
	
	double sigma = sigma_v/10;
	auto min=svel.begin(), max = svel.end();
	double delta = (*max - *min)/sigma;
	int countv;
	double vmin, vmax;

	for(int ii=0; ii <= delta*2;ii++){
		countv = 0;
		vmin = *min+(sigma*ii); vmax = *min+(sigma*(ii+1));
		for(auto it=svel.begin(); it != svel.end(); it++){
			if(*it >= vmin && *it < vmax){countv += 1;}
				}
		HistoV << (vmin+vmax)/2 << "\t" << countv << endl;
	}
	HistoV.close();
	Grafica();
  return 0;
}
