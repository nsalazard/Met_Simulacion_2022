// Simular el movimiento de N moleculas en un gas 2D
#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
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
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t), "; // ","
  cout <<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";
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
		if (Molecula[i].Getx() > Lx -R0 ){
		 h = Molecula[i].Getx()-(Lx -R0);
		Faux.load((-1)*K*pow(abs(h),1.5),0,0);
		Molecula[i].AdicioneFuerza(Faux); 
		}
		if (Molecula[i].Getx() < R0 ){
		h = R0- Molecula[i].Getx() ;
		Faux.load(K*pow(abs(h),1.5),0,0);
		Molecula[i].AdicioneFuerza(Faux); 
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
  F=(12*Ep/(d21*d21))*((pow((r_0/(d21)),12))-(pow(r_0/(d21),6)));
	F1 = F*n;
  Molecula2.AdicioneFuerza(F1);   Molecula1.AdicioneFuerza(F1*(-1));
  }   


//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
   cout<<"set terminal gif animate"<<endl; 
  cout<<"set output '5e_k05.gif'"<<endl;
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

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Molecula[N]; //N+4
  Colisionador Hertz;
  Crandom ran64(1);
  double m0=1.0, R0=2.5, kT=0.5, V0=sqrt(2*kT/m0);
  int i,ix,iy;
  double t,tdibujo,dt=1e-4,tmax=200.0,tcuadro=tmax/1000;
  double dx=Lx/(Nx+1), dy=Lx/(Ny+1);
  double Theta, OmegaMax=1.0;
  
  InicieAnimacion(); //Dibujar

  //Inicializar las paredes
  double Rpared=100*Lx, Mpared=100*m0;
  //------------------(  x0,       y0,Vx0,Vy0,theta0,omega0,m0,R0)
 /* Molecula[N+0].Inicie(Lx/2,Ly+Rpared,  0,  0,     0,     0,Mpared,Rpared); //Pared de arriba
  Molecula[N+1].Inicie(Lx/2,  -Rpared,  0,  0,     0,     0,Mpared,Rpared); //Pared de abajo
  Molecula[N+2].Inicie(Lx+Rpared,Ly/2,  0,  0,     0,     0,Mpared,Rpared); //Pared derecha
  Molecula[N+3].Inicie(  -Rpared,Ly/2,  0,  0,     0,     0,Mpared,Rpared); //Pared izquierda*/
  //Inicializar las mol??culas
  for(ix=0;ix<Nx;ix++)
    for(iy=0;iy<Ny;iy++){
      Theta=2*M_PI*ran64.r();
      //-----------------------(x0,y0,Vx0,Vy0, theta0,omega0  ,m0,R0)
      Molecula[Nx*iy+ix].Inicie((ix+1)*dx,(iy+1)*dy, V0,0,Theta,0,m0,R0);//OJO
    }
  for(t=0,tdibujo=0 ; t<tmax ; t+=dt,tdibujo+=dt){
    //Dibujar
    if(tdibujo>tcuadro){
      
      InicieCuadro();
      for(i=0;i<N;i++) Molecula[i].Dibujese();
      TermineCuadro();
      
      tdibujo=0;
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
  return 0;
}
