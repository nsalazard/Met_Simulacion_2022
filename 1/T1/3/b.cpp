#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

//-----Constantes globales-------------
const int N=2;
const double G=1.0;

//-----Constantes de PEFRL--------------
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

//------Declaración de las clases-------
class Cuerpo;
class Colisionador;

//---------- Clase Cuerpo --------------
class Cuerpo{
private:
  vector3D r,V,F;  double m,R;
public:
  void Inicie(double x0,double y0,
	      double Vx0,double Vy0,double m0,double R0);
  void BorreFuerza(void){F.load(0,0,0);};
  void SumeFuerza(vector3D F0){F+=F0;};
  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  double Getx(void){return r.x();};  
	double Getxrot(double Omega, double t);//inline
  double Gety(void){return r.y();}; //inline
	double Getyrot(double omega, double t); //inline
	void DibujeseRot(double omega, double t);
	void Dibujese(void);
	void PrintRot(double omega, double t);

  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,
	      double Vx0,double Vy0,double m0,double R0){
  r.load(x0,y0,0);  V.load(Vx0,Vy0,0); m=m0; R=R0;
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(dt*coeficiente);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(dt*coeficiente/m);
}

void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}


//---------- Punto B, coordenadas rotadas -------------

double Cuerpo::Getxrot(double omega, double t){
  return r.x()*std::cos(omega*t)+r.y()*std::sin(omega*t);
}

double Cuerpo::Getyrot(double omega, double t){
  return -r.x()*std::sin(omega*t)+r.y()*std::cos(omega*t);
}

void Cuerpo::DibujeseRot(double omega, double t){
  cout<<" , "<<Getxrot(omega,t)<<"+"<<R<<"*cos(t),"<<Getyrot(omega,t)<<"+"<<R<<"*sin(t)";
}

void Cuerpo::PrintRot(double omega, double t){
  cout <<Getxrot(omega,t)<<","<<Getyrot(omega,t)<<",";
}
//---------- Clase Colisionador --------------
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Planeta);
  void CalculeFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2);    
};
void Colisionador::CalculeFuerzas(Cuerpo * Planeta){
  int i,j;
  //Borrar fuerzas
  for(i=0;i<N;i++)
    Planeta[i].BorreFuerza();
  //Calcular las fuerzas entre todas las parejas de planetas
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++)
      CalculeFuerzaEntre(Planeta[i],Planeta[j]);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2){
  vector3D r21,n,F1; double d21,F;
  r21=Planeta2.r-Planeta1.r; d21=r21.norm(); n=r21/d21;
  F=G*Planeta1.m*Planeta2.m*pow(d21,-2.0);
  F1=F*n; Planeta1.SumeFuerza(F1); Planeta2.SumeFuerza(F1*(-1));
}

//-------------------------- Funciones de Animacion -------------------

void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'Nb.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-1100:1100]"<<endl;
  cout<<"set yrange[-1100:1100]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
  cout<<"plot 0,0 ";
}
void TermineCuadro(void){
  cout<<endl;
}

//-----------------------  Programa Principal ------------------------

int main(void){
  Cuerpo Planeta[N];
  Colisionador Newton;
  double m0=1047, m1=1, r=1000;
  double M=m0+m1, Mu=m0*m1/M;
  double x1=r*m0/M, x0=-r*m1/M;
  double omega=std::sqrt(G*M*std::pow(r,-3.0));
  double V1=omega*x1,V0=omega*x0, T=2*M_PI/omega;
  double R0=100, R1=40;
  double t, tdibujo, tmax=20*T, tcuadro=T/200, dt=0.1;
	int i;
  //---------------(x0, y0,Vx0,   Vy0, m0,R0)
  Planeta[0].Inicie(x0, 0, 0 ,   V0 , m0, R0);
  Planeta[1].Inicie(x1, 0, 0 ,   V1 , m1, R1);
  
  //InicieAnimacion(); 
  //Dibujar animacion
  for(t=0, tdibujo=0; t<tmax; t+=dt, tdibujo+=dt){
    if(tdibujo>tcuadro){   
      //InicieCuadro();
	for(int i=0; i<N; i++) {Planeta[i].PrintRot(omega,t);}
      TermineCuadro();
      tdibujo=0;
    }
    //PEFRL
    for(i=0; i<N; i++)Planeta[i].Mueva_r(dt,Zeta);
    Newton.CalculeFuerzas(Planeta);  for(i=0; i<N; i++)Planeta[i].Mueva_V(dt,Coeficiente1);
    for(i=0; i<N; i++)Planeta[i].Mueva_r(dt,Chi);
    Newton.CalculeFuerzas(Planeta);  for(i=0; i<N; i++)Planeta[i].Mueva_V(dt,Lambda);
    for(i=0; i<N; i++)Planeta[i].Mueva_r(dt,Coeficiente2);
    Newton.CalculeFuerzas(Planeta);  for(i=0; i<N; i++)Planeta[i].Mueva_V(dt,Lambda);
    for(i=0; i<N; i++)Planeta[i].Mueva_r(dt,Chi);
    Newton.CalculeFuerzas(Planeta);  for(i=0; i<N; i++)Planeta[i].Mueva_V(dt,Coeficiente1);
    for(i=0; i<N; i++)Planeta[i].Mueva_r(dt,Zeta);
    
  }   
  return 0;
}
