#include <iostream>
#include <cmath>
#include "vector.h"
#include <fstream>
using namespace std;

//Constantes globales
const int N=3;
const double g=980;
const double L = 12.0;

//constantes de PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

//Declaración de las clases
class Cuerpo;
class Colisionador;

//---------- Clase Cuerpo --------------
class Cuerpo{
private:
	double Theta, Omega,Tau;  double m,R, I, L, x0;
public:
  void Inicie(double Theta0,double Omega0,double m0, double R0, double L0, double x00);
  void BorreTorque(void){Tau=0;};
  void SumeTorque(double Tau0){Tau+=Tau0;};
  void Mueva_Theta(double dt,double coeficiente);
  void Mueva_Omega(double dt,double coeficiente);
  void Dibujese(void);
  double GetX(void){return x0 +L*sin(Theta);};
  double GetY(void){return -L*cos(Theta);};
  double GetTheta(void){return Theta;};
  double GetTau(void){return Tau;}; //Inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double Theta0,double Omega0,double m0, double R0, double L0, double x00){
  Theta = Theta0; Omega = Omega0; m=m0; R= R0; L= L0; I =m*L*L; x0 = x00;
}
void Cuerpo::Mueva_Theta(double dt,double coeficiente){
  Theta += Omega*( coeficiente*dt);
}
void Cuerpo::Mueva_Omega(double dt,double coeficiente){
  Omega += Tau*( coeficiente*dt/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<GetX()<<"+"<<R<<"*cos(t),"<<GetY()<<"+"<<R<<"*sin(t)";
	 cout<<" , "<<x0<<"+"<<L/7<<"*t*sin("<<Theta<<"),-"<<L/7<<"*t*cos("<<Theta<<")";
}
//---------- Clase Colisionador --------------
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Pendulo, double K);
  void CalculeFuerzaEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2, double K);    
};
void Colisionador::CalculeFuerzas(Cuerpo * Pendulo, double K){
  int i,j;
  //Borrar todas las fuerzas
  for(i=0;i<N;i++)
    Pendulo[i].BorreTorque();
  //Adicionar fuerzas individuales (que no son de interaccion)
  for(i=0;i<N;i++)
    Pendulo[i].SumeTorque(-Pendulo[i].L*Pendulo[i].m*g*sin(Pendulo[i].Theta));
  //Calcular las fuerzas entre todas las parejas de Pendulos
  for(i=0;i<N-1;i++)
    for(j=i+1;j<N;j++)
  	CalculeFuerzaEntre(Pendulo[i],Pendulo[j], K);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2, double K){
  double S = (Pendulo1.GetX() + Pendulo1.R) -(Pendulo2.GetX() - Pendulo2.R);
  if(S>0){
   double F = K* pow(S, 1.5);
   double T2= F*Pendulo2.L;
   Pendulo1.SumeTorque(-T2); Pendulo2.SumeTorque(T2);
   }
}

//----------- Funciones Globales -----------

void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'CunadeNewton.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-14:34]"<<endl;
  cout<<"set yrange[-18:0]"<<endl;
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

  int main(){  
  double Kc[7]= {0.1e9, 0.2e9, 0.5e9, 1e9, 2e9, 5e9, 10e9};
  string name[7] = {"01","02","05","1","2","5","10"};
  for(int ii = 0; ii <7; ii++){
  double K = Kc[ii];
  ofstream data;
  data.open ("D"+ name[ii] +".txt");
		
  Cuerpo Pendulo[N];
  Colisionador Newton;
  double m0=100, R0=1.5, L0=12;
  double T=2*M_PI*sqrt(L0/g);
  double t,tdibujo,tcuadro=T/50;
  double tmax=0.5*T,dt=0.0001;
  int i;
  //---------------(Theta0,Omega0,m0,R0,L0,x00)
  Pendulo[0].Inicie(-M_PI/12,0,m0,R0,L0,R0);
	for(i=1; i <N; i++)
  Pendulo[i].Inicie(0,0,m0,R0,L0,(2*i+1)*R0);
 
  //InicieAnimacion();
  for(t=0,tdibujo=0; t<tmax; t+=dt,tdibujo+=dt){
    //Dibujar
    /*if(tdibujo>tcuadro){
      InicieCuadro();
      for(i=0;i<N;i++) 	Pendulo[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }       */
    //Imprimir Tau pendulo del medio
    data<<t<<" "<<Pendulo[1].GetTau()<<endl;
    // Mover por PEFRL
    for(i=0;i<N;i++) Pendulo[i].Mueva_Theta(dt,Zeta);
    Newton.CalculeFuerzas(Pendulo,  K);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Omega(dt,Coeficiente1);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Theta(dt,Chi);
    Newton.CalculeFuerzas(Pendulo,  K);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Omega(dt,Lambda);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Theta(dt,Coeficiente2);
    Newton.CalculeFuerzas(Pendulo, K);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Omega(dt,Lambda);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Theta(dt,Chi);
    Newton.CalculeFuerzas(Pendulo, K);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Omega(dt,Coeficiente1);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Theta(dt,Zeta);   
  }
   data.close();
  }
  
  return 0;
}
