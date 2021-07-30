#include <bits/stdc++.h>
using namespace std;

/*
  IMB1
  two-objective
  recommended: n=10;
  PF
   0 <= f1 <= 0.2
   f2= 1- sqrt(f1)
    0 <= xi <= 1
  
*/
void imb1(vector<double> &x, vector<double> &f)
{
  int n = x.size();
  double gx = 0.0;
  if(x[0] > 0.2)
  {
     for(int i = 1; i < n; i++)
     {
	double ti = x[i] - sin(0.5*M_PI*x[0]);
	gx += 0.5*(-0.9*ti*ti + pow(fabs(ti), 0.6)); 
     }
  }
  f[0] = (1.0+gx)*x[0];
  f[1] = (1.0+gx)*(1.0-sqrt(x[0]));
}
/*
  IMB2
  two-objective
  recommended: n=10;
  PF
   0.4 <= f1 <= 0.6
   f2= 1- sqrt(f1)
    0 <= xi <= 1
  PS: x_j = sin(0.5*M_PI*x_1)
  
*/
void imb2(vector<double> &x, vector<double> &f)
{
  int n = x.size();
  double gx = 0.0;
  if( x[0]< 0.4 || x[0]>0.6)
  {
     for(int i = 1; i < n; i++)
     {
	double ti = x[i] - sin(0.5*M_PI*x[0]);
	gx += 0.5*(-0.9*ti*ti + pow(fabs(ti), 0.6)); 
     }
  }
  f[0] = (1.0+gx)*x[0];
  f[1] = (1.0+gx)*(1.0-x[0]);
}
/*
  IMB3
  two-objective
  recommended: n=10;
  PF
   f1^2 + f2^2 = 1
    0 <= xi <= 1
  PS: x_j = sin(0.5*M_PI*x_1)
  
*/
void imb3(vector<double> &x, vector<double> &f)
{
  int n = x.size();
  double gx = 0.0;
  if(x[0]< 0.8)
  {
     for(int i = 1; i < n; i++)
     {
	double ti = x[i] - sin(0.5*M_PI*x[0]);
	gx += 0.5*(-0.9*ti*ti + pow(fabs(ti), 0.6)); 
     }
  }
  f[0] = (1.0+gx)*cos(M_PI*x[0]*0.5);
  f[1] = (1.0+gx)*sin(M_PI*x[0]*0.5);
}
/*
  IMB4
  three-objective
  recommended: n=10;
  PF
   f1^2 + f2^2 + f3^2 = 1
    0 <= xi <= 1
  PS: x_j = (x1+x2)/2
  
*/
void imb4(vector<double> &x, vector<double> &f)
{
  int n = x.size();
  double gx = 0.0;
  if(x[0]<(2.0/3.0))
  {
     for(int i = 2; i < n; i++)
     {
	double ti = x[i] - ((x[0]+x[1])*0.5);
	gx += (-0.9*ti*ti + pow(fabs(ti), 0.6)); 
     }
     gx *= 2.0*cos(M_PI*x[0]*0.5);
  }
  f[0] = (1.0+gx)*x[0]*x[1];
  f[1] = (1.0+gx)*x[0]*(1.0-x[1]);
  f[2] = (1.0+gx)*(1.0-x[0]);
}
/*
  IMB5
  three-objective
  recommended: n=10;
  PF
   f1^2 + f2^2 + f3^2 = 1
    0 <= xi <= 1
  PS: x_j = (x1+x2)/2
  
*/
void imb5(vector<double> &x, vector<double> &f)
{
  int n = x.size();
  double gx = 0.0;
  if(x[0]>(0.5))
  {
     for(int i = 2; i < n; i++)
     {
	double ti = x[i] - ((x[0]+x[1])*0.5);
	gx += (-0.9*ti*ti + pow(fabs(ti), 0.6)); 
     }
     gx *= 2.0*cos(M_PI*x[0]*0.5);
  }
  f[0] = (1.0+gx)*cos(M_PI*x[0]*0.5)*cos(M_PI*x[1]*0.5);
  f[1] = (1.0+gx)*cos(M_PI*x[0]*0.5)*sin(M_PI*x[1]*0.5);
  f[2] = (1.0+gx)*sin(M_PI*x[0]*0.5);
}
/*
  IMB6
  three-objective
  recommended: n=10;
  PF
   f1^2 + f2^2 + f3^2 = 1
    0 <= xi <= 1
  PS: x_j = (x1+x2)/2
  
*/
void imb6(vector<double> &x, vector<double> &f)
{
  int n = x.size();
  double gx = 0.0;
  if(x[0]>(0.75))
  {
     for(int i = 2; i < n; i++)
     {
	double ti = x[i] - ((x[0]+x[1])*0.5);
	gx += (-0.9*ti*ti + pow(fabs(ti), 0.6)); 
     }
     gx *= 2.0*cos(M_PI*x[0]*0.5);
  }
  f[0] = (1.0+gx)*x[0]*x[1];
  f[1] = (1.0+gx)*x[0]*(1.0-x[1]);
  f[2] = (1.0+gx)*(1.0 - x[0]);
}
/*
  IMB7
  two-objective
  recommended: n=10;
  PF:  f2=1-sqrt(f1)
    0 <= xi <= 1
  PS: x_j = 0.5
  
*/
void imb7(vector<double> &x, vector<double> &f)
{
  int n = x.size();
  double gx = 0.0;
  for(int i = 1; i < n; i++)
  {
     double ti = x[i] - 0.5, si=x[i]-sin(0.5*M_PI*x[0]);
     if(x[0]>=0.5 && x[0]<=0.8)
       gx += (-0.9*si*si + pow(fabs(si), 0.6)); 
     else 
       gx += pow(fabs(ti), 0.6);
  }
  f[0] = (1.0+gx)*x[0];
  f[1] = (1.0+gx)*(1.0-sqrt(x[0]));
}
/*
  IMB8
  two-objective
  recommended: n=10;
  PF:  f2=1-f1
    0 <= xi <= 1
  PS: x_j = 0.5
  
*/
void imb8(vector<double> &x, vector<double> &f)
{
  int n = x.size();
  double gx = 0.0;
     for(int i = 1; i < n; i++)
     {
	double ti = x[i] - 0.5, si=x[i]-sin(0.5*M_PI*x[0]);
	if(x[0]>=0.5 && x[0]<=0.8)
	  gx += (-0.9*si*si + pow(fabs(si), 0.6)); 
	else 
	  gx += pow(fabs(ti), 0.6);
     }
  f[0] = (1.0+gx)*x[0];
  f[1] = (1.0+gx)*(1.0-x[0]);
}
/*
  IMB9
  two-objective
  recommended: n=10;
  PF:  f2=sqrt(1-f1)
    0 <= xi <= 1
  PS: x_j = 0.5
  
*/
void imb9(vector<double> &x, vector<double> &f)
{
  int n = x.size();
  double gx = 0.0;
  for(int i = 1; i < n; i++)
  {
     double ti = x[i] - 0.5, si = x[i]-sin(0.5*M_PI*x[0]);
     if(x[0]>=0.5 && x[0]<=0.8)
       gx += (-0.9*si*si + pow(fabs(si), 0.6)); 
     else 
       gx += pow(fabs(ti), 0.6);
  }
  f[0] = (1.0+gx)*cos(M_PI*x[0]*0.5);
  f[1] = (1.0+gx)*sin(M_PI*x[0]*0.5);
}
/*
  IMB10
  three-objective
  recommended: n=10;
  PF:  f1+f2+f3=1
    0 <= xi <= 1
  PS: x_j = sin(0.5*M_PI*x1)
  
*/
void imb10(vector<double> &x, vector<double> &f)
{
  int n = x.size();
  double gx = 0.0;
  for(int i = 2; i < n; i++)
  {
     double ti = x[i] - x[0]*x[1], si=x[i]-((x[0]+x[1])/2.0);
     if((x[0]>=0.2 && x[0]<=0.8) && ( x[1] >=0.2 && x[1]<=0.8))
       gx += 2.0*(-0.9*si*si + pow(fabs(si), 0.6)); 
     else 
       gx += pow(fabs(ti), 0.6);
  }
  f[0] = (1.0+gx)*x[0]*x[1];
  f[1] = (1.0+gx)*x[0]*(1.0-x[1]);
  f[2] = (1.0+gx)*(1.0-x[0]);
}
int main()
{
  srand(time(0));
  vector<double>f(3), x(10); 
  for(int i = 0 ; i <100000; i++)
  {
    for(int j = 0; j < x.size(); j++) //optimal imb10
    {
	x[j]=rand()/((double)RAND_MAX+1);
	if(j==0||j==1)continue;
	if((x[0]>=0.2 && x[0]<=0.8) && ( x[1]>=0.2 && x[1]<=0.8))x[j]=(x[0]+x[1])/2.0;
	else  
x[j]=x[0]*x[1];
    }
    imb10(x, f);
    for(auto &m:f)cout <<m<<" ";
    cout<<endl;
  }
  return 0;
}
