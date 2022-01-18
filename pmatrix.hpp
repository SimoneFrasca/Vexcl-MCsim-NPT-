#ifndef _PMATRIX_
#define _PMATRIX_ 

#include "./randnumgen.hpp"
#include "./pvector.hpp"
#include <iostream>
using namespace std;

template <class ntype, int N>
class pmatrix
{
   pvector<ntype,3> m[N];
   int idx;
public:
   pmatrix()
     {
     }
   ~pmatrix()
     {}
     
//----------------------------------------------------------------------     
   // Comma inizialization
   pmatrix<ntype,N>&  operator<<(const ntype& param)
     {
       idx=0;
       m[0][0] = param;
       return (*this);
     }
     
   pmatrix<ntype,N>& operator,(const ntype& param)
     {
       idx++;
       if (idx >= N*N)
         {
           cout << "[ERROR] Too many elements provided in comma initialization\n";
           exit(1);
         }
       m[idx/N][idx%N]=param;
       return (*this);
     }

//----------------------------------------------------------------------     

   // M[i][j]=...
   ntype* operator[](int i)
     {
       return &(m[i][0]);
     }

//----------------------------------------------------------------------    
   // addition 
   // note that omitting the const modifier in the argument, m2 can not be a temporary object like pmatrix<double,3>()
   pmatrix<ntype,N> operator+(pmatrix<ntype,N> &m2)
     {
       // mt = m1+m2
       pmatrix<ntype,N> mt;
       for (auto i=0; i < N; i++)
         {
           for (auto j=0; j < N; j++)
             {
               mt[i][j] = (*this).m[i][j]+m2[i][j];
             }
         }
        return mt;
     }
     
    pmatrix<ntype,N> operator+=(pmatrix<ntype,N> &m2)
     {
       for (auto i=0; i < N; i++)
         {
           for (auto j=0; j < N; j++)
             {
               (*this).m[i][j] += m2[i][j];
             }
         }
        return (*this);
     }
     
//----------------------------------------------------------------------
   // subtraction
   pmatrix<ntype,N> operator-(pmatrix<ntype,N> &m2)
     {
       // mt = m1+m2
       pmatrix<ntype,N> mt;
       for (auto i=0; i < N; i++)
         {
           for (auto j=0; j < N; j++)
             {
               mt[i][j] = (*this).m[i][j]-m2[i][j];
             }
         }
        return mt;
     }
   
        
    pmatrix<ntype,N> operator-=(pmatrix<ntype,N> &m2)
     {
       for (auto i=0; i < N; i++)
         {
           for (auto j=0; j < N; j++)
             {
               (*this).m[i][j] -= m2[i][j];
             }
         }
        return (*this);
     }
     
//----------------------------------------------------------------------  
    //matrix devided by a scalar
    
    pmatrix<ntype,N> operator/(const ntype& s)
   {
      pmatrix<ntype,N> mt;
       for (auto i=0; i < N; i++)
         {
           for (auto j=0; j < N; j++)
             {
               mt[i][j] = (*this).m[i][j]/s;
             }
         }
        return mt;
   } 
   
       
   pmatrix<ntype,N> operator/=(const ntype& s)
   {
       for (auto i=0; i < N; i++)
         {
           for (auto j=0; j < N; j++)
             {
              (*this).m[i][j] /= s;
             }
         }
        return (*this);
   } 
        
//----------------------------------------------------------------------  
    //matrix times a scalar
    
    pmatrix<ntype,N> operator*(const ntype& s)
   {
      pmatrix<ntype,N> mt;
       for (auto i=0; i < N; i++)
         {
           for (auto j=0; j < N; j++)
             {
               mt[i][j] = (*this).m[i][j]*s;
             }
         }
        return mt;
   } 
   
       
   pmatrix<ntype,N> operator*=(const ntype& s)
   {
       for (auto i=0; i < N; i++)
         {
           for (auto j=0; j < N; j++)
             {
              (*this).m[i][j] *= s;
             }
         }
        return (*this);
   } 
    
   // scalar times matrix    
   friend pmatrix<ntype,N> operator*(ntype s, pmatrix<ntype,3> m)
     {
       pmatrix<ntype,N> mt;
       for(auto i=0; i < N; i++)
         {
           for (auto j=0; j < N; j++)
             {
               mt[i][j] = m[i][j]*s;
             }
         }
       return mt;
     }   
   
//----------------------------------------------------------------------
   // matrix times vector
   // pvector<double,3>(1,2,3)
   pvector<ntype,N> operator*(pvector<ntype,N> &v2) 
     {
       // mt = (*this)*v2
       pvector<ntype,N> vt;
       for (auto i=0; i < N; i++)
         {
           vt[i] = 0.0; 
           for (auto j=0; j < N; j++)
             {
               vt[i] += m[i][j]*v2[j];
             }
         }
        return vt;
     }
     
   // vector times matrix
   friend pvector<ntype,N> operator*(pvector<ntype>& v1, pmatrix<ntype,N>& m2)
     {
       // vt = m1*v2
       pvector<ntype,N> vt;
       for (auto i=0; i < N; i++)
         {
           vt[i] = 0.0; 
           for (auto j=0; j < N; j++)
             {
               vt[i] += v1[j]*m2[j][i];
             }
         }
        return vt;
     }
          
//----------------------------------------------------------------------  
   //matrix times matrix
   pmatrix<ntype,N> operator*(pmatrix<ntype,N> &m2)
     {
       // mt = m1*m2
       pmatrix<ntype,N> mt;
       for (auto i=0; i < N; i++)
         {
           for (auto j=0; j < N; j++)
             {
               mt[i][j] = 0.0;
               for (auto l=0; l < N; l++)
                 {
                   mt[i][j] += (*this).m[i][l]*m2[l][j];
                 }
             }
         }
        return mt;
     }
          
//----------------------------------------------------------------------  
   // transpose
   pmatrix<ntype,N> transpose(pmatrix<ntype,N> &m2)
     {
       // mt = m1*m2
       pmatrix<ntype,N> mt;
       for (auto i=0; i < N; i++)
         {
           for (auto j=0; j < N; j++)
             {
               mt[i][j] = m2[j][i];
             }
         }
        return mt;
     }
     
//----------------------------------------------------------------------  
pvector<ntype,N> get_col(int i)
     {
       pvector<ntype,N> vt;
       for (auto j=0; j < N; j++)
         {
           vt[j] = m[j][i];
         }
      return vt;
     }
     
//----------------------------------------------------------------------  
   pvector<ntype,N> get_row(int i)
     {
       pvector<ntype,N> vt;
       for (auto j=0; j < N; j++)
         {
           vt[j] = m[i][j];
         }
       return vt;
     }
     
//----------------------------------------------------------------------  
   void set_col(int i, const pvector<ntype,N>& v)
     {
       for (auto j=0; j < N; j++)
         {
           m[j][i]=v[j];
         }
     }
     
//----------------------------------------------------------------------  
   void set_row(int i, pvector<ntype,N>& v)
     {
       for (auto j=0; j < N; j++)
         {
           m[i][j]=v[j];
         }
     }
     
//----------------------------------------------------------------------  
   void random()
     {
       if constexpr (N == 3)
         { 
           pvector<ntype,N> v1, v2, v3;
           v1.random_orient();
           v2 << 1, 0, 0;
           if (v2==v1)
             {
               v2 << 0,1,0;
             }
           v2 = v2 - (v2*v1)*v1; // GS orthonormalization
           v2 = v2 * (1.0/v2.norm());
           v3 = v1 ^ v2;
           ntype theta=rng.ranf()*2.0*M_PI;
           v2 = cos(theta)*v2 + sin(theta)*v3;
           v3 = v1^v2;
           set_row(0,v1);
           set_row(1,v2);
           set_row(2,v3);
         }
     }
     
     
//----------------------------------------------------------------------  
     void orientation() {
		double sp, cp, ct, st, csi;
		pvector<ntype,3> x,y,z;
		x.random_orient();
		
		sp = x[1];
		cp = sqrt(1-x[1]*x[1]);
		ct = x[2]/cp;
		st = x[0]/cp;
		csi = rng.ranf()*M_PI;
		 
		y << ct*cos(csi)-sp*st*sin(csi), sin(csi)*cp, -st*cos(csi)-sin(csi)*sp*ct;
		z = (y^x)*randsgn();
		
        set_row(0,y);
        set_row(1,x);
        set_row(2,z);
	}
//----------------------------------------------------------------------  
   void show(const char* str="")
    {
      int i, j;
      if (str!=NULL)
        cout << str;
      cout << "{";
      for (i=0; i < N; i++)
        {
          cout << "{";
          for (j=0; j < N; j++)
            { 
              cout << m[i][j];
              if (j < N-1)
                cout << ",";
            }
          cout << "}";
          if (i < N-1)
            cout << ",\n";
        }
      cout << "}\n";
    }
    
//----------------------------------------------------------------------  
   	double randsgn() {
		if (rng.ranf() > 0.5) {
			return 1;
		}
		else {
			return -1;
		}
	}
};

#endif
