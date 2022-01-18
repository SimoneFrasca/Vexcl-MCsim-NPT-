#ifndef _PVECTOR_
#define _PVECTOR_
// TODO
// implement:
// - scalar times vector
// - cross product
// - norm (abs(v))
// rint(v)

#include<iostream>
#include<cstdlib>
#include<cstdarg>
#include "./randnumgen.hpp"
#define VARIADIC_ARGS_CPP_STYLE
using namespace std;
template <class ntype,int NT=3>
class pvector {
 ntype v[NT]; // array of NT objects of type ntype
public:
 pvector() // constructor
   {
     for (auto i=0; i < NT; i++)
       v[i] = 0;
   }
 ~pvector() // destructor
   {
   }
#ifdef VARIADIC_ARGS_CPP_STYLE

 template <class...Args>
   void construct(ntype value, Args...args) // this means you can have several argument with different types.
   // The compiler deduces all the types associated to the args!
   // Value is the first argument and args are all the remaining arguments 
     {
       if (idx >= NT)
         { 
           cout << "[ERROR] too many arguments in constructor\n";
           exit(0);
         }
       v[idx++] = (ntype) value; 
       construct(args...); // call construct with all arguments except T (recursion)
                           // args are all arguments except T (see also above)
     }
 // stop recursion since construct is called as last recursion step without arguments
 void construct()
   {}

 // constructor with variable number of arguments with variable type
 template <class...Args>
   pvector(Args...args) 
     {
       idx=0;
       construct(args...);
     }
#else
 // this solution works but it is rather unpleasant...
 pvector(ntype f,...)
   {
     va_list va;
     va_start(va, f);
     v[0]=f;
     for (auto i=1; i < NT; i++)
       { 
         v[i] = va_arg(va, ntype);  
       }
     va_end(va);
   }
#endif
 pvector(ntype a, ntype b=0, ntype c=0)
   {
     if (NT >= 1)
       {
         v[0]= a;
       }
     if (NT >= 2)
       {
         v[1] = b; 
       }
     if (NT >= 3)
       {
         v[2] = c;
       }
   } 

 void set(int i, ntype a)
   {
     v[i] = a;
   }

 ntype get(int i)
   {
     return v[i];
   }
 int idx;
 pvector<ntype,NT>& operator<<(const ntype& a)
   {
     v[0] = a;
     idx=0;
     return (*this);
   }
 pvector<ntype,NT>& operator,(const ntype& a)
    {
      idx++;
      if (idx >= NT)
        {
          cout << "[ERROR] vector has only" << NT << " elements\n";
          exit(1);
        }
      v[idx] = a;
      return (*this);
    }
    
 ntype& operator[](const int i)
   {
     return v[i];
   }

 bool operator==(const pvector<ntype,NT>& B)
    {
      for (auto i = 0; i < NT; i++) {
		  if (v[i] != B.v[i]) {
				return false;
			  }
		  }
		return true;
    }

   
 // sum is *this + B
 pvector<ntype,NT> operator+(const pvector<ntype,NT>& B)
   {
     pvector<ntype,NT> vt;
     for (auto i=0; i < NT; i++)
       {
         vt.v[i] = v[i] + B.v[i];
       }
     return vt;
   } 
   
 // diff. is *this - B
 pvector<ntype,NT> operator-(const pvector<ntype,NT>& B)
   {
     pvector<ntype,NT> vt;
     for (auto i=0; i < NT; i++)
       {
         vt.v[i] = v[i] - B.v[i];
       }
     return vt;
   } 
   
 // diff. is *this - B
 pvector<ntype,NT> operator-=(const pvector<ntype,NT>& B)
   {
     pvector<ntype,NT> vt;
     for (auto i=0; i < NT; i++)
       {
         vt.v[i] -= B.v[i];
       }
     return vt;
   } 
   
 // scalar product
 ntype operator*(const pvector<ntype,NT>& v1)
   {
     ntype s=0;
     for (auto i=0; i < NT; i++)
       {
         s+=v[i]*v1.v[i];
       }
     return s;
   }
  
 // product for a scalar (vector times a scalar), 1st way
 pvector<ntype,NT> operator*(const ntype& s)
   {
     pvector<ntype,NT> vt;
     for (auto i=0; i < NT; i++)
       {
         vt.v[i] = v[i]*s;
       }
     return vt;
   } 
   
 // product for a scalar, 2nd way
 pvector<ntype,NT>& operator*=(const ntype &s)
   {
     for (auto i=0; i < NT; i++)
       {
        v[i] *= s;
       }
     return (*this);
   }  
   
 // product for a scalar (vector times a scalar), 1st way
 pvector<ntype,NT> operator/(const ntype& s)
   {
     pvector<ntype,NT> vt;
     for (auto i=0; i < NT; i++)
       {
         vt.v[i] = v[i]/s;
       }
     return vt;
   } 
   
 // product for a scalar, 2nd way
 pvector<ntype,NT>& operator/=(const ntype &s)
   {
     for (auto i=0; i < NT; i++)
       {
        v[i] /= s;
       }
     return (*this);
   }
   
 // scalar times vector  s*v, 3rd way 
 friend pvector<ntype,NT> operator*(const ntype& s, const pvector<ntype, NT>& v2)
   {
     pvector<ntype, NT> vt;
     for (int i=0; i < NT; i++)
       vt.v[i] = s*v2.v[i];
     return vt;
   }
 
 // cross product
 pvector<ntype, NT> operator^(const pvector<ntype,NT>& v2)
   {
     pvector<ntype, NT> vt;
     if constexpr (NT==3)
       {
         vt.v[0] = v[1]*v2.v[2] - v[2]*v2.v[1];
         vt.v[1] = v[2]*v2.v[0] - v[0]*v2.v[2];
         vt.v[2] = v[0]*v2.v[1] - v[1]*v2.v[0];
       }
     else
       {
         cout << "[ERROR] cross product is just for 3d vector\n";
         exit(0);
       }
     return vt;
   }
   
 ntype norm()
   {
     return sqrt((*this)*(*this));
   }
   
 // \Delta r_{ij} = r_i - r_j
 // \Delta r'_{ij} = \Delta r'_{ij} - L*\rint(\Delta r_{ij}/L)
pvector<ntype,NT> rint(pvector<ntype,3> L)
   {
     pvector <ntype, NT> vt;
     for (int i=0; i < NT; i++)
       {
         vt.v[i]=L.v[i]*std::rint(v[i]/L.v[i]);
       }
     return vt;
   }

 // method which generates random number in the interval [0,1]
 ntype ranf()
   {
     return rng.ranf();
   }

 // generate a random vector of the uniform
 void random_orient()
   {
     if (NT==3)
       {
         // x1 in [-1, 1]
         // Marsaglia algorithm 
         // Tha annals of Mathematical statistics 43, 645 (1972)
         ntype x1, x2, xsq, xi; 
         do
           {
             x1 = 2.0*ranf()-1.0;
             x2 = 2.0*ranf()-1.0;
             xsq= x1*x1+x2*x2; 
           }
         while (xsq > 1.0);
         xi = sqrt(1.0 - xsq);
         v[0] = 2.0*x1*xi;
         v[1] = 2.0*x2*xi;
         v[2] = 1.0 - 2.0*xsq;
       }
     else
       {
         cout << "[ERROR] random unit vector implemented only for 3d vectors\n";
         exit(1); 
       }
   }

 // TODO
 void show(void)
   {
     cout << "(";
     for (auto i=0; i < NT; i++)
       {
         cout << v[i] ;
         if (i < NT-1)
           cout << ",";
         else
           cout << ")";
       }
   }
};
#endif
