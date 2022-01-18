#include "./pmatrix.hpp"
#include "./check_overlap.hpp"
#include "./randnumgen.hpp"
#include <fstream>

template <class ntype>

class Particle
{
	ntype ranf(void)
	{
		return rng.ranf();
	}
public:

	pvector<ntype,3> r, rold, shift; 
  
	// copia la corrente posiziona in un vettore
	void store (void) {
	rold = r;
	}
	
	// ritorna la posizione precedente in caso di rifiuto della mossa
	void restore(void)
	{
		r=rold;
	}
    
	// estrae una posizione random nella box LxLxL
	void random_insert(pvector<ntype,3>& L)
	{
		for (auto k=0; k < 3; k++)
		r[k] = (ranf()-0.5)*L[k];
	}
    
	// esegue una mossa di traslazione
	void tra_move(ntype delta)
		{
		  for (auto k=0; k < 3; k++)
			r[k] += 2.0*delta*(ranf()-0.5);
		}
		
	  // inizializzazione posizione
	  Particle()
		{
		  r << 0.0,0.0,0.0;
		}
	};

template <class ntype>
class Molecule: public Particle<ntype>
{
	public:
	pvector<ntype,3> sax; // assi della hard box
	using pt = Particle<ntype>;
	pmatrix<ntype,3> R, Rold; // orientazione molecola
	  
	// estrae una posizione ed una orientazione random della particella
	void random_insert(pvector<ntype,3>& L)
		{
		pt::random_insert(L);
		R.orientation();
	}

	// salva la corrente orientazione della particella
	void store(void) {
		pt::store();
		Rold = R;
	}
  
	// ritorna l'orientazione precedente in caso di rifiuto della mossa 
	void restore(void)
		{
		pt::restore();
		R=Rold;
	}

	// esegue una rotazione random
	void rot_move(ntype delth)
	{
		pvector<ntype,3> u;
		pmatrix<ntype,3> omega, M, N, S;
		// generazione di un'orientazione random sulla sfera unitaria
		u.random_orient();
		omega << 0, -u[2], u[1], u[2], 0, -u[0], -u[1], u[0], 0;  
		ntype dtheta = 2.0*(rng.ranf()-0.5)*delth;
		S = (1-cos(dtheta))*(omega*omega);
		M = -sin(dtheta)*omega + S;
		N = R * M;
		R += N;
	}

	// inizialiazzazione orientazione
	Molecule()
	{
		R << 1, 0, 0, 0, 1, 0, 0, 0, 1; // identity matrix
	}
};

template <class ntype>
class HardBox: public Molecule<ntype>
	{
	using mt = Molecule<ntype>;
	public:
	using mt::sax;
	using mt::shift;
	  
	  /*bool overlap(HardBox<ntype> p2)
		{
		  double saxA[3], rA[3], RA[3][3], saxB[3], rB[3], RB[3][3];
		  for (auto k=0; k < 3; k++)
			{
			  saxA[k] = sax[k]/2;
			  saxB[k] = p2.sax[k]/2;
			  rA[k] = mt::r[k];
			  rB[k] = p2.r[k] + shift[k]; // rendo conto dell'immagine minima
			  for (auto k2= 0; k2 < 3; k2++)
			{
			  RA[k][k2] = mt::R[k][k2];
			  RB[k][k2] = p2.R[k][k2];
			}
		}
	  return (check_overlap_boxes(saxA, rA, RA, saxB, rB, RB)<0.0)?true:false;
	}*/
   
   // la funzione analizza l'overlap di due hard boxes
   bool my_overlap(HardBox<ntype> p2)
	{       
		my_check_overlap_boxes<int> overlap;
		return overlap.check_overlap(sax, mt::r, p2.r+mt::shift, mt::R, p2.R);	
	}
   
  // inizializzazione dei semiassi   
  void set_semiaxes(pvector<ntype>& s)
	{
	  sax=s;
	}
	
  HardBox()
	{
	}
};

