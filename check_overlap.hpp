/* Il seguente programma valuta la distanza di due hb misurando la loro
 * distanza nel seguente modo:
 * 1) si misura la distanza tra i centri 
 *   1a) se è maggiore della lunghezza del diametro non si sovrappongono
 *   1b) altrimenti si prosegue
 * 2) si divide in cubi e si valuta quali sono i più vicini
 * 3) la segmentazione continua fino a che la distanza minima tra due 
 * cubi è minore di un lato (per cui si sovrappongono) oppure si arresta
 * se maggiore del diametro del cubo, restituendo nessuna sovrapposizione*/
 
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>


using namespace std;
template  <class ntype, int NT = 3>

class my_check_overlap_boxes {
	int over = 0;
	
	public:	
	
	pvector<double,NT> v1, v2, v1t, v2t;
	pvector<double,NT> w1, w2, z, zt, h1, h2;
	pvector<double,NT> v1b, v2b, sax;
	pmatrix<double,NT> R1, R2;
	double d = 0, m = 1, d1, d2, min, max, x0;
	
	//possibili spostamenti dei punti
	int a[8] = {1,1,1,1,-1,-1,-1,-1};
	int b[8] = {-1,-1,1,1,-1,-1,1,1};
	int c[8] = {1,-1,1,-1,1,-1,1,-1};
	int overtemp; 

	// funzioni che tengono il controllo  dell'overlap
	// se esiste possibilità di overlap restituisce 1
	// se l'overlap è certo restiuisce 2
	// se non c'è minima possibilità di overlap restituisce 0
	void check(int a) {
		if (a == 1 and over != 2) {
			over = 1;
			}
		if (a == 2) {
			over = 2;
			}
		}
		
	void overlap(double d, double dt) {
		overtemp = 0; 
		if (d <= min) {
			overtemp = 2;
			}
		if (d < max and d > min) {
			overtemp = 1;
			if (d < dt) {
				zt = z;
				v1b = w1;
				v2b = w2;
			} 
		}
		check(overtemp);
	}


	// funzione che verifica la distanza dei centri dei parallelepipedi
	int firststep() {
		over = 0;
		min = 2.*sax[2];
		max = 2.*sax.norm();
		z = v1-v2;
		double d = sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);
		if (d <= min) {
			over = 2;
			}
		if (d < max and d > min) {
			over = 1;
		}
		return over;
	}

	// si valuta la distanza tra le segmentazioni dei cubi
	int otherstep() {
		double d, dt;
		int m = 1;
		
		zt << 2.*sax[0], 2.*sax[0], 2.*sax[0];
		over = 0;
		
		min = 2.*sax[2];
		max = 2.*sqrt(1.+sax[1]*sax[1]+sax[2]*sax[2]);
		// prima segmentazione (si divide il parallepipedo in cubi)
		for (auto j = 0; j < (int)(sax[0]); j++) {
			for (auto l = 0; l < (int)(sax[0]); l++) {
				d1 = l*2-(int)(sax[0])+1;
				d2 = j*2-(int)(sax[0])+1;
				w1 << R1[0][0]*d1, R1[0][1]*d1, R1[0][2]*d1;
				w2 << R2[0][0]*d2, R2[0][1]*d2, R2[0][2]*d2;
				w1 = w1+v1;
				w2 = w2+v2;
				z = w1-w2;
				d = sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);
				dt = sqrt(zt[0]*zt[0]+zt[1]*zt[1]+zt[2]*zt[2]);
				overlap(d,dt);
			}
		}
		
		
		if (over == 0 or over == 2) { 
			return over;
		}
		
		if (rint(sax[1]) != 1) {
			v1t = v1b;
			v2t = v2b; 
			over = 0;
			min = 2.*sax[2];
			max = 2.*sqrt(2.+sax[2]*sax[2]);
			// prima segmentazione (si divide il parallepipedo in cubi)
			for (auto j = 0; j < (int)(sax[1]); j++) {
				for (auto l = 0; l < (int)(sax[1]); l++) {
					d1 = l*2-(int)(sax[1])+1;
					d2 = j*2-(int)(sax[1])+1;
					w1 << R1[1][0]*d1, R1[1][1]*d1, R1[1][2]*d1;
					w2 << R2[1][0]*d2, R2[1][1]*d2, R2[1][2]*d2;
					w1 = w1+v1t;
					w2 = w2+v2t;
					z = w1-w2;
					d = sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);
					dt = sqrt(zt[0]*zt[0]+zt[1]*zt[1]+zt[2]*zt[2]);
					overlap(d,dt);
				}
			}
		}
		
		
		if (over == 0 or over == 2) { 
			return over;
		}
		
		if (rint(sax[2]) != 1) {
			v1t = v1b;
			v2t = v2b; 
			over = 0;
			min = 2.;
			max = 2.*sqrt(3.);
			// prima segmentazione (si divide il parallepipedo in cubi)
			for (auto j = 0; j < (int)(sax[2]); j++) {
				for (auto l = 0; l < (int)(sax[2]); l++) {
					d1 = l*2-(int)(sax[2])+1;
					d2 = j*2-(int)(sax[2])+1;
					w1 << R1[2][0]*d1, R1[2][1]*d1, R1[2][2]*d1;
					w2 << R2[2][0]*d2, R2[2][1]*d2, R2[2][2]*d2;
					w1 = w1+v1t;
					w2 = w2+v2t;
					z = w1-w2;
					d = sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);
					dt = sqrt(zt[0]*zt[0]+zt[1]*zt[1]+zt[2]*zt[2]);
					overlap(d,dt);
				}
			}
		}
		
		
		if (over == 0 or over == 2) { 
			return over;
		}
		
		// successive segmentazioni dei cubi
		while (over == 1 and m < 7) {
			v1t = v1b;
			v2t = v2b; 
			over = 0;
			min = 2./pow(2,m);
			max = 2.*sqrt(3)/pow(2,m);
			for (auto j = 0; j < 8; j++) {
				for (auto l = 0; l < 8; l++) {
					h1 << a[l]*R1[0][0]+b[l]*R1[1][0]+c[l]*R1[2][0], a[l]*R1[0][1]+b[l]*R1[1][1]+c[l]*R1[2][1], a[l]*R1[0][2]+b[l]*R1[1][2]+c[l]*R1[2][2];
					h2 << a[j]*R2[0][0]+b[j]*R2[1][0]+c[j]*R2[2][0], a[j]*R2[0][1]+b[j]*R2[1][1]+c[j]*R2[2][1], a[j]*R2[0][2]+b[j]*R2[1][2]+c[j]*R2[2][2];
					w1 = (h1*2./pow(2,m+1))+v1t;
					w2 = (h2*2./pow(2,m+1))+v2t;
					z = w1-w2;
					d = sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);
					dt = sqrt(zt[0]*zt[0]+zt[1]*zt[1]+zt[2]*zt[2]);
					overlap(d,dt);
				}
			}
			m++;
		}
		
		// se dopo 7 segmentazioni del cubo non ho un risultato positivo 
		// si assume nessun overlap 
		if (m == 7) over = 2;
		
		return over;
	}

	
	bool check_overlap(pvector<double,NT> s, pvector<double,NT> v, pvector<double,NT> w, pmatrix<double,NT> M, pmatrix<double,NT> N)
    {
		sax = s; 
		v1 = v; // posizione prima hb
		v2 = w; // posizione seconda hb
		R1 = M; // orientazione prima hb
		R2 = N; // orientazione seconda hb
		if (firststep() == 2) {
			return true;
		}
		
		if (firststep() == 1) {
			if (otherstep() == 2) {
				return true;
			}
		}
      return false;
    }
};
