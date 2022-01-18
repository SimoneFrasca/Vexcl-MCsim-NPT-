#include "./mcsim.hpp"

int main(int argc, char** argv)
{
	vexcl<double> V; 
	double d;
	pvector<double, 3> L, sax;
	  
	rng.rseed();
	if (argc != 4) {
		cout << "Errore, numero errato di argomenti\nInserire nell'ordine:\n1) numero di simulazioni\n2) rate di output\n3) lunghezza semiasse maggiore\n";
		exit(0);  
	}
		
	V.set_steps(atoll(argv[1]));    
	V.set_num_parts(2);
	V.set_outsteps(atoi(argv[2]));
	sax << atof(argv[3]), 1., 1.;
	
	V.set_X0(sax);
	d = 2.*sax.norm();
	L << 2.0*d, 2.0*d, 2.0*d;
	V.set_box(L); 
	V.calcvexcl();

	return 0;
}
