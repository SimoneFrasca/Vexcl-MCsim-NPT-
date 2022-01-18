#include "./mcsim.hpp"
int main(int argc, char** argv)
{
	mcsim<double> mcbox; 
	rng.rseed();
	int nx, ny, nz;
	pvector<double,3> s;
	
	if (argc != 13) {
		cout << "Errore, numero errato di argomenti\nInserire nell'ordine:  \n1)Step simulazione\n2)rate di output\n3)lunghezza blocchi\n4)nx ny nz\n5)lunghezza semiasse maggiore\n6)delta tra\n7)delta rot\n8)Temperatura e pressione\n9)volume escluso\n";
		exit(0);  
	}
	
	// Step di simulazione 
	mcbox.set_steps(atoi(argv[1]));
	// rate di output
	mcbox.set_outsteps(atoi(argv[2]));
	// lunghezza blocchi
	mcbox.set_len(atoi(argv[3]));
	
	// numero di particelle per asse
	nx = atoi(argv[4]);
	ny = atoi(argv[5]);
	nz = atoi(argv[6]);
	mcbox.set_xyz(nx,ny,nz);
	
	// Lunghezza semiassi
	s << atoi(argv[7]), 1., 1.;
	if (s[0]*nx != s[1]*ny or s[0]*nx != s[2]*nz) {
		cout << "ATTENZIONE, il volume di simulazione deve essere cubico\n";
		exit(0);
	}
	mcbox.set_X0(s);
	
	// delta di traslazione e delta di rotazione
	mcbox.set_delta_tra(atof(argv[8]));
	mcbox.set_delta_rot(atof(argv[9]));
	
	// Temperatura
	mcbox.set_T(atof(argv[10]));
	// Pressione
	mcbox.set_P(atof(argv[11]));
	
	// volume escluso
	mcbox.set_v_excl(atof(argv[12]));
	
	// Simulazione MC
	mcbox.runMC();
	
	return 0;
}
