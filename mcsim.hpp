// Questa parte del programma ha lo scopo sia di calcolare il volume escluso 
// di una hard box, che eseguire una simulazione MC in un ensamble NPT
// di un sistema di N hard box a temperatura T e pressione P
#include "./particles.hpp"
#include <vector>
#include <fstream>
#include <time.h>

template <class ntype> 


//######################################################################
//                     STATO DELLA SIMULAZIONE
//######################################################################
// si tiene conto dello stato della simulazione per l'ottimizzazione 
// dei parametri 
class state_of_sim 
{
	public:
	long int tra_rej, rot_rej, tot_rot, tot_tra;
	// mosse di traslazione/rotazione totali/rifiutate
	long int vol_rej, tot_vol;
	// mosse di volume totali/rifiutate
	state_of_sim()
		{
		tra_rej = rot_rej = tot_rot = tot_tra = vol_rej = tot_vol = 0;
		}
	};

//######################################################################
//                  PARAMETRI DI SIMULAZIONE
//######################################################################
template <class ntype> 
// Parametri del sistema
class parsC
{
	public:
	
	long int Nsteps, Np, outsteps; // numero di step della simulazione MC, numero di particelle, numeri di passi per l'output
	pvector<ntype,3> L;		// lati della box di simulazione
	int nx , ny, nz;		// numero di particelle per coordinata
	ntype phi0;				// volume d'interazione
	pvector<ntype,3> sax;	// semiassi dei parallelepipedi
	ntype T, P;				// temperatura e pressione
	ntype delta_tra, delta_rot, vmax; // parametri di traslazione/rotazione e variazione volumere
	ntype MSD;				// mean square displacement
	int tau, len, nblock;	// parametri per eseguire la media a blocchi 
							// tau = primo tempo per cui si inizia a misurare
							// len = numero di elementi per blocco
							// numero di blocchi
	vector<ntype> stat = {0,0};		// media e incertezza
	ntype vhb, v_excl;
	
	parsC()
	{
		Nsteps = Np = outsteps = 0;
		len = 1;
		MSD = 0;
		nx = ny = nz = 0;
		sax << 0, 0, 0;
		phi0= 0; 
		delta_tra = delta_rot = 0;
		vmax=0;
		T = P = 0; 
		vhb = v_excl = 0;
		}
	};


//######################################################################
//                  LINKED LIST CELL 
//######################################################################
template <class ntype> 
class Linked_List_Cell{
	public:
	int M;				// numero di celle per lato
	int NCELL;			// numero totale di celle
	int MAPSIZE;		// dimesione della mappa dei primi vicini
	ntype cell_size;	// lato della cella
	int* HEAD;			// vettore HEAD
	int *MAP;			// vettore MAP
	int *LIST;			// vettore LIST
	
	Linked_List_Cell() {
		M = 0;
		MAPSIZE = 0;
		NCELL = 0;
		cell_size = 0;
	}

	// ritorna l'indice della cella
	int ICELL(int x, int y, int z) {
		return ((x+M)%M)+((y+M)%M)*M + ((z+M)%M)*M*M;
	}

	// inizializza i primi vicimi pe rogni cella
	void MAPS() {
		int x, y, z, IMAP;
		
		for (x = 0; x < M; x++) {
			for (y = 0; y < M; y++) {
				for (z = 0; z < M; z++) {
					IMAP = ICELL(x,y,z)*26;
					MAP[IMAP  ] = ICELL(x+1,y  ,z  );
					MAP[IMAP+1] = ICELL(x+1,y+1,z  );
					MAP[IMAP+2] = ICELL(x  ,y+1,z  );
					MAP[IMAP+3] = ICELL(x-1,y+1,z  );
					MAP[IMAP+4] = ICELL(x-1,y  ,z  );
					MAP[IMAP+5] = ICELL(x-1,y-1,z  );
					MAP[IMAP+6] = ICELL(x  ,y-1,z  );
					MAP[IMAP+7] = ICELL(x+1,y-1,z  );
					
					MAP[IMAP+8 ] = ICELL(x  ,y  ,z+1);
					MAP[IMAP+9 ] = ICELL(x+1,y  ,z+1);
					MAP[IMAP+10] = ICELL(x+1,y+1,z+1);
					MAP[IMAP+11] = ICELL(x  ,y+1,z+1);
					MAP[IMAP+12] = ICELL(x-1,y+1,z+1);
					MAP[IMAP+13] = ICELL(x-1,y  ,z+1);
					MAP[IMAP+14] = ICELL(x-1,y-1,z+1);
					MAP[IMAP+15] = ICELL(x  ,y-1,z+1);
					MAP[IMAP+16] = ICELL(x+1,y-1,z+1);
					
					MAP[IMAP+17] = ICELL(x  ,y  ,z-1);
					MAP[IMAP+18] = ICELL(x+1,y  ,z-1);
					MAP[IMAP+19] = ICELL(x+1,y+1,z-1);
					MAP[IMAP+20] = ICELL(x  ,y+1,z-1);
					MAP[IMAP+21] = ICELL(x-1,y+1,z-1);
					MAP[IMAP+22] = ICELL(x-1,y  ,z-1);
					MAP[IMAP+23] = ICELL(x-1,y-1,z-1);
					MAP[IMAP+24] = ICELL(x  ,y-1,z-1);
					MAP[IMAP+25] = ICELL(x+1,y-1,z-1);
				}
			}
		}
	}

	// posiziona la le particelle nelle corrette celle
	void LINKS(int N, vector<HardBox<ntype>> parts, pvector<ntype,3> L)
	{
		int icell, i;
		double invcell_size, halfL;
		pvector<ntype,3> r;
		 
		invcell_size = 1./cell_size;
		halfL = L[0]/2.;
		
		for (icell = 0; icell< NCELL; icell++) {
			HEAD[icell] = -1;
		}
			
		for (i=0; i < N; i++) {
			r = parts[i].r-parts[i].r.rint(L);
			icell = (int)((r[0]+halfL)*invcell_size) + (int)((r[1]+halfL)*invcell_size)*M + (int)((r[2]+halfL)*invcell_size)*M*M;
			LIST[i] = HEAD[icell];
			HEAD[icell]=i;
		}
	}

};

//######################################################################
//                   CALCOLO VOLUME ESCLUSO
//######################################################################
template <class ntype> 
class vexcl
{
	fstream of; // stream di output
	public:

	vector<HardBox<ntype>> parts;
	parsC<ntype> pars;
	int outsteps = 0;
  
	void calcvexcl(void)
	{
		long int over = 0; // conteggio degli overlap
		// La prima particella è sempre centrata in 0 con orientazione 
		// lungo gli assi
		parts[0].r << 0.0, 0.0, 0.0;
		parts[0].R << 1, 0, 0, 0, 1, 0, 0, 0, 1;
		ntype Vol = pars.L[0]*pars.L[1]*pars.L[2]; // volume box
      
		of.open("vexcl.dat", ios::out|ios::trunc);
		for (auto i=0; i < pars.Nsteps; i++)
		{
			// generazione della seconda particella
			parts[1].random_insert(pars.L);
			if (parts[0].my_overlap(parts[1]))
			{
				over++;
			}
			if (i % pars.outsteps == 0)
			{
				cout << "step = " << i << " V_excl = " << Vol*((ntype)over)/(i+1)  << "\n";        
				of << i << " " << Vol*((ntype)over)/(i+1) << "\n"; 
			}
		}
		of.close();
		cout << "V_excl = " <<  Vol*((ntype)over)/pars.Nsteps << "\n";
	}
	
	// inizializzazione dei parametri
	
	// inizializzazione box di simulazione
	void set_box(pvector<ntype> Lt)
	{
		pars.L = Lt;
	}
    
	// numero di particelle da confrontare per l'overlap
	void set_num_parts(int Np)
	{
		parts.resize(Np);
	}
	  
	// inizializzazione dei semiassi
	void set_X0(pvector<ntype,3> s)
	{
		pars.sax = s;
		for (auto& p: parts)
		{
			p.set_semiaxes(s);
		}
	}
	  
	// numero di step di simulazione
	void set_steps(long int Nt)
		{
		pars.Nsteps = Nt;
		}
	  
	// rate di visualizzazione
	void set_outsteps(long int Nt)
	{
		pars.outsteps = Nt;
	}
};
  

//######################################################################
//                  SIMULAZIONE MONTE CARLO IN NPT 
//######################################################################
template <class ntype> 
class mcsim
{

public:

	vector<HardBox<ntype>> parts, in_parts, no_bound_parts;
	// parts = contiene le informazioni sulla particelle, sotto le condizioni
	// delle pbc
	// in_parts = particelle nella configurazione iniziale
	// no_bound_parts = particelle senza le pbc
	parsC<ntype> pars;
	state_of_sim<ntype> state;
	Linked_List_Cell<ntype> llc;
	ntype *block;		// blocco per eseguire la media a blocchi
	
	// Preparazione configurazione iniziale su un simple cubic
	void prepare_initial_conf(void)
	{
		int i = 0; 
		pars.Np = int(pars.nx * pars.ny * pars.nz);
		pvector<ntype,3> dr;
		// si definiscono i tre assi dei parallelepipedi
		dr << 2.*pars.sax[0], 2.*pars.sax[1], 2.*pars.sax[2];
		parts.resize(pars.Np);
		pars.vhb = dr[0]*dr[1]*dr[2];
		
		// inizializzazione della phi0 invertendo l'EOS secondo
		// l'espansione del viriale
		pars.phi0 = sqrt((pars.vhb*pars.vhb/(pars.v_excl*pars.v_excl))+2.*pars.P/pars.T*pars.vhb*pars.vhb/pars.v_excl)-pars.vhb/pars.v_excl;
		
		// alfa dice quanto sono vicine i parallelepipedi nella configurazione iniziale
		// dato il volume di interazione
		ntype alpha = cbrt(1./pars.phi0);
		
		// di conseguenza alla scelta di alpha viene costruita la box
		pars.L << pars.nx*dr[0]*alpha, pars.ny*dr[1]*alpha, pars.nz*dr[2]*alpha;
		// si inizializza il vmax per la mossa del volume come una sua millesima parte
		pars.vmax = pars.L[0]*pars.L[1]*pars.L[2]/1000;
		
		// eseguo l'inizializzazione per le Np = nx*ny*nz particelle
		int ix, iy, iz;  
		  for (ix = 0; ix < pars.nx; ix++) { //cicli sulle celle primitive
			for (iy = 0; iy < pars.ny; iy++) { //cicli sulle celle primitive
				for (iz = 0; iz < pars.nz; iz++) { //cicli sulle celle primitive
					// semiassi
					parts[i].sax = dr/2.;
					// si posizionano le particelle su un reticolo equispaziate
					parts[i].r << dr[0]*(ix+0.5)*alpha, dr[1]*(iy+0.5)*alpha, dr[2]*(iz+0.5)*alpha;
					// le particelle sono orientate tutte nella stessa direzione
					parts[i].R << 1, 0, 0, 0, 1, 0, 0, 0, 1;
					i++;
			}
		  }
		}
		
		// si esegue una traslazione della box in modo che il centro sia in 0
		for (i = 0; i < pars.Np; i++)
			{
			  parts[i].r = parts[i].r - pars.L*0.5;
			}
			
		in_parts = parts;
		no_bound_parts = parts;
	}   
  	
	void runMC()
	{
		int np, movetype, j = 0, k = 0, l = 0, ok = 0;
		// inizializzazione delle particelle nella box
		prepare_initial_conf();
		int ntot = pars.Np+1;
		// inizializzazione della LLC
		
		llc.cell_size = pars.L[0]/(int)(pars.L[0]/(2.*pars.sax.norm()));
		llc.M = (int)(pars.L[0]/llc.cell_size);
		llc.NCELL = llc.M*llc.M*llc.M;
		llc.MAPSIZE = 26*llc.NCELL;
		llc.LIST = (int*)calloc(pars.Np, sizeof(double));
		llc.HEAD = (int*)calloc(llc.NCELL, sizeof(double));
		llc.MAP = (int*)calloc(llc.MAPSIZE, sizeof(double));
		
		llc.MAPS();
		llc.LINKS(pars.Np,parts,pars.L);
		
		// file che restituiscono l'evoluzione temporale di phi e della MSD
		// in funzione del numero dei passi
		ofstream g;
		g.open("phi_evolution.dat");
		ofstream f;
		f.open("MSD.dat");
		for (long int i = 1; i <= pars.Nsteps; i++) {
			for (auto n = 0; n < pars.Np; n++) {
				bool rejected = false;
				// scelta della particella
				np = (int)(rng.ranf()*ntot);
				
				// mossa di volume con probabilità 1/Np
				if (np == pars.Np) {
					if (move_box()) state.vol_rej++;
				}
				else
				{ 
					// si salva una copia della configurazione attuale
					store_part(np);
					// si esegue una mossa
					movetype=mcmotion(np);
					no_bound_parts[np].r = no_bound_parts[np].r + parts[np].r-parts[np].rold;
					no_bound_parts[np].R = no_bound_parts[np].R + parts[np].R-parts[np].Rold;
					// si impongono per parts le pbc
					pbc(np);
					// si aggiornano gli elementi nella LLC
					llc.LINKS(pars.Np,parts,pars.L);
					// se si overlappano ritorno alla configurazione precedente
					if (check_all_overlaps(np))
					{
						restore_part(np);
						llc.LINKS(pars.Np,parts,pars.L);
						rejected=true;
					}
					if (movetype==0) // translation
					{
						state.tot_tra++;
						if (rejected) state.tra_rej++;
					}
					else
						{
						state.tot_rot++;
						if (rejected) state.rot_rej++;
						}
					}
				}
			// calcolo del MSD
			calc_MSD();
			 
			// condizione per verificare l'equilibrio del sistema
			if (pars.MSD/pars.Np/(i+1) > pars.L[0]*pars.L[0]) ok = 1;
			
			if (ok == 1 && l==0) {
				if (k == 0) {
					pars.tau = i;
					pars.nblock = (int)(floor((pars.Nsteps-pars.tau)/pars.len));
					block = (ntype*)calloc(pars.nblock, sizeof(double));
					k = 1;
				}
				block[j] += pars.vhb*pars.Np/(pars.L[0]*pars.L[1]*pars.L[2]);
				if (i-pars.tau > pars.len*(j+1)) {
					j++;
					if (j == pars.nblock) l = 1;
				}
			}
			
			// ottimizzazione dei parametri per le mosse
			if (i % 1000==0 && i < 20000) calc_acceptance_and_adjust();
			
			// visualizzazione su terminale
			if (i % pars.outsteps == 0)
			{
				cout << "step=" << i << " l=" << pars.L[1]<< " \tphi=" << pars.vhb*pars.Np/(pars.L[0]*pars.L[1]*pars.L[2]) << " \tV=" << (pars.L[0]*pars.L[1]*pars.L[2]) << " \tvmax=" << pars.vmax << " \tM=" << llc.M << " \tcell=" << llc.cell_size << " \ta_r=" << (state.tot_rot-state.rot_rej)*1./state.tot_rot << " \ta_t=" << (state.tot_tra-state.tra_rej)*1./state.tot_tra << " \ta_v=" << (state.tot_vol-state.vol_rej)*1./state.tot_vol << " " << state.tot_vol << " " << state.vol_rej << "\n";
			}
			
			f << pars.MSD/pars.Np/i << "\n";
			g << pars.vhb*pars.Np/(pars.L[0]*pars.L[1]*pars.L[2]) << " \t" << pars.L[0]*pars.L[1]*pars.L[2] << "\n";
		}
		
		f.close();
		g.close();
		save_mgl();
		save_mgl_tot();
		
		// calcolo di phi medio ed errore
		calc_mean();
		calc_std();
		cout << "FINE SIMULAZIONE:\n";
		cout << "N = " << pars.Np << ", P = " << pars.P << ", T = " << pars.T << ", phi = " << pars.stat[0] << "+/-" << pars.stat[1] << "\n";
		
		free(block);
		free(llc.HEAD);
		free(llc.LIST);
		free(llc.MAP);
	}

	// calcolo della media
	void calc_mean() 
	{
		double m = 0;
		for (auto i = 0; i < pars.nblock; i++) m += block[i]/pars.len;
		pars.stat[0] = m*1./pars.nblock;
		}

	// calcolo dell'errore
	void calc_std() 
	{
		double s = 0;
		for (auto i = 0; i < pars.nblock; i++) s += (block[i]/pars.len-pars.stat[0])*(block[i]/pars.len-pars.stat[0]);
		pars.stat[1] = sqrt(s*1./pars.nblock);
	}

	// calcolo del MSD
	void calc_MSD() 
	{
		for (auto i=0; i < pars.Np; i++) 
		{
			pars.MSD += (in_parts[i].r-no_bound_parts[i].r).norm()*(in_parts[i].r-no_bound_parts[i].r).norm();
		}
	}

	// mossa Monte Carlo per la particella
	int mcmotion(int i)
	{	
		// si eseguoe con la stessa probabilità o una traslazione o una rotazione
		if (rng.ranf()<0.5)
		{
			// tralazione
			parts[i].tra_move(pars.delta_tra);
			return 0;
		}
		else 
		{
			// rotazione
			parts[i].rot_move(pars.delta_rot);
			return 1;
		}	
	} 

	// aggiornamento delle pbc  
	void pbc(int np)
	{
		pvector<ntype,3> delr;
		delr = parts[np].r;
		auto shift = delr.rint(pars.L);
		parts[np].r = parts[np].r - shift;
	}

	// definizione dello shifting della box per il calcolo dell'overlap
	// al bordo zona
	pvector<ntype,3> shifting(pvector<ntype,3> a,pvector<ntype,3> b) 
	{
		pvector<ntype,3> s;
		auto r = a-b;
		s[0] = std::floor(std::abs(r[0])/(pars.L[0]-2.*pars.sax.norm()))*sgn(r[0]);
		s[1] = std::floor(std::abs(r[1])/(pars.L[0]-2.*pars.sax.norm()))*sgn(r[1]);
		s[2] = std::floor(std::abs(r[2])/(pars.L[0]-2.*pars.sax.norm()))*sgn(r[2]);
		if (s[0]!=0) s[0] = s[0]/std::abs(s[0]);
		if (s[1]!=0) s[1] = s[1]/std::abs(s[1]);
		if (s[2]!=0) s[2] = s[2]/std::abs(s[2]);;
		return s*pars.L[0];
	}
		

	// controllo dell'overlap tra particelle
	bool check_all_overlaps(int np)
	{
		int j, icell, jcell, jcell0, nabor;
		double halfL, invcell_size;
		pvector<ntype,3> r;

		invcell_size = 1./llc.cell_size;
		halfL = pars.L[0]/2.;
			
		r = parts[np].r-parts[np].r.rint(pars.L);
		icell = (int)((r[0]+halfL)*invcell_size) + (int)((r[1]+halfL)*invcell_size)*llc.M + (int)((r[2]+halfL)*invcell_size)*llc.M*llc.M;
		j = llc.HEAD[icell];
		
		// controllo dell'overlap nella stessa cella 
		while (j>-1) {
			if (np!=j) 
			{	
				parts[np].shift = shifting(parts[np].r,parts[j].r);
				if (parts[np].my_overlap(parts[j]))  return true;
			}
			j = llc.LIST[j];
		}
		
		// controllo dell'overlap nella 26 celle vicine
		jcell0 = 26*icell;
		for (nabor = 0; nabor < 26; nabor++) {
			jcell = llc.MAP[jcell0 + nabor];
			j = llc.HEAD[jcell];
			while (j>-1) {
				if (np!=j) {
					parts[np].shift = shifting(parts[np].r,parts[j].r);
					if (parts[np].my_overlap(parts[j]))  return true;
					}
				j = llc.LIST[j];
			}
		}
		
		return false;
	}
	  
		// mossa di volume
	bool move_box(void)
	{
		state.tot_vol++;
		bool rejected = false;
		ntype Vo, Vn, fact;
		store_all_parts();
		Vo = pars.L[0]*pars.L[1]*pars.L[1];
		
		// nuovo volume
		Vn = Vo + 2.*pars.vmax*(rng.ranf()-0.5); 
		// fattore di dilatazione
		fact = pow(Vn/Vo,1.0/3.0);
		// nuovo lunghezza della lato della box
		pars.L = pars.L*fact;
		// ricollocazione delle particelle a seguito della mossa di volume
		for (int i=0; i < pars.Np; i++)
		{
			parts[i].r = parts[i].r*fact;
			no_bound_parts[i].r = no_bound_parts[i].r*fact;
		}
		// viene ridefinito il lato della cella
		llc.cell_size = pars.L[0]/llc.M;
		
		// controllo di overlap
		for (int i=0; i < pars.Np; i++)
		{
			if (check_all_overlaps(i)==true)
			{
				rejected=true;
				break;
			}
		}

		// se non ci sono overlap si procede a valutare la mossa MC
		if (!rejected)
		{
			ntype beta=1.0/pars.T, DG;
			// variazione energia libera di Gibbs
			DG = -(beta*pars.P*(Vn-Vo)-pars.Np*log(Vn/Vo));
			if (rng.ranf() > exp(DG))
			{
				rejected=true;
			}
		}
		// se la mossa è rifiutata si ripristina il vecchio volume e le vecchie posizioni
		if (rejected)
		{
			restore_all_parts();
			pars.L = pars.L * (1.0/fact);  
			llc.cell_size = pars.L[0]/llc.M;
		}
		
		if (llc.cell_size < 2.*pars.sax.norm()) {
			cout << "ERRORE, cella troppo piccola\n";
			exit(0);
		} 

		return rejected;
	}	
			
	// ottimizzazione dell'accettanza
	void calc_acceptance_and_adjust(void)
	{
		ntype r_tra = 0, r_rot = 0, r_vol = 0;
		// traslazione
		if (state.tot_tra > 0) {
			r_tra = (ntype)((state.tot_tra - state.tra_rej))/ state.tot_tra;
		}
		if (r_tra > 0.5) 
		{
			pars.delta_tra *= 1.1;
		}
		else
		{
			pars.delta_tra /= 1.1;
		}
		// rotazione
		if (state.tot_rot > 0) {
			r_rot = (ntype)((state.tot_rot - state.rot_rej))/ state.tot_rot;
		}
		if (r_rot > 0.5) 
		{
			pars.delta_rot *= 1.1;
		}
		else
		{
			pars.delta_rot /= 1.1;
		}
		// volume
		if (state.tot_vol > 0) {
			r_vol = (ntype)((state.tot_vol - state.vol_rej))/ state.tot_vol;
		}		
		if (r_vol > 0.5) 
		{
			pars.vmax *= 1.1;
		}
		else
		{
			pars.vmax /= 1.1;
		}
	}
	 
	// viene salvata la configurazione finale
	void save_mgl()
	{
		ofstream f;
		f.open("cnf-final.mgl");
		for (auto i = 0; i < pars.Np; i++) {
			for (auto k = 0; k < 3; k++) {
			f << parts[i].r[k] << " ";
			}
			for (auto k1 = 0; k1 < 3; k1++) {
				for (auto k2 = 0; k2 < 3; k2++) {
					f << parts[i].R[k1][k2] << " ";
				}
			}
			f << " B " << 2.*pars.sax[0] << " " <<  2.*pars.sax[1] << " " <<  2.*pars.sax[2] << " C[red]\n";
		}
		f.close();
	}
	
	// viene salvata la configurazione iniziale e finale
	void save_mgl_tot(){
		 ofstream f;
		 f.open("cnf-tot.mgl");
		 for (auto i = 0; i < pars.Np; i++) {
			for (auto k = 0; k < 3; k++) {
				f << parts[i].r[k] << " ";
			}
			for (auto k1 = 0; k1 < 3; k1++) {
				for (auto k2 = 0; k2 < 3; k2++) {
					f << parts[i].R[k1][k2] << " ";
				}
			}
			f << " B " << 2.*pars.sax[0] << " " <<  2.*pars.sax[1] << " " <<  2.*pars.sax[2] << " C[red]\n";
			for (auto k = 0; k < 3; k++) {
				f << in_parts[i].r[k] << " ";
			}
			for (auto k1 = 0; k1 < 3; k1++) {
				for (auto k2 = 0; k2 < 3; k2++) {
					f << in_parts[i].R[k1][k2] << " ";
				}
			}
			 f << " B " << 2.*pars.sax[0] << " " <<  2.*pars.sax[1] << " " <<  2.*pars.sax[2] << " C[blue]\n";
		}
		f.close();
	}
	
	// si salva la configurazione attuale per la singola particella
	void store_part(int ip)
	{
		parts[ip].store();
		no_bound_parts[ip].store();
	}
	
	// si ripristina la configurazione della particella  
	void restore_part(int ip)
	  {
		  parts[ip].restore();
		  no_bound_parts[ip].store();
		  }
	  
	// si salva la configurazione attuale di tutte le particelle
	void store_all_parts()
	{
		for (int i=0;  i < pars.Np; i++) {
				store_part(i);
			}
		} 

	// si ripristina la configurazione di tutte le particelle  
	void restore_all_parts(void)
	{
		for (int i=0; i < pars.Np; i++)
		{
			restore_part(i);
		}
	}
    
	// inizializzazione box di simulazione
	void set_box(pvector<ntype> Lt)
	{
		pars.L = Lt;
	}
    
	// numero di particelle nella box
	void set_num_parts(int Np)
	{
		parts.resize(Np);
		in_parts.resize(Np);
		no_bound_parts.resize(Np);
	}
  
	// inizializzazione dei semiassi
	void set_X0(pvector<ntype,3> s)
	{
		pars.sax = s; 
		for (auto& p: parts)
		{
			p.set_semiaxes(s);
		}
	}
  
	// numero di step di simulazione
	void set_steps(long int Nt)
	{
		pars.Nsteps = Nt;
	}
	  
	// numero di particelle per lato nella configurazione iniziale
	void set_xyz(int x, int y, int z)
	{
		pars.nx = x;
		pars.ny = y;
		pars.nz = z;
	}
	
	// inizializzazione delta mossa di traslazione
	void set_delta_tra(ntype dt)
	{
		pars.delta_tra = dt;
	}
	
	// inizializzazione delta mossa di rotazione
	void set_delta_rot(ntype dr)
	{
		pars.delta_rot = dr;
	}
	
	// inizializzazione temperatura
	void set_T(ntype T)
	{
		pars.T = T;
	}
	
	// inizializzazione pressione	
	void set_P(ntype P)
	{
		pars.P = P;
	}
	
	// inizializzazione volume escluso
	void set_v_excl(ntype V)
	{
		pars.v_excl = V;
	}
	
	// rate di visualizzazione
	void set_outsteps(long int Nt)
	{
		pars.outsteps = Nt;
	}
	
	// numero di passi per bloccp
	void set_len(long int N)
	{
		pars.len = N;
	}
		 
	int sgn(ntype x) {
		if (x<0) return -1;
		if (x>0) return 1;
		return 0;
	}
};

