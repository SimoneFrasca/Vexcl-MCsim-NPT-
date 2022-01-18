Benvenuto: Questo programma è stato progettato con lo scopo di:
1) calcolare il volume escluso di una hard box
2) ricavare l'equazione di stato di un sistema di N hard box

#1 Esecuzione del programma

 a) andare nella cartella dove sono contenuti tutti i file del programma

 b) per prima cosa da terminale eseguire il comando "make clean"

 c) poi eseguire il comando "make"

 d) Nel caso si voglia calcolare il volume escluso eseguire
	"./vexcl.exe N_STEP N_ouptut maj_sax"
	Dove N_step è il numero di simulazioni, N_output il rate di output
	che si vuole visualizzare e maj_sax la lunghezza del semiasse maggiore.
	Alla fine della simulazione in output si avrà il valore del volume escluso
	
 e) Nel caso si voglia eseguire una simulazione MC di un sistema di hard box
	in un ensamble NPT
	"./mcsim.exe N_step N_output Len_BLOCK nx ny nz maj_sax delta_tra delta_rot T P v_excl"
	Dove:
	* N_step è il numero di step della simulazione MC
	* N_output è il rate di visualizzazione dell'output
	* Len_BLOCK è il numero di elemnti per blocco per fare la media a blocchi
	* nx, ny, nz è il numero di hb per asse
	* maj_sax è la lunghezza del semiasse maggiore
	* delta_tra e delta_rot gli incrementi massimi per mosse di traslazione
	  e rotazione
	* T e P sono la temperatura e la pressione
	* v_excl è il volume escluso
 Alla fine della simulazione si avrà in output un breve resoconto delle quantità 
 fissate all'inizio e del valore della frazione di volume

NOTA: l'algoritmo nella versione attuale permette di simulre solo e esclusivamente
particelle a forma di parallelepipedo con base di semiassi di lunghezza 
1x1 e il semiasse maggiore di lunghezza pari ad un numero intero maggiore di 1
(eventuali aggiornamenti andranno ad eliminare questa limitazione)

#2 Consigli 
 a) Calcolo del volume escluso: per avere una stima il più possibile corretta 
    è consigliato eseguire almento 10^5 simulazioni
 
 b) Simulazione MC: per la simulazione MC si consiglia:
   1) fare attenzione che il volume sia cubico, per cui la quantità
      ni x saxi, dove ni è il numero di particelle lungo la direzione i
      nella inizializzazione e saxi è la lunghezza del semiasse nella stessa sirezione,
      sia la stessa per tutte e tre le dimensioni 
   2) la quantità minima di passi della simulazione è 10^4, mentre un
      blocco deve contenere circa 10^3 elementi

#3 File prodotti e visualizzazione
  Il programma genererà i seguenti file:
  1) Un file "vexcl.dat" che registra ad ogni passo di simulazione il volume escluso
  2) Un file "phi_evolution.dat" che contiene ad ogni passo della simulazione MC il valore 
     della frazione di volume
  3) Un file "MSD.dat" che contiene il valore del MSD per ogni passo di simulazione MC

  Per i file che permettono la visualizzazione abbiamo:
  1) "cnf-final.mgl" che contiene la configurazione finale delle hard box
  2) "cnf-tot.mgl" che contiene la configurazione finale e iniziale delle hard box
 
  Questi salvano la posizione, l'orientazione e la lunghezza dei semiassi 
  di ogni singola hard box
