# Tesi NWSSSP
### Author: _Daniele Cabassi_

Algoritmo (1): NWSSSP (Negative-Weight Single-Source Shortesh Path)  
Algoritmo (2): MFMCF (Maximum Flow and Minimum-Cost Flow)

(1) E' un algoritmo randomizzato che calcola i percorsi minimi da un'unica sorgente in O(m log<sup>8</sup>(n) log W) in un grafo con archi con pesi interi e che possono anche essere negativi.

(2) E' un algoritmo che calcola esattamente il flusso massimo ed il flusso di costo minimo su grafi direzionati con costi, capacita' minima e massima interi e limitati polinomialmente in tempo O(m<sup>1+o(1)</sup>)

**Near Linear Time (1)**:
_Una funzione f:N->N e' "near linear" se f(n) appartiene agli O(n<sup>1+eps</sup>) per ogni eps>0.
Quindi n log<sup>k</sup> (n) e' near linear._

**Almost Linear Time (2)**:
_Praticamente la stessa cosa_

Per NWSSSP, i limiti precedenti erano di O((m+n<sup>1.5</sup>)log W) su grafi moderatamente densi e O(m<sup>4/3+o(1)</sup> log W) su grafi sparsi, da algoritmi trovati entrambi nel 2020.

Prima di questo (1) algoritmo, esistevano gia' (dal 2001) degli algoritmi near-linear ma solo per il caso particolare di grafi planari direzionati, in O(n log<sup>2</sup> (n)/(log log n)).

**Grafo planare**:
_Nella teoria dei grafi si definisce grafo planare un grafo che può essere raffigurato in un piano in modo che non si abbiano archi che si intersecano
(Curva di Jordan?)_

I recenti algoritmi precedenti si basavano su metodi sofisticati di ottimizzazione continua. Questo algoritmo usa invece la "graph decomposition" e strumenti di combinatoria. E' il primo algoritmo combinatorio per NWSSSP che va sotto il limite O(m sqrt(n) log W) di Gabow e Tarjan del 1989.

**Algoritmo combinatorio**:
_Algoritmo che risolve un problema combinatorio, che consiste nel trovare un oggetto ottimo da un insieme finito di oggetti, dove l'insieme di soluzioni possibile e' discreto o puo' essere ridotto ad un insieme discreto._

Notazione:

**m**: numero di archi  
**n**: numero di vertici  
**G(V, E, w)**: grafo con insieme V di vertici, E di archi e pesi interi w per ogni arco e in E.  
**s**: vertice di partenza appartenente a V  
**dist<sub>G</sub>(s,v)**: distanza minima nel grafo da s a v con s e v appartenenti a V  
**W**: minimo intero >= 2 tale che w(e) >= -W per ogni e in E. Sostanzialmente il valore assoluto del peso minore tra gli archi.  

Problema: trovare la distanza minima da s a v per ogni v in V.  
Dijkstra funziona solo per archi con pesi non negativi.   
Bellman-Ford fornisce una soluzione per questo problema: se c'e' un ciclo negativo lo trova, altrimenti ritorna dist<sub>G</sub>(s,v) per ogni vertice in V. Tuttavia, esegue in O(mn).

Contesto:  
dagli anni '50 ci sono stati diversi miglioramenti:  
* Shimbel, Ford, Bellman negli anni '50, O(mn)  
* Gabow, Tarhan, Goldberg negli anni '80 e '90, O(m sqrt(n) logW)  
* Negli ultimi anni ci son stati miglioramenti in algoritmi di ottimizzazione continua che han portato ad algoritmi piu' veloci per i problemi di trasbordo (trasferimento di passeggeri o carichi in cui il percorso puo' includere tappe intermedie) ed i problemi di flusso con costo minimo. Questo ha portato ai limiti di tempo descritti sopra anche per negative-weights SSSP (1).  

La ricerca per (1) ha fondamentalmente seguito due domande:  
1-1. Si puo' ottenere un tempo near-linear per la risoluzione del Problema per tutti i tipi di grafi?  
1-2. Possiamo ottenere algoritmi efficienti senza strutture dati complesse?  

La seconda domanda riguarda il fatto che ci sia un modo per risolvere il Problema senza utilizzare metodi sofisticati di ottimizzazione continua e una serie di algoritmi complessi algebrici e sui grafi.

La tesi punta a risolvere entrambe le domande, presentando un algoritmo combinatorio che riduce il tempo di esecuzione ad essere near-linear.

**Teorema 1.1)** _Esiste un algoritmo randomizzato (Las Vegas) che esegue in O(m log<sup>8</sup>(n) log(W)) con alta probabilita' per un grafo G<sub>in</sub> con m archi e nodo di partenza s<sub>in</sub>. Tale algoritmo, o ritorna l'albero dei percorsi minimi da s<sub>in</sub>, oppure ritorna un ciclo negativo._

L'algoritmo utilizzato di decomposizione del grafo e' chiamato "Low Diameter Decomposition" e viene studiato fin dagli anni '80. Tale algoritmo funziona con grafi con pesi _non_ negativi, ma puo' essere utilizzato per sviluppare un algoritmo di riduzione scalare per SSSP anche con pesi negativi.

Per il futuro: La domanda 1-2) riguarda un'ampia gamma di problemi. Il panorama attuale degli algoritmi sui grafi e' che per i problemi fondamentali, compresi quelli che si insegnano in universita' e quelli che si utilizzano di piu' in pratica, lo state-of-the-art delle soluzioni e' un algoritmo complesso per un problema piu' generale di flusso di costo minimo (ad esempio: SSSP con anche archi negativi, bipartite matching, il problema dell'assegnazione, taglio s-t, flusso massimo...)  
Questo suggerisce l'idea di trovare degli algoritmi semplici piu' specifici per tali problemi fondamentali. Il lavoro della tesi su SSSP con anche archi negativi e' un passo verso questa direzione.

Indipendentemente dall'algoritmo in questione (1), e' stato trovato un algoritmo di ottimizzazione continua nello stesso anno (2022), che e' il (2) e pareggia i limiti di tempo per SSSP con anche pesi negativi come un caso particolare del problema piu' generale che risolve. I due algoritmi sono completamente diversi, e per quanto si sappia non c'e' alcuna sovrapposizione nelle tecniche utilizzate.

**Weak Diameter**:
_In una decomposizione, il diametro della stessa e' la distanza massima tra due vertici della decomposizione. Se e' possibile scegliere anche archi non appartenenti alla decomposizione, allora si sta parlando di "weak diameter", altrimenti di "strong diameter"_

Una delle funzioni principali e' "LowDiamDecomposition":  
**Lemma 1.2**: _Esiste un algoritmo LowDiamDecomposition(G, D) tale che:_  
* in INPUT si abbia un grafo con m-archi, n-vertici, una funzione di pesi non-negativi degli archi e un intero positivo D.  
* in OUTPUT si abbia un insieme di archi E<sup>rem</sup> con le seguenti caratteristiche:  
    * ogni SCC di G \ E<sup>rem</sup> ha "weak diameter" <= D  
    * per ogni arco e in E, la probabilita' che e appartenga ad E<sup>rem</sup> e' O((w(e)*log<sup>2</sup>(n)/D)+n<sup>-10</sup>)  
* L'algoritmo esegua in O(m log<sup>2</sup>(n) + n log<sup>3</sup>(n))  

Algoritmo: https://raw.githubusercontent.com/CawaAlreadyTaken/Tesi_SSSP_negative_weights/main/LowDiameterDecomposition.png

Per comodita', si sceglie di concentrarsi su un algoritmo che ritorna correttamente l'albero delle distanze minime se non c'e' alcun ciclo negativo, altrimenti non garantisce nulla.  
In seguito si trattera' il caso in cui c'e' un ciclo negativo.  

Il paper non ha puntato ad ottimizzare i fattori di log, ma si potrebbe ad esempio trovare il ciclo negativo direttamente, invece che utilizzare direttamente la riduzione sopra, che in caso di ciclo negativo incontra un fattore O(log<sup>2</sup> n) ulteriore.  

Definiamo:
* V(G) = V  
* E(G) = E
* E<sup>neg</sup>(G) := {e in E | w(e) < 0}.  
* W<sub>G</sub> := max(2, -min<sub>e in E</sub>{w(e)})

Dato un qualsiasi insieme di archi S (sottoinsieme di E), definiamo:  
* $w(S) = \sum_{e \in S} w(e)$

Diciamo che un ciclo C nel grafo G e' di peso negativo se w(C) < 0.  
Se c'e' almeno un ciclo di peso negativo per qualche percorso da u a v, diciamo che dist<sub>G</sub>(u, v) = -inf

**Lemma 2.2**: 
_Supponiamo che esista un algoritmo per il Teorema 1.1 per il caso speciale in cui, nel grafo G<sub>in</sub>:_
* $w(e) >= -1 \forall e \in E$
* _Ogni vertice in G<sub>in</sub> abbia out-degree costante_  

_Con tempo di esecuzione T(m, n)_  
_Allora, esiste un algoritmo come nel Teorema 1.1  per ogni grafo in input G<sub>in</sub> con pesi interi con tempo di esecuzione O(T(m, m) log(W<sub>G<sub>in</sub></sub>))_

**Dimostrazione**
Si consideri il primo punto, ovvero che w<sub>in</sub> >= -1. Lo "scaling framework" di Goldberg dimostra che l'esistenza di un algoritmo per questo caso implica l'esistenza di un algoritmo per ogni grafo con pesi interi al costo di un fattore log(W<sub>G</sub>) extra.  
Per il secondo punto, creiamo Theta(out-degree(v)) copie di ogni vertice v, cosi' che ogni copia abbia out-degree costante.  
Il grafo finale ha O(E) vertici e O(E) archi.

**Definizione 2.3, (G<sub>s</sub>, w<sub>s</sub>, G<sup>B</sup>, w<sup>B</sup>, G<sup>B</sup><sub>s</sub>, w<sup>B</sup><sub>s</sub>)**:
_Dato un grafo qualsiasi G = (V, E, w), sia G<sub>s</sub> = (V U {s}, E U {s, v}<sub>$v \in V$</sub>, w<sub>s</sub>) il grafo con una finta sorgente aggiunta, in cui c'e' un arco di peso nullo da s ad ogni v in V e nessun arco entrante in s. G<sub>s</sub> ha un ciclo di peso negativo se e solo se G ce l'ha, e dist<sub>G<sub>s</sub></sub>(s, v) = min<sub>$u \in V$</sub>dist<sub>G</sub>(u, v)._  
_Per un qualsiasi intero B, sia G<sup>B</sup> = (V, E, w<sup>B</sup>) il grafo ottenuto aggiungendo B a tutti gli archi con pesi negativi in G._  
_Formalmente, w<sup>B</sup>(e) = w(e)+B per ogni e in E<sup>neg</sup>(G) e w<sup>B</sup>(e) = w(e) per ogni e in E \ E<sup>neg</sup>(G)._

**Definizione 2.4,** $(\eta <sub>G</sub>(v), P<sub>G</sub>(v) )$:
__

