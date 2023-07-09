import matplotlib.pyplot as plt
import numpy as np

# Definisci le tue funzioni di complessità qui. 
# Per esempio, ho definito due funzioni: O(n) e O(n^2).
def func1(m):
    return (m*m/2)/33000000

def func2(m):
    return (m*np.log2(m/2)**5)/33000000

# Crea un array di valori per n
n = np.array(range(1, 15000000))

# Calcola i valori y per ciascuna funzione
y1 = func1(n)
y2 = func2(n)

# Crea il grafico
plt.figure(figsize=(10, 6))

# Disegna le funzioni
plt.plot(n, y1, label='Bellman-Ford: O(mn)')
plt.plot(n, y2, label='SPmain: O(mlog^5(n))')

# Aggiungi titolo e etichette
plt.title('Confronto complessità Bellman-Ford - SPmain')
plt.xlabel('Numero di archi (m)')
plt.ylabel('Tempo di esecuzione in secondi')

plt.axvline(x=1200000, color='r', linestyle='--', label='Limite di archi per i test')

# Mostra la legenda
plt.legend()

#plt.savefig('pltSparseFunctionsFar.pdf', format='pdf')
plt.savefig('pltSparseFunctionsFar.png')

# Mostra il grafico
#plt.show()
