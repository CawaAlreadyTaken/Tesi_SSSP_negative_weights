import matplotlib.pyplot as plt

# Definisci i tuoi punti qui. Per esempio, ho definito due serie di punti.
x1 = [20, 100, 200, 600, 1200, 2000, 5000, 10000, 20000, 60000, 120000, 200000, 300000, 400000, 600000, 800000, 1200000]
y1 = [0.004, 0.005, 0.004, 0.007, 0.013, 0.042, 0.173, 0.714, 2.62, 24.5, 121.4, 387.1, 958.9, 1800.0, 4781, 12624, 41220]

x2 = [20, 100, 200, 600, 1200, 2000, 5000, 10000, 20000, 60000, 120000, 200000, 300000, 400000, 600000, 800000, 1200000]
y2 = [0.016, 0.075, 0.203, 0.863, 2.343, 4.545, 16.498, 48.778, 130.16, 787.8, 2654.6, 6440.3, 11973.2, 24641.5, 54332, 90746, 175843]

# Crea il grafico
plt.figure(figsize=(10, 6))

# Disegna i punti uniti da segmenti
plt.plot(x1, y1, '-o', label='Bellman-Ford')
plt.plot(x2, y2, '-o', label='SPmain')

# Aggiungi titolo e etichette
plt.title('Confronto implementazione Bellman-Ford - SPmain')
plt.xlabel('Numero di archi (m)')
plt.ylabel('Tempo di esecuzione in secondi')

# Mostra la legenda
plt.legend()

plt.savefig('pltSparsePunti.png')

# Mostra il grafico
#plt.show()
