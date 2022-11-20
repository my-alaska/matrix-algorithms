import numpy as np
import matplotlib.pyplot as plt

sizes = []
times = []
operations = []

oper='det'

with open(f'{oper}.txt', 'r') as file:
    for line in file.readlines():
        line = line.strip()
        if line == '':
            continue
        items = line.split(' ')
        items = [float(item) for item in items]
        sizes.append(items[0])
        times.append(items[1])
        operations.append(items[2])



fig = plt.figure(figsize=(10, 5))

ax = fig.add_subplot(121)
ax.scatter(sizes, times, marker='+', color='red', label='Pomiary')
x = np.arange(0, 129)
y = x**2.81*0.0000106
ax.plot(x, y, label='x^2.81 * 0.0000106')
ax.set_title("Czas obliczenia wyznacznika macierzy\nw zależności od rozmiaru macierzy")
ax.set_xlabel('Rozmiar macierzy')
ax.set_ylabel('Czas [s]')
ax.legend()
# plt.show()

ax = fig.add_subplot(122)
ax.scatter(sizes, operations, marker='+', color='red', label='Pomiary')
x = np.arange(0, 129)
y = x**2.81*9
ax.plot(x, y, label='x^2.81 * 9')
ax.set_title("Operacje potrzebne do obliczenia wyznacznika\nw zależności od rozmiaru macierzy")
ax.set_xlabel('Rozmiar macierzy')
ax.set_ylabel('Liczba operacji')
ax.legend()

plt.show()
fig.savefig('det.png')
