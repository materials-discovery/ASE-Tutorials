import matplotlib.pyplot as plt 

x = [1,2,3,4,5,6,7,8]
y = []

plt.figure(figsize=(8, 6))
plt.plot(x, y, color='blue', marker='o')
plt.xlabel('Kgrid')
plt.ylabel('energy / eV')
plt.grid()
plt.savefig('K-grid.png')
