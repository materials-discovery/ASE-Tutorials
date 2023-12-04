import matplotlib.pyplot as plt 

x = [1,2,3,4,5,6,7,8]
y = [-198.53484512933608, -213.10376807721624, -215.065959198655, -215.50332981238003, -215.62957893262117, -215.6711784717626, -215.686171127931, -215.691957628]

plt.figure(figsize=(8, 6))
plt.plot(x, y, color='blue', marker='o')
plt.xlabel('Kgrid')
plt.ylabel('energy / eV')
plt.grid()
plt.savefig('K-grid.png')
