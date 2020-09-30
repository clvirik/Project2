import matplotlib.pyplot as plt
import glob
import numpy as np
plt.rcParams["text.usetex"] =True
plt.rcParams.update({'font.size': 14})
#plt.rcParams["mathtext.fontset"] = "cm"

infile = open("Time.txt", 'r')
N = []; iterations = []; jacobi_time = []; arma_time = [];
for line in infile:
    words = line.split()
    N.append(float(words[0])*100)
    iterations.append(float(words[1]))
    jacobi_time.append(float(words[2]))
    arma_time.append(float(words[3]))

poly = np.polyfit(np.asarray(N), np.asarray(iterations), 2)

regression = np.zeros_like(iterations)
for i in range(len(poly)):
    regression += np.asarray(N)**(len(poly)-1-i)*poly[i]

plt.figure(1)
plt.plot(N, iterations, label='Iterations required')
plt.plot(N, regression, '--', label=f'Quadratic regression ({poly[0]:.2f}' +  r'$N^2$)')
plt.xlabel('N'); plt.ylabel('Number of iterations')
plt.legend()
plt.tight_layout()
plt.savefig("iteraitons.png")


plt.figure(2)
plt.subplot(1,2,1)
plt.title('Jacobi-method')
plt.plot(N, jacobi_time)
plt.tight_layout()
plt.xlabel('N'); plt.ylabel('Time [ms]')

plt.subplot(1,2,2)
plt.title('Armadillo')
plt.plot(N, arma_time)
plt.xlabel('N'); plt.ylabel('Time [ms]')
plt.tight_layout()
plt.savefig('Times.png')

infile = open("OneElectron.txt", 'r')
N = []; rho = []; error = []
for line in infile:
    words = line.split()
    error.append(np.log10(float(words[2])))

error_mat = np.zeros((5, 20))
for i in range(5):
    for j in range(20):
        error_mat[i][j] = error[20*i + j]

error_mat = np.asarray(error_mat)

rho_plot, N_plot = np.meshgrid(np.linspace(2, 5.8, 20), np.linspace(100, 300, 5))

plt.figure(3)
plt.contourf(N_plot, rho_plot, error_mat, levels=100)
plt.colorbar(label=r'$\log_{10}$(error)')
plt.xlabel('N'); plt.ylabel(r'$\rho_{max}$')
plt.tight_layout()
plt.savefig('OneElectron.png')