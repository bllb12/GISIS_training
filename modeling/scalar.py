import numpy as np
import matplotlib.pyplot as plt

class Wavefield_1D():
    
    def __init__(self):
        
        self._type = "1D wave propagation in constant density acoustic isotropic media"

        # parameterFile = 'C:\\Users\\Gamer\\OneDrive\\Área de Trabalho\\GISIS_training\\modeling\\data_test.txt' #Caminho que identifica o arquivo no qual os dados estão contidos.
        # print('Loading parameters from disk %s'%parameterFile)

        # with open(parameterFile,'r') as archive:
        #     linhas = [linha.strip() for linha in archive.readlines()]

        self.nt = 1001 #  - número de amostras
        self.dt = 0.001 #  - intervalo de tempo
        self.fmax = 30.0 #  - frequência máxima

        self.nx = 2000 #  - número de pontos no eixo x
        self.dz = 1.0 #  - espaçamento entre os pontos no eixo z
        self.num_steps = 1000 # - número de passos no tempo
        self.rho = 1000 # - densindade do meio

        self.model = np.zeros(self.nx)
        self.p = np.zeros(self.nx)
        self.v = np.zeros(self.nx)


        self.prof = [0, 1250, 1500, 1750, 2000] #  - profundidade das interfaces
        self.velocities = [1500, 1800, 2100, 2400, 2700] #  - velocidades das camadas

        self.source_idxs = list([int(100 * 0.6), int(100 * 0.7), int(100 * 0.8)])
        self.source_times = list([0.01, 0.04, 0.07])
        self.source_frequencies = list([15.0, 20.0, 25.0])
        self.receiver_idxs = list([int(self.nx * 0.7), int(self.nx * 0.8), int(self.nx * 0.9)])
        self.receiver_depths = list([1250, 1500, 1750, 2000])
    
    def get_type(self):
        print(self._type)

    def set_ricker(self): # Wavelet Ricker

        t0 = 2.0*np.pi/self.fmax
        fc = self.fmax/(3.0*np.sqrt(np.pi))
        td = np.arange(self.nt)*self.dt - t0

        arg = np.pi*(np.pi*fc*td)**2.0

        self.wavelet = (1.0 - 2.0*arg)*np.exp(-arg)

    def plot_ricker(self):

        t = np.arange(self.nt)*self.dt

        fig, ax = plt.subplots(figsize = (10, 5), clear = True)

        ax.plot(t, self.wavelet)
        ax.set_title("Wavelet", fontsize = 18)
        ax.set_xlabel("Time [s]", fontsize = 15)
        ax.set_ylabel("Amplitude", fontsize = 15) 
        
        ax.set_xlim([0, np.max(t)])
        
        fig.tight_layout()
        plt.show()

    def set_model(self):
        interfaces = []
        for i in range(len(self.prof)):
            start_depth = int(self.prof[i - 1] / self.dz)
            end_depth = int(self.prof[i] / self.dz)
            vel = self.velocities[i - 1]
            interfaces.append((start_depth, end_depth, vel))

        for j in range(len(interfaces)):
            start, end, vel = interfaces[j]
            start_depth = int(start)
            end_depth = int(end)
            self.model[start_depth:end_depth] = vel
            
    def propagate(self):
        for n in range(1, self.num_steps + 1):
            # Fontes pontuais na superfície
            for idx, time, frequency in zip(self.source_idxs, self.source_times, self.source_frequencies):
                if n * self.dt < time:
                    self.p[idx] += np.sin(2.0 * np.pi * frequency * n * self.dt)

            # Atualização dos campos de pressão e velocidade
            self.v[1:self.nx-1] += ((self.model[1:self.nx-1] * self.dt**2 / self.dz**2) * (self.p[2:self.nx] - 2 * self.p[1:self.nx-1] + self.p[0:self.nx-2]))
            self.p[1:self.nx-1] += self.v[1:self.nx-1] * self.dt

    def plot_model(self):
        fig, ax = plt.subplots(figsize=(8, 10), clear=True)

        # Adicionando cada camada com uma cor distinta
        colors = ['green', 'orange', 'blue', 'red']
        for i in range(len(self.prof) - 1):
            ax.axhspan(self.prof[i], self.prof[i + 1], facecolor=colors[i], alpha=0.3)

        ax.plot(self.model, np.arange(self.nx) * self.dz, label='Modelo de Velocidade', color='black')

        # Adicionando fontes
        for idx, time, frequency in zip(self.source_idxs, self.source_times, self.source_frequencies):
            ax.plot(self.model[idx], idx * self.dz, '*', label=f'Fonte {frequency} Hz', markersize=10)

        # Adicionando receptores
        for idx, depth in zip(self.receiver_idxs, self.receiver_depths):
            ax.plot(self.model[idx], depth, 'v', label=f'Receptor a {depth}m', markersize=10)

        ax.set_title("Modelo de Velocidade com Propagação de Onda", fontsize=18)
        ax.set_xlabel("Velocidade [m/s]", fontsize=15)
        ax.set_ylabel("Profundidade [m]", fontsize=15)
        ax.set_ylim(max(np.arange(self.nx) * self.dz), min(np.arange(self.nx) * self.dz))
        ax.legend()

        fig.tight_layout()
        plt.show()
        

class Wavefield_2D(Wavefield_1D):
    
    def __init__(self):
        super().__init__()
        
        self._type = "2D wave propagation in constant density acoustic isotropic media"

class Wavefield_3D(Wavefield_2D):
    
    def __init__(self):
        super().__init__()
        
        self._type = "3D wave propagation in constant density acoustic isotropic media"    
