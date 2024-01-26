import numpy as np
import matplotlib.pyplot as plt

class Wavefield_1D():
    
    def __init__(self):
        
        self._type = "1D wave propagation in constant density acoustic isotropic media"

        # TODO: read parameters from a file

        # self.nt = 1001
        # self.dt = 1e-3
        # self.fmax = 30.0
        parameterFile = 'C:\\Users\\Gamer\\OneDrive\\Área de Trabalho\\GISIS_training\\modeling\\data.txt'
        print('Loading parameters from disk %s'%parameterFile)
        
        with open(parameterFile,'r') as archive:
            linhas = [linha.strip() for linha in archive.readlines()]

        self.nt = int(linhas[0])
        self.dt = float(linhas[1])
        self.fmax = float(linhas[2])

    def get_type(self):
        print(self._type)

    def set_wavelet(self):
    
        t0 = 2.0*np.pi/self.fmax
        fc = self.fmax/(3.0*np.sqrt(np.pi))
        td = np.arange(self.nt)*self.dt - t0

        arg = np.pi*(np.pi*fc*td)**2.0

        self.wavelet = (1.0 - 2.0*arg)*np.exp(-arg)

    def plot_wavelet(self):
       
        t = np.arange(self.nt)*self.dt

        fig, ax = plt.subplots(figsize = (10, 5), clear = True)

        ax.plot(t, self.wavelet)
        ax.set_title("Wavelet", fontsize = 18)
        ax.set_xlabel("Time [s]", fontsize = 15)
        ax.set_ylabel("Amplitude", fontsize = 15) 
        
        ax.set_xlim([0, np.max(t)])
        
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
