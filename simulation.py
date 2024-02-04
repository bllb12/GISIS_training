from modeling import scalar

def simulation():

    id = 0

    myWave = [scalar.Wavefield_1D(), 
              scalar.Wavefield_2D(),
              scalar.Wavefield_3D()] 

    # print(myWave[id]._type)
    myWave[id].get_type()

    myWave[id].set_ricker()
    myWave[id].plot_ricker()
    myWave[id].set_model()
    myWave[id].propagate()
    myWave[id].plot_model()

if __name__ == "__main__":
    simulation()
