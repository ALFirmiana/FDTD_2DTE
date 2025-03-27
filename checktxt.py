import numpy as np
import matplotlib.pyplot as plt

file_path = "build/"
# |%%--%%| <1hQv8kKawq|sfh9GANB74>

Bz = []
for i in range(0, 1001, 10):
    file = file_path + f"Bz_{i}.txt"
    Bz.append(np.loadtxt(file))
np.shape(Bz)

# |%%--%%| <sfh9GANB74|kVRQavqLlQ>

for i in range(0, 102, 10):
    # plt.imshow(Bz[i], vmax=1, vmin=-1)
    plt.imshow(Bz[i], vmax=0.1, vmin=-0.1)
    plt.colorbar()
    plt.show()

# |%%--%%| <kVRQavqLlQ|NuGrbXSxiU>

Ex = []
for i in range(0, 1001, 10):
    file = file_path + f"Ex_{i}.txt"
    Ex.append(np.loadtxt(file))

for i in range(0, 102, 10):
    plt.imshow(Ex[i], vmax=1, vmin=-1)
    plt.colorbar()
    plt.show()

# |%%--%%| <NuGrbXSxiU|gVDk4Ghvix>


Ey = []
for i in range(0, 1001, 10):
    file = file_path + f"Ey_{i}.txt"
    Ey.append(np.loadtxt(file))

for i in range(0, 102, 10):
    plt.imshow(Ey[i], vmax=1, vmin=-1)
    plt.colorbar()
    plt.show()

