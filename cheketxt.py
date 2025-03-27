import numpy as np
import matplotlib.pyplot as plt

file_path = "build/"
# |%%--%%| <1hQv8kKawq|sfh9GANB74>

Bz = []
for i in range(101):
    file = file_path + f"Bz_{i}.txt"
    Bz.append(np.loadtxt(file))
np.shape(Bz)

# |%%--%%| <sfh9GANB74|kVRQavqLlQ>

for i in range(20):
    plt.imshow(Bz[i], vmax=1, vmin=-1)
    plt.colorbar()
    plt.show()

# |%%--%%| <kVRQavqLlQ|NuGrbXSxiU>

Ex = []


