import numpy as np
import matplotlib.pyplot as plt

m = 1.0 #9.1*10**(-31)  (on utilise pas car perte précision)# masse de la particule
hbar = 1.0 # 1.054571628*10**(-34) (on utilise pas car perte précision)

#fonction permettant de créer un paquet d'onde de gauss
def GaussWP(k0, a, ux, t):
    x = ux #-1/k0*2
    #j'ai séparé le diviseur
    diver = m * a ** 2 + 2j * hbar * t

    # facteur d'amplitude
    amp = (1 / (8 * np.pi ** 3)) ** 0.25 * np.sqrt((4 * np.pi * m * a) / diver)

    # Terme exponentiel
    partexp = (m / 4) * ((a ** 2 * k0 + 2j * x) ** 2 / diver) - (a ** 2 * k0 ** 2) / 4

    return amp * np.exp(partexp)

# Params
k0 = 0.821
x_values = np.linspace(-2/k0, 2/k0, 500)
#x_values = np.linspace(0, 10, 500)

psi = GaussWP(0.821,1, ux=x_values, t=0.0)

# Extract composantes
partie_reelle = np.real(psi)
partie_imaginaire = np.imag(psi)
densite_proba = np.abs(psi)**2

# Visualisation
fig, axs = plt.subplots(2, 1, figsize=(8, 8))

axs[0].plot(x_values, densite_proba, label="Densité de probabilité", color="red")
axs[0].set_title("Densité de probabilité")
axs[0].grid(True)
axs[0].legend()

axs[1].plot(x_values, partie_reelle, label="Partie réelle", color="orange")
axs[1].plot(x_values, partie_imaginaire, label="Partie imaginaire", color="blue")
axs[1].set_title("Paquet d'onde")
axs[1].set_xlabel("Position x")
axs[1].grid(True)
axs[1].legend()

plt.tight_layout()
plt.show()
