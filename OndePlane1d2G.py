from numpy import pi, exp, sqrt, real, imag, zeros, linspace
import matplotlib.pyplot as plt
PI = pi

# --- Question 2.a ---
# En Python, le nombre imaginaire 'i' s'écrit '1j'
def PlaneWave(amp, k, omega, x, t):
    return amp * exp(1j * (k * x - omega * t))

#def PlaneWaveButBetter(amp, k, omega, x, t): # je ne sais pas ce que j'avais tenté de faire ici ..
    return amp * exp(1j*1j * (k * x - omega * t))

# --- Question 2.b ---

# 1. On crée notre axe des abscisses x
dk = 0.8
x_values = linspace(-pi/dk, pi/dk, 500)

# 2. On calcule les fonctions d'onde à des instants t (ex: t=0)

psi =  PlaneWave(amp=1.0, k=0.8, omega=0.5, x=x_values, t=0.0)

psi2 = PlaneWave(amp=0.5, k=0.8-dk, omega=0.5, x=x_values, t=0.0)

psi3 =PlaneWave(amp=0.5, k=0.8+dk, omega=0.5, x=x_values, t=0.0)

#la somme de fonction pour un paquet d'onde
psi4 = (PlaneWave(amp=0.5, k=0.8-dk, omega=0.5, x=x_values, t=0.0)+
    PlaneWave(amp=1.0, k=0.8, omega=0.5, x=x_values, t=0.0)+
    PlaneWave(amp=0.5, k=0.8+dk, omega=0.5, x=x_values, t=0.0))
# 3. On extrait la partie réelle et la partie imaginaire
partie_reelle = real(psi)
partie_imaginaire = imag(psi)

# 4. On crée le graphique exactement comme demandé dans l'énoncé
fig, axs  = plt.subplots(2)

# On trace les deux courbes sur le même graphique 2D
axs[0].plot(x_values, real(psi4), label="Somme", color="blue")
axs[0].plot(x_values, real(psi2), label="-dk", color="orange",linestyle="--")
axs[0].plot(x_values, real(psi3), label="+dk", color="green",linestyle="--")
axs[0].plot(x_values, real(psi), label="0", color="red",linestyle="--")


axs[1].plot(x_values, real(psi), label="Partie réelle", color="blue")
axs[1].plot(x_values, imag(psi), label="Partie imag", color="red",linestyle="--")

#axs[0].plot(x_values, real(psi), label="partie Réelle", color="blue")
#axs[0].plot(x_values, imag(psi), label="partie Imaginaire", color="red",linestyle="--")
#axs[4].plot(x_values, real(psi), label="Partie Imaginaire", color="red", linestyle="--")

# Mise en forme
axs[0].set_title("Somme des parties réelles de 3 ondes planes")
axs[0].set_xlabel("Position x")
axs[0].set_ylabel("Amplitude")
axs[0].legend()
axs[0].grid(True)

axs[1].set_title("Parties réelle et imaginaire d'une onde plane")
axs[1].set_xlabel("Position x")
axs[1].set_ylabel("Amplitude")
axs[1].legend()
axs[1].grid(True)

# Affichage
plt.show()
