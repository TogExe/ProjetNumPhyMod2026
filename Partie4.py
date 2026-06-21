import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# --- Constantes ---
m = 1.0
hbar = 1.0
# --- Variables de mon experience ---

k0 = 2.67
a_param = 3.0  # largeur initiale


V0 = 3.15  # Hauteur de la barrière
longueur_barriere= 3# longueur de la barrière
distance_barriere = 40


# --- Fonctions ---
def thomas_algorithm(a, b, c, d):
    n = len(d)
    cp = c.copy()
    dp = d.copy()
    for i in range(1, n):
        w = a[i - 1] / b[i - 1]
        b[i] = b[i] - w * cp[i - 1]
        dp[i] = dp[i] - w * dp[i - 1]
    x = np.zeros(n, dtype=complex)
    x[-1] = dp[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (dp[i] - cp[i] * x[i + 1]) / b[i]
    return x
def multBpsi(diag, low, up, psi):
    d_vec = np.zeros(len(psi), dtype=complex)
    d_vec[0] = diag[0] * psi[0] + up[0] * psi[1]
    d_vec[-1] = low[-1] * psi[-2] + diag[-1] * psi[-1]

    for i in range(1, len(psi) - 1):
        d_vec[i] = low[i - 1] * psi[i - 1] + diag[i] * psi[i] + up[i] * psi[i + 1]
    return d_vec

def GaussWP(k0, a, ux, t):
    x = ux
    diver = m * a ** 2 + 2j * hbar * t
    amp = (1 / (8 * np.pi ** 3)) ** 0.25 * np.sqrt((4 * np.pi * m * a) / diver)
    partexp = (m / 4) * ((a ** 2 * k0 + 2j * x) ** 2 / diver) - (a ** 2 * k0 ** 2) / 4
    return amp * np.exp(partexp)


# --- Paramètres du maillage ---
# On augmente la taille pour bien voir le mouvement
startx, endx = -60, 150
nx = 3000
x_values = np.linspace(startx, endx, nx)
dx = x_values[1] - x_values[0]

nt = 400
anime_t = 60.0
dt = anime_t / nt
# --- Création du Potentiel V ---

V = np.zeros(nx)
#V0 = 0.0  # Hauteur de la barrière
# On place la barrière par exemple entre x=40 et x=45
for i in range(nx):
    if distance_barriere <= x_values[i] <= distance_barriere+longueur_barriere:
        V[i] = V0

# --- Initialisation (t=0 uniquement) ---
psi_matrice = np.zeros((nt, nx), dtype=complex)
psi_matrice[0, :] = GaussWP(k0, a_param, ux=x_values, t=0.0)
# --- Boucle de Crank-Nicolson ---
r = 1j * hbar * dt / (2 * m * dx ** 2)

# Les diagonales incluent maintenant le potentiel V !
diag_A = 1 +  r + 1j * V * dt / (2 * hbar)
lower_A = -r/2 * np.ones(nx - 1, dtype=complex)
upper_A = -r/2 * np.ones(nx - 1, dtype=complex)

diag_B = 1 -  r - 1j * V * dt / (2 * hbar)
lower_B = r /2* np.ones(nx - 1, dtype=complex)
upper_B = r /2* np.ones(nx - 1, dtype=complex)

psi_actuelle = psi_matrice[0, :]

print("Calcul de l'évolution en cours...")
indice_max_initial = np.argmax(np.abs(psi_matrice[0, :])**2)
position_initiale = x_values[indice_max_initial]
tau_0_num = None  # Variable pour stocker le temps trouvé


distance_b = distance_barriere + longueur_barriere
print("valeurs : V0="+str(V0)+" k="+str(k0)+" énergie e : "+str(hbar**2*k0**2/2*m))
print("distance de la barriere : ", distance_barriere)
print("longueur de la barriere : ", longueur_barriere)
for n in range(1, nt):

    #On résoud A*psy(n+1) = B*psy(n)
    B_psi = multBpsi(diag_B, lower_B, upper_B, psi_actuelle)
    psi_future = thomas_algorithm(lower_A, diag_A.copy(), upper_A, B_psi)
    psi_matrice[n, :] = psi_future
    psi_actuelle = psi_future

    # calcul pour trouver la position
    densite_n = np.abs(psi_actuelle) ** 2
    indice_max_actuel = np.argmax(densite_n)
    position_actuelle = x_values[indice_max_actuel]

    # on regarde la distance parcourue (position_actuelle -position_initiale)
    distance_parcourue = position_actuelle - position_initiale
    if distance_parcourue >= distance_barriere:
        # Le temps correspondant est le numéro de l'itération multiplié par le pas de temps
        tau_0_num = n * dt
        print(f"Distance parcourue : {distance_parcourue:.2f} m")
        print(f"Temps de parcours numérique (τ_0,num) = {tau_0_num:.4f} s")
        distance_barriere = 1000 + distance_barriere
    if distance_parcourue >= distance_b:
        # Le temps correspondant est le numéro de l'itération multiplié par le pas de temps
        tau_0_num = n * dt
        print(f"Distance parcourue : {distance_parcourue:.2f} m")
        print(f"Temps de parcours numérique (τ_0,num) = {tau_0_num:.4f} s")
        distance_b = 1000+ distance_b




print("Calcul terminé ! Lancement de l'animation.")

# --- Animation ---
fig_anim, ax_anim = plt.subplots(figsize=(10, 5))
densite_initiale = np.abs(psi_matrice[0, :]) ** 2

ligne_onde, = ax_anim.plot(x_values, densite_initiale, color="red", label="|ψ(x,t)|²")
# On dessine la barrière pour bien la voir
ax_anim.fill_between(x_values, 0, V/10, color="gray", alpha=0.3, label="Barrière de potentiel V(x)")

ax_anim.set_xlim(-20, 100)
ax_anim.set_ylim(0, np.max(densite_initiale) * 3.2)
ax_anim.set_xlabel("Position x")
ax_anim.set_ylabel("Densité de probabilité")
ax_anim.legend()


def update(frame):
    ligne_onde.set_ydata(np.abs(psi_matrice[frame, :]) ** 2)
    return ligne_onde,


ani = animation.FuncAnimation(fig_anim, update, frames=nt, interval=20, blit=True)
plt.show()
