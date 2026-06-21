"""
Partie 3 — Résolution numérique de l'équation de Schrödinger
SchrodingerNumerique.py

Contenu :
  3.1 — Algorithme de dérivation (première et seconde)
  3.2 — Résolution numérique de l'équation de Schrödinger (différences finies)
        Schéma explicite : dérivée temporelle d'ordre 1, dérivée spatiale d'ordre 2
  - Validation : comparaison avec le paquet d'ondes gaussien analytique (partie 2)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Import du paquet gaussien analytique (partie 2)
# On redéfinit ici les constantes pour l'autonomie du fichier
hbar = 1.0545718e-34   # J·s
m    = 9.10938e-31     # kg


# ═══════════════════════════════════════════════════════════════
# 3.1 — ALGORITHME DE DÉRIVATION
# ═══════════════════════════════════════════════════════════════

def derivee_premiere(f, dx):
    """
    Dérivée première par différences finies centrées (ordre 2).

    f'(x_i) ≈ [f(x_{i+1}) - f(x_{i-1})] / (2·dx)

    Aux bords, on utilise les différences avant/arrière (ordre 1).

    Paramètres
    ----------
    f  : tableau numpy 1D (valeurs de la fonction)
    dx : pas spatial (scalaire)

    Retour
    ------
    df : tableau numpy 1D de même taille que f
    """
    n = len(f)
    df = np.zeros(n, dtype=f.dtype)

    # Bord gauche (différence avant)
    df[0] = (f[1] - f[0]) / dx

    # Points intérieurs (différences centrées)
    df[1:-1] = (f[2:] - f[:-2]) / (2 * dx)

    # Bord droit (différence arrière)
    df[-1] = (f[-1] - f[-2]) / dx

    return df


def derivee_seconde(f, dx):
    """
    Dérivée seconde par différences finies centrées (ordre 2).

    f''(x_i) ≈ [f(x_{i+1}) - 2·f(x_i) + f(x_{i-1})] / dx²

    Conditions aux bords : Dirichlet (f = 0 aux extrémités).

    Paramètres
    ----------
    f  : tableau numpy 1D
    dx : pas spatial

    Retour
    ------
    d2f : tableau numpy 1D de même taille
    """
    n = len(f)
    d2f = np.zeros(n, dtype=f.dtype)

    # Points intérieurs
    d2f[1:-1] = (f[2:] - 2*f[1:-1] + f[:-2]) / dx**2

    # Bords : conditions de Dirichlet (Ψ=0 aux bords de la boîte)
    d2f[0]  = 0.0
    d2f[-1] = 0.0

    return d2f


def test_derivees():
    """
    Valide les algorithmes de dérivation sur f(x) = x².
    Dérivée théorique : f'(x) = 2x,  f''(x) = 2.
    """
    x  = np.linspace(0, 2, 500)
    dx = x[1] - x[0]
    f  = x**2

    df_num  = derivee_premiere(f, dx)
    d2f_num = derivee_seconde(f, dx)

    df_th  = 2 * x          # théorique
    d2f_th = np.full_like(x, 2.0)

    # Erreurs relatives (en excluant les bords)
    idx = slice(5, -5)
    err1 = np.max(np.abs(df_num[idx]  - df_th[idx])  / np.abs(df_th[idx]))
    err2 = np.max(np.abs(d2f_num[idx] - d2f_th[idx]) / np.abs(d2f_th[idx]))

    print(f"Erreur relative max — dérivée 1re : {err1:.2e}")
    print(f"Erreur relative max — dérivée 2de : {err2:.2e}")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    axes[0].plot(x, df_th,  label='2x (théorique)',  color='royalblue')
    axes[0].plot(x, df_num, label='f\'(x) numérique', color='tomato', linestyle='--')
    axes[0].set_title("Dérivée première de $x^2$")
    axes[0].legend(); axes[0].grid(alpha=0.3)

    axes[1].plot(x, d2f_th,  label='2 (théorique)',    color='royalblue')
    axes[1].plot(x, d2f_num, label='f\'\'(x) numérique', color='tomato', linestyle='--')
    axes[1].set_title("Dérivée seconde de $x^2$")
    axes[1].legend(); axes[1].grid(alpha=0.3)

    plt.tight_layout()
    plt.show()


# ═══════════════════════════════════════════════════════════════
# 3.2 — RÉSOLUTION NUMÉRIQUE DE L'ÉQUATION DE SCHRÖDINGER
# ═══════════════════════════════════════════════════════════════
#
# Équation de Schrödinger 1D :
#   iℏ ∂Ψ/∂t = -ℏ²/(2m) ∂²Ψ/∂x²  +  V(x)·Ψ
#
# Schéma aux différences finies (Euler explicite en temps) :
#   Ψ(x, t+dt) = Ψ(x, t) + dt/(iℏ) × [-ℏ²/(2m) ∂²Ψ/∂x² + V·Ψ]
#              = Ψ(x, t) - i·dt/ℏ  × [-ℏ²/(2m) ∂²Ψ/∂x² + V·Ψ]
#
# ⚠️  Le schéma d'Euler explicite est instable pour de grands dt.
#     Condition de stabilité (CFL) :  dt < ℏ·dx² / (π·ℏ²/m) = m·dx²/(π·ℏ)
# ═══════════════════════════════════════════════════════════════

def GaussWP_analytique(k0, a, x, t):
    """Paquet d'ondes gaussien analytique (formule 6)."""
    prefactor = (1.0 / (8.0 * np.pi**3))**0.25
    denom     = m * a**2 + 2j * hbar * t
    sqrt_term = np.sqrt(4.0 * np.pi * m * a / denom)
    exponent  = m / 4.0 * (a**2 * k0 + 2j * x)**2 / denom  -  a**2 * k0**2 / 4.0
    return prefactor * sqrt_term * np.exp(exponent)


def resoudre_schrodinger(k0, a, x, t_array, V=None, verbose=True):
    """
    Résout numériquement l'équation de Schrödinger 1D par différences finies.

    Paramètres
    ----------
    k0      : vecteur d'onde central du paquet initial (m⁻¹)
    a       : largeur initiale du paquet (m)
    x       : grille spatiale 1D (numpy array, nx points)
    t_array : grille temporelle 1D (numpy array, nt points)
    V       : potentiel V(x) (numpy array nx, ou None pour V=0)

    Retour
    ------
    psi_2d : tableau 2D de forme (nx, nt) — Ψ(x, t)
    """
    nx = len(x)
    nt = len(t_array)
    dx = x[1] - x[0]
    dt = t_array[1] - t_array[0]

    if V is None:
        V = np.zeros(nx)

    # Critère de stabilité CFL
    dt_max = m * dx**2 / (np.pi * hbar)
    if dt > dt_max and verbose:
        print(f"⚠️  dt={dt:.2e} > dt_max={dt_max:.2e} : risque d'instabilité !")
        print(   "   Réduire dt ou augmenter dx.")

    # Tableau 2D : psi_2d[i, n] = Ψ(x_i, t_n)
    psi_2d = np.zeros((nx, nt), dtype=complex)

    # Condition initiale : paquet gaussien à t=t_array[0]
    psi_2d[:, 0] = GaussWP_analytique(k0, a, x, t_array[0])

    # Boucle temporelle (Euler explicite)
    for n in range(nt - 1):
        psi_n    = psi_2d[:, n]
        d2psi_dx2 = derivee_seconde(psi_n, dx)

        # Hamiltonien appliqué à Ψ :  Ĥ Ψ = -ℏ²/(2m) Ψ'' + V·Ψ
        H_psi = -hbar**2 / (2 * m) * d2psi_dx2 + V * psi_n

        # Évolution temporelle :  Ψ(t+dt) = Ψ(t) - i·dt/ℏ · ĤΨ
        psi_2d[:, n+1] = psi_n - 1j * dt / hbar * H_psi

        # Conditions aux bords : Dirichlet (boîte infinie)
        psi_2d[0,  n+1] = 0.0
        psi_2d[-1, n+1] = 0.0

    return psi_2d


# ─────────────────────────────────────────────
# Validation : comparaison numérique vs analytique
# ─────────────────────────────────────────────

def validation(k0=1e10, a=2e-9, n_temps_check=5):
    """
    Compare la solution numérique au paquet gaussien analytique.
    Affiche |Ψ_num|² et |Ψ_ana|² à plusieurs instants.
    """
    # Paramètres de grille (unités atomiques pratiques : on travaille en nm/fs)
    tau   = m * a**2 / (2 * hbar)               # temps caractéristique
    v_g   = hbar * k0 / m                        # vitesse de groupe

    # Grille spatiale : fenêtre large centrée sur x=0 à t=0
    sigma0 = a / np.sqrt(2)
    L      = max(30 * sigma0, v_g * 3 * tau)
    nx     = 2000
    x      = np.linspace(-L/2, L/2, nx)
    dx     = x[1] - x[0]

    # Grille temporelle
    dt_stable = 0.4 * m * dx**2 / (np.pi * hbar)   # facteur 0.4 pour la stabilité
    t_end     = 2 * tau
    nt        = int(t_end / dt_stable) + 1
    t_array   = np.linspace(0, t_end, nt)

    print(f"Grille : nx={nx}, nt={nt}, dx={dx:.2e} m, dt={dt_stable:.2e} s")
    print(f"Temps caractéristique τ = {tau:.2e} s")
    print("Calcul en cours...")

    psi_2d = resoudre_schrodinger(k0, a, x, t_array, verbose=True)

    print("Calcul terminé. Affichage...")

    # Instants à comparer
    idx_times = np.linspace(0, nt-1, n_temps_check, dtype=int)
    colors = plt.cm.plasma(np.linspace(0, 0.85, n_temps_check))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    for i, n in enumerate(idx_times):
        t_n   = t_array[n]
        label = f't = {t_n/tau:.2f}τ'

        psi_num = psi_2d[:, n]
        psi_ana = GaussWP_analytique(k0, a, x, t_n)

        axes[0].plot(x*1e9, np.abs(psi_num)**2, color=colors[i],
                     label=label, linewidth=1.3)
        axes[1].plot(x*1e9, np.abs(psi_ana)**2, color=colors[i],
                     label=label, linewidth=1.3)

    for ax, titre in zip(axes, ['Numérique $|\\Psi_\\mathrm{num}|^2$',
                                  'Analytique $|\\Psi_\\mathrm{ana}|^2$']):
        ax.set_xlabel('x (nm)')
        ax.set_ylabel('Densité de probabilité (m⁻¹)')
        ax.set_title(titre)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    plt.suptitle('Validation : résolution numérique vs paquet gaussien analytique', y=1.01)
    plt.tight_layout()
    plt.show()

    # Erreur L2 au dernier instant
    psi_fin_num = psi_2d[:, -1]
    psi_fin_ana = GaussWP_analytique(k0, a, x, t_array[-1])
    err_L2 = np.sqrt(np.sum(np.abs(psi_fin_num - psi_fin_ana)**2) * dx)
    print(f"\nErreur L2 à t_final : {err_L2:.4e}  (0 = parfait)")

    return psi_2d, x, t_array


# ─────────────────────────────────────────────
# Point d'entrée
# ─────────────────────────────────────────────

if __name__ == '__main__':
    print("=== 3.1 Test des algorithmes de dérivation ===")
    test_derivees()

    print("\n=== 3.2 Validation numérique vs analytique (V=0) ===")
    print("(peut prendre quelques secondes...)\n")
    psi_2d, x, t_arr = validation(k0=1e10, a=2e-9, n_temps_check=5)
