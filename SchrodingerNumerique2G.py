import numpy as np
import matplotlib.pyplot as plt

m    = 1.0
hbar = 1.0

# 3.1 — ALGORITHME DE DÉRIVATION

# 3.1.1a — Définition de la dérivée

# La dérivée d'une fonction réelle f en un point x est définie par :
#
#   f'(x) = lim_{h→0} [f(x+h) - f(x)] / h
#
# Sur un tableau discret de pas dx, on approche cette limite par :
#   f'(x_i) ≈ [f(x_{i+1}) - f(x_i)] / dx       (différence avant, ordre 1)
#
# Pour une meilleure précision, on utilise les différences centrées (ordre 2) :
#   f'(x_i) ≈ [f(x_{i+1}) - f(x_{i-1})] / (2·dx)

# 3.1.1b — Algorithme sur tableau 1D (dérivée première)

def derivee_premiere(f, dx):
    """
    Calcule la dérivée première d'un tableau 1D par différences finies centrées.

    f'(x_i) ≈ [f(x_{i+1}) - f(x_{i-1})] / (2·dx)

    Aux bords : différences avant/arrière (ordre 1).

    Paramètres
    ----------
    f  : tableau numpy 1D
    dx : pas spatial (scalaire)

    Retour
    ------
    df : tableau numpy 1D de même taille
    """
    n  = len(f)
    df = np.zeros(n, dtype=f.dtype)

    # Bord gauche (différence avant)
    df[0] = (f[1] - f[0]) / dx

    # Points intérieurs (différences centrées)
    df[1:-1] = (f[2:] - f[:-2]) / (2 * dx)

    # Bord droit (différence arrière)
    df[-1] = (f[-1] - f[-2]) / dx

    return df

# 3.1.1c — Fonctions x² et 2x

def carre(x):
    """Retourne x²."""
    return x**2

def derivee_carre(x):
    """Retourne la dérivée analytique de x², soit 2x."""
    return 2 * x

# 3.1.1d — Validation : dérivée numérique de x² vs 2x

def valider_derivee_premiere():
    x  = np.linspace(0, 2, 500)
    dx = x[1] - x[0]
    f  = carre(x)

    df_num = derivee_premiere(f, dx)
    df_th  = derivee_carre(x)

    # Erreur relative (on exclut les bords et x≈0 pour éviter div/0)
    masque    = np.abs(df_th) > 1e-6
    err_rel   = np.abs(df_num[masque] - df_th[masque]) / np.abs(df_th[masque])
    err_max   = np.max(err_rel)

    print("=== 3.1.1d Validation dérivée première ===")
    print(f"  Erreur relative max sur f'(x²) vs 2x : {err_max:.2e}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    axes[0].plot(x, df_th,  label="2x  (analytique)",    color="royalblue")
    axes[0].plot(x, df_num, label="f'(x²) numérique",    color="tomato", linestyle="--")
    axes[0].set_title("Dérivée première de $x^2$")
    axes[0].set_xlabel("x")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    axes[1].plot(x[masque], err_rel, color="seagreen")
    axes[1].set_title("Erreur relative $|f'_{num} - 2x| / |2x|$")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("Erreur relative")
    axes[1].grid(alpha=0.3)

    plt.suptitle("Validation dérivée première (différences finies centrées)")
    plt.tight_layout()
    plt.show()

# 3.1.2 — Dérivée seconde
#
# Définition discrète :
#   f''(x_i) ≈ [f(x_{i+1}) - 2·f(x_i) + f(x_{i-1})] / dx²
#
# C'est l'approximation d'ordre 2 de la dérivée seconde.

def derivee_seconde(f, dx):
    """
    Calcule la dérivée seconde d'un tableau 1D par différences finies centrées.

    f''(x_i) ≈ [f(x_{i+1}) - 2·f(x_i) + f(x_{i-1})] / dx²

    Conditions aux bords : Dirichlet (f=0), cohérent avec une boîte infinie.

    Paramètres
    ----------
    f  : tableau numpy 1D
    dx : pas spatial

    Retour
    ------
    d2f : tableau numpy 1D de même taille
    """
    n   = len(f)
    d2f = np.zeros(n, dtype=f.dtype)

    # Points intérieurs
    d2f[1:-1] = (f[2:] - 2*f[1:-1] + f[:-2]) / dx**2

    # Bords : Dirichlet (Ψ=0 aux parois de la boîte)
    d2f[0]  = 0.0
    d2f[-1] = 0.0

    return d2f


def derivee_seconde_carre(x):
    """Dérivée seconde analytique de x² : constante = 2."""
    return np.full_like(x, 2.0)


def valider_derivee_seconde():
    x  = np.linspace(0, 2, 500)
    dx = x[1] - x[0]
    f  = carre(x)

    d2f_num = derivee_seconde(f, dx)
    d2f_th  = derivee_seconde_carre(x)

    # Erreur (on exclut les bords où la condition Dirichlet est imposée)
    idx    = slice(5, -5)
    err_max = np.max(np.abs(d2f_num[idx] - d2f_th[idx]))

    print("\n=== 3.1.2 Validation dérivée seconde ===")
    print(f"  Erreur absolue max sur f''(x²) vs 2 : {err_max:.2e}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    axes[0].plot(x, d2f_th,  label="2  (analytique)",      color="royalblue")
    axes[0].plot(x, d2f_num, label="f''(x²) numérique",    color="tomato", linestyle="--")
    axes[0].set_title("Dérivée seconde de $x^2$")
    axes[0].set_xlabel("x")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    axes[1].plot(x[idx], np.abs(d2f_num[idx] - d2f_th[idx]), color="seagreen")
    axes[1].set_title("Erreur absolue $|f''_{num} - 2|$")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("Erreur absolue")
    axes[1].grid(alpha=0.3)

    plt.suptitle("Validation dérivée seconde (différences finies centrées)")
    plt.tight_layout()
    plt.show()


# 3.2 — ALGORITHME POUR L'ÉQUATION DE SCHRÖDINGER

# 3.2.1 — Rappel : équation de Schrödinger 1D avec potentiel constant V0
#
# iℏ ∂Ψ/∂t = -ℏ²/(2m) ∂²Ψ/∂x²  +  V0·Ψ
#
# Schéma d'Euler explicite en temps :
#   Ψ(x, t+dt) = Ψ(x, t)  -  i·dt/ℏ · [ -ℏ²/(2m) Ψ''(x,t)  +  V0·Ψ(x,t) ]
#
# Condition de stabilité (CFL) :
#   dt  <  m·dx² / (π·ℏ)


# 3.2.2–3 — Paquet d'ondes gaussien (condition initiale) + grilles

def GaussWP(k0, a, x, t):
    """Paquet d'ondes gaussien analytique (formule 6 du projet)."""
    diver   = m * a**2 + 2j * hbar * t
    amp     = (1 / (8 * np.pi**3))**0.25 * np.sqrt((4 * np.pi * m * a) / diver)
    partexp = (m / 4) * ((a**2 * k0 + 2j * x)**2 / diver) - (a**2 * k0**2) / 4
    return amp * np.exp(partexp)

# 3.2.4 — Algorithme d'évolution (Euler explicite)

def euler_schrodinger(psi0, x, t_array, V0=0.0):
    """
    Résout iℏ ∂Ψ/∂t = [-ℏ²/(2m) ∂²/∂x² + V0] Ψ  par Euler explicite.

    Schéma :
      Ψ^{n+1}_i = Ψ^n_i  -  i·dt/ℏ · [ -ℏ²/(2m) · (Ψ^n_{i+1} - 2Ψ^n_i + Ψ^n_{i-1})/dx²
                                         + V0 · Ψ^n_i ]

    Paramètres
    ----------
    psi0    : condition initiale (tableau complexe de taille nx)
    x       : grille spatiale (nx points)
    t_array : grille temporelle (nt points)
    V0      : hauteur du potentiel constant (défaut 0 = particule libre)

    Retour
    ------
    psi_2d[ix, it] = Ψ(x_ix, t_it)   tableau 2D (nx × nt)
    """
    nx  = len(x)
    nt  = len(t_array)
    dx  = x[1] - x[0]
    dt  = t_array[1] - t_array[0]

    # Vérification stabilité CFL
    dt_max = m * dx**2 / (np.pi * hbar)
    if dt > dt_max:
        print(f"  ⚠️  dt={dt:.2e} > dt_max={dt_max:.2e} : instabilité possible !")

    # Tableau 2D : psi_2d[ix, it]
    psi_2d       = np.zeros((nx, nt), dtype=complex)
    psi_2d[:, 0] = psi0
    psi            = psi0.copy()

    for n in range(nt - 1):
        d2psi    = derivee_seconde(psi, dx)
        H_psi    = -hbar**2 / (2*m) * d2psi + V0 * psi
        psi      = psi - 1j * dt / hbar * H_psi
        psi[0]   = 0.0    # Dirichlet
        psi[-1]  = 0.0
        psi_2d[:, n+1] = psi

    return psi_2d


# 3.2.5 — Confrontation numérique vs analytique (V0=0)

def confrontation(k0=2.0, a=3.0, n_instants=5):
    """
    Compare la solution d'Euler explicite (V0=0) au paquet gaussien analytique.
    Affiche |Ψ_num|² et |Ψ_ana|² à n_instants instants.
    Calcule l'erreur L2 à chaque instant.
    """
    # ── Grille ────────────────────────────────────────────────────────────────
    v_g    = hbar * k0 / m
    tau    = m * a**2 / (2*hbar)       # temps d'étalement caractéristique

    L   = max(40 * a, v_g * 2 * tau * 3)
    nx  = 2000
    x   = np.linspace(-L/2, L/2, nx)
    dx  = x[1] - x[0]

    # Pas de temps stable (facteur 0.35 pour la marge)
    dt    = 0.35 * m * dx**2 / (np.pi * hbar)
    t_end = 2 * tau
    nt    = int(t_end / dt) + 1
    t_arr = np.linspace(0, t_end, nt)

    print(f"\n=== 3.2.5 Confrontation numérique vs analytique ===")
    print(f"  k0={k0}, a={a}, nx={nx}, nt={nt}")
    print(f"  dx={dx:.4f}, dt={dt:.2e}, τ={tau:.4f}")
    print("  Calcul Euler explicite (V0=0)...")

    psi0   = GaussWP(k0, a, x, 0.0)
    psi_2d = euler_schrodinger(psi0, x, t_arr, V0=0.0)

    print("  Calcul terminé.")

    # ── Graphique ─────────────────────────────────────────────────────────────
    idx_list = np.linspace(0, nt-1, n_instants, dtype=int)
    colors   = plt.cm.plasma(np.linspace(0, 0.85, n_instants))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    erreurs_L2 = []

    for i, n in enumerate(idx_list):
        t_n      = t_arr[n]
        psi_num  = psi_2d[:, n]
        psi_ana  = GaussWP(k0, a, x, t_n)
        label    = f"t = {t_n:.2f} = {t_n/tau:.2f}τ"

        axes[0].plot(x, np.abs(psi_num)**2, color=colors[i], label=label, linewidth=1.2)
        axes[1].plot(x, np.abs(psi_ana)**2, color=colors[i], label=label, linewidth=1.2)

        err = np.sqrt(np.sum(np.abs(psi_num - psi_ana)**2) * dx)
        erreurs_L2.append((t_n, err))

    for ax, titre in zip(axes, ["Numérique $|\\Psi_{num}|^2$  (Euler explicite)",
                                  "Analytique $|\\Psi_{ana}|^2$  (formule exacte)"]):
        ax.set_xlabel("x")
        ax.set_ylabel("Densité de probabilité")
        ax.set_title(titre)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    plt.suptitle(f"Confrontation numérique / analytique — V₀=0, k₀={k0}, a={a}")
    plt.tight_layout()
    plt.show()

    # ── Erreur L2 ─────────────────────────────────────────────────────────────
    print("\n  Erreurs L2 à chaque instant affiché :")
    for t_n, err in erreurs_L2:
        print(f"    t = {t_n:.4f} ({t_n/tau:.2f}τ)  →  erreur L2 = {err:.4e}")

    fig2, ax2 = plt.subplots(figsize=(7, 4))
    ts   = [e[0] for e in erreurs_L2]
    errs = [e[1] for e in erreurs_L2]
    ax2.plot(ts, errs, "o-", color="tomato")
    ax2.set_xlabel("t")
    ax2.set_ylabel("Erreur L2")
    ax2.set_title("Erreur L2 entre solution numérique et analytique")
    ax2.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

    return psi_2d, x, t_arr


# POINT D'ENTRÉE

if __name__ == "__main__":

    print("=" * 55)
    print("  PARTIE 3 — Résolution numérique de Schrödinger")
    print("=" * 55)

    # ── 3.1.1 Dérivée première
    valider_derivee_premiere()

    # ── 3.1.2 Dérivée seconde
    valider_derivee_seconde()

    # ── 3.2 Évolution Schrödinger 
    print("\n=== 3.2 Paramètres de la simulation ===")
    k0 = 2.0
    a  = 3.0
    print(f"  k0={k0}, a={a}, m={m}, hbar={hbar}")
    print(f"  V0=0 (particule libre, pour valider l'algo)")

    # Illustration du tableau 2D et des grilles (3.2.2–3)
    tau   = m * a**2 / (2*hbar)
    v_g   = hbar * k0 / m
    L     = 40 * a
    nx    = 500   # réduit pour l'illustration
    nt    = 200
    x_ex  = np.linspace(-L/2, L/2, nx)
    t_ex  = np.linspace(0, 2*tau, nt)
    dx_ex = x_ex[1] - x_ex[0]
    dt_ex = t_ex[1] - t_ex[0]

    print(f"\n=== 3.2.2–3 Tableau 2D et grilles ===")
    print(f"  x : {nx} points sur [{-L/2:.1f}, {L/2:.1f}],  dx={dx_ex:.4f}")
    print(f"  t : {nt} points sur [0, {2*tau:.4f}],        dt={dt_ex:.4f}")
    print(f"  Tableau psi_2d : shape ({nx}, {nt})  = {nx*nt} complexes")

    psi_2d_ex = np.zeros((nx, nt), dtype=complex)
    psi_2d_ex[:, 0] = GaussWP(k0, a, x_ex, 0.0)
    print(f"  psi_2d[:,0] initialisé avec GaussWP(k0={k0}, a={a}, t=0)")
    print(f"  Reste du tableau : zéros (sera rempli par l'algo)")

    # ── 3.2.5 Confrontation (grille fine) ───────────────────────────────────
    psi_2d, x, t_arr = confrontation(k0=k0, a=a, n_instants=5)
