"""
Partie 4 — Effet tunnel : temps de traversée d'une barrière de potentiel
EffetTunnel.py

Structure :
  4.1 Aspects numériques
      a. Propagation du paquet d'ondes avec barrière rectangulaire
      b. Temps τ0_num  (particule libre, V=0)  pour parcourir une distance a
      c. Temps τt_num  (particule tunnelisant, V=V0)
      d. Influence de la largeur a de la barrière
      e. Influence de la hauteur V0

  4.2 Comparaison au cas classique
      - E > V0  : particule passe (vitesse réduite)
      - 0<E<V0  : particule bloquée classiquement, mais tunnel QM

  4.3 Aspects analytiques
      a. États stationnaires + coefficients de transmission/réflexion
      b. Vitesse de phase et de groupe du paquet gaussien
      c. Temps τ0_th  (libre, gaussien)
      d. Influence de a sur τ0_th
      e. Paquet transmis
      f. Temps τt_th  (Büttiker–Landauer / temps de groupe)
      g. Influence de a sur τt_th
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle

# ─── Constantes physiques (SI) ───────────────────────────────────────────────
hbar = 1.0545718e-34   # J·s
m    = 9.10938e-31     # kg (électron)
eV   = 1.60218e-19     # J

# ══════════════════════════════════════════════════════════════════════════════
# UTILITAIRES COMMUNS
# ══════════════════════════════════════════════════════════════════════════════

def GaussWP(k0, a, x, t):
    """Paquet d'ondes gaussien analytique (formule 6)."""
    prefactor = (1.0 / (8.0 * np.pi**3))**0.25
    denom     = m * a**2 + 2j * hbar * t
    sqrt_term = np.sqrt(4.0 * np.pi * m * a / denom)
    exponent  = m / 4.0 * (a**2 * k0 + 2j * x)**2 / denom - a**2 * k0**2 / 4.0
    return prefactor * sqrt_term * np.exp(exponent)


def derivee_seconde(f, dx):
    """Dérivée seconde par différences finies centrées (ordre 2, Dirichlet aux bords)."""
    d2f = np.zeros_like(f)
    d2f[1:-1] = (f[2:] - 2*f[1:-1] + f[:-2]) / dx**2
    return d2f


def barriere_potentiel(x, x1, x2, V0):
    """Potentiel rectangulaire : V0 pour x1 ≤ x ≤ x2, 0 ailleurs."""
    V = np.zeros_like(x, dtype=float)
    V[(x >= x1) & (x <= x2)] = V0
    return V


def resoudre_schrodinger(psi0, x, t_array, V):
    """
    Résolution numérique de l'équation de Schrödinger 1D.
    Schéma d'Euler explicite (différences finies).

    iℏ ∂Ψ/∂t = [-ℏ²/(2m) ∂²/∂x² + V(x)] Ψ
    → Ψ(t+dt) = Ψ(t) - i·dt/ℏ · ĤΨ(t)

    Retourne psi_2d[ix, it] = Ψ(x_ix, t_it).
    Stocke uniquement si store_all=True pour économiser la RAM.
    """
    nx = len(x)
    nt = len(t_array)
    dx = x[1] - x[0]
    dt = t_array[1] - t_array[0]

    psi_2d = np.zeros((nx, nt), dtype=complex)
    psi_2d[:, 0] = psi0

    psi = psi0.copy()
    for n in range(nt - 1):
        H_psi = -hbar**2 / (2*m) * derivee_seconde(psi, dx) + V * psi
        psi   = psi - 1j * dt / hbar * H_psi
        psi[0] = psi[-1] = 0.0          # Dirichlet
        psi_2d[:, n+1] = psi

    return psi_2d


def grille(k0, a, a_bar, V0, facteur_L=40, nx=3000, n_periodes=800):
    """
    Construit une grille spatiale et temporelle adaptée au problème.

    Le paquet part à gauche de la barrière, traverse, et on attend
    qu'il soit complètement passé.

    Retourne x, t_array, x1_bar, x2_bar, V, dt
    """
    sigma0  = a / np.sqrt(2)          # largeur initiale
    v_g     = hbar * k0 / m           # vitesse de groupe
    tau_tau = m * a**2 / (2*hbar)     # temps caractéristique d'étalement

    # Boîte : le paquet démarre à -L/3, la barrière est centrée en 0
    L    = max(facteur_L * sigma0, 6 * a_bar + 2 * v_g * 3*tau_tau)
    x    = np.linspace(-L/2, L/2, nx)
    dx   = x[1] - x[0]

    # Position initiale du paquet : à gauche de la barrière
    x_start = -L/4

    # Barrière centrée en 0
    x1_bar = -a_bar / 2
    x2_bar =  a_bar / 2

    # Potentiel
    V = barriere_potentiel(x, x1_bar, x2_bar, V0)

    # Pas de temps stable
    dt = 0.3 * m * dx**2 / (np.pi * hbar)

    # Durée : temps pour traverser L/2 à v_g, avec marge
    t_end   = (L/2 + abs(x_start)) / v_g * 3
    nt      = max(int(t_end / dt), n_periodes)
    t_array = np.linspace(0, t_end, nt)

    return x, t_array, x_start, x1_bar, x2_bar, V, dt


# ══════════════════════════════════════════════════════════════════════════════
# 4.1a  PROPAGATION AVEC BARRIÈRE
# ══════════════════════════════════════════════════════════════════════════════

def simulation_barriere(k0, a, a_bar, V0, nx=2000, n_snapshots=6, titre=None):
    """
    Simule la propagation du paquet d'ondes gaussien rencontrant
    une barrière rectangulaire de hauteur V0 et de largeur a_bar.

    Affiche n_snapshots instantanés de |Ψ(x,t)|².
    """
    x, t_array, x_start, x1_bar, x2_bar, V, dt = grille(
        k0, a, a_bar, V0, nx=nx)

    E_moy = (hbar * k0)**2 / (2*m)   # énergie cinétique moyenne

    print(f"  Énergie moyenne  E = {E_moy/eV:.3f} eV")
    print(f"  Hauteur barrière V0= {V0/eV:.3f} eV")
    print(f"  Ratio E/V0       = {E_moy/V0:.3f}")
    print(f"  Grille : nx={nx}, nt={len(t_array)}, dx={x[1]-x[0]:.2e} m")

    # Condition initiale centrée en x_start
    psi0 = GaussWP(k0, a, x - x_start, 0.0)

    psi_2d = resoudre_schrodinger(psi0, x, t_array, V)

    # ── Graphique ────────────────────────────────────────────────────────────
    idx_snaps = np.linspace(0, len(t_array)-1, n_snapshots, dtype=int)
    colors    = plt.cm.viridis(np.linspace(0, 0.9, n_snapshots))

    fig, ax = plt.subplots(figsize=(11, 5))
    for i, n in enumerate(idx_snaps):
        t_n = t_array[n]
        ax.plot(x*1e9, np.abs(psi_2d[:, n])**2 / np.max(np.abs(psi0)**2),
                color=colors[i], label=f't = {t_n:.2e} s', linewidth=1.1)

    # Barrière
    V_norm = V / V0 * 0.5 if V0 > 0 else V
    ax.fill_between(x*1e9, V_norm, alpha=0.25, color='red', label=f'Barrière V₀={V0/eV:.2f} eV')

    ax.set_xlabel('x (nm)')
    ax.set_ylabel(r'$|\Psi|^2$ (normalisé)')
    ax.set_title(titre or f'Effet tunnel — $k_0$={k0:.2e} m⁻¹, $a_{{bar}}$={a_bar*1e9:.1f} nm, $V_0$={V0/eV:.2f} eV')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

    return psi_2d, x, t_array, x_start, x1_bar, x2_bar


# ══════════════════════════════════════════════════════════════════════════════
# 4.1b–c  MESURE DES TEMPS τ0_num ET τt_num
# ══════════════════════════════════════════════════════════════════════════════

def position_moyenne(psi_col, x):
    """Calcule ⟨x⟩ = ∫ x |Ψ|² dx à un instant donné."""
    dx  = x[1] - x[0]
    rho = np.abs(psi_col)**2
    return np.sum(x * rho) * dx


def mesurer_temps(k0, a, a_bar, V0, nx=2000):
    """
    Mesure numérique de :
      τ0_num : temps que met le centre du paquet LIBRE pour parcourir a_bar
      τt_num : temps que met le centre du paquet TUNNEL pour traverser la barrière

    Méthode : on suit ⟨x(t)⟩ et on repère les instants où le centre
    du paquet entre (x = x1_bar) et sort (x = x2_bar) de la barrière.
    """
    results = {}

    for label, V0_val in [('libre', 0.0), ('tunnel', V0)]:
        x, t_array, x_start, x1_bar, x2_bar, V, _ = grille(
            k0, a, a_bar, V0_val, nx=nx)
        psi0   = GaussWP(k0, a, x - x_start, 0.0)
        psi_2d = resoudre_schrodinger(psi0, x, t_array, V)

        # Suivre ⟨x(t)⟩
        x_moy = np.array([position_moyenne(psi_2d[:, n], x)
                          for n in range(len(t_array))])

        # Instant d'entrée : ⟨x⟩ = x1_bar
        # Instant de sortie : ⟨x⟩ = x2_bar
        idx_in  = np.argmin(np.abs(x_moy - x1_bar))
        idx_out = np.argmin(np.abs(x_moy - x2_bar))

        tau = t_array[idx_out] - t_array[idx_in]
        results[label] = {'tau': tau, 'x_moy': x_moy, 't': t_array}
        print(f"  τ_{label[:3]}_num = {tau:.4e} s")

    return results


# ══════════════════════════════════════════════════════════════════════════════
# 4.1d–e  INFLUENCE DE a ET V0
# ══════════════════════════════════════════════════════════════════════════════

def influence_largeur(k0, a, V0, a_bar_values, nx=1500):
    """
    Trace τ0_num et τt_num en fonction de la largeur a_bar de la barrière.
    """
    tau0_list = []
    taut_list = []

    for a_bar in a_bar_values:
        print(f"  a_bar = {a_bar*1e9:.1f} nm ...")
        res = mesurer_temps(k0, a, a_bar, V0, nx=nx)
        tau0_list.append(res['libre']['tau'])
        taut_list.append(res['tunnel']['tau'])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(a_bar_values*1e9, np.array(tau0_list)*1e15,
            'o-', color='royalblue', label=r'$\tau_{0,\mathrm{num}}$ (libre)')
    ax.plot(a_bar_values*1e9, np.array(taut_list)*1e15,
            's--', color='tomato',     label=r'$\tau_{t,\mathrm{num}}$ (tunnel)')
    ax.set_xlabel('Largeur de la barrière $a$ (nm)')
    ax.set_ylabel('Temps (fs)')
    ax.set_title(f'Influence de la largeur — $V_0$={V0/eV:.2f} eV')
    ax.legend(); ax.grid(alpha=0.3)
    plt.tight_layout(); plt.show()

    return np.array(tau0_list), np.array(taut_list)


def influence_V0(k0, a, a_bar, V0_values, nx=1500):
    """
    Trace τt_num en fonction de la hauteur V0 de la barrière.
    """
    taut_list = []
    E_moy = (hbar * k0)**2 / (2*m)

    for V0 in V0_values:
        if V0 >= E_moy:   # régime tunnel pur
            print(f"  V0 = {V0/eV:.2f} eV (> E={E_moy/eV:.2f} eV, tunnel) ...")
        else:
            print(f"  V0 = {V0/eV:.2f} eV (< E, passage classique) ...")
        res = mesurer_temps(k0, a, a_bar, V0, nx=nx)
        taut_list.append(res['tunnel']['tau'])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.axvline(E_moy/eV, color='gray', linestyle=':', label=f'E = {E_moy/eV:.2f} eV')
    ax.plot(np.array(V0_values)/eV, np.array(taut_list)*1e15,
            's-', color='seagreen', label=r'$\tau_{t,\mathrm{num}}$')
    ax.set_xlabel('Hauteur de la barrière $V_0$ (eV)')
    ax.set_ylabel('Temps (fs)')
    ax.set_title(f'Influence de V₀ — $a_{{bar}}$={a_bar*1e9:.1f} nm')
    ax.legend(); ax.grid(alpha=0.3)
    plt.tight_layout(); plt.show()

    return np.array(taut_list)


# ══════════════════════════════════════════════════════════════════════════════
# 4.2  COMPARAISON AU CAS CLASSIQUE
# ══════════════════════════════════════════════════════════════════════════════

def comparaison_classique(k0, a_bar, V0_values):
    """
    Compare τ classique et τ quantique (numérique) pour différents V0.

    Cas classique :
      - E > V0 : v_in = ℏk0/m,  v_bar = ℏ√(k0²-2mV0/ℏ²)/m
                 τ_cl = a_bar / v_bar
      - E < V0 : la particule est réfléchie → τ_cl = +∞
    """
    E_moy = (hbar * k0)**2 / (2*m)
    v_g   = hbar * k0 / m

    print("\n=== Comparaison classique / quantique ===")
    print(f"{'V0 (eV)':>10} {'E/V0':>8} {'τ_cl (fs)':>12} {'τ0_num (fs)':>13} {'τt_num (fs)':>13}")
    print("-"*60)

    for V0 in V0_values:
        if V0 < E_moy:
            k_bar  = np.sqrt(k0**2 - 2*m*V0/hbar**2)
            v_bar  = hbar * k_bar / m
            tau_cl = a_bar / v_bar * 1e15
        else:
            tau_cl = float('inf')   # bloquée classiquement

        tau0 = a_bar / v_g * 1e15   # libre

        print(f"{V0/eV:>10.3f} {E_moy/V0:>8.3f} {tau_cl:>12.3f} {tau0:>13.3f}  (numérique non calculé ici)")

    print("\nNote : τt_num > τ0_num dans le régime tunnel — le paquet est ralenti.")
    print("       τt_num peut être < τ_cl (E>V0) car le paquet gaussien est étalé en k.")


# ══════════════════════════════════════════════════════════════════════════════
# 4.3a  ÉTATS STATIONNAIRES — COEFFICIENTS DE TRANSMISSION ET RÉFLEXION
# ══════════════════════════════════════════════════════════════════════════════

def coeff_transmission(k0_val, a_bar, V0):
    """
    Calcule analytiquement le coefficient de transmission T et de réflexion R
    pour une barrière rectangulaire de largeur a_bar et de hauteur V0,
    pour une onde plane de vecteur d'onde k0_val.

    Régions :
      I  (x < 0)       : Ψ = e^{ikx} + r·e^{-ikx},    k  = sqrt(2mE)/ℏ
      II (0 < x < a_bar): Ψ = A·e^{iqx} + B·e^{-iqx}, q  = sqrt(2m(E-V0))/ℏ
                           (imaginaire pur si E < V0 → atténuation exponentielle)
      III(x > a_bar)   : Ψ = t·e^{ikx}

    Conditions aux limites en x=0 et x=a_bar → système linéaire → t, r.

    T = |t|²,  R = |r|²,  T + R = 1.
    """
    E  = (hbar * k0_val)**2 / (2*m)
    k  = k0_val

    if E > V0:
        q = np.sqrt(2*m*(E - V0)) / hbar   # réel
    else:
        q = 1j * np.sqrt(2*m*(V0 - E)) / hbar   # imaginaire → kappa

    # Matrice de transfert (formule standard)
    M11 = np.cos(q * a_bar) - 1j/2 * (k/q + q/k) * np.sin(q * a_bar)
    M21 = -1j/2 * (q/k - k/q) * np.sin(q * a_bar)

    t = 1.0 / M11 * np.exp(-1j * k * a_bar)
    r = M21 / M11

    T = np.abs(t)**2
    R = np.abs(r)**2
    return T, R, t, r


def plot_transmission(k0, a, a_bar_vals, V0):
    """
    Trace T et R en fonction de la largeur a_bar pour une énergie fixée.
    """
    T_list, R_list = [], []
    for a_b in a_bar_vals:
        T, R, _, _ = coeff_transmission(k0, a_b, V0)
        T_list.append(T)
        R_list.append(R)

    E_moy = (hbar*k0)**2/(2*m)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(a_bar_vals*1e9, T_list, 'o-', color='seagreen', label='T (transmission)')
    ax.plot(a_bar_vals*1e9, R_list, 's--', color='tomato',  label='R (réflexion)')
    ax.axhline(1, color='gray', linestyle=':', linewidth=0.8)
    ax.set_xlabel('Largeur de la barrière $a$ (nm)')
    ax.set_ylabel('Coefficient')
    ax.set_title(f'T et R analytiques — $V_0$={V0/eV:.2f} eV, $E$={E_moy/eV:.2f} eV')
    ax.legend(); ax.grid(alpha=0.3)
    ax.set_ylim(-0.05, 1.1)
    plt.tight_layout(); plt.show()


# ══════════════════════════════════════════════════════════════════════════════
# 4.3b  VITESSE DE PHASE ET DE GROUPE DU PAQUET GAUSSIEN
# ══════════════════════════════════════════════════════════════════════════════

def vitesses(k0):
    """
    Relation de dispersion particule libre : ω(k) = ℏk²/(2m)
    Vitesse de phase   : v_φ = ω/k   = ℏk/(2m)
    Vitesse de groupe  : v_g = dω/dk = ℏk/m
    """
    omega0 = hbar * k0**2 / (2*m)
    v_phi  = hbar * k0 / (2*m)
    v_g    = hbar * k0 / m
    E      = hbar * omega0

    print(f"\n=== Vitesses (k0 = {k0:.3e} m⁻¹) ===")
    print(f"  ω0        = {omega0:.3e} rad/s")
    print(f"  E         = {E/eV:.4f} eV")
    print(f"  v_phase   = ℏk₀/(2m) = {v_phi:.3e} m/s")
    print(f"  v_groupe  = ℏk₀/m    = {v_g:.3e} m/s  (= v_classe)")
    print(f"  v_g / v_φ = {v_g/v_phi:.1f}  (toujours 2 pour particule libre)")
    return v_phi, v_g


# ══════════════════════════════════════════════════════════════════════════════
# 4.3c–d  TEMPS τ0_th (PARTICULE LIBRE, GAUSSIEN)
# ══════════════════════════════════════════════════════════════════════════════

def tau0_theorique(k0, a, a_bar):
    """
    Temps théorique pour que le centre du paquet gaussien libre
    parcoure la distance a_bar.

    Le centre du paquet se déplace à la vitesse de groupe :
        ⟨x(t)⟩ = ⟨x(0)⟩ + v_g · t

    Donc : τ0_th = a_bar / v_g = a_bar · m / (ℏ k0)
    """
    v_g   = hbar * k0 / m
    tau0  = a_bar / v_g
    print(f"\n  τ0_th = a_bar / v_g = {a_bar*1e9:.2f} nm / {v_g:.3e} m/s = {tau0:.4e} s")
    return tau0


def influence_a_tau0(k0, a, a_bar_values):
    """
    τ0_th = a_bar · m / (ℏ k0)  est linéaire en a_bar (trivial, mais illustratif).
    """
    tau0_vals = np.array([tau0_theorique(k0, a, ab) for ab in a_bar_values])

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(a_bar_values*1e9, tau0_vals*1e15, 'o-', color='royalblue')
    ax.set_xlabel('$a$ (nm)')
    ax.set_ylabel('$\\tau_{0,th}$ (fs)')
    ax.set_title('Temps de traversée libre $\\tau_{0,th}$ vs largeur $a$\n(linéaire : $\\tau_0 = ma/(\\hbar k_0)$)')
    ax.grid(alpha=0.3)
    plt.tight_layout(); plt.show()

    return tau0_vals


# ══════════════════════════════════════════════════════════════════════════════
# 4.3e  PAQUET TRANSMIS
# ══════════════════════════════════════════════════════════════════════════════

def paquet_transmis_expression(k0, a, a_bar, V0):
    """
    Le paquet transmis s'exprime comme la transformée de Fourier inverse
    de g(k) · t(k), où t(k) est l'amplitude de transmission de l'onde plane k.

    Ψ_trans(x, t) = 1/√(2π) ∫ g(k) · t(k, a_bar, V0) · e^{i(kx - ω(k)t)} dk

    Pour un paquet gaussien étroit (Δk ≪ k0), on peut approximer :
        t(k) ≈ t(k0)  (approximation de phase stationnaire / paquet étroit)

    Alors :
        Ψ_trans(x, t) ≈ t(k0) · Ψ_libre(x, t)

    Le paquet transmis est donc le paquet gaussien libre multiplié par t(k0),
    décalé en phase et en amplitude.
    """
    T, R, t_amp, r_amp = coeff_transmission(k0, a_bar, V0)
    E = (hbar*k0)**2 / (2*m)

    print(f"\n=== Paquet transmis (approximation paquet étroit) ===")
    print(f"  t(k0)   = {t_amp:.4f}  (amplitude complexe de transmission)")
    print(f"  T = |t|² = {T:.4f}")
    print(f"  R = |r|² = {R:.4f}")
    print(f"  T + R     = {T+R:.6f}  (doit être = 1)")
    print(f"\n  Ψ_trans(x,t) ≈ t(k0) · Ψ_gaussien_libre(x - a_bar, t)")
    print(f"  Le paquet transmis est gaussien, de même largeur, mais d'amplitude réduite de √T.")

    return T, R, t_amp


# ══════════════════════════════════════════════════════════════════════════════
# 4.3f–g  TEMPS τt_th (BÜTTIKER–LANDAUER) ET INFLUENCE DE a
# ══════════════════════════════════════════════════════════════════════════════

def tau_tunnel_theorique(k0, a_bar, V0):
    """
    Temps de traversée tunnel par la méthode du temps de groupe (Wigner).

    Le temps de groupe (phase stationnaire) est défini par :
        τ_t = dφ_t / dω  = (ℏ / v_g) · dφ_t/dk

    où φ_t = arg[t(k)] est la phase de l'amplitude de transmission.

    On calcule dφ_t/dk numériquement par différences finies.

    Référence : Büttiker & Landauer (1982), Hauge & Støvneng (1989).
    """
    dk = k0 * 1e-5   # pas infinitésimal pour la dérivée numérique

    _, _, t_m, _ = coeff_transmission(k0 - dk, a_bar, V0)
    _, _, t_p, _ = coeff_transmission(k0 + dk, a_bar, V0)

    phi_m = np.angle(t_m)
    phi_p = np.angle(t_p)

    dphi_dk = (phi_p - phi_m) / (2 * dk)

    # τ_t = ℏ · dφ/dE = ℏ · dφ/dk · dk/dE = (m/ℏk0) · dφ/dk
    tau_t = (m / (hbar * k0)) * dphi_dk

    E_moy = (hbar*k0)**2 / (2*m)
    T, _, _, _ = coeff_transmission(k0, a_bar, V0)

    print(f"  τt_th (Wigner) = {tau_t:.4e} s")
    print(f"  T = {T:.4f},  E = {E_moy/eV:.3f} eV,  V0 = {V0/eV:.3f} eV")

    return tau_t


def influence_a_tunnel(k0, a_bar_values, V0):
    """
    Trace τ0_th et τt_th en fonction de a_bar.
    Met en évidence le 'Hartman effect' : τt devient indépendant de a_bar
    pour de grandes barrières en régime tunnel pur.
    """
    v_g = hbar * k0 / m
    tau0_vals = a_bar_values / v_g
    taut_vals = np.array([tau_tunnel_theorique(k0, ab, V0) for ab in a_bar_values])

    E_moy = (hbar*k0)**2/(2*m)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(a_bar_values*1e9, tau0_vals*1e15,
            'o-', color='royalblue', label=r'$\tau_{0,th} = a/v_g$ (libre)')
    ax.plot(a_bar_values*1e9, taut_vals*1e15,
            's--', color='tomato', label=r'$\tau_{t,th}$ (Wigner/tunnel)')

    ax.set_xlabel('Largeur de la barrière $a$ (nm)')
    ax.set_ylabel('Temps (fs)')
    ax.set_title(f'Temps théoriques vs largeur\n'
                 f'$V_0$={V0/eV:.2f} eV, $E$={E_moy/eV:.2f} eV\n'
                 f'Effet Hartman : τ_t sature quand a → ∞ en régime tunnel')
    ax.legend(); ax.grid(alpha=0.3)
    plt.tight_layout(); plt.show()

    return tau0_vals, taut_vals


# ══════════════════════════════════════════════════════════════════════════════
# BILAN FINAL — TABLEAU COMPARATIF
# ══════════════════════════════════════════════════════════════════════════════

def bilan(k0, a, a_bar, V0):
    """
    Affiche un tableau comparatif : classique / numérique / analytique.
    """
    E_moy = (hbar*k0)**2 / (2*m)
    v_g   = hbar * k0 / m

    tau0_th = a_bar / v_g
    tau_t_th = tau_tunnel_theorique(k0, a_bar, V0)
    T, R, t_amp, _ = coeff_transmission(k0, a_bar, V0)

    if E_moy > V0:
        k_bar  = np.sqrt(2*m*(E_moy - V0)) / hbar
        v_bar  = hbar * k_bar / m
        tau_cl = a_bar / v_bar
    else:
        tau_cl = float('inf')

    print("\n" + "═"*65)
    print("  BILAN COMPARATIF")
    print("═"*65)
    print(f"  Paramètres : k0={k0:.2e} m⁻¹, a_bar={a_bar*1e9:.1f} nm, V0={V0/eV:.2f} eV")
    print(f"  Énergie moyenne E = {E_moy/eV:.3f} eV  ({'E > V0' if E_moy>V0 else 'E < V0 → tunnel'})")
    print(f"  Coefficient de transmission T = {T:.4f}")
    print("-"*65)
    print(f"  τ0_th  (libre, analytique)  = {tau0_th:.4e} s  = {tau0_th*1e15:.2f} fs")
    print(f"  τt_th  (tunnel, Wigner)     = {tau_t_th:.4e} s = {tau_t_th*1e15:.2f} fs")
    print(f"  τ_cl   (classique)          = {'∞ (bloqué)' if np.isinf(tau_cl) else f'{tau_cl:.4e} s = {tau_cl*1e15:.2f} fs'}")
    print("═"*65)


# ══════════════════════════════════════════════════════════════════════════════
# POINT D'ENTRÉE
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':

    # ── Paramètres de référence ─────────────────────────────────────────────
    k0    = 1.0e10          # m⁻¹  (~électron de 3.8 eV)
    a     = 2.0e-9          # m    (largeur initiale du paquet)
    a_bar = 3.0e-9          # m    (largeur de la barrière)
    V0    = 2.0 * eV        # J    (hauteur barrière, > E pour effet tunnel)
    E_moy = (hbar*k0)**2 / (2*m)
    print(f"Énergie cinétique moyenne : E = {E_moy/eV:.3f} eV")
    print(f"Hauteur barrière          : V0= {V0/eV:.3f} eV")
    print(f"Régime                    : {'tunnel (E < V0)' if E_moy < V0 else 'passage classique (E > V0)'}\n")

    # ── 4.1a : Propagation avec barrière ────────────────────────────────────
    print("=== 4.1a Simulation propagation avec barrière ===")
    simulation_barriere(k0, a, a_bar, V0, nx=1500, n_snapshots=6)

    # ── 4.1b–c : Mesure des temps numériques ────────────────────────────────
    print("\n=== 4.1b–c Mesure τ0_num et τt_num ===")
    res = mesurer_temps(k0, a, a_bar, V0, nx=1500)

    # ── 4.1d : Influence de la largeur ──────────────────────────────────────
    print("\n=== 4.1d Influence de la largeur a_bar ===")
    a_bar_vals = np.linspace(1e-9, 6e-9, 6)
    influence_largeur(k0, a, V0, a_bar_vals, nx=1200)

    # ── 4.1e : Influence de V0 ──────────────────────────────────────────────
    print("\n=== 4.1e Influence de V0 ===")
    V0_vals = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 4.0]) * eV
    influence_V0(k0, a, a_bar, V0_vals, nx=1200)

    # ── 4.2 : Comparaison classique ─────────────────────────────────────────
    print("\n=== 4.2 Comparaison au cas classique ===")
    comparaison_classique(k0, a_bar, V0_vals)

    # ── 4.3a : Coefficients de transmission ─────────────────────────────────
    print("\n=== 4.3a Coefficients T et R ===")
    a_bar_vals2 = np.linspace(0.5e-9, 8e-9, 50)
    plot_transmission(k0, a, a_bar_vals2, V0)

    # ── 4.3b : Vitesses ─────────────────────────────────────────────────────
    print("\n=== 4.3b Vitesses ===")
    vitesses(k0)

    # ── 4.3c–d : τ0_th ──────────────────────────────────────────────────────
    print("\n=== 4.3c–d τ0_th et influence de a ===")
    tau0_theorique(k0, a, a_bar)
    influence_a_tau0(k0, a, a_bar_vals)

    # ── 4.3e : Paquet transmis ──────────────────────────────────────────────
    print("\n=== 4.3e Paquet transmis ===")
    paquet_transmis_expression(k0, a, a_bar, V0)

    # ── 4.3f–g : τt_th et influence de a ────────────────────────────────────
    print("\n=== 4.3f–g τt_th (Wigner) et influence de a ===")
    tau_tunnel_theorique(k0, a_bar, V0)
    influence_a_tunnel(k0, a_bar_vals2[:20], V0)

    # ── Bilan final ──────────────────────────────────────────────────────────
    bilan(k0, a, a_bar, V0)
