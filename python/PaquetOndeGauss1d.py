"""
Partie 2 — Paquet d'ondes gaussien
PaquetOndeGauss1d.py

Contenu :
  - Constantes physiques (hbar, masse électron)
  - Fonction GaussWP(k0, a, x, t)  — formule analytique exacte (éq. 6)
  - Visualisation Re[Ψ] et Im[Ψ] à t = 0
  - Discussion de la difficulté numérique et solution (normalisation / fenêtre)
"""

from numpy import pi, exp, sqrt, real, imag, linspace, abs as npabs
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────
# Constantes physiques (SI)
# ─────────────────────────────────────────────

hbar = 1.0545718e-34   # Constante de Planck réduite  (J·s)
m    = 9.10938e-31     # Masse de l'électron           (kg)


# ─────────────────────────────────────────────
# Fonction paquet d'ondes gaussien
# ─────────────────────────────────────────────

def GaussWP(k0, a, x, t):
    """
    Paquet d'ondes gaussien — solution analytique exacte de l'éq. de Schrödinger
    pour une particule libre (équation 6 du projet).

    Paramètres
    ----------
    k0 : vecteur d'onde central (m⁻¹)
    a  : largeur initiale du paquet dans l'espace réel (m)
         À t=0, la largeur spatiale est Δx ≈ a/√2.
    x  : tableau de positions (m)
    t  : instant (s)

    Retour
    ------
    Ψ(x, t) complexe (normalisé)

    Formule
    -------
    Ψ(x,t) = (1/(8π³))^{1/4}
              × sqrt(4π m a / (m a² + 2iℏt))
              × exp[ m/4 × (a²k₀ + 2ix)² / (ma² + 2iℏt)  −  a²k₀²/4 ]
    """
    prefactor  = (1.0 / (8.0 * pi**3))**0.25
    denom      = m * a**2 + 2j * hbar * t
    sqrt_term  = sqrt(4.0 * pi * m * a / denom)
    numerator  = (a**2 * k0 + 2j * x)**2
    exponent   = m / 4.0 * numerator / denom  -  a**2 * k0**2 / 4.0

    return prefactor * sqrt_term * exp(exponent)


# ─────────────────────────────────────────────
# Vérification de la normalisation
# ─────────────────────────────────────────────

def check_normalisation(k0, a, t=0.0, npts=10000):
    """
    Vérifie numériquement ∫|Ψ|² dx ≈ 1
    en intégrant sur une fenêtre large autour du paquet.
    """
    # Fenêtre : ±5 largeurs de paquet autour du centre classique
    x_center = hbar * k0 / m * t              # position classique à t
    sigma_t   = a / sqrt(2) * sqrt(1 + (2*hbar*t / (m*a**2))**2)
    x = linspace(x_center - 10*sigma_t, x_center + 10*sigma_t, npts)
    dx = x[1] - x[0]

    psi = GaussWP(k0, a, x, t)
    norm = (npabs(psi)**2 * dx).sum()
    return norm, x


# ─────────────────────────────────────────────
# Visualisations
# ─────────────────────────────────────────────

def plot_gauss_wp(k0=1e10, a=1e-9, t=0.0, npts=2000):
    """
    Trace Re[Ψ] et Im[Ψ] du paquet d'ondes gaussien en fonction de x.

    Difficulté numérique (question d) :
    ------------------------------------
    Avec les unités SI (hbar ~ 1e-34, m ~ 1e-31), les exposants de l'exponentielle
    deviennent numériquement très grands ou très petits, ce qui provoque des
    dépassements (overflow / underflow) ou une perte de précision.

    Solution / astuce (question e) :
    ---------------------------------
    Deux approches :
      1. Utiliser des unités atomiques (ℏ = m_e = 1), ce qui rend tous les
         nombres proches de l'unité et évite les problèmes numériques.
      2. Choisir des paramètres (k0, a) cohérents avec l'échelle atomique :
         a ~ quelques nm, k0 ~ 1/a, puis centrer la fenêtre d'affichage
         automatiquement sur le paquet.
    Ici on utilise la stratégie 2 : la fenêtre est adaptée à la largeur du paquet.
    """
    norm, x = check_normalisation(k0, a, t, npts)
    psi = GaussWP(k0, a, x, t)

    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)

    ax = axes[0]
    ax.plot(x * 1e9, real(psi), color='royalblue', linewidth=1.2,
            label=r'$\mathrm{Re}[\Psi(x,t)]$')
    ax.plot(x * 1e9, imag(psi), color='tomato', linestyle='--', linewidth=1.2,
            label=r'$\mathrm{Im}[\Psi(x,t)]$')
    ax.set_ylabel(r'$\Psi$ (m$^{-1/2}$)')
    ax.set_title(f'Paquet d\'ondes gaussien  —  $t = {t:.2e}$ s\n'
                 f'$k_0 = {k0:.2e}$ m⁻¹,  $a = {a*1e9:.1f}$ nm  '
                 f'  [norme numérique = {norm:.6f}]')
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(x * 1e9, npabs(psi)**2, color='seagreen', linewidth=1.5,
            label=r'$|\Psi(x,t)|^2$')
    ax.fill_between(x * 1e9, npabs(psi)**2, alpha=0.2, color='seagreen')
    ax.set_xlabel('x (nm)')
    ax.set_ylabel(r'$|\Psi|^2$ (m$^{-1}$)')
    ax.set_title('Densité de probabilité de présence')
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.show()


def plot_gauss_evolution(k0=1e10, a=1e-9, n_times=4, npts=2000):
    """
    Trace l'évolution temporelle de |Ψ(x,t)|² pour plusieurs instants.
    Illustre l'étalement du paquet gaussien au cours du temps.
    """
    # Temps caractéristique : τ = m a² / (2ℏ)
    tau = m * a**2 / (2 * hbar)
    times = linspace(0, 3 * tau, n_times)

    fig, ax = plt.subplots(figsize=(9, 5))
    colors = plt.cm.viridis(linspace(0, 1, n_times))

    for i, t in enumerate(times):
        _, x = check_normalisation(k0, a, t, npts)
        psi = GaussWP(k0, a, x, t)
        label = f't = {t/tau:.1f} τ'
        ax.plot(x * 1e9, npabs(psi)**2, color=colors[i], label=label, linewidth=1.3)

    ax.set_xlabel('x (nm)')
    ax.set_ylabel(r'$|\Psi(x,t)|^2$ (m$^{-1}$)')
    ax.set_title(f'Évolution temporelle du paquet d\'ondes gaussien\n'
                 f'($τ = ma²/2ℏ = {tau:.2e}$ s)')
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()


# ─────────────────────────────────────────────
# Point d'entrée
# ─────────────────────────────────────────────

if __name__ == '__main__':
    # Paramètres physiques typiques (électron, échelle nanométrique)
    k0 = 1e10    # m⁻¹  (~  électron de 3.8 eV)
    a  = 2e-9    # m    (paquet initial de largeur ~2 nm)

    print("=== Vérification normalisation à t=0 ===")
    norm, _ = check_normalisation(k0, a, t=0.0)
    print(f"  ∫|Ψ|² dx = {norm:.8f}  (attendu : 1.0)")

    print("\n=== Visualisation à t = 0 ===")
    plot_gauss_wp(k0=k0, a=a, t=0.0)

    print("\n=== Évolution temporelle ===")
    plot_gauss_evolution(k0=k0, a=a, n_times=5)
