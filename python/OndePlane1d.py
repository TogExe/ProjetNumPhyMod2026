"""
Partie 1 — Ondes planes 1D
OndePlane1d.py

Contenu :
  - Fonction PlaneWave(amp, k, omega, x, t)
  - Visualisation partie réelle et imaginaire d'une onde plane
  - Superposition de 3 ondes planes + enveloppe (section 1.2)
"""

from numpy import pi, exp, sqrt, real, imag, zeros, linspace, cos
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────
# 1.1 — Fonction onde plane
# ─────────────────────────────────────────────

def PlaneWave(amp, k, omega, x, t):
    """
    Retourne une onde plane à 1 dimension.

    Paramètres
    ----------
    amp   : amplitude (en m^{-1/2})
    k     : vecteur d'onde (en m^{-1})
    omega : pulsation (en rad/s)
    x     : position ou tableau de positions (en m)
    t     : instant (en s)

    Retour
    ------
    Ψ(x, t) = amp * exp(i*(k*x - omega*t))   (complexe)
    """
    return amp * exp(1j * (k * x - omega * t))


# ─────────────────────────────────────────────
# 1.2 — Superposition de 3 ondes planes
# ─────────────────────────────────────────────

def SumThreeWaves(amp, k0, dk, x, t=0):
    """
    Somme de trois ondes planes à t=0 (par défaut) :
      - Ψ1 : amplitude amp,   vecteur d'onde k0
      - Ψ2 : amplitude amp/2, vecteur d'onde k0 - dk/2
      - Ψ3 : amplitude amp/2, vecteur d'onde k0 + dk/2

    Résultat analytique :
      Ψ = amp * [1 + cos(dk*x/2)] * exp(i*k0*x)

    Retour
    ------
    psi_sum : onde résultante (complexe)
    envelope_pos : enveloppe positive  +amp*(1+cos(dk*x/2))
    envelope_neg : enveloppe négative  -amp*(1+cos(dk*x/2))
    """
    psi1 = PlaneWave(amp,     k0,        0, x, t)
    psi2 = PlaneWave(amp / 2, k0 - dk/2, 0, x, t)
    psi3 = PlaneWave(amp / 2, k0 + dk/2, 0, x, t)
    psi_sum = psi1 + psi2 + psi3

    envelope_pos =  amp * (1 + cos(dk * x / 2))
    envelope_neg = -amp * (1 + cos(dk * x / 2))

    return psi_sum, envelope_pos, envelope_neg


# ─────────────────────────────────────────────
# Visualisations
# ─────────────────────────────────────────────

def plot_plane_wave(amp=1.0, k=2*pi, omega=1.0, t=0.0, npts=500):
    """Trace Re[Ψ] et Im[Ψ] d'une onde plane en fonction de x."""
    x = linspace(-2, 2, npts)           # intervalle en mètres (arbitraire)
    psi = PlaneWave(amp, k, omega, x, t)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x, real(psi), label=r'$\mathrm{Re}[\Psi]$', color='royalblue')
    ax.plot(x, imag(psi), label=r'$\mathrm{Im}[\Psi]$', color='tomato', linestyle='--')
    ax.set_xlabel('x (m)')
    ax.set_ylabel(r'$\Psi(x, t)$')
    ax.set_title(f'Onde plane à 1D  —  $t = {t}$ s,  $k = {k:.2f}$ m⁻¹')
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_superposition(amp=1.0, k0=20*pi, dk=4*pi, npts=1000):
    """
    Trace les parties réelles des 3 ondes individuelles,
    de leur somme, et l'enveloppe — sur [-π/Δk, π/Δk].
    """
    x_max = pi / dk
    x = linspace(-x_max, x_max, npts)

    psi1 = PlaneWave(amp,     k0,        0, x, 0)
    psi2 = PlaneWave(amp / 2, k0 - dk/2, 0, x, 0)
    psi3 = PlaneWave(amp / 2, k0 + dk/2, 0, x, 0)
    psi_sum, env_pos, env_neg = SumThreeWaves(amp, k0, dk, x)

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    # --- Graphique du haut : 3 ondes individuelles ---
    ax = axes[0]
    ax.plot(x, real(psi1), label=r'$\mathrm{Re}[\Psi_1]$ ($k_0$)',          color='royalblue')
    ax.plot(x, real(psi2), label=r'$\mathrm{Re}[\Psi_2]$ ($k_0-\Delta k/2$)', color='tomato',     linestyle='--')
    ax.plot(x, real(psi3), label=r'$\mathrm{Re}[\Psi_3]$ ($k_0+\Delta k/2$)', color='seagreen',   linestyle=':')
    ax.set_ylabel('Amplitude')
    ax.set_title('Parties réelles des 3 ondes planes individuelles')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # --- Graphique du bas : somme + enveloppe ---
    ax = axes[1]
    ax.plot(x, real(psi_sum), label=r'$\mathrm{Re}[\Psi_\mathrm{somme}]$', color='purple', linewidth=1.2)
    ax.plot(x, env_pos, label='Enveloppe $+$', color='darkorange', linestyle='--', linewidth=1.5)
    ax.plot(x, env_neg, label='Enveloppe $-$', color='darkorange', linestyle='--', linewidth=1.5)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('Amplitude')
    ax.set_title(r'Somme des 3 ondes et enveloppe — intervalle $[-\pi/\Delta k,\, \pi/\Delta k]$')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.show()


# ─────────────────────────────────────────────
# Point d'entrée
# ─────────────────────────────────────────────

if __name__ == '__main__':
    print("=== 1.1 Onde plane ===")
    plot_plane_wave(amp=1.0, k=2*pi, omega=1.0, t=0.0)

    print("=== 1.2 Superposition de 3 ondes planes ===")
    plot_superposition(amp=1.0, k0=20*pi, dk=4*pi)
