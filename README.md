# ProjetNumPhyMod2026

Étude de l'**effet tunnel quantique** : temps de traversée d'une barrière rectangulaire de potentiel par un paquet d'ondes gaussien.

Groupe 2G — CY Tech

---

## Structure du projet

| Fichier | Partie | Contenu |
|---------|--------|---------|
| `OndePlane1d2G.py` | 1 | Onde plane 1D, superposition de 3 ondes, enveloppe |
| `PaquetOndeGauss1d2G.py` | 2 | Paquet d'ondes gaussien, densité de probabilité |
| `SchrodingerNumerique2G.py` | 3 | Dérivées numériques, résolution d'Euler explicite, validation |
| `Partie4.py` | 4 | Effet tunnel complet : numérique (Crank-Nicolson) + analytique |

---

## Partie 1 — Ondes planes

`OndePlane1d2G.py`

- Fonction `PlaneWave(amp, k, omega, x, t)` → onde plane complexe à 1D
- Tracé des parties réelle et imaginaire
- Superposition de 3 ondes planes (vecteurs d'onde $k_0$, $k_0 \pm \Delta k/2$) avec enveloppe analytique :

$$\Psi(x,0) = A\left[1 + \cos\!\left(\frac{\Delta k\, x}{2}\right)\right] e^{ik_0 x}$$

---

## Partie 2 — Paquet d'ondes gaussien

`PaquetOndeGauss1d2G.py`

- Constantes : `hbar = m = 1` (unités naturelles pour éviter les dépassements numériques)
- Fonction `GaussWP(k0, a, x, t)` → formule analytique exacte (équation 6 du projet) :

$$\Psi(x,t) = \left(\frac{1}{8\pi^3}\right)^{1/4} \sqrt{\frac{4\pi m a}{ma^2+2i\hbar t}} \exp\!\left[\frac{m}{4}\frac{(a^2k_0+2ix)^2}{ma^2+2i\hbar t} - \frac{a^2k_0^2}{4}\right]$$

- **Question 2d — Difficulté numérique :** En unités SI ($\hbar \sim 10^{-34}$, $m \sim 10^{-31}$), les exposants de l'exponentielle deviennent extrêmement grands → dépassement numérique (*overflow*).
- **Question 2e — Solution :** Utiliser les **unités naturelles** $\hbar = m = 1$, ce qui ramène tous les nombres à des valeurs proches de 1 et évite tout problème numérique.

---

## Partie 3 — Résolution numérique de Schrödinger

`SchrodingerNumerique2G.py`

### 3.1 Dérivées numériques (différences finies centrées)

- **Dérivée première :** $f'(x_i) \approx \dfrac{f(x_{i+1}) - f(x_{i-1})}{2\,dx}$
- **Dérivée seconde :** $f''(x_i) \approx \dfrac{f(x_{i+1}) - 2f(x_i) + f(x_{i-1})}{dx^2}$
- Validation sur $f(x) = x^2$ : erreur relative < $10^{-4}$

### 3.2 Schéma d'Euler explicite

$$i\hbar\,\frac{\partial\Psi}{\partial t} = -\frac{\hbar^2}{2m}\frac{\partial^2\Psi}{\partial x^2} + V_0\,\Psi$$

$$\Psi^{n+1}_i = \Psi^n_i - \frac{i\,dt}{\hbar}\left[-\frac{\hbar^2}{2m}\frac{\Psi^n_{i+1} - 2\Psi^n_i + \Psi^n_{i-1}}{dx^2} + V_0\,\Psi^n_i\right]$$

Condition de stabilité (CFL) : $dt < \dfrac{m\,dx^2}{\pi\,\hbar}$

Confrontation avec le paquet gaussien analytique ($V_0 = 0$) : erreur L2 < $10^{-3}$.

---

## Partie 4 — Effet tunnel

`Partie4.py`

### Aspects numériques (Crank-Nicolson)

- Animation de la propagation du paquet rencontrant la barrière
- Mesure de $\tau_{0,\text{num}}$ (particule libre) et $\tau_{t,\text{num}}$ (tunnel) par suivi du maximum de $|\Psi|^2$
- Influence de la largeur $a$ et de la hauteur $V_0$ sur les temps de traversée

### Aspects analytiques

- Coefficients de transmission/réflexion (états stationnaires) :

$$T = \frac{1}{1 + \dfrac{(k^2+\kappa^2)^2}{4k^2\kappa^2}\sinh^2(\kappa a)} \qquad (E < V_0)$$

- Temps de traversée libre : $\tau_{0,\text{th}} = a / v_g = ma / (\hbar k_0)$
- Temps tunnel (méthode de Wigner) : $\tau_{t,\text{th}} = \dfrac{m}{\hbar k_0}\dfrac{d\phi_t}{dk}$
- Mise en évidence de l'**effet Hartman** : $\tau_{t,\text{th}}$ sature quand $a \to \infty$

### Comparaison

| Quantité | Libre (analytique) | Tunnel (Wigner) | Classique |
|----------|-------------------|-----------------|-----------|
| $\tau$ | $a/v_g$ | $\frac{m}{\hbar k_0}\frac{d\phi_t}{dk}$ | $\infty$ si $E < V_0$ |

---

## Dépendances

```
numpy
matplotlib
```

Installation : `pip install numpy matplotlib`

## Utilisation

```bash
python OndePlane1d2G.py          # Partie 1
python PaquetOndeGauss1d2G.py    # Partie 2
python SchrodingerNumerique2G.py # Partie 3
python Partie4.py                # Partie 4 (quelques minutes)
```
