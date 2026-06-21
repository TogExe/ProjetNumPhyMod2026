# Projet Numérique en Physique Moderne — Réponses théoriques

> **CY Tech — PréIng 2 MIP**  
> Effet tunnel : étude d'une particule quantique traversant une barrière de potentiel rectangulaire.

---

## Table des matières

1. [Partie 1 — Ondes planes](#partie-1--ondes-planes)
   - [1.1 Généralités](#11-généralités)
   - [1.2 Superposition d'ondes planes](#12-superposition-dondes-planes)
2. [Partie 2 — Paquets d'ondes](#partie-2--paquets-dondes)
   - [2.1 Paquet d'ondes gaussien](#21-paquet-dondes-gaussien)
3. [Partie 3 — Résolution numérique](#partie-3--résolution-numérique)
   - [3.1 Algorithme de dérivation](#31-algorithme-de-dérivation)
   - [3.2 Algorithme pour l'équation de Schrödinger](#32-algorithme-pour-léquation-de-schrödinger)

---

## Partie 1 — Ondes planes

### 1.1 Généralités

#### 1.1.1a — Onde plane à 3 dimensions

L'expression d'une onde plane à trois dimensions est :

$$\Psi(\vec{r}, t) = A \cdot \exp\!\left[i\!\left(\vec{k} \cdot \vec{r} - \omega t\right)\right]$$

| Grandeur | Signification | Dimension | Unité SI |
|----------|--------------|-----------|----------|
| $\vec{k}$ | Vecteur d'onde — direction de propagation, $\|\vec{k}\| = 2\pi/\lambda$ | $[\text{longueur}]^{-1}$ | $\text{m}^{-1}$ |
| $\omega$ | Pulsation (fréquence angulaire), $\omega = 2\pi\nu$ | $[\text{temps}]^{-1}$ | $\text{rad·s}^{-1}$ |

---

#### 1.1.1b — Réduction à une dimension

En prenant la propagation selon l'axe $(Ox)$ :

$$\boxed{\Psi(x, t) = A \cdot e^{i(kx - \omega t)}}$$

- **Partie réelle :** $\text{Re}(\Psi) = A \cos(kx - \omega t)$
- **Partie imaginaire :** $\text{Im}(\Psi) = A \sin(kx - \omega t)$

---

#### 1.1.1c — Dimension de l'amplitude $A$

La condition de normalisation impose :

$$\int_{-\infty}^{+\infty} |\Psi(x,t)|^2 \, dx = |A|^2 \int_{-\infty}^{+\infty} dx = 1$$

Cette intégrale diverge pour une onde plane infinie, mais formellement $[A]^2 \cdot [\text{longueur}] = 1$, donc :

$$\boxed{[A] = \text{m}^{-1/2}}$$

---

#### 1.1.1d — L'onde plane est solution de l'équation de Schrödinger libre

L'équation de Schrödinger à 1D pour une particule libre ($V = 0$) :

$$i\hbar \frac{\partial \Psi}{\partial t} = -\frac{\hbar^2}{2m} \frac{\partial^2 \Psi}{\partial x^2}$$

**Calcul du membre gauche :**

$$i\hbar \frac{\partial}{\partial t} A e^{i(kx-\omega t)} = i\hbar \cdot (-i\omega) \cdot \Psi = \hbar\omega\, \Psi$$

**Calcul du membre droit :**

$$-\frac{\hbar^2}{2m} \frac{\partial^2}{\partial x^2} A e^{i(kx-\omega t)} = -\frac{\hbar^2}{2m}(ik)^2 \Psi = \frac{\hbar^2 k^2}{2m}\,\Psi$$

**Condition de compatibilité :**

$$\boxed{\hbar\omega = \frac{\hbar^2 k^2}{2m}} \quad \Longleftrightarrow \quad E = \frac{p^2}{2m}$$

Cette relation est la relation de dispersion : elle est satisfaite, donc l'onde plane est bien solution. ✓

---

#### 1.1.1e — Relation de dispersion, vitesse de phase et vitesse de groupe

$$\omega(k) = \frac{\hbar k^2}{2m}$$

| Grandeur | Expression | Valeur |
|----------|-----------|--------|
| Vitesse de phase | $v_\varphi = \dfrac{\omega}{k}$ | $\dfrac{\hbar k}{2m}$ |
| Vitesse de groupe | $v_g = \dfrac{d\omega}{dk}$ | $\dfrac{\hbar k}{m}$ |

On obtient la relation remarquable $v_g = 2\, v_\varphi$.

---

#### 1.1.1f — Comparaison avec la vitesse classique $v = p/m$

Avec $p = \hbar k$ (relation de de Broglie) :

$$v_{\text{class}} = \frac{p}{m} = \frac{\hbar k}{m}$$

- La **vitesse de groupe** $v_g = \hbar k / m$ **coïncide** avec la vitesse classique de la particule.
- La **vitesse de phase** $v_\varphi = \hbar k / 2m = v_{\text{class}}/2$ n'a pas de sens physique direct.

> La vitesse de groupe est donc bien la vitesse physique de propagation du paquet d'ondes.

---

#### 1.1.1g — Condition de normalisation

| Système | Bornes | Condition |
|---------|--------|-----------|
| Particule libre | $]-\infty, +\infty[$ | $\displaystyle\int_{-\infty}^{+\infty} \|\Psi(x,t)\|^2\,dx = 1$ |
| Puits infini (largeur $L$) | $[0, L]$ | $\displaystyle\int_{0}^{L} \|\Psi(x,t)\|^2\,dx = 1$ |

---

#### 1.1.1h — L'onde plane n'est pas physiquement acceptable

Pour une onde plane :

$$|\Psi(x,t)|^2 = |A|^2 = \text{constante}$$

L'intégrale de normalisation diverge :

$$\int_{-\infty}^{+\infty} |A|^2\, dx \longrightarrow +\infty$$

L'onde plane **ne peut pas être normalisée** sur $\mathbb{R}$ : elle ne représente pas un état physique réalisable. Elle reste un outil mathématique fondamental (base de décomposition pour les paquets d'ondes), mais ne correspond à aucun état propre physique d'une particule libre.

---

### 1.2 Superposition d'ondes planes

#### 1.2.1a — La superposition reste solution

L'équation de Schrödinger est **linéaire** : si $\Psi_1$ et $\Psi_2$ sont solutions pour le même potentiel, alors pour tout $\alpha, \beta \in \mathbb{C}$ :

$$\alpha\Psi_1 + \beta\Psi_2 \text{ est aussi solution.}$$

Toute superposition (finie ou continue) d'ondes planes est donc solution — c'est le **principe de superposition**.

---

#### 1.2.1b — Somme des trois ondes planes à $t = 0$

Les trois ondes (à $t = 0$) :

$$\Psi_1 = A\,e^{ik_0 x}, \quad \Psi_2 = \frac{A}{2}\,e^{i(k_0 - \Delta k/2)x}, \quad \Psi_3 = \frac{A}{2}\,e^{i(k_0 + \Delta k/2)x}$$

Leur somme :

$$\Psi(x,0) = A\,e^{ik_0 x}\!\left[1 + \frac{1}{2}e^{-i\Delta k\, x/2} + \frac{1}{2}e^{+i\Delta k\, x/2}\right]$$

En utilisant $e^{i\theta} + e^{-i\theta} = 2\cos\theta$ :

$$\boxed{\Psi(x,0) = A\,e^{ik_0 x}\!\left[1 + \cos\!\left(\frac{\Delta k}{2}\,x\right)\right]}$$

---

#### 1.2.1c — Densité de probabilité de présence

$$\boxed{|\Psi(x,0)|^2 = |A|^2 \left[1 + \cos\!\left(\frac{\Delta k}{2}\,x\right)\right]^2}$$

La densité n'est plus uniforme : elle présente des **maxima** espacés de $4\pi/\Delta k$ et des **nœuds** (zéros) là où $\cos(\Delta k\, x/2) = -1$. La localisation spatiale augmente avec $\Delta k$.

---

## Partie 2 — Paquets d'ondes

### 2.1 Paquet d'ondes gaussien

#### 2.1.1a — Expression générale à 1D sans $\omega$

Pour une particule libre se déplaçant selon $(Ox)$, on substitue $\omega = \hbar k^2 / 2m$ :

$$\boxed{\Psi(x,t) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{+\infty} g(k)\, \exp\!\left[i\!\left(kx - \frac{\hbar k^2}{2m}\,t\right)\right] dk}$$

---

#### 2.1.1b — Paquet d'ondes gaussien : expression

La fonction $g$ est une gaussienne centrée en $k_0$ :

$$g(k) = \sqrt{a}\,(2\pi)^{-1/4}\exp\!\left[-\frac{a^2(k-k_0)^2}{4}\right]$$

Le paquet d'ondes s'écrit alors :

$$\Psi(x,t) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{+\infty} \sqrt{a}\,(2\pi)^{-1/4}\exp\!\left[-\frac{a^2(k-k_0)^2}{4}\right] \exp\!\left[i\!\left(kx - \frac{\hbar k^2 t}{2m}\right)\right] dk$$

---

#### 2.1.1c — Résultat de l'intégrale

L'intégrale gaussienne donne le résultat analytique suivant (fourni dans le sujet) :

$$\Psi(x,t) = \left(\frac{1}{8\pi^3}\right)^{1/4} \sqrt{\frac{4\pi m a}{ma^2 + 2i\hbar t}} \exp\!\left[\frac{m}{4}\,\frac{\left(a^2 k_0 + 2ix\right)^2}{ma^2 + 2i\hbar t} - \frac{a^2 k_0^2}{4}\right]$$

Le paramètre $a$ est la **largeur initiale** du paquet dans l'espace réel : à $t = 0$, l'enveloppe gaussienne a un écart-type $\sigma_0 = a/\sqrt{2}$.

---

#### 2.1.1d — Vérification de la normalisation

À $t = 0$, le module carré du paquet gaussien s'écrit :

$$|\Psi(x,0)|^2 = \frac{1}{\sigma_0\sqrt{2\pi}}\,\exp\!\left[-\frac{x^2}{2\sigma_0^2}\right], \quad \sigma_0 = \frac{a}{\sqrt{2}}$$

C'est une **distribution gaussienne normalisée** :

$$\int_{-\infty}^{+\infty} |\Psi(x,0)|^2\,dx = 1 \quad \checkmark$$

Pour $t > 0$, l'équation de Schrödinger conserve la norme : le paquet s'élargit mais $\int |\Psi|^2\,dx = 1$ à tout instant.

---

#### 2.1.1e — Relation entre $g(k)$ et $\Psi(x, 0)$

Le paquet d'ondes est une **transformée de Fourier** :

$$\Psi(x,0) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{+\infty} g(k)\, e^{ikx}\, dk$$

Par **inversion de la transformée de Fourier** :

$$\boxed{g(k) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{+\infty} \Psi(x,0)\, e^{-ikx}\, dx}$$

$g(k)$ est la transformée de Fourier de $\Psi(x,0)$ — relation qui traduit la **dualité position-impulsion** (inégalité de Heisenberg : $\Delta x \cdot \Delta k \geq 1/2$).

---

## Partie 3 — Résolution numérique

### 3.1 Algorithme de dérivation

#### 3.1.1a — Définition de la dérivée

$$\boxed{f'(x) = \lim_{h \to 0} \frac{f(x+h) - f(x)}{h}}$$

En pratique numérique, on préfère le **schéma centré** (erreur d'ordre $h^2$ au lieu de $h$) :

$$f'(x) \approx \frac{f(x+h) - f(x-h)}{2h}$$

---

#### 3.1.1b — Algorithme : dérivée première d'un tableau 1D

```
Entrée : tableau f[0..npts-1], pas h = x[1]-x[0]
Sortie : tableau df[0..npts-1]

# Points intérieurs (schéma centré, ordre 2)
Pour i de 1 à npts-2 :
    df[i] = (f[i+1] - f[i-1]) / (2*h)

# Points de bord (schéma décentré, ordre 1)
df[0]      = (f[1] - f[0]) / h
df[npts-1] = (f[npts-1] - f[npts-2]) / h
```

> **Erreur locale :** $\mathcal{O}(h^2)$ pour les points intérieurs, $\mathcal{O}(h)$ aux bords.

---

#### 3.1.2 — Dérivée seconde

**Définition :**

$$f''(x) = \lim_{h \to 0} \frac{f(x+h) - 2f(x) + f(x-h)}{h^2}$$

**Algorithme :**

```
Entrée : tableau f[0..npts-1], pas h
Sortie : tableau d2f[0..npts-1]

# Points intérieurs
Pour i de 1 à npts-2 :
    d2f[i] = (f[i+1] - 2*f[i] + f[i-1]) / h²

# Bords : conditions aux limites imposées par le problème physique
# (ex. Ψ = 0 aux bords pour un puits infini)
d2f[0]      = 0
d2f[npts-1] = 0
```

> **Erreur locale :** $\mathcal{O}(h^2)$ — ce schéma est utilisé pour discrétiser le terme $\partial^2\Psi/\partial x^2$ dans l'équation de Schrödinger.

---

### 3.2 Algorithme pour l'équation de Schrödinger

#### 3.2.1 — Équation de Schrödinger à 1D, potentiel constant $V_0$

$$\boxed{i\hbar\,\frac{\partial \Psi(x,t)}{\partial t} = -\frac{\hbar^2}{2m}\,\frac{\partial^2 \Psi(x,t)}{\partial x^2} + V_0\,\Psi(x,t)}$$

**Discrétisation (schéma explicite — méthode d'Euler)** :

En notant $\Psi_i^n = \Psi(x_i, t_n)$ :

$$\Psi_i^{n+1} = \Psi_i^n + \frac{\Delta t}{i\hbar}\left[-\frac{\hbar^2}{2m}\,\frac{\Psi_{i+1}^n - 2\Psi_i^n + \Psi_{i-1}^n}{\Delta x^2} + V_0\,\Psi_i^n\right]$$

> **Stabilité :** le schéma explicite impose une condition CFL sur les pas : $\Delta t \ll m\,\Delta x^2 / \hbar$.  
> Pour une meilleure stabilité, on peut utiliser le **schéma de Crank-Nicolson** (implicite, inconditionnellement stable).

---

## Récapitulatif des formules clés

| Grandeur | Formule |
|----------|---------|
| Onde plane 1D | $\Psi = A\,e^{i(kx-\omega t)}$ |
| Relation de dispersion (libre) | $\omega = \hbar k^2 / 2m$ |
| Vitesse de groupe | $v_g = \hbar k/m = p/m$ |
| Vitesse de phase | $v_\varphi = \hbar k/2m = v_g/2$ |
| Paquet d'ondes | $\Psi(x,t) = \frac{1}{\sqrt{2\pi}}\int g(k)\,e^{i(kx - \hbar k^2 t/2m)}\,dk$ |
| Transformée de Fourier inverse | $g(k) = \frac{1}{\sqrt{2\pi}}\int\Psi(x,0)\,e^{-ikx}\,dx$ |
| Dérivée seconde (numérique) | $f''_i \approx (f_{i+1} - 2f_i + f_{i-1})/\Delta x^2$ |
| Équation de Schrödinger 1D | $i\hbar\,\partial_t\Psi = -\frac{\hbar^2}{2m}\partial_{xx}\Psi + V\Psi$ |

---

*Rédigé dans le cadre du projet numérique de physique moderne — CY Tech.*
