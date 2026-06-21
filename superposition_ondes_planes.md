# 1.2 Superposition d'ondes planes

## 1. Notions physiques

---

## a. La superposition reste solution de l'équation de Schrödinger

L'équation de Schrödinger pour une particule libre en 1D est **linéaire** :

$$i\hbar \frac{\partial \Psi}{\partial t} = -\frac{\hbar^2}{2m} \frac{\partial^2 \Psi}{\partial x^2}$$

Si $\Psi_1$ et $\Psi_2$ sont deux solutions, alors pour tous scalaires $\lambda_1, \lambda_2 \in \mathbb{C}$ :

$$\Psi = \lambda_1 \Psi_1 + \lambda_2 \Psi_2$$

est aussi solution, car les opérateurs de dérivation sont linéaires :

$$i\hbar \frac{\partial \Psi}{\partial t} = \lambda_1 \underbrace{i\hbar \frac{\partial \Psi_1}{\partial t}}_{= \hat{H}\Psi_1} + \lambda_2 \underbrace{i\hbar \frac{\partial \Psi_2}{\partial t}}_{= \hat{H}\Psi_2} = \hat{H}(\lambda_1\Psi_1 + \lambda_2\Psi_2) = \hat{H}\Psi$$

**Conclusion :** Par le **principe de superposition** (hérité de la linéarité de l'équation), toute combinaison linéaire de solutions est encore solution.

---

## b. Expression de l'onde résultante à $t = 0$

### Données

Les trois ondes planes à $t = 0$ ont pour vecteurs d'onde et amplitudes :

| Onde | Vecteur d'onde | Amplitude |
|------|---------------|-----------|
| $\Psi_1$ | $k_0$ | $A$ |
| $\Psi_2$ | $k_0 - \Delta k/2$ | $A/2$ |
| $\Psi_3$ | $k_0 + \Delta k/2$ | $A/2$ |

À $t = 0$ :

$$\Psi_j(x, 0) = A_j \, e^{i k_j x}$$

### Somme des trois ondes

$$\Psi(x, 0) = A e^{ik_0 x} + \frac{A}{2} e^{i(k_0 - \Delta k/2)x} + \frac{A}{2} e^{i(k_0 + \Delta k/2)x}$$

On factorise $A e^{ik_0 x}$ :

$$\Psi(x, 0) = A e^{ik_0 x} \left[ 1 + \frac{1}{2} e^{-i\Delta k \, x/2} + \frac{1}{2} e^{+i\Delta k \, x/2} \right]$$

On reconnaît la formule d'Euler : $\cos\theta = \frac{e^{i\theta} + e^{-i\theta}}{2}$

$$\Psi(x, 0) = A e^{ik_0 x} \left[ 1 + \cos\!\left(\frac{\Delta k \, x}{2}\right) \right]$$

$$\boxed{\Psi(x, 0) = A \left[1 + \cos\!\left(\frac{\Delta k \, x}{2}\right)\right] e^{ik_0 x}}$$

---

## c. Densité de probabilité de présence à $t = 0$

La densité de probabilité est $\rho(x) = |\Psi(x, 0)|^2$.

Puisque $|e^{ik_0 x}|^2 = 1$ et que le facteur entre crochets est **réel** :

$$\rho(x) = |\Psi(x, 0)|^2 = A^2 \left[1 + \cos\!\left(\frac{\Delta k \, x}{2}\right)\right]^2$$

En développant avec $1 + \cos\theta = 2\cos^2(\theta/2)$ :

$$1 + \cos\!\left(\frac{\Delta k \, x}{2}\right) = 2\cos^2\!\left(\frac{\Delta k \, x}{4}\right)$$

$$\boxed{\rho(x) = 4A^2 \cos^4\!\left(\frac{\Delta k \, x}{4}\right)}$$

**Propriétés :**
- Maximum en $x = 0$ : $\rho(0) = 4A^2$
- Premiers zéros en $x = \pm \dfrac{2\pi}{\Delta k}$
- L'enveloppe est $4A^2 \cos^4\!\left(\dfrac{\Delta k \, x}{4}\right)$, modulée par la porteuse $e^{ik_0 x}$

---

## d. Représentation graphique

> Le graphique représente, sur l'intervalle $\left[-\dfrac{\pi}{\Delta k},\, \dfrac{\pi}{\Delta k}\right]$ :
> - La partie réelle de chaque onde plane ($\Psi_1$, $\Psi_2$, $\Psi_3$)
> - La partie réelle de leur somme $\Psi$
> - L'enveloppe $\pm A\left[1 + \cos\!\left(\dfrac{\Delta k \, x}{2}\right)\right]$

### Parties réelles des trois ondes individuelles

$$\text{Re}[\Psi_1] = A \cos(k_0 x)$$

$$\text{Re}[\Psi_2] = \frac{A}{2} \cos\!\left[\left(k_0 - \frac{\Delta k}{2}\right)x\right]$$

$$\text{Re}[\Psi_3] = \frac{A}{2} \cos\!\left[\left(k_0 + \frac{\Delta k}{2}\right)x\right]$$

### Partie réelle de la somme

$$\text{Re}[\Psi(x,0)] = A \left[1 + \cos\!\left(\frac{\Delta k \, x}{2}\right)\right] \cos(k_0 x)$$

### Enveloppe

L'enveloppe de la partie réelle est donnée par le facteur de modulation (toujours positif sur l'intervalle considéré) :

$$\mathcal{E}(x) = \pm A \left[1 + \cos\!\left(\frac{\Delta k \, x}{2}\right)\right]$$

Elle atteint son maximum $2A$ en $x = 0$ et s'annule aux bords $x = \pm \dfrac{\pi}{\Delta k}$ :

$$\mathcal{E}\!\left(\pm\frac{\pi}{\Delta k}\right) = \pm A\left[1 + \cos\!\left(\frac{\pi}{2}\right)\right] = \pm A \cdot 1 = \pm A$$

> ⚠️ Sur l'intervalle $\left[-\dfrac{\pi}{\Delta k}, \dfrac{\pi}{\Delta k}\right]$, l'enveloppe ne s'annule pas encore (les zéros se trouvent en $\pm \dfrac{2\pi}{\Delta k}$). On observe donc **un seul paquet central** bien localisé, correspondant à une particule dont la position est approximativement connue avec une incertitude $\Delta x \sim \dfrac{1}{\Delta k}$, cohérent avec la relation d'Heisenberg $\Delta x \cdot \Delta k \gtrsim 1$.

---

### Schéma qualitatif

```
  Amplitude
     ^
  2A |          *  ← somme Ψ(x,0)
     |        *   *
   A |      *       *         ← enveloppe +A[1+cos(...)]
     |    *           *
   0 +--*-----------------*---> x
     | -π/Δk           +π/Δk
  -A |    *           *
     |      *       *         ← enveloppe -A[1+cos(...)]
 -2A |        *   *
     |          *
     |
     (oscillations rapides à k₀ visibles à l'intérieur)
```

**Interprétation physique :** La superposition de trois ondes planes proches (en $k$) crée un **paquet d'ondes** localisé. Plus $\Delta k$ est grand, plus le paquet est étroit en espace (meilleure localisation), illustrant la relation d'incertitude de Heisenberg : $\Delta x \cdot \Delta p \gtrsim \hbar/2$.
