# 1. Notions physiques — Onde plane et particule libre

---

## a. Onde plane à trois dimensions d'espace

L'expression générale d'une **onde plane à 3D** est :

$$\Psi(\vec{r}, t) = A \, e^{i(\vec{k} \cdot \vec{r} - \omega t)}$$

### Signification, dimension physique et unité de $\vec{k}$ et $\omega$

| Grandeur | Signification physique | Dimension | Unité SI |
|----------|----------------------|-----------|----------|
| $\vec{k}$ | Vecteur d'onde : indique la direction et le nombre de périodes spatiales par unité de longueur | $[\text{L}^{-1}]$ | $\text{m}^{-1}$ |
| $\omega$ | Pulsation (fréquence angulaire) : nombre de périodes temporelles par unité de temps (en radians) | $[\text{T}^{-1}]$ | $\text{rad} \cdot \text{s}^{-1}$ |

---

## b. Expression à une dimension — Parties réelle et imaginaire

En se restreignant à **une dimension d'espace** (axe $x$), le vecteur d'onde se réduit à $k$ scalaire, et l'onde plane devient :

$$\Psi(x, t) = A \, e^{i(kx - \omega t)}$$

En développant à l'aide de la formule d'Euler $e^{i\theta} = \cos\theta + i\sin\theta$ :

$$\boxed{\Psi(x, t) = A\cos(kx - \omega t) + i \, A\sin(kx - \omega t)}$$

- **Partie réelle :** $\text{Re}[\Psi] = A\cos(kx - \omega t)$
- **Partie imaginaire :** $\text{Im}[\Psi] = A\sin(kx - \omega t)$

> 🔵 **Convention :** On travaille désormais uniquement en une dimension d'espace.

---

## c. Dimension et unité de l'amplitude $A$

La fonction d'onde $\Psi(x, t)$ doit être **normalisable**. La condition de normalisation (voir question g) impose :

$$\int_{-\infty}^{+\infty} |\Psi(x,t)|^2 \, dx = 1$$

L'intégrale $|\Psi|^2 \, dx$ est sans dimension (probabilité). Donc $|\Psi|^2$ a la dimension de l'inverse d'une longueur :

$$[|\Psi|^2] = \text{L}^{-1} \implies [\Psi] = \text{L}^{-1/2}$$

L'amplitude $A$ devant l'exponentielle a donc :

$$\boxed{[A] = \text{L}^{-1/2}, \quad \text{unité SI : m}^{-1/2}}$$

---

## d. Vérification : solution de l'équation de Schrödinger pour une particule libre

L'équation de Schrödinger pour une **particule libre** (potentiel $V = 0$) en 1D est :

$$i\hbar \frac{\partial \Psi}{\partial t} = -\frac{\hbar^2}{2m} \frac{\partial^2 \Psi}{\partial x^2}$$

**Calcul du membre de gauche :**

$$\frac{\partial \Psi}{\partial t} = -i\omega \, A e^{i(kx - \omega t)}$$

$$i\hbar \frac{\partial \Psi}{\partial t} = i\hbar (-i\omega) \Psi = \hbar\omega \, \Psi$$

**Calcul du membre de droite :**

$$\frac{\partial^2 \Psi}{\partial x^2} = (ik)^2 A e^{i(kx - \omega t)} = -k^2 \Psi$$

$$-\frac{\hbar^2}{2m} \frac{\partial^2 \Psi}{\partial x^2} = \frac{\hbar^2 k^2}{2m} \Psi$$

**Condition de cohérence :**

$$\hbar\omega = \frac{\hbar^2 k^2}{2m}$$

$$\boxed{\omega = \frac{\hbar k^2}{2m}}$$

C'est exactement la **relation de dispersion** de la particule libre. L'onde plane est donc bien solution de l'équation de Schrödinger pour une particule libre. ✓

---

## e. Relation de dispersion, vitesse de phase et vitesse de groupe

### Relation de dispersion

$$\boxed{\omega(k) = \frac{\hbar k^2}{2m}}$$

### Vitesse de phase $v_\varphi$

La vitesse de phase est la vitesse à laquelle se déplace un plan de phase constant :

$$v_\varphi = \frac{\omega}{k} = \frac{\hbar k^2 / 2m}{k} = \frac{\hbar k}{2m}$$

$$\boxed{v_\varphi = \frac{\hbar k}{2m}}$$

### Vitesse de groupe $v_g$

La vitesse de groupe est la vitesse à laquelle se propage l'enveloppe d'un paquet d'ondes :

$$v_g = \frac{d\omega}{dk} = \frac{d}{dk}\left(\frac{\hbar k^2}{2m}\right) = \frac{2\hbar k}{2m} = \frac{\hbar k}{m}$$

$$\boxed{v_g = \frac{\hbar k}{m}}$$

> On remarque que $v_g = 2 v_\varphi$.

---

## f. Comparaison avec $v = p/m$

En mécanique quantique, la relation de de Broglie donne l'impulsion :

$$p = \hbar k$$

Donc la vitesse classique d'une particule d'impulsion $p$ est :

$$v = \frac{p}{m} = \frac{\hbar k}{m}$$

**Comparaison :**

| Vitesse | Expression | Valeur |
|---------|-----------|--------|
| Phase $v_\varphi$ | $\hbar k / 2m$ | $v/2$ |
| Groupe $v_g$ | $\hbar k / m$ | $v$ |
| Classique $v = p/m$ | $\hbar k / m$ | $v$ |

$$\boxed{v_g = \frac{p}{m} = v \quad \text{et} \quad v_\varphi = \frac{v}{2}}$$

**Conclusion :** La vitesse de groupe coïncide avec la vitesse classique de la particule. C'est la vitesse de groupe qui porte l'information physique, et non la vitesse de phase.

---

## g. Condition de normalisation et bornes d'intégration

### Particule libre (étendue à tout l'espace)

La condition de normalisation s'écrit :

$$\int_{-\infty}^{+\infty} |\Psi(x,t)|^2 \, dx = 1$$

Pour une onde plane $\Psi = A e^{i(kx-\omega t)}$, on a $|\Psi|^2 = |A|^2 = \text{constante}$, donc :

$$\int_{-\infty}^{+\infty} |A|^2 \, dx = |A|^2 \cdot \infty \to \infty$$

Cette intégrale **diverge** : la particule libre ne peut pas être normalisée sur $\mathbb{R}$ entier.

### Particule dans un puits de profondeur infinie

Pour une particule confinée dans un puits infini de largeur $L$ (de $x = 0$ à $x = L$), la fonction d'onde est nulle à l'extérieur. La condition de normalisation devient :

$$\boxed{\int_{0}^{L} |\Psi(x,t)|^2 \, dx = 1}$$

Les bornes sont donc $0$ et $L$ (ou $-L/2$ et $+L/2$ selon la convention choisie pour le puits).

---

## h. Les ondes planes ne sont pas des solutions physiquement acceptables

Une solution physiquement acceptable (état quantique physique) doit être **de carré sommable** :

$$\int_{-\infty}^{+\infty} |\Psi(x,t)|^2 \, dx < +\infty$$

Or pour l'onde plane :

$$|\Psi(x,t)|^2 = |A|^2 = \text{constante} \neq 0$$

L'intégrale sur $\mathbb{R}$ est donc :

$$\int_{-\infty}^{+\infty} |A|^2 \, dx = |A|^2 \cdot (+\infty) = +\infty$$

**L'onde plane n'est pas normalisable** : elle est étendue sur tout l'espace avec une densité de probabilité uniforme, ce qui signifierait que la particule est équiprobable en tout point de l'espace, y compris à l'infini — physiquement impossible.

$$\boxed{\Psi(x,t) = A e^{i(kx - \omega t)} \notin L^2(\mathbb{R}) \implies \text{non physique}}$$

> 💡 **Remède :** On construit des **paquets d'ondes** (superpositions d'ondes planes avec une distribution $\phi(k)$ en vecteur d'onde) qui sont, eux, normalisables et représentent des états physiques réels :
>
> $$\Psi(x,t) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{+\infty} \phi(k) \, e^{i(kx - \omega(k)t)} \, dk$$
