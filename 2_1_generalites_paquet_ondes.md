# 2.1 Généralités — Paquet d'ondes

L'équation de Schrödinger étant **linéaire**, si deux fonctions d'ondes sont solutions (pour un même potentiel), leur somme reste solution et constitue un état possible (superposition d'états — c'est l'expérience du chat de Schrödinger). À trois dimensions d'espace, l'expression la plus générale d'une telle superposition est un **paquet d'ondes** :

$$\Psi(\vec{r}, t) = [2\pi]^{-3/2} \iiint g(\vec{k}) \exp\!\left[i\vec{r}\cdot\vec{k} - i\omega t\right] d^3\vec{k} \tag{1}$$

où $\vec{k}$ est le vecteur d'onde et $\vec{r}$ le vecteur position (ou rayon vecteur).

---

### Analyse de l'expression

- Le **pré-facteur** $[2\pi]^{-3/2}$ est ajouté de manière à ce que $\Psi$ soit normalisable.
- La **fonction** $g$ est une fonction de trois variables, à savoir les trois composantes de $\vec{k}$.
- Le **terme en exponentiel** correspond à l'onde plane.
- L'expression générale fait intervenir une **intégrale triple** sur les trois coordonnées du vecteur $\vec{k}$.

---

### Expression en base cartésienne

Dans une base cartésienne $(\vec{e}_x, \vec{e}_y, \vec{e}_z)$, on écrit :

$$\vec{r} = x\vec{e}_x + y\vec{e}_y + z\vec{e}_z \qquad \text{et} \qquad \vec{k} = k_x\vec{e}_x + k_y\vec{e}_y + k_z\vec{e}_z$$

Le paquet d'ondes s'écrit alors :

$$\Psi(x, y, z, t) = [2\pi]^{-3/2} \iiint g(k_x, k_y, k_z) \exp\!\left[i\vec{r}\cdot\vec{k} - i\omega t\right] dk_x\, dk_y\, dk_z \tag{2}$$

avec $\vec{r}\cdot\vec{k} = xk_x + yk_y + zk_z$.

---

### Remarque sur les notations

En physique, on rencontre fréquemment la notation condensée :

$$\Psi(\vec{r}, t) = [2\pi]^{-3/2} \int g(\vec{k}) \exp\!\left[i\vec{r}\cdot\vec{k} - i\omega t\right] d^3k \tag{3}$$

ou encore, si on note $k = |\vec{k}|$ :

$$\Psi(r, t) = [2\pi]^{-3/2} \int g(k) \exp\!\left[ir\cdot k - i\omega t\right] dk \tag{4}$$

> **Note :** Les expressions (3) et (4) sont des notations abrégées de (2) ; $d^3k = dk_x\,dk_y\,dk_z$ désigne l'élément de volume dans l'espace des vecteurs d'onde.
