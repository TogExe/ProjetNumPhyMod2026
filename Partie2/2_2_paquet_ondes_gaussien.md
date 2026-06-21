# 2.2 Paquets d'ondes gaussien

Un cas particulier de paquet d'ondes est celui pour lequel la fonction $g$ est une **gaussienne** (loi normale). Ce cas particulier est assez important notamment parce que les calculs sont plus simples et parce que ce type de paquet d'ondes est réalisable en laboratoire.

## En dimension une

$$g(k) = \sqrt{a}\,[2\pi]^{-1/4} \exp\!\left[-\frac{a^2}{2}\frac{(k-k_0)^2}{4}\right] \tag{5}$$

où $a$ est une grandeur qui sera interprétée par la suite.

---

## 1. Notions physiques

### a. Expression générale d'un paquet d'ondes (1D, particule libre, sans $\omega$)

Pour une particule libre se déplaçant selon l'axe $(Ox)$, en faisant abstraction de $\omega$ (i.e. à $t = 0$, ou en notation formelle sans préciser la relation de dispersion), le paquet d'ondes en 1D s'écrit comme la transformée de Fourier de $g(k)$ :

$$\boxed{\Psi(x, t) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{+\infty} g(k)\, e^{i(kx - \omega t)}\, dk}$$

Sans faire intervenir $\omega$ (expression formelle) :

$$\Psi(x) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{+\infty} g(k)\, e^{ikx}\, dk$$

---

### b. Expression générale du paquet d'ondes gaussien

En substituant l'expression (5) dans la formule générale :

$$\Psi(x, t) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{+\infty} \sqrt{a}\,[2\pi]^{-1/4} \exp\!\left[-\frac{a^2(k - k_0)^2}{4}\right] e^{i(kx - \omega t)}\, dk$$

$$\boxed{\Psi(x, t) = \frac{\sqrt{a}}{(2\pi)^{3/4}} \int_{-\infty}^{+\infty} \exp\!\left[-\frac{a^2(k-k_0)^2}{4} + i(kx - \omega t)\right] dk}$$

---

### c. Calcul de l'intégrale — Paquet d'ondes gaussien à l'instant $t$

Le calcul de l'intégrale (par complétion du carré dans l'exponentielle, en utilisant l'intégrale gaussienne $\int e^{-\alpha u^2} du = \sqrt{\pi/\alpha}$) donne :

$$\boxed{\Psi(x, t) = \left(\frac{1}{8\pi^3}\right)^{1/4} \sqrt{\frac{4\pi ma}{ma^2 + 2i\hbar t}} \exp\!\left[\frac{m}{4} \cdot \frac{\left(a^2 k_0 + 2ix\right)^2}{ma^2 + 2i\hbar t} - \frac{a^2 k_0^2}{4}\right]} \tag{6}$$

> **Note sur le calcul :** On effectue le changement de variable $u = k - k_0$, puis on complète le carré dans l'argument de l'exponentielle. La relation de dispersion de la particule libre $\omega(k) = \hbar k^2 / 2m$ est utilisée pour exprimer $\omega t$ en fonction de $k$.

---

### d. Vérification de la normalisation

Pour montrer que ce paquet d'ondes est normalisé, on calcule $\int_{-\infty}^{+\infty} |\Psi(x,t)|^2\, dx = 1$.

**Méthode :** On utilise le **théorème de Parseval**, qui stipule que la transformée de Fourier conserve la norme :

$$\int_{-\infty}^{+\infty} |\Psi(x,t)|^2\, dx = \int_{-\infty}^{+\infty} |g(k)|^2\, dk$$

Il suffit donc de vérifier que $g(k)$ est normalisée :

$$\int_{-\infty}^{+\infty} |g(k)|^2\, dk = \int_{-\infty}^{+\infty} a\,[2\pi]^{-1/2} \exp\!\left[-\frac{a^2(k-k_0)^2}{2}\right] dk$$

On reconnaît une gaussienne d'écart-type $\sigma = 1/a$. En utilisant $\int_{-\infty}^{+\infty} e^{-\alpha u^2} du = \sqrt{\pi/\alpha}$ avec $\alpha = a^2/2$ :

$$= a\,[2\pi]^{-1/2} \cdot \sqrt{\frac{2\pi}{a^2}} = a \cdot \frac{1}{\sqrt{2\pi}} \cdot \frac{\sqrt{2\pi}}{a} = 1 \checkmark$$

$$\boxed{\int_{-\infty}^{+\infty} |\Psi(x,t)|^2\, dx = 1}$$

Le paquet d'ondes gaussien est bien normalisé pour tout $t$.

---

### e. Expression de $g(k)$ en fonction de $\Psi(x, t=0)$

Par la **transformée de Fourier inverse**, si :

$$\Psi(x, 0) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{+\infty} g(k)\, e^{ikx}\, dk$$

alors :

$$\boxed{g(k) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{+\infty} \Psi(x, 0)\, e^{-ikx}\, dx}$$

**Interprétation :** $g(k)$ est la transformée de Fourier de $\Psi(x, 0)$. Le paquet d'ondes et sa distribution en vecteur d'onde forment une paire de Fourier. Cette dualité est à l'origine de la relation d'incertitude de Heisenberg :

$$\Delta x \cdot \Delta k \geq \frac{1}{2} \implies \Delta x \cdot \Delta p \geq \frac{\hbar}{2}$$

Dans le cas gaussien, cette inégalité est **saturée** (le paquet gaussien réalise le minimum d'incertitude) :

$$\Delta x \cdot \Delta p = \frac{\hbar}{2}$$
