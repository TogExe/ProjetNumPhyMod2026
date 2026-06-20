1. Notions physiques
a. Onde plane à 3D : expression, signification de k⃗\vec{k}
k et ω\omega
ω
Une onde plane progressive harmonique à trois dimensions s'écrit :
Ψ(r⃗,t)=Aexp⁡[i(k⃗⋅r⃗−ωt)]\Psi(\vec{r}, t) = A \exp\left[i(\vec{k}\cdot\vec{r} - \omega t)\right]Ψ(r,t)=Aexp[i(k⋅r−ωt)]
où AA
A est l'amplitude (complexe en général).

k⃗\vec{k}
k est le vecteur d'onde. Sa direction donne la direction de propagation de l'onde, et sa norme k=∥k⃗∥k = \|\vec{k}\|
k=∥k∥ est reliée à la longueur d'onde par k=2π/λk = 2\pi/\lambda
k=2π/λ. En mécanique quantique, k⃗\vec{k}
k est relié à la quantité de mouvement par p⃗=ℏk⃗\vec{p} = \hbar \vec{k}
p​=ℏk.

Dimension : [k⃗]=L−1[\vec{k}] = L^{-1}
[k]=L−1 (inverse d'une longueur)
Unité SI : m−1\text{m}^{-1}
m−1 (rad/m si on veut être explicite sur l'angle)


ω\omega
ω est la pulsation (fréquence angulaire). Elle est reliée à la fréquence ν\nu
ν par ω=2πν\omega = 2\pi\nu
ω=2πν, et en mécanique quantique à l'énergie par E=ℏωE = \hbar\omega
E=ℏω.

Dimension : [ω]=T−1[\omega] = T^{-1}
[ω]=T−1 (inverse d'un temps)
Unité SI : rad⋅s−1\text{rad}\cdot\text{s}^{-1}
rad⋅s−1 (ou s−1\text{s}^{-1}
s−1)




b. Expression à 1D, parties réelle et imaginaire
À une dimension d'espace, k⃗⋅r⃗→kx\vec{k}\cdot\vec{r} \to kx
k⋅r→kx, donc :
Ψ(x,t)=Aexp⁡[i(kx−ωt)]\Psi(x,t) = A\exp[i(kx - \omega t)]Ψ(x,t)=Aexp[i(kx−ωt)]
En notant A=A0eiφ0A = A_0 e^{i\varphi_0}
A=A0​eiφ0​ (amplitude complexe générale), et en utilisant la formule d'Euler eiθ=cos⁡θ+isin⁡θe^{i\theta} = \cos\theta + i\sin\theta
eiθ=cosθ+isinθ :
Ψ(x,t)=A0cos⁡(kx−ωt+φ0)+iA0sin⁡(kx−ωt+φ0)\Psi(x,t) = A_0 \cos(kx - \omega t + \varphi_0) + iA_0\sin(kx-\omega t + \varphi_0)Ψ(x,t)=A0​cos(kx−ωt+φ0​)+iA0​sin(kx−ωt+φ0​)
Donc, si on prend par simplicité AA
A réelle (φ₀ = 0) :
Re[Ψ(x,t)]=Acos⁡(kx−ωt),Im[Ψ(x,t)]=Asin⁡(kx−ωt)\text{Re}[\Psi(x,t)] = A\cos(kx - \omega t), \qquad \text{Im}[\Psi(x,t)] = A\sin(kx-\omega t)Re[Ψ(x,t)]=Acos(kx−ωt),Im[Ψ(x,t)]=Asin(kx−ωt)

c. Dimension et unité de l'amplitude AA
A
La fonction d'onde Ψ(x,t)\Psi(x,t)
Ψ(x,t) doit vérifier la condition de normalisation (cf. question g) :
∫∣Ψ(x,t)∣2 dx=1\int |\Psi(x,t)|^2\,dx = 1∫∣Ψ(x,t)∣2dx=1
Cette intégrale doit être sans dimension (c'est une probabilité totale égale à 1). Donc ∣Ψ∣2|\Psi|^2
∣Ψ∣2 doit avoir la dimension de l'inverse d'une longueur (pour qu'intégrée sur dxdx
dx, le résultat soit sans dimension) :
[∣Ψ∣2]=L−1⇒[Ψ]=L−1/2[|\Psi|^2] = L^{-1} \quad \Rightarrow \quad [\Psi] = L^{-1/2}[∣Ψ∣2]=L−1⇒[Ψ]=L−1/2
L'amplitude AA
A a donc :

Dimension : L−1/2L^{-1/2}
L−1/2
Unité SI : m−1/2\text{m}^{-1/2}
m−1/2

*(Remarque : c'est cohérent avec le fait qu'une onde plane pure n'est en réalité pas normalisable sur R\mathbb{R}
R tout entier — voir question h — mais formellement, l'analyse dimensionnelle reste celle-ci.)*

d. Vérifier que Ψ(x,t)\Psi(x,t)
Ψ(x,t) est solution de l'équation de Schrödinger (particule libre)
L'équation de Schrödinger pour une particule libre (potentiel V=0V=0
V=0) à 1D s'écrit :
iℏ∂Ψ∂t=−ℏ22m∂2Ψ∂x2i\hbar \frac{\partial \Psi}{\partial t} = -\frac{\hbar^2}{2m}\frac{\partial^2 \Psi}{\partial x^2}iℏ∂t∂Ψ​=−2mℏ2​∂x2∂2Ψ​
On pose Ψ(x,t)=Aexp⁡[i(kx−ωt)]\Psi(x,t) = A\exp[i(kx-\omega t)]
Ψ(x,t)=Aexp[i(kx−ωt)]. Calculons chaque membre.
Membre de gauche :
∂Ψ∂t=−iωΨ⇒iℏ∂Ψ∂t=iℏ(−iω)Ψ=ℏω Ψ\frac{\partial \Psi}{\partial t} = -i\omega \Psi \quad\Rightarrow\quad i\hbar\frac{\partial \Psi}{\partial t} = i\hbar(-i\omega)\Psi = \hbar\omega\,\Psi∂t∂Ψ​=−iωΨ⇒iℏ∂t∂Ψ​=iℏ(−iω)Ψ=ℏωΨ
Membre de droite :
∂Ψ∂x=ikΨ⇒∂2Ψ∂x2=(ik)2Ψ=−k2Ψ\frac{\partial \Psi}{\partial x} = ik\Psi \quad\Rightarrow\quad \frac{\partial^2 \Psi}{\partial x^2} = (ik)^2\Psi = -k^2\Psi∂x∂Ψ​=ikΨ⇒∂x2∂2Ψ​=(ik)2Ψ=−k2Ψ
⇒−ℏ22m∂2Ψ∂x2=−ℏ22m(−k2Ψ)=ℏ2k22mΨ\Rightarrow\quad -\frac{\hbar^2}{2m}\frac{\partial^2\Psi}{\partial x^2} = -\frac{\hbar^2}{2m}(-k^2\Psi) = \frac{\hbar^2 k^2}{2m}\Psi⇒−2mℏ2​∂x2∂2Ψ​=−2mℏ2​(−k2Ψ)=2mℏ2k2​Ψ
Égalité des deux membres : l'équation est vérifiée si et seulement si
ℏω=ℏ2k22m\hbar\omega = \frac{\hbar^2 k^2}{2m}ℏω=2mℏ2k2​
C'est précisément la relation de dispersion de la particule libre (cf. question e). Donc Ψ\Psi
Ψ est bien solution, à condition que ω\omega
ω et kk
k vérifient cette relation.

e. Relation de dispersion, vitesse de phase, vitesse de groupe
Relation de dispersion (obtenue ci-dessus) :
ω(k)=ℏk22m\omega(k) = \frac{\hbar k^2}{2m}ω(k)=2mℏk2​
Vitesse de phase vφv_\varphi
vφ​ : vitesse de déplacement des plans de phase constante (où kx−ωt=cstekx-\omega t = \text{cste}
kx−ωt=cste) :
vφ=ωk=ℏk2mv_\varphi = \frac{\omega}{k} = \frac{\hbar k}{2m}vφ​=kω​=2mℏk​
Vitesse de groupe vgv_g
vg​ : vitesse de déplacement de l'enveloppe d'un paquet d'ondes, donnée par :
vg=dωdk=ℏkmv_g = \frac{d\omega}{dk} = \frac{\hbar k}{m}vg​=dkdω​=mℏk​
On remarque que vg=2vφv_g = 2v_\varphi
vg​=2vφ​ pour cette relation de dispersion quadratique.

f. Comparaison avec v=p/mv = p/m
v=p/m
En mécanique quantique, p=ℏkp = \hbar k
p=ℏk (relation de de Broglie). Donc :
v=pm=ℏkmv = \frac{p}{m} = \frac{\hbar k}{m}v=mp​=mℏk​
On constate que :
vg=v=ℏkmtandis quevφ=v2≠v\boxed{v_g = v = \frac{\hbar k}{m}} \qquad \text{tandis que} \qquad v_\varphi = \frac{v}{2} \neq vvg​=v=mℏk​​tandis quevφ​=2v​=v
Interprétation physique : c'est la vitesse de groupe qui correspond à la vitesse classique de la particule (vitesse de déplacement de l'enveloppe du paquet d'ondes, donc de la probabilité de présence), et non la vitesse de phase. C'est un résultat important : la vitesse de phase n'a pas de sens physique direct en termes de "vitesse de la particule".

g. Condition de normalisation
Pour une particule libre sur R\mathbb{R}
R (espace infini), la condition de normalisation s'écrit :
∫−∞+∞∣Ψ(x,t)∣2 dx=1\int_{-\infty}^{+\infty} |\Psi(x,t)|^2\, dx = 1∫−∞+∞​∣Ψ(x,t)∣2dx=1
Cas d'un puits de profondeur infinie : si la particule est confinée dans un puits [0,L][0, L]
[0,L] (ou [−L/2,L/2][-L/2, L/2]
[−L/2,L/2] selon convention), la fonction d'onde est nulle en dehors du puits, donc les bornes de l'intégrale deviennent :
∫0L∣Ψ(x,t)∣2 dx=1\int_{0}^{L} |\Psi(x,t)|^2\, dx = 1∫0L​∣Ψ(x,t)∣2dx=1
(ou ∫−L/2L/2\int_{-L/2}^{L/2}
∫−L/2L/2​ selon le repère choisi), au lieu de −∞-\infty
−∞ à +∞+\infty
+∞.

h. Pourquoi l'onde plane libre n'est pas physiquement acceptable
On calcule ∣Ψ(x,t)∣2|\Psi(x,t)|^2
∣Ψ(x,t)∣2 pour l'onde plane Ψ(x,t)=Aexp⁡[i(kx−ωt)]\Psi(x,t) = A\exp[i(kx-\omega t)]
Ψ(x,t)=Aexp[i(kx−ωt)] :
∣Ψ(x,t)∣2=ΨΨ∗=Aei(kx−ωt)⋅A∗e−i(kx−ωt)=∣A∣2|\Psi(x,t)|^2 = \Psi\Psi^* = A e^{i(kx-\omega t)} \cdot A^* e^{-i(kx-\omega t)} = |A|^2∣Ψ(x,t)∣2=ΨΨ∗=Aei(kx−ωt)⋅A∗e−i(kx−ωt)=∣A∣2
Ce module au carré est constant, indépendant de xx
x et de tt
t. Dès lors, l'intégrale de normalisation diverge :
∫−∞+∞∣A∣2 dx→∞(sauf si A=0)\int_{-\infty}^{+\infty} |A|^2\, dx \to \infty \quad (\text{sauf si } A=0)∫−∞+∞​∣A∣2dx→∞(sauf si A=0)
Conséquences physiques :

Il est impossible de normaliser Ψ\Psi
Ψ à 1 sur R\mathbb{R}
R (sauf cas trivial nul), donc on ne peut pas interpréter ∣Ψ∣2|\Psi|^2
∣Ψ∣2 comme une véritable densité de probabilité au sens strict.
Le fait que ∣Ψ∣2|\Psi|^2
∣Ψ∣2 soit uniforme sur tout l'espace signifie que la particule aurait une probabilité de présence identique en tout point de l'espace, de −∞-\infty
−∞ à +∞+\infty
+∞ : la particule serait totalement délocalisée, ce qui ne correspond à aucune particule physique réelle (qui doit être localisée dans une région finie, même grande).
