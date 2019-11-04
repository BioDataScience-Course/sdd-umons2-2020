# Modèle linéaire {#mod-lineaire}




##### Objectifs {-}

- Comprendre le modèle linéaire (ANOVA et régression linéaire tout en un)

- Appréhender la logique des matrices de contraste

- Découvrir l'ANCOVA

- Comprendre le mécanisme du modèle linéaire généralisé


##### Prérequis {-}

- L'ANOVA (modules 10 & 11 du cours [SDD 1](http://biodatascience-course.sciviews.org/sdd-umons/)), ainsi que la régression linéaires (modules 1 et 2 du présent cours) doivent être maitrisés avant d'aborder cette matière.


## Variables numériques ou facteurs

L'ANOVA analyse une **variable dépendante numérique** en fonction d'une ou plusieurs **variables indépendantes qualitatives**. Ces variables sont dites "facteurs" non ordonnés (objets de classe `factor`), ou "facteurs" ordonnés (objets de classe `ordered`) dans R. 

La régression linéaire analyse une **variable dépendante numérique** en fonction d'une ou plusieurs **variables indépendantes numérique** (quantitatives) également. Ce sont des objets de classe `numeric` (ou éventuellement `integer`, mais assimilé à `numeric` concrètement) dans R.

Donc, la principale différence entre ANOVA et régression linéaire telles que nous les avnos abordés jusqu'ici réside dans la **nature** de la ou des variables indépendantes, c'est-à-dire, leur type. Pour rappel, il existe deux grandes catégories de variables\ : quantitatives et qualitatives, et deux sous-catégories pour chacune d'elle. Cela donne quatyre types principaux de variables, formant plus de 90% des cas rencontrés\ :

- variables quantitatives continues représentables par des nombres réels (`numeric` dans R),

- variables quantitatives discrètes pour des dénombrements d'événements finis par exemple, et représentables par des nombres entiers (`integer` dans R),

- variables qualitatives ordonnées pour des variables prenant un petit nombre de valeurs, mais pouvant être ordonnées de la plus petite à la plus grande (`ordered` dans R),

- variables qualitatives non ordonnées prenant également un petit nombre de valeurs possibles, mais sans ordre particulier (`factor` dans R).

\BeginKnitrBlock{warning}<div class="warning">
Par la suite, un encodage correct des variables sera *indispensable* afin de distinguer correctement ces différentes situations. En effet, R considèrera automatiquement comment mener l'analyse en fonction de la classe des variables fournies. Donc, si la classe est incorrecte, l'analyse le sera aussi\ ! Si vous avez des doutes concernant les types de variables, relisez la section [type de variables](http://biodatascience-course.sciviews.org/sdd-umons/types-de-variables.html) avant de continuer ici. 
</div>\EndKnitrBlock{warning}


## ANOVA et régression linéaire

Avez-vous remarqué une ressemblance particulière entre la régression linéaire que nous avons réalisé précédement et l'analyse de variance (ANOVA)\ ? Les plus observateurs auront mis en avant que la fonction de base dans R est la même dans les deux cas\ : `lm()`. Cette fonction est donc capable de traiter aussi bien des variables réponses qualitatives que quantitatives, et effectue alors une ANOVA dans un cas ou une régression linéaire dans l'autre.

Par ailleurs, nous avons vu que l’ANOVA et la régression linéaire se représentent par des modèles semblables\ :

- $y  = \mu + \tau_i + \epsilon$ pour l’ANOVA et

- $y = \beta_1 + \beta_2 x + \epsilon$ pour la régression linéaire, avec

- $\epsilon \sim \mathcal{N}(0, \sigma)$ dans les deux cas.

Donc, nous retrouvons bien au niveau du modèle mathématique sous-jacent la différence principale entre les deux qui réside dans le type de variable indépendante (ou explicative)\ :

- Variable **qualitative** pour l’ANOVA,
- Variable **quantitative** pour la régression linéaire.

Le calcul est, en réalité, identique en interne. Il est donc possible de généraliser ces deux approches en une seule appelée **modèle linéaire**, mais à condition d'utiliser une astuce pour modifier nos modèles afin qu'ils soient intercompatibles.


### Modèle linéaire commun

Le nœud du problème revient donc à transformer nos modèles mathématiques pour qu'ils puissent être fusionnés en un seul. Comment homogénéiser ces deux modèles\ ?

- $y = \mu + \tau_i + \epsilon$ pour l’ANOVA et

- $y = \beta_1 + \beta_2 x + \epsilon$ pour la régression linéaire.

Avant de poursuivre, réfléchisser un peu par vous-même. Quelles sont les différences qu'il faut contourner\ ? Est-il possible d'effectuer une ou plusieurs transformations des variables pour qu'elles se comportent de manière similaire dans les deux cas\ ?


### Réencodage des variables de l'ANOVA

Considérons dans un premier temps, un cas très simple\ : une ANOVA à un facteur avec une variable indépendante qualitative (`factor`) à deux niveaux^[Concrètement, un cas aussi simple se traite habituellement à l'aide d'un test *t* de Student, mais pour notre démonstration, nous allons considérer ici utiliser une ANOVA à un facteur plutôt.]. Nous pouvons écrire\ :

$$
y = \mu + \tau_1 I_1 + \tau_2 I_2 + \epsilon
$$

avec $I_i$, une variable dite **indicatrice** créée de toute pièce qui prend la valeur 1 lorsque le niveau
correspond à _i_, et 0 dans tous les autres cas. Vous pouvez vérifier par vous-même que l'équation ci-dessus fonctionnera exactement de la même manière que le modèle utilisé jusqu'ici pour l'ANOVA. En effet, poiur un individu de la population 1, $I_1$ vaut 1 et $\tau_1$ est utilisé, alors que comme $I_2$ vaut 0, $\tau_2$ est annulé dans l'équation car $\tau_2 I_2$ vaut également 0. Et c'est exactement l'inverse qui se produit pour un individu de la population 2, de sorte que c'est $\tau_2$ qui est utilisé cette fois-ci.

Notez que notre nouvelle formulation, à l'aide de variables indicatrices ressemble fortement à la régression linéaire. La seule différence par rapport à cette dernière est que nos variables $I_i$ ne peuvent prendre que des valeurs 0 ou 1 (en tous cas, pour l'instant), alors que les $x_i$ dans la régression linéaire multiple sont des variables quantitatives qui peuvent prendre une infinité de valeurs différentes (nombres réels).

Nouys pouvons encore réécrire notre équation comme suit pour qu'elle se rapproche encore plus de celle de la régression linéaire simple. Passons par l'introduction de deux termes identiques $\tau_1 I_2$ additionné et soustrait, ce qui revient au même qu'en leur absence\ :

$$
y = \mu + \tau_1 I_1 + \tau_1 I_2 - \tau_1 I_2 + \tau_2 I_2 + \epsilon
$$

- En considérant $\beta_2 = \tau_2 - \tau_1$, cela donne\ :

$$
y = \mu + \tau_1 I_1 + \tau_1 I_2 + \beta_2 I_2 + \epsilon
$$

- En considérant $\beta_1 = \mu + \tau_1 = \mu + \tau_1 I_1 + \tau_1 I_2$ (car quelle que soit la population à laquelle notre individu appartient, il n'y a jamais qu'une seule des deux valeurs $I_1$ ou $I_2$ non nulle et dans tous les cas le résultat est donc égal à $\tau_1$), on obtient\ :

$$
y = \beta_1 + \beta_2 I_2 + \epsilon
$$

Cette dernière formulation est strictement équivalente au modèle de la régression linéaire simple dans laquelle la variable $x$ a simplement été remplacée par notre variable indicatrice $I_2$. Ceci se généralise pour une variable indépendante à $k$ niveaux, avec $k - 1$ variables indicatrices au final.

\BeginKnitrBlock{note}<div class="note">
En prenant soin de réencoder le modèle de l'ANOVA relatif aux variables indépendantes qualitatives, nous pouvons à présent mélanger les termes des deux modèles en un seul\ : notre fameux modèle linéaire. Nous aurons donc, quelque chose du genre (avec les $x_i$ correspondant aux variables quantitatives et les $I_j$ des variables indicatrices pour les différents niveaux des variables qualitatives)\ :

$$
y = \beta_1 + \beta_2 x_1 + \beta_3 x_2 + ... + \beta_n I_1 + \beta_{n+1} I_2 ... + \epsilon  
$$
</div>\EndKnitrBlock{note}


## Matrice de contraste

Donc, pour les variables qualitatives, nous considérons un ensemble de variables indicatrices (dans le cas précédent, la moyenne correspondant au premier niveau est considérée comme valeur de référence pour toutes les autres, et les variables indicatrices pour toutes les autres prennent la valeur de 1 séparément à chaque fois que le niveau correspondant est rencontré (_k_ niveaux)\ :

$$
y = \beta_1 + \beta_2 I_2 + \beta_3 I_3 + ... + \beta_k I_k + \epsilon
$$

Il s’agit de ce qu’on appelle une **matrice de contrastes** de type traitement (voir l'instruction `contr.treatment(4)` pour générer la matrice à 4 niveaux *-en ligne les niveaux, en colonne les valeurs que prennent les $I_{k-1}$ variables indicatrices-*).

Les contrastes doivent être de préférence **orthogonaux par rapport à l’ordonnée à l’origine**, ce qui signifie que la somme de leurs pondérations doit être nulle pour tous les contrastes définis (*donc, en colonnes*). Or, les contrastes de type traitement ne sont pas orthogonaux.


### Autres matrices de contrastes courantes

- Somme à zéro\ : `contr.sum(4)`

- Matrice de contrastes de Helmert\ : chaque niveau est comparé à la
moyenne des niveaux précédents\ : `contr.helmert(4)`

- Matrice de contrastes polynomiaux\ : adapté aux facteurs ordonnés pour
lesquels on s’attend à une certaine évolution du modèle du niveau le
plus petit au plus grand\ : `contr.poly(4)`


```r
plot(contr.poly(10)[, 1], type = "b")
plot(contr.poly(10)[, 2], type = "b")
plot(contr.poly(10)[, 3], type = "b")
```

R utilise par défaut des **contrastes de traitement pour les facteurs non ordonnés** et des **contrastes polynomiaux pour des facteurs ordonnés**. voir\ : `getOption("contrasts")`.


## ANCOVA

Avant l'apparition du modèle linéaire, une version particulière d'un mélange de régression linéaire et d'une ANOVA avec une variable indépendante quantitative et une autre variable indépendante qualitative s'appelait une ANCOVA (ANalyse de la COVaraiance). Un tel modèle d'ANCOVA peut naturellement également se résoudre à l'aide de la fonction `lm()` qui, en outre, peut faire bien plus.


### Exemple

Masse de nouveaux nés en fonction du poids de la mère et du fait qu’elle fume ou non.


```r
SciViews::R
babies <- read("babies", package = "UsingR")
# wt = masse du bébé à la naissance en onces et 999 = valeur manquante
# wt1 = masse de la mère à la naissance en livres et 999 = valeur manquante
# smoke = 0 (non), = 1 (oui), = 2 (jusqu'à grossesse),
#       = 3 (plus depuis un certain temps) and = 9 (inconnu)
babies %>.% select(., wt, wt1, smoke) %>.%
  filter(., wt1 < 999, wt < 999, smoke < 9) %>.%
  mutate(., wt = wt * 0.02835) %>.%
  mutate(., wt1 = wt1 * 0.4536) -> Babies
Babies$smoke <- as.factor(Babies$smoke)
# Descriptions graphiques
boxplot(data = Babies, wt ~ smoke)
boxplot(data = Babies, wt1 ~ smoke)
# ANCOVA
anova(lm(data = Babies, wt ~ smoke * wt1))
```


## Modèle linéaire généralisé

Le modèle linéaire nous a permis de *\alert{**généraliser la régression linéaire multiple** (applicable seulement sur des variables quantitatives) à des variables réponses qualitatives grâce aux variables indicatrices $I_i$. Le **modèle linéaire généralisé** reprend cette idée, mais permet en plus d’avoir d’autres variables dépendantes (ou réponses) que quantitatives, ou avec des distributions des résidus différentes. Dans R, c’est la fonction `glm()` qui se charge de calculer un modèle linéaire généralisé.

Nous rajoutons une **fonction de lien** *f*(*y*) qui transforme la variable initiale en une variable quantitative dont la relation avec les variables explicatives est linéarisée\ :

$$
f(y) = \beta_1 + \beta_2 I_2 + \beta_3 I_3 + ... + \beta_k I_k + \beta_l x_1 + \beta_m x_2 + ... + \epsilon
$$

Par exemple, pour une variable réponse binaire (distribution binomiale), avec une réponse de type logistique

$$
y = 1/(1 + e^{- \beta x})
$$

la transformation **logit** est une bonne fonction de lien\ : $\ln(y / (1 - y)) = \beta x$.


### Exemple

Recherche d’effet de variables qualitatives et quantitatives sur une réponse binaire\ :


```r
SciViews:: R
babies <- read("babies", package = "UsingR")
# Transformation : garder aussi la variable 'gestation' en jours
# et 'ht', la taille de la mère en pouces à convertir en m (/ 39.37)
babies %>.% select(., gestation, smoke, wt1, ht) %>.%
  filter(., gestation < 999, smoke < 9, wt1 < 999, ht < 999) %>.%
  # Transformer wt1 en kg et ht en cm
  mutate(., wt1 = wt1 * 0.4536) %>.%
  mutate(., ht = ht / 39.37) -> Babies_prem
Babies_prem$smoke <- as.factor(Babies_prem$smoke)
# Déterminer quels sont les enfants prématurés (nés avant 37 semaines)
Babies_prem$premat <- as.factor(as.numeric(Babies_prem$gestation < 7*37))
# BMI peut être plus parlant que la masse pour la mère?
Babies_prem %>.%
  mutate(bmi = wt1 / ht^2) -> Babies_prem

# Modèle linéaire généralisé avec fonction de lien de type logit
summary(glm(data = Babies_prem, premat ~ smoke + bmi,
  family = binomial(link = logit)))
```
