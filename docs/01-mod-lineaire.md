# (PART) Cours II: analyse et modélisation {-}

# Modèle linéaire {#lm}



##### Objectifs {-}

- TODO

##### Prérequis {-}

Avant de poursuivre, vous allez réaliser une séance d'exercice couvrant les points essentiels des notions abordées dans le livre [science des données biologiques partie 1](http://biodatascience-course.sciviews.org/sdd-umons/). 

\BeginKnitrBlock{bdd}<div class="bdd">Démarrez la SciViews Box et RStudio. Dans la fenêtre **Console** de RStudio, entrez l'instruction suivante suivie de la touche `Entrée` pour ouvrir le tutoriel concernant les bases de R :

    BioDataScience2::run("01a_rappel")

N’oubliez pas d’appuyer sur la touche ESC pour reprendre la main dans R à la fin d’un tutoriel dans la console R</div>\EndKnitrBlock{bdd}


## La régression linéaire

Nous allons découvrir les base du modèle linéaire de façon intuitive. Nous utilisons le jeu de données `trees` qui rassemble la mesure du diamètre, de la hauteur et du volume de bois de cerisiers noirs. 


```r
# importation des données
trees <- read("trees", package = "datasets", lang = "fr")
```

Rapellons nous que dans le [chapitre 12 du livre science des données 1](http://biodatascience-course.sciviews.org/sdd-umons/association-de-deux-variables.html), nous avons étudié l'association de deux variables numériques. Nous utilisons donc une matrice de corrélation afin de mettre en évidence la corrélation entre nos 3 variables qui composent le jeu de donnée `trees`.

La fonction `correlation()` nous renvoie un tableau de la matrice de correlation avec l'indice de Pearson.


```r
(trees_corr <- correlation(trees))
```

```
# Matrix of Pearson's product-moment correlation:
# (calculation uses everything)
#          diameter height volume
# diameter 1.000    0.519  0.967 
# height   0.519    1.000  0.597 
# volume   0.967    0.597  1.000
```

Nous pouvons également observer cette matrice sous la forme d'un graphique plus convivial.


```r
plot(trees_corr, type = "lower")
```

<img src="01-mod-lineaire_files/figure-html/unnamed-chunk-4-1.png" width="672" style="display: block; margin: auto;" />

Cependant, n'oubliez pas qu'il est indispensable de visualiser les nuages de points pour ne pas tobmer dans le piège mis en avant par le [jeu de données artificiel appelé “quartet d’Anscombe”](http://biodatascience-course.sciviews.org/sdd-umons/association-de-deux-variables.html#importance-des-graphiques) qui montre très bien comment des données très différentes peuvent avoir même moyenne, même variance et même coefficient de corrélation.


```r
GGally::ggscatmat(as.data.frame(trees), 1:3)
```

<img src="01-mod-lineaire_files/figure-html/unnamed-chunk-5-1.png" width="672" style="display: block; margin: auto;" />

Nous observons une forte corrélation linéaire entre le volume et la hauteur des cerisiers noirs.  Interessons nous à cette association.


```r
chart(trees, volume ~ diameter) +
  geom_point()
```

<img src="01-mod-lineaire_files/figure-html/unnamed-chunk-6-1.png" width="672" style="display: block; margin: auto;" />

Si vous deviez ajouter une droite permettant de représenter au mieux les données, où est ce que vous la placeriez ? 

Une droite respecte l'équation mathématique suivante: $$ y = ax + b$$ 

dont  `a` est la pente (*slope* en anglais) et `b` est l'ordonnée à l'origine (*intercept* en anglais).


```r
# Sélection de pente et d'ordonnée à l'origine
models <- tibble(
  model = paste("mod", 1:4, sep = "-"),
  slope = c(5, 5.6, 6, 0),
  intercept = c(-0.5, -1, -1.5, 1)
)

chart(trees, volume ~ diameter) +
  geom_point() +
  geom_abline(data = models, aes(slope = slope, intercept = intercept, color = model)) + 
  labs( color = "Modèle")
```

<img src="01-mod-lineaire_files/figure-html/unnamed-chunk-7-1.png" width="672" style="display: block; margin: auto;" />

Nous avons 4 droites qui veulent représenter au mieux les observations. Quel est le meilleur modèle selon vous ? 

Nous voulons identifier le meilleur modèle, c'est à dire le modèle le plus proche de nos données. Nous avons besoin d'une règle qui va nous permettre de quantifier la distance de nos observations à notre modèle afin d'obtenir le modèle avec la plus faible distance possible de l'ensemble de nos observations. 

Décomposons le problème étape par étape et intéressons nous au `mod-1`

- Connaitre les valeurs de y prédite par le modèle 



```r
# Calculer la valeur de y pour chaque valeur de x suivant le model souhaité
## Création de notre fonction 
model <- function(slope, intercept, x) {
  prediction <- intercept + slope * x
  attributes(prediction) <- NULL
  prediction
}
## Application de notre fonction 
mod1 <- model(slope = 5, intercept = -0.5, x = trees$diameter)
## Affichage des résultats
mod1
```

```
#  [1] 0.555 0.590 0.620 0.835 0.860 0.870 0.895 0.895 0.910 0.920 0.935
# [12] 0.950 0.950 0.985 1.025 1.140 1.140 1.190 1.240 1.255 1.280 1.305
# [23] 1.340 1.530 1.570 1.695 1.720 1.775 1.785 1.785 2.115
```

- Connaitre la distance entre les observations mesurées en y et les observations prédites en y par le modèle

Nous pouvons premièrement visualiser cette distance graphiquement :


```r
chart(trees, volume ~ diameter) +
  geom_point() +
  geom_abline(slope = 5, intercept = -0.5) +
  geom_segment(
    aes(x = diameter, y = volume, 
      xend = diameter, yend = mod1))
```

<img src="01-mod-lineaire_files/figure-html/unnamed-chunk-9-1.png" width="672" style="display: block; margin: auto;" />


Nous pouvons ensuite facilement calculer cette distance comme ci-dessous : 


```r
# Calculer la distance entre y observé et y prédit
## création de notre fonction
distance <- function(observation, prediction) {
  diff <- observation - prediction
  attributes(diff) <- NULL
  diff
}
## Application de la fonction
dist1 <- distance(observation = trees$volume, prediction = mod1)
## affichage des résultats
dist1
```

```
#  [1] -0.263 -0.298 -0.331 -0.371 -0.328 -0.312 -0.453 -0.380 -0.270 -0.357
# [11] -0.250 -0.355 -0.344 -0.382 -0.484 -0.511 -0.183 -0.414 -0.512 -0.550
# [21] -0.303 -0.407 -0.312 -0.445 -0.364 -0.126 -0.143 -0.124 -0.327 -0.341
# [31]  0.065
```

- Defenir une règle pour obtenir une valeur unique de la distance de nos observations mesurées en y par rapport aux observations prédites en y par le modèle

Une première idée serait de sommer l'ensemble de nos distances comme ci-dessous :


```r
sum(dist1)
```

```
# [1] -10.175
```

Appliquons la suite d'étapes ci-dessus pour nos 4 modèles afin de les comparer






\BeginKnitrBlock{bdd}<div class="bdd">Essayez par vous même de trouver le meilleur modèle dans l'application shiny suivante :

Démarrez la SciViews Box et RStudio. Dans la fenêtre **Console** de RStudio, entrez l'instruction suivante suivie de la touche `Entrée` pour ouvrir le tutoriel concernant les bases de R :

    BioDataScience2::app("01a_modele_lineaire") TODO

N’oubliez pas d’appuyer sur la touche ESC pour reprendre la main dans R à la fin d’un tutoriel dans la console R</div>\EndKnitrBlock{bdd}


```r
summary(lm(data = trees, formula = volume ~ diameter))
```

```
# 
# Call:
# lm(formula = volume ~ diameter, data = trees)
# 
# Residuals:
#       Min        1Q    Median        3Q       Max 
# -0.231211 -0.087021  0.003533  0.100594  0.271725 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.04748    0.09553  -10.96 7.85e-12 ***
# diameter     5.65154    0.27649   20.44  < 2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1206 on 29 degrees of freedom
# Multiple R-squared:  0.9351,	Adjusted R-squared:  0.9329 
# F-statistic: 417.8 on 1 and 29 DF,  p-value: < 2.2e-16
```

