# Régression linéaire II {#lm2}




##### Objectifs {-}

- Appréhender les différentes formes de régressions linéaires par les moindres carrés.

- Choisir sa régression linéaire de manière judicieuse.

- Savoir utiliser les outils de diagnostic de la régression linéaire, en particulier l'analyse des résidus.


##### Prérequis {-}

- Le module précédent est une entrée en matière indispensable qui est complétée par le contenu du présent module.


### Résumé avec `summary()`(suite)

- Residual standard error\ :

Il s'agit de l'écart-type résiduel, considérant que les degrés de liberté du modèle est le nombre d'observations $n$ soustrait du nombre de paramètres à estimer (ici 2).

$$\sqrt{\frac{\sum(y_i - ŷ_i)^2}{n-2}}$$

- Multiple R-squared\ :

Il s'agit de la valeur du coefficient $R^2$ qui exprime la fraction de variance exprimé par le modèle. Souvenons nous que la variance totale respecte la propiété d'additivité. La variance conditionnelle $s^2_{y\left|x\right|}$ peut être décomposée comme une somme de carrés ($SC$) divisés par des degrés de liberté associés, avec\ :

$$SC(total) = SC(reg) + SC(résidus)$$

$$SC(total) = \sum_{i=0}^n(y_i - \bar y_i)^2$$

$$SC(reg) = \sum_{i=0}^n(ŷ_i - \bar y_i)^2$$

$$SC(résidus) = \sum_{i=0}^n(y_i - ŷ_i)^2$$

A partir de la décomposition de ces sommes de carrés, le coefficient $R^2$ se définit comme\ :

$$R^2 = \frac{SC(reg)}{SC(total)}$$

La valeur du $R^2$ est comprise entre 0 (lorsque le modèle est très mauvais et n'explique rien) et 1 (lorsque le modèle est parfait et "capture" toute la variance des données\ ; dans ce cas, tous les résidus valent zéro). Donc, **plus le coefficient $R^2$ se rapproche de un, plus le modèle explique bien les données**.

- Adjusted R-squared\:

La valeur du coefficient $R^2$ ajustée. Le calcul de cette valeur sera abordé dans la suite de ce livre.

- F-statistic\ :

Tout comme pour l'ANOVA, le test de la significativité de la régression car  $MS(reg)/MS(résidus)$ suit une distribution F à respectivement 1 et $n-2$ degré de liberté, avec $MS$ les carrés moyens, c'est-à-dire les sommes des carrés $SC$ divisés par leurs degrés de liberté respectifs.

- p-value\ : 

Il s'agit de la valeur p associé à la statistique de F, donc à l'ANOVA associée à la régression linéaire.

### Comparaison de régressions

Vous pouvez à présent comparer ces résultats avec un tableau et les six graphiques d'analyse des résidus sans la valeur supérieur à 0.5m de diamètre. **Attention, On ne peut supprimer une valeur sans raison valable.** La suppression de pointsd aberrants doit en principe être faite avant de débuter l'analyse. La raison de la suppression de ce point est lié au fait qu'il soit seul et unique point supérieur à 0.5m de diamètre. Nous le faisons ici à titre de comparaison.


```r
trees <- read("trees", package = "datasets", lang = "fr")
lm. <- lm(data = trees, volume ~ diameter)

trees_red <- filter(trees, diameter < 0.5)
lm1 <- lm(data = trees_red, volume ~ diameter)

chart(trees, volume ~ diameter) +
  geom_point() + 
  geom_abline(
    aes(intercept = lm.$coefficients[1], slope = lm.$coefficients[2]), 
    color = "red", size = 1.5) +
  labs( color = "Modèle")  +
  scale_color_viridis_c(direction = -1) +
  geom_abline(
    aes(intercept = lm1$coefficients[1], slope = lm1$coefficients[2]), 
    color = "blue", size = 1.5)
```

<img src="02-reg-lineaire-2_files/figure-html/unnamed-chunk-1-1.png" width="672" style="display: block; margin: auto;" />

Tentez d'analyser le tableau de notre régression.


```r
summary(lm1)
```

```
# 
# Call:
# lm(formula = volume ~ diameter, data = trees_red)
# 
# Residuals:
#       Min        1Q    Median        3Q       Max 
# -0.215129 -0.068502 -0.001149  0.070522  0.181398 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.94445    0.09309  -10.15 6.98e-11 ***
# diameter     5.31219    0.27540   19.29  < 2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1082 on 28 degrees of freedom
# Multiple R-squared:   0.93,	Adjusted R-squared:  0.9275 
# F-statistic: 372.1 on 1 and 28 DF,  p-value: < 2.2e-16
```

Tentez d'analyser les graphiques d'analyse des résidus ci-dessous.


```r
#plot(lm1, which = 1)
lm1 %>.%
  chart(broom::augment(.), .resid ~ .fitted) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
  labs(x = "Fitted values", y = "Residuals") +
  ggtitle("Residuals vs Fitted") 
```

<img src="02-reg-lineaire-2_files/figure-html/unnamed-chunk-3-1.png" width="672" style="display: block; margin: auto;" />

```r
#plot(lm1, which = 2)
lm1 %>.%
  chart(broom::augment(.), aes(sample = .std.resid)) +
  geom_qq() +
  geom_qq_line(colour = "darkgray") +
  labs(x = "Theoretical quantiles", y = "Standardized residuals") +
  ggtitle("Normal Q-Q") 
```

<img src="02-reg-lineaire-2_files/figure-html/unnamed-chunk-3-2.png" width="672" style="display: block; margin: auto;" />

```r
#plot(lm1, which = 3)
lm1 %>.%
  chart(broom::augment(.), sqrt(abs(.std.resid)) ~ .fitted) +
  geom_point() +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
  labs(x = "Fitted values",
    y = expression(bold(sqrt(abs("Standardized residuals"))))) +
  ggtitle("Scale-Location") 
```

<img src="02-reg-lineaire-2_files/figure-html/unnamed-chunk-3-3.png" width="672" style="display: block; margin: auto;" />

```r
#plot(lm1, which = 4)
lm1 %>.%
  chart(broom::augment(.), .cooksd ~ seq_along(.cooksd)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = seq(0, 0.1, by = 0.05), colour = "darkgray") +
  labs(x = "Obs. number", y = "Cook's distance") +
  ggtitle("Cook's distance") 
```

<img src="02-reg-lineaire-2_files/figure-html/unnamed-chunk-3-4.png" width="672" style="display: block; margin: auto;" />

```r
#plot(lm1, which = 5)
lm1 %>.%
  chart(broom::augment(.), .std.resid ~ .hat %size=% .cooksd) +
  geom_point() +
  geom_smooth(se = FALSE, size = 0.5, method = "loess", formula = y ~ x) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs(x = "Leverage", y = "Standardized residuals") +
  ggtitle("Residuals vs Leverage") -> a

#plot(lm1, which = 6)
lm1 %>.%
  chart(broom::augment(.), .cooksd ~ .hat %size=% .cooksd) +
  geom_point() +
  geom_vline(xintercept = 0, colour = NA) +
  geom_abline(slope = seq(0, 3, by = 0.5), colour = "darkgray") +
  geom_smooth(se = FALSE, size = 0.5, method = "loess", formula = y ~ x) +
  labs(x = expression("Leverage h"[ii]), y = "Cook's distance") +
  ggtitle(expression("Cook's dist vs Leverage h"[ii]/(1-h[ii]))) -> b
```

### Critère d'Akaike

A faire...


## Régression linéaire multiple 

TODO 

## Régression linéaire polynomiale

TODO
