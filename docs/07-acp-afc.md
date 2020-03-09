# ACP & AFC {#acp-afc}




##### Objectifs {-}

- Apprendre à réaliser une ordination de données quantitatives à l'aide de l'ACP.

- Savoir ordiner des variables qualitatives sous forme de tavbleaux cas par variables ou de tables de contingences à double entrée à l'aide de l'AFC.

##### Prérequis {-}

- Le module 6, et en particulier la partie sur le MDS doivent être assimilés avant d'attaquer le présent module.


## Analyse en composantes principales

Notre première approche d'ordination avec le MDS dans le précédent module nous a permis de comprendre l'intérêt de représenter des données multivariées sur des **cartes**. Malheureusement, les techniques itératives et basées sur les matrices de distances du MDS rendent cette technique peu propice pour analyser des gros jeux de données. En effet, le temps de calcul et le besoin en mémoire vive grandissent de manière exponentielle avec la taille des jeux de données. Heureusement, il existe aussi des techniques d'ordination qui se calculent plus facilement et plus rapidement sur de très gros jeux de données. L'**Analyse en Composantes Principales** ou ACP ("Principal Component Analysis" ou PCA en anglais) est une méthode de base qu'il est indispensable de connaitre et de comprendre. La plupart des autres techniques d'ordination plus sophistiquées sont des variante de l'ACP.

- Des **relations linéaires** sont suspectées entres les variables  (si elles ne sont pas linéaires, penser à transformer les données auparavant pour les linéairser).

- Ces relations conduisent à une répartition des individus (le nuage de points) qui forme une **structure que l’on cherchera à interpréter**.

- Pour **visualiser** cette structure, les données sont simplifiées (réduites) de **N variables à n (n < N et n = 2 ou 3 généralement)**. La représentation sous forme d'un nuage de points s'appelle une **carte**.

- La réduction des dimensions se fait avec une perte minimale d'information au sens de la variance des données.


### ACP dans SciViews::R

L'ACP est facilitée dans `SciViews::R`, mais au stade actuel, tout le code nécessaire (en particulier pour réaliser les graphiques avec `chart()`) n'est pas encore complètement intégré dans les packages. Ainsi, vous pouvez copier-coller le code du chunk suivant au début de vos scripts ou dans un chunk de `setup` dans vos documenbts R Markdown/Notebook.


```r
SciViews::R()
library(broom)

# broom implements only methods for prcomp objects, not princomp, while pcomp
# is compatible with princomp... but prcomp is simpler. So, conversion is easy
as.prcomp <- function(x, ...)
  UseMethod("as.prcomp")

as.prcomp.default <- function(x, ...)
  stop("No method to convert this object into a 'prcomp'")

as.prcomp.prcomp <- function(x, ...)
  x

as.prcomp.princomp <- function(x, ...)
  structure(list(sdev = as.numeric(x$sdev), rotation = unclass(x$loadings),
    center = x$center, scale = x$scale, x = as.matrix(x$scores)),
    class = "prcomp")

# Comparison of pcomp() -> as.prcomp() with prcomp() directly
# Almost the same, only no rownames for x (is it important?)
#iris_prcomp_pcomp <- as.prcomp(pcomp(iris[, -5], scale = TRUE))
#iris_prcomp <- prcomp(iris[, -5], scale = TRUE)

# Now, broom methods can be defined simply by converting into prcomp objects
augment.princomp <- function(x, data = NULL, newdata, ...)
  if (missing(newdata)) {
  augment(as.prcomp(x), data = data, ...)
  } else {
    augment(as.prcomp(x), data = data, newdata = newdata, ...)
  }

tidy.princomp <- function(x, matrix = "u", ...)
  tidy(as.prcomp(x), matrix = matrix, ...)

# There is no glance.prcomp() method

# There is a problem with pcomp() that returns a data.frame in scores,
# while it is a matrix in the original princomp object. pca() corrects this
pca <- function(x, ...) {
  res <- SciViews::pcomp(x, ...)
  # Change scores into a matrix
  res$scores <- as.matrix(res$scores)
  res
}

scale_axes <- function(data, aspect.ratio = 1) {
  range_x <- range(data[, 1])
  span_x <- max(range_x) - min(range_x)
  range_y <- range(data[, 2])
  span_y <- max(range_y) - min(range_y)
  if (span_y * aspect.ratio < span_x) {
    # Adjust range_x
    span_x_2 <- span_y / aspect.ratio / 2
    range_x_mid <- sum(range_x) / 2
    range_x <- c(range_x_mid - span_x_2, range_x_mid + span_x_2)
  } else {
    # Adjust range_y
    span_y_2 <- span_x * aspect.ratio / 2
    range_y_mid <- sum(range_y) / 2
    range_y <- c(range_y_mid - span_y_2, range_y_mid + span_y_2)
  }
  list(x = range_x, y = range_y)
}

autoplot.pcomp <- function(object,
type = c("screeplot", "altscreeplot", "loadings", "correlations", "scores", "biplot"),
choices = 1L:2L, name = deparse(substitute(object)), ar.length = 0.1,
circle.col = "gray", col = "black", fill = "gray", scale = 1, aspect.ratio = 1,
repel = FALSE, labels, title, xlab, ylab, ...) {
  type = match.arg(type)

  if (missing(title))
    title <- paste(name, type, sep = " - ")

  contribs <- paste0(names(object$sdev), " (",
    round((object$sdev^2/object$totdev^2) * 100, digits = 1), "%)")[choices]

  scores <- as.data.frame(object$scores[, choices])
  names(scores) <- c("x", "y")
  if (!missing(labels)) {
    if (length(labels) != nrow(scores))
      stop("You must provide a character vector of length ", nrow(scores),
        " for 'labels'")
    scores$labels <- labels
  } else {# Default labels are row numbers
    scores$labels <- 1:nrow(scores)
  }

  lims <- scale_axes(scores, aspect.ratio = aspect.ratio)

  if (!missing(col)) {
    if (length(col) != nrow(scores))
      stop("You must provide a vector of length ", nrow(scores), " for 'col'")
    scores$color <- col
    scores_formula <- y ~ x %col=% color %label=% labels
  } else {
    if (missing(labels)) {
      scores_formula <- y ~ x %label=% labels
    } else {
      scores_formula <- y ~ x %col=% labels %label=% labels
    }
  }

  res <- switch(type,
    screeplot = object %>.% # Classical screeplot
      tidy(., "pcs") %>.%
      chart(data = ., std.dev^2 ~ PC) +
      geom_col(col = col, fill = fill) +
      labs(y = "Variances", title = title),

    altscreeplot = object %>.% # screeplot represented by dots and lines
      tidy(., "pcs") %>.%
      chart(data = ., std.dev^2 ~ PC) +
      geom_line(col = col) +
      geom_point(col = "white", fill = col, size = 2, shape = 21, stroke = 3) +
      labs(y = "Variances", title = title),

    loadings = object %>.% # Plots of the variables
      tidy(., "variables") %>.%
      spread(., key = PC, value = value) %>.%
      #rename_if(., is.numeric, function(x) paste0("PC", x)) %>.%
      select(., c(1, choices + 1)) %>.%
      set_names(., c("labels", "x", "y")) %>.%
      chart(data = ., y ~ x %xend=% 0 %yend=% 0 %label=% labels) +
        annotate("path", col = circle.col,
          x = cos(seq(0, 2*pi, length.out = 100)),
          y = sin(seq(0, 2*pi, length.out = 100))) +
        geom_hline(yintercept = 0, col = circle.col) +
        geom_vline(xintercept = 0, col = circle.col) +
        geom_segment(arrow = arrow(length = unit(ar.length, "inches"),
          ends = "first")) +
        ggrepel::geom_text_repel(hjust = "outward", vjust = "outward") +
        coord_fixed(ratio = 1) +
        labs(x = contribs[1], y = contribs[2], title = title),

    correlations = object %>.% # Correlations plot
      Correlation(.) %>.%
      as_tibble(., rownames = "labels") %>.%
      select(., c(1, choices + 1)) %>.%
      set_names(., c("labels", "x", "y")) %>.%
      chart(data = ., y ~ x %xend=% 0 %yend=% 0 %label=% labels) +
      annotate("path", col = circle.col,
        x = cos(seq(0, 2*pi, length.out = 100)),
        y = sin(seq(0, 2*pi, length.out = 100))) +
      geom_hline(yintercept = 0, col = circle.col) +
      geom_vline(xintercept = 0, col = circle.col) +
      geom_segment(arrow = arrow(length = unit(ar.length, "inches"),
        ends = "first")) +
      ggrepel::geom_text_repel(hjust = "outward", vjust = "outward") +
      coord_fixed(ratio = 1) +
      labs(x = contribs[1], y = contribs[2], title = title),

    scores = scores %>.% # Plot of the individuals
      chart(data = ., scores_formula) +
      geom_hline(yintercept = 0, col = circle.col) +
      geom_vline(xintercept = 0, col = circle.col) +
      coord_fixed(ratio = 1, xlim = lims$x, ylim = lims$y, expand = TRUE) +
      labs(x = contribs[1], y = contribs[2], title = title) +
      theme(legend.position = "none"),

    biplot = object %>.% # Biplot using ggfortify function
      as.prcomp(.) %>.%
      ggfortify:::autoplot.prcomp(., x = choices[1], y = choices[2],
        scale = scale, size = -1, label = TRUE, loadings = TRUE,
        loadings.label = TRUE) +
      geom_hline(yintercept = 0, col = circle.col) +
      geom_vline(xintercept = 0, col = circle.col) +
      theme_sciviews() +
      labs(x = contribs[1], y = contribs[2], title = title),

    stop("Unrecognized type, must be 'screeplot', 'altscreeplot', loadings', 'correlations', 'scores' or 'biplot'")
  )

  if (type == "scores") {
    if (isTRUE(repel)) {
      res <- res + geom_point() + ggrepel::geom_text_repel()
    } else {# Use text
      res <- res + geom_text()
    }
  }

  if (!missing(xlab))
    res <- res + xlab(xlab)
  if (!missing(ylab))
    res <- res + ylab(ylab)
  res
}

chart.pcomp <- function(data, choices = 1L:2L, name = deparse(substitute(data)),
..., type = NULL, env = parent.frame())
  autoplot.pcomp(data, choices = choices, name = name, ..., type = type, env = env)
class(chart.pcomp) <- c("function", "subsettable_type")
```


### Indiens diabétiques

Les indiens Pimas sont des amérindiens originaires du nord du Mexique qui sont connus pour compter le plus haut pourcentage d'obèses et de diabétiques de toutes les éthnies. Ils ont fait l'objet de plusieurs études scientifiques d'autant plus que les Pimas en Arizona développent principalement cette obésité et ce diabète, alors que les Pimas mexicains les ont plus rarement. Il est supposé que leur mode de vie différent aux Etats_Units pourrait en être la raison. Voici un jeu de données qui permet d'explorer un peu ceci\ :


```r
pima <- read("PimaIndiansDiabetes2", package = "mlbench")
pima
```

```
# # A tibble: 768 x 9
#    pregnant glucose pressure triceps insulin  mass pedigree   age diabetes
#       <dbl>   <dbl>    <dbl>   <dbl>   <dbl> <dbl>    <dbl> <dbl> <fct>   
#  1        6     148       72      35      NA  33.6    0.627    50 pos     
#  2        1      85       66      29      NA  26.6    0.351    31 neg     
#  3        8     183       64      NA      NA  23.3    0.672    32 pos     
#  4        1      89       66      23      94  28.1    0.167    21 neg     
#  5        0     137       40      35     168  43.1    2.29     33 pos     
#  6        5     116       74      NA      NA  25.6    0.201    30 neg     
#  7        3      78       50      32      88  31      0.248    26 pos     
#  8       10     115       NA      NA      NA  35.3    0.134    29 neg     
#  9        2     197       70      45     543  30.5    0.158    53 pos     
# 10        8     125       96      NA      NA  NA      0.232    54 pos     
# # … with 758 more rows
```

Ce jeu de données contient des vlaeurs manquantes. Le graphique suivant permet de [visualiser l'importance des "dégâts"](https://cran.r-project.org/web/packages/naniar/vignettes/naniar-visualisation.html)\ :



```r
naniar::vis_miss(pima)
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-3-1.png" width="672" style="display: block; margin: auto;" />

Moins de 10% des données sont manquantes, et c'est principalement dans les variables `insulin` et `triceps`. Si nous souhaitons un tableau sans variables manquantes, nous pouvons décider d'éliminer des lignes et ou des colonnes (variables), mais ici nous souhaitons garder toutes les variables et réduisons donc uniquement le nombre de lignes avec la fonction `drop_na()`.



```r
pima <- drop_na(pima)
pima
```

```
# # A tibble: 392 x 9
#    pregnant glucose pressure triceps insulin  mass pedigree   age diabetes
#       <dbl>   <dbl>    <dbl>   <dbl>   <dbl> <dbl>    <dbl> <dbl> <fct>   
#  1        1      89       66      23      94  28.1    0.167    21 neg     
#  2        0     137       40      35     168  43.1    2.29     33 pos     
#  3        3      78       50      32      88  31      0.248    26 pos     
#  4        2     197       70      45     543  30.5    0.158    53 pos     
#  5        1     189       60      23     846  30.1    0.398    59 pos     
#  6        5     166       72      19     175  25.8    0.587    51 pos     
#  7        0     118       84      47     230  45.8    0.551    31 pos     
#  8        1     103       30      38      83  43.3    0.183    33 neg     
#  9        1     115       70      30      96  34.6    0.529    32 pos     
# 10        3     126       88      41     235  39.3    0.704    27 neg     
# # … with 382 more rows
```

Notre tableau est presque amputé de la moitié, mais il nous reste tout de même encore 392 cas, soit assez pour notre analyse. Avant de nous lancer dans une ACP, nous devons décrire les données, repérer les variables quantitatives d'intérêt, et synthétiser les corrélations linéaires (coefficients de corrélation de Pearson) entre ces variables.


```r
skimr::skim(pima)
```

```
# Skim summary statistics
#  n obs: 392 
#  n variables: 9 
# 
# ── Variable type:factor ───────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  variable missing complete   n n_unique                top_counts ordered
#  diabetes       0      392 392        2 neg: 262, pos: 130, NA: 0   FALSE
# 
# ── Variable type:numeric ──────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  variable missing complete   n   mean     sd     p0   p25    p50    p75
#       age       0      392 392  30.86  10.2  21     23     27     36   
#   glucose       0      392 392 122.63  30.86 56     99    119    143   
#   insulin       0      392 392 156.06 118.84 14     76.75 125.5  190   
#      mass       0      392 392  33.09   7.03 18.2   28.4   33.2   37.1 
#  pedigree       0      392 392   0.52   0.35  0.085  0.27   0.45   0.69
#  pregnant       0      392 392   3.3    3.21  0      1      2      5   
#  pressure       0      392 392  70.66  12.5  24     62     70     78   
#   triceps       0      392 392  29.15  10.52  7     21     29     37   
#    p100     hist
#   81    ▇▂▂▁▁▁▁▁
#  198    ▁▅▇▆▅▃▂▂
#  846    ▇▆▂▁▁▁▁▁
#   67.1  ▂▆▇▅▂▁▁▁
#    2.42 ▇▆▃▁▁▁▁▁
#   17    ▇▃▂▁▁▁▁▁
#  110    ▁▁▂▆▇▆▁▁
#   63    ▃▆▇▇▆▃▁▁
```

Nous avons une variable facteur `diabetes` à exclure de l'analyse, mais la variable `pregnant`, est une variable numérique discrète (nombre d'enfants portés). Nous l'éliminerons aussi de l'analyse.

La fonction `correlation()` du package `SciViews` nous permet d'inspecter les corrélations entre les variables choisies (donc toutes à l'exception de `pregnant` et `diabetes` qui ne sont pas quantitatives continues)\ :


```r
pima_cor <- correlation(pima[, 2:8])
knitr::kable(pima_cor, digits = 2)
```

            glucose   pressure   triceps   insulin   mass   pedigree    age
---------  --------  ---------  --------  --------  -----  ---------  -----
glucose        1.00       0.21      0.20      0.58   0.21       0.14   0.34
pressure       0.21       1.00      0.23      0.10   0.30      -0.02   0.30
triceps        0.20       0.23      1.00      0.18   0.66       0.16   0.17
insulin        0.58       0.10      0.18      1.00   0.23       0.14   0.22
mass           0.21       0.30      0.66      0.23   1.00       0.16   0.07
pedigree       0.14      -0.02      0.16      0.14   0.16       1.00   0.09
age            0.34       0.30      0.17      0.22   0.07       0.09   1.00

```r
plot(pima_cor)
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-6-1.png" width="672" style="display: block; margin: auto;" />

Quelques corrélations positives d'intensités moyennes se dégagent ici, notamment entre `mass` et `triceps` (épaisseur du pli cutané au niveau du triceps), ainsi qu'entre `glucose` (taux de glucose dans le sang), `insulin` (taux d'insuline dans le sang) et `age`. Par contre, la pression artérielle (`pressure`) et le `pedigree` (variable qui quantifie la susceptibilité au diabète en fonction de la parenté) semblent peu corrélés avec les autres variables.

L'ACP est en fait équivalente à une Analyse en Coordonnées Principales sur une matrice de distances euclidiennes (MDS métrique), mais en plus efficace en terme de calculs. Nous pouvons donc nous lancer dans l'analyse et en comprendre les résultats en gardant ceci à l'esprit.

Nous utiliserons la fonction `pca()` qui prend un argument `data =` et une formule du type `~ var1 + var2 + .... + varn`, ou plus simplement, directement un tableau contenant uniquement les variables à analyser comme argument unique. Comme les différentes variables sont mesurées dans des unités différentes, nous devons les standardiser (écart type ramené à un pour toutes). Ceci est réalisé par la fonction `pca()` en lui indiquant `scale = TRUE`. Donc\ :


```r
pima_pca <- pca(data = pima, ~ glucose + pressure + triceps + insulin + mass +
  pedigree + age, scale = TRUE)
```

Ou alors, nous sélectionnons les variables d'intérêt avec `select()` et appliquons `pca()` directement sur ce tableau, ce qui donnera le même résultat.


```r
pima %>.%
  select(., glucose:age) %>.%
  pca(., scale = TRUE) -> pima_pca
```

Le nuage de points dans l'espace initial à sept dimensions a été centré (origine ramenée au centre de gravité du nuage de points = moyenne des variables) par l'ACP. Ensuite une rotation des axes a été réalisée pour orienter son plus grand axe selon un premier **axe principal 1** ou **PC1** . Ensuite **PC2** est construit orthogonal au premier et dans la seconde direction de plus grande variabilité du nuage de points, et ainsi de suite pour les autres axes. Ainsi les axes PC1, PC2, PC3, ... représentent une **part de variance** de plus en plus faible par rapport à la variance totale du jeu de données. Ceci est présenté dans le résumé\ :


```r
summary(pima_pca)
```

```
# Importance of components (eigenvalues):
#                          PC1   PC2   PC3   PC4    PC5   PC6    PC7
# Variance               2.412 1.288 1.074 0.878 0.6389 0.399 0.3098
# Proportion of Variance 0.345 0.184 0.153 0.126 0.0913 0.057 0.0443
# Cumulative Proportion  0.345 0.529 0.682 0.807 0.8988 0.956 1.0000
# 
# Loadings (eigenvectors, rotation matrix):
#          PC1    PC2    PC3    PC4    PC5    PC6    PC7   
# glucose   0.441 -0.455        -0.198         0.736       
# pressure  0.329  0.101 -0.613  0.206  0.654        -0.171
# triceps   0.439  0.488               -0.367        -0.644
# insulin   0.402 -0.418  0.263 -0.388  0.123 -0.642 -0.129
# mass      0.446  0.506        -0.181                0.711
# pedigree  0.198         0.625  0.711  0.251              
# age       0.325 -0.337 -0.384  0.471 -0.592 -0.168  0.179
```

Le premier tableau `Importance of components (eigenvalues):` montre la part de variance présentée sur chacun des sept axes de l'ACP (PC1, PC2, ..., PC7). Le fait qu'il s'agit de *valeurs propres* (*eigenvalues* en anglais) apparaitra plus clair lorsque vous aurez lu les explications détaillées plus bas. Ces parts de variance s'additionnent pour donner la variance totale du nuage de points dans les sept dimensions (propriété d'additivité des variances). Pour facilité la lecture, la `Proportion de Variance` en %  est reprise également, ainsi que les proportions cumulées. Ainsi, les deux premiers axes de l'ACP capturent ici 53% de la variance totale. Et il faudrait considérer les cinq premiers axes pour capturer 90% de la variance totale. Cependant, les trois premiers axes cumulent tout de même plus des 2/3 de la variance. Nous pouvons restreindre notre analyse à ces trois axes-là.

Le second tableau `Loadings (eigenvectors, rotation matrix):` est la matrice de transformation des coordonnées initiales sur les lignes en coordonnées PC1 à PC7 en colonnes. Nous pouvons y lire l'**importante** des variables initiales sur les axes de l'ACP. Par exemple, l'axe PC3 contraste essentiellement `pressure` et `pedigree`.

Le **graphique des éboulis** sert à visualiser la "chute" de la variance d'un axe principal à l'autre, et aide à choisir le nombre d'axes à conserver (espace à dimensions réduites avec perte minimale d'information). Deux variantes en diagramme en barres versticales `chart$screeplot()` ou `chart$scree()` ou sous forme d'une ligne brisée `chart$altscree()` sont disponibles\ :


```r
chart$scree(pima_pca, fill = "cornsilk")
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-10-1.png" width="672" style="display: block; margin: auto;" />


```r
chart$altscree(pima_pca)
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-11-1.png" width="672" style="display: block; margin: auto;" />

La diminution est importante entre le premier et le second axe, mais plus progressive ensuite. Ceci traduit une structure plus complexe dans les données qui ne se réduit pas facilement à un très petit nombre d'axes. Nous pouvons visualiser le **premier plan principal** constitué par PC1 et PC2, tout en gardant à l'esprit que seulement 53% de la variance totale y est capturée. Donc, nous pouvons nous attendre à des déformations non négligeables des données dans ce plan, et d'autres aspects qui n'y sont pas (correctement) représentés. Nous verrons qu'il est porteur, toutefois, d'information utile.

Deux types de représentations peuvent être réalisées à partir d'ici\ : la représentation dans **l'espace des variables**, et la représentation complémentaire dans **l'espace des individus**. Ces deux représentations sont complémentaires et s'analysent conjointement. L'espace des variables représente les axes initiaux projettés comme des ombres dans le plan choisi de l'ACP (rappelez-vous l'analogie avec les ombres chinoises). Il se réalise à l'aide de `chart$loadings()`. Par exemple pour PC1 et PC2 nous indiquons `choices = c(1, 2)` (ou rien du tout, puisque ce sont les valeurs par défaut))\ :


```r
chart$loadings(pima_pca, choices = c(1, 2))
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-12-1.png" width="672" style="display: block; margin: auto;" />

Ce graphique s'interpète comme suit\ :

- Plus la norme (longueur) du vecteur qui représente une variable est grande et se rapporche de un (matérialisé par le cer cle gris), plus la variable est bien représentée dans le plan choisi. On évitera d'interpréter ici les variables qui ont des normes petites, comme `pedigree` ou `pressure`.

- Des vecteurs qui pointent dans la même direction représentent des variables **directement corrélés** entre elles. C'est le cas de `glucose`, `insulin` et `age`d'une part, et par ailleurs aussi de `mass` et `triceps`.

- Des vecteurs qui pointent en directions opposées représentent des variables **inversément proportionnelles**. Il n'y en a pas ici.

- Des vecteurs orthogonaux représentent des variables **non corrélées** entre elles. ainsi le groupoe `glucose`/`insulin`/`age` n'est pas corrélé avec le groupe `mass`/`triceps`

- Les PCs sont orientés en fonction des variables initiales, ou à défaut, les zones du graphique sont orientés. Ici, les gros sont dans le haut à droite du graphique, alors que ceux qui sont agés, et ont beaucoup de sucre et d'insuline dans le sang sont en bas à droite. A l'opposé, on trouve les plus maigres en bas à gauche et les jeunes ayant moins de glucose et d'insuline dans le sang en haut à gauche du graphique.

Cela donne déjà une vision synthétique des différentes corrélations entre la variables. Naturellement, on peut très bien choisir d'autres axes, pour peu qu'ils représentent une part de variance relativement importante. Par exemple, ici, nous pouvons représente le plan constitué par PC1 et PC3, puisque nous avons décidé de retenir les 3 premiers axes\ :


```r
chart$loadings(pima_pca, choices = c(1, 3))
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-13-1.png" width="672" style="display: block; margin: auto;" />

Nous voyons que `pedigree` et `pressure` (inversément proportionnels) sont bien mieux représentés le long de PC3. Ici l'axe PC3 est plus facile à orienter\ : en haut les pédigrées élevés et les pressions qartérielles basses, et en bas le contraire. Nous avons déjà lu cette informatioin dans le tableau des vecteurs propres de `summary()`.

Le graphice entre PC2 et PC3 complète l'analyse, mais n'apportant rien de plus, il peut être typiquement éliminé de votre rapport.


```r
chart$loadings(pima_pca, choices = c(2, 3))
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-14-1.png" width="672" style="display: block; margin: auto;" />

La seconde représentation se fait dans **l'espace des individus**. Ici, nous allons projeter les points relatifs à chaque individu dans le plan de l'ACP choisi. Cela se réalise à l'aide de `chart$scores()` (l'aspect ratio est le rapport hauteur/largeur peut s'adapter)\ :


```r
chart$scores(pima_pca, choices = c(1, 2), aspect.ratio = 3/5)
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-15-1.png" width="672" style="display: block; margin: auto;" />

Ce graphique est peu lisible tel quel. Généralement, nous représentons d'autres informations utiles sous forme de labels et ou de couleurs différentes. Nous pouvons ainsi contraster les individus qui ont le diabète de ceux qui ne l'ont pas sur ce graphique et aussi ajouter des ellipses de confiance à 95% autour des deux groupes pour aider à la cerner à l'aide de `stat_ellipse()`\ :


```r
chart$scores(pima_pca, choices = c(1, 2),
  labels = pima$diabetes) +
  stat_ellipse()
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-16-1.png" width="672" style="display: block; margin: auto;" />

Ce graphique est nettement plus intéressant. Il s'interprète comme suit\ :

- Nous savons que les individus plus âgés et ayant plus de glucose et d'insuline dans le sang sont dans le bas à droite du graphique. Or le groupe des diabétique, s'il ne se détache pas complètement tend à s'étaler plus dans cette région.

- A l'inverse, le groupe des non diabétiques s'étale vers la gauche, c'est-à-dire dans une région reprenant les individus les plus jeunes et les moins gros.

Le graphique entre PC1 et PC3 (analyse du troisième axe) donne ceci\ :


```r
chart$scores(pima_pca, choices = c(1, 3),
  labels = pima$diabetes) +
  stat_ellipse()
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-17-1.png" width="672" style="display: block; margin: auto;" />

Ici, la séparation se fait essentiellement sur l'axe horizontal (PC1). Donc, les différentes de pédigrée (élevé dans le haut du graphique) et de pression artérielle (élevée dans le bas du graphique) semblent être moins liés au diabète. Le graphique PC3 _versus_ PC2 peut aussi être réalisé, mais il n'apporte rien de plus (et en pratique, nous l'éliminerions d'un rapport).


```r
chart$scores(pima_pca, choices = c(2, 3),
  labels = pima$diabetes) +
  stat_ellipse()
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-18-1.png" width="672" style="display: block; margin: auto;" />

Etant donné que les deux graphiques (variables et individus) s'interprètent conjointement, nous pourrions être tentés de les superposer, cela s'appelle un **biplot**. Mais se pose alors un problème\ : celui de mettre à l'échelle les deux représentations pour qu'elles soient cohérentes entre elles. Ceci n'est pas facile et différentes représentations coexistent. L'argument `scale =` de la fonction `chart$biplot()` permet d'utiliser différentes mises à l'échelle. Enfin, ce type de graphique tend à être souvent bien trop encombré. Il est donc plus difficile à lire que les deux graphiques des variables et individus séparés. Voici ce que cela donne pour notre jeu de données exemple\ :


```r
chart$biplot(pima_pca)
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-19-1.png" width="672" style="display: block; margin: auto;" />

Bien moins lisible, en effet\ !


### Biométrie d'oursin

Analysons à présent un autre jeu de données qui nous montrera l'importance de la transformation (linéarisation), du choix de réduire ou non (argument `scale =`), et l'effet d'un **effet saturant**, et comment s'en débarrasser. Il s'agit de la biométrie effectuée sur deux populations de l'oursin violet *Paracentrotus lividus*, une en élevage et une oautre provenant du milieu naturel. Nous avons abondamment utilisé ce jeu de données en SDD I dans la section visualisation. Nous le connaissons bien, mais reprenons certains éléments essentiels ici...


```r
urchin <- read("urchin_bio", package = "data.io", lang = "FR")
urchin
```

```
# # A tibble: 421 x 19
#    origin diameter1 diameter2 height buoyant_weight weight solid_parts
#    <fct>      <dbl>     <dbl>  <dbl>          <dbl>  <dbl>       <dbl>
#  1 Pêche…       9.9      10.2    5               NA  0.522       0.478
#  2 Pêche…      10.5      10.6    5.7             NA  0.642       0.589
#  3 Pêche…      10.8      10.8    5.2             NA  0.734       0.677
#  4 Pêche…       9.6       9.3    4.6             NA  0.370       0.344
#  5 Pêche…      10.4      10.7    4.8             NA  0.610       0.559
#  6 Pêche…      10.5      11.1    5               NA  0.610       0.551
#  7 Pêche…      11        11      5.2             NA  0.672       0.605
#  8 Pêche…      11.1      11.2    5.7             NA  0.703       0.628
#  9 Pêche…       9.4       9.2    4.6             NA  0.413       0.375
# 10 Pêche…      10.1       9.5    4.7             NA  0.449       0.398
# # … with 411 more rows, and 12 more variables: integuments <dbl>,
# #   dry_integuments <dbl>, digestive_tract <dbl>,
# #   dry_digestive_tract <dbl>, gonads <dbl>, dry_gonads <dbl>,
# #   skeleton <dbl>, lantern <dbl>, test <dbl>, spines <dbl>,
# #   maturity <int>, sex <fct>
```

Ici aussi nous avons des valeurs manquantes\ :


```r
naniar::vis_miss(urchin)
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-21-1.png" width="672" style="display: block; margin: auto;" />

Ces valeurs manquantes sont rassemblées essentiellement dans les variables `buoyant_weight`, `dry_integuments`, les mesures relatives au squelette (`skeleton`, `lantern`, `test` et `spines`), et surtout au niveau de `sex` (impossible de déterminer le sexe des individus les plus jeunes). Si nous éliminons purement et simplement les lignes qui ont au moins une valeur manquante, nous perdons tous les individus jeunes, et c'est dommage. Nous allons donc d'abord éliminer les variables `sex`, ainsi que les quatres variables liées au squelette. Dans un second temps, nous appliquerons `drop_na()` sur ce qui reste\ :


```r
urchin %>.%
  select(., -(skeleton:spines), -sex) %>.%
  drop_na(.) -> urchin2
urchin2
```

```
# # A tibble: 319 x 14
#    origin diameter1 diameter2 height buoyant_weight weight solid_parts
#    <fct>      <dbl>     <dbl>  <dbl>          <dbl>  <dbl>       <dbl>
#  1 Pêche…      16.7      16.8    8.4          0.588   2.58        2.04
#  2 Pêche…      19.9      20      9.2          1.10    4.26        3.66
#  3 Pêche…      19.9      19.2    8.5          0.629   2.93        2.43
#  4 Pêche…      19.3      19.8   10.2          0.781   3.71        3.09
#  5 Pêche…      18.8      20      9.3          0.761   3.59        2.99
#  6 Pêche…      21.5      20.9    9.6          1.13    4.98        4.42
#  7 Pêche…      17.4      16.5    7.8          0.477   2.33        1.97
#  8 Pêche…      21        21.2   10.8          1.23    5.4         4.55
#  9 Pêche…      17.8      18.8    8.6          0.548   2.58        2.07
# 10 Pêche…      19.7      19.6    9.7          0.862   3.59        3.08
# # … with 309 more rows, and 7 more variables: integuments <dbl>,
# #   dry_integuments <dbl>, digestive_tract <dbl>,
# #   dry_digestive_tract <dbl>, gonads <dbl>, dry_gonads <dbl>,
# #   maturity <int>
```

Il nous reste 319 lignes des 421 initiales. Nous n'avons perdu qu'un quart des données, tout en nous privant seulement de quatres variables quantitatives liées au squelette (`sex`étant une variable qualitative, elle ne peut de toutes façons pas être introduite dans l'analyse, mais elle aurait pu servir pour colorer les individus).


```r
skimr::skim(urchin2)
```

```
# Skim summary statistics
#  n obs: 319 
#  n variables: 14 
# 
# ── Variable type:factor ───────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  variable missing complete   n n_unique                top_counts ordered
#    origin       0      319 319        2 Cul: 188, Pêc: 131, NA: 0   FALSE
# 
# ── Variable type:integer ──────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  variable missing complete   n mean   sd p0 p25 p50 p75 p100     hist
#  maturity       0      319 319 0.37 0.71  0   0   0   0    2 ▇▁▁▁▁▁▁▂
# 
# ── Variable type:numeric ──────────────────────────────────────────────────────────────────────────────────────────────────────────────
#             variable missing complete   n  mean    sd     p0    p25   p50
#       buoyant_weight       0      319 319  4.27  3.84  0.31   1.35   3.18
#            diameter1       0      319 319 32.78 11.71 14.6   23.25  31.5 
#            diameter2       0      319 319 32.71 11.67 15     23.45  31.6 
#      digestive_tract       0      319 319  1.9   2.03  0.034  0.45   1.21
#  dry_digestive_tract       0      319 319  0.23  0.21  0.015  0.075  0.17
#           dry_gonads       0      319 319  0.51  0.82  0      0.029  0.17
#      dry_integuments       0      319 319  7.16  6.3   0.58   2.22   5.42
#               gonads       0      319 319  1.72  2.65  0      0.1    0.63
#               height       0      319 319 16.78  6.25  7.3   11.1   16.2 
#          integuments       0      319 319 12.32 10.64  1.09   4      9.4 
#          solid_parts       0      319 319 16.52 15.27  1.46   4.96  11.73
#               weight       0      319 319 21.8  21.37  1.61   6.08  15.25
#    p75   p100     hist
#   5.67  17.73 ▇▃▂▂▁▁▁▁
#  39.65  65.6  ▆▇▆▆▃▃▂▁
#  39.6   65.6  ▇▇▆▆▃▃▂▁
#   2.54  10.37 ▇▃▂▁▁▁▁▁
#   0.31   1.02 ▇▅▂▁▁▁▁▁
#   0.64   5    ▇▂▁▁▁▁▁▁
#   9.42  28.8  ▇▃▂▂▁▁▁▁
#   2.2   15.93 ▇▂▁▁▁▁▁▁
#  21.5   32.2  ▇▆▆▆▅▃▂▁
#  16.07  47.22 ▇▃▃▁▁▁▁▁
#  21.69  73.14 ▇▅▂▂▁▁▁▁
#  28.14 100.51 ▇▅▁▁▁▁▁▁
```

Nous avons 12 variables quatitatives continues. Notez la distribution très asymétrique et similaire (voir colonne `hist`) de toutes ces variables. les variables `origin` et `maturity` ne pourront pas être utilisées, mais seront éventuellement utiles pour colorer les points dans nos graphiques. Qu'en est-il des corrélations entre les 12 variables\ ?


```r
urchin2_cor <- correlation(urchin2[, 2:13])
knitr::kable(urchin2_cor, digits = 2)
```

                       diameter1   diameter2   height   buoyant_weight   weight   solid_parts   integuments   dry_integuments   digestive_tract   dry_digestive_tract   gonads   dry_gonads
--------------------  ----------  ----------  -------  ---------------  -------  ------------  ------------  ----------------  ----------------  --------------------  -------  -----------
diameter1                   1.00        1.00     0.98             0.95     0.96          0.96          0.97              0.96              0.91                  0.93     0.80         0.79
diameter2                   1.00        1.00     0.97             0.95     0.96          0.96          0.97              0.96              0.91                  0.93     0.80         0.79
height                      0.98        0.97     1.00             0.93     0.92          0.93          0.94              0.93              0.88                  0.91     0.76         0.75
buoyant_weight              0.95        0.95     0.93             1.00     0.99          0.99          0.99              1.00              0.92                  0.94     0.88         0.87
weight                      0.96        0.96     0.92             0.99     1.00          0.99          0.99              0.99              0.95                  0.96     0.88         0.87
solid_parts                 0.96        0.96     0.93             0.99     0.99          1.00          0.99              0.99              0.95                  0.95     0.91         0.90
integuments                 0.97        0.97     0.94             0.99     0.99          0.99          1.00              1.00              0.93                  0.95     0.87         0.85
dry_integuments             0.96        0.96     0.93             1.00     0.99          0.99          1.00              1.00              0.92                  0.94     0.87         0.86
digestive_tract             0.91        0.91     0.88             0.92     0.95          0.95          0.93              0.92              1.00                  0.98     0.81         0.81
dry_digestive_tract         0.93        0.93     0.91             0.94     0.96          0.95          0.95              0.94              0.98                  1.00     0.82         0.82
gonads                      0.80        0.80     0.76             0.88     0.88          0.91          0.87              0.87              0.81                  0.82     1.00         0.99
dry_gonads                  0.79        0.79     0.75             0.87     0.87          0.90          0.85              0.86              0.81                  0.82     0.99         1.00

```r
plot(urchin2_cor)
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-24-1.png" width="672" style="display: block; margin: auto;" />

Toutes les corrélations sont positives, et certaines sont très élevées. Cela indique que plusieurs variables sont (pratiquement complètement)  redondantes, par exemple, `diameter1` et `diameter2`. Un effet princiapl semble dominer.

Si nous refaisons quelques graphiques, nous nous rappelons que les relations *ne sont pas* linéaires, par exemple, entre `diameter1` et `weight`\ :


```r
chart(data = urchin2, weight ~ diameter1) +
  geom_point()
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-25-1.png" width="672" style="display: block; margin: auto;" />

Ce type de relation, dite allométrique se linéarise très bien en effectuant une transformation double-log, comme nous pouvons le constater sur le graphique suivant\ :


```r
chart(data = urchin2, log(weight) ~ log(diameter1)) +
  geom_point()
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-26-1.png" width="672" style="display: block; margin: auto;" />

\BeginKnitrBlock{warning}<div class="warning">Il est crucial de bien nettoyer son jeu de données avant une ACP, et aussi, de vérifier que les relations sont linéaires, sinon, de transformer les données de manière appropriée. Rappelez-vous que l'ACP s'intéresse aux corrélations **linéaires** entre vos variables. </div>\EndKnitrBlock{warning}

Attention toutefois à la transformation logarithmique appliquée sur des données qui peuvent contenir des zéros (par exemple, `gonads` ou `dry_gonads`). Dans ce cas, la transformartion logarithme(x + 1) réalisée avec la fonction `log1p()` est plus indiquée. Nous allons ici trtansformer **toutes** les variables en log(x + 1). C'est assez fastidieux à faire avec `mutate()`, mais nous pouvons l'utiliser directemnt sur le tableau entier réduit aux variables quantitatives continues seules lors de l'appel à `pca()` comme suit\ :


```r
urchin2 %>.%
  select(., -origin, -maturity) %>.% # Elimine les variables non quantitatives
  log1p(.) %>.% # Transforme toutes les autres en log(x + 1)
  pca(., scale = TRUE) -> urchin2_pca # Effectue l'ACP après standardisation
```

Nous avons standardisé les données puisqu'elles sont mesurées dans des unités différentes (longueurs en mm, masses en g). Voici ce que donne notre ACP\ :


```r
summary(urchin2_pca)
```

```
# Importance of components (eigenvalues):
#                           PC1    PC2    PC3     PC4     PC5     PC6
# Variance               11.219 0.5010 0.1813 0.03862 0.02601 0.01657
# Proportion of Variance  0.935 0.0418 0.0151 0.00322 0.00217 0.00138
# Cumulative Proportion   0.935 0.9767 0.9918 0.99503 0.99720 0.99858
#                            PC7     PC8     PC9    PC10    PC11    PC12
# Variance               0.00931 0.00336 0.00210 0.00108 0.00082 0.00034
# Proportion of Variance 0.00078 0.00028 0.00017 0.00009 0.00007 0.00003
# Cumulative Proportion  0.99936 0.99964 0.99981 0.99990 0.99997 1.00000
# 
# Loadings (eigenvectors, rotation matrix):
#                     PC1    PC2    PC3    PC4    PC5    PC6    PC7   
# diameter1            0.295 -0.162        -0.441  0.237 -0.174  0.177
# diameter2            0.295 -0.166        -0.449  0.249 -0.178  0.157
# height               0.291 -0.218  0.154 -0.106 -0.902              
# buoyant_weight       0.296         0.120  0.509  0.124              
# weight               0.296 -0.149                                   
# solid_parts          0.297 -0.106  0.127                0.234       
# integuments          0.296 -0.159  0.157  0.153  0.125              
# dry_integuments      0.296 -0.115  0.160  0.455  0.114              
# digestive_tract      0.288        -0.571 -0.187         0.485 -0.519
# dry_digestive_tract  0.283        -0.702  0.217        -0.278  0.513
# gonads               0.271  0.568  0.226                0.575  0.430
# dry_gonads           0.259  0.697        -0.104        -0.465 -0.453
#                     PC8    PC9    PC10   PC11   PC12  
# diameter1            0.242 -0.706                     
# diameter2            0.124  0.688  0.263              
# height                                                
# buoyant_weight       0.530                0.265 -0.504
# weight              -0.145  0.116 -0.912              
# solid_parts         -0.594         0.216  0.638       
# integuments         -0.396         0.148 -0.702 -0.371
# dry_integuments      0.105         0.128 -0.133  0.774
# digestive_tract      0.201                            
# dry_digestive_tract -0.161                            
# gonads               0.143                            
# dry_gonads          -0.111
```

Whaaa\ ! Plus de 93% de la variance représentée sur le premier axe. Ça parait parfait\ ! Voici le graphique des éboulis\ :


```r
chart$scree(urchin2_pca)
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-30-1.png" width="672" style="display: block; margin: auto;" />

Ne vous réjousissez pas trop vite. Nous avons ici un **effet saturant** lié au fait que toutes les variables sont positivement corrélées entre elles. Cet effet est évident. Ici, c'est la taille. Nous allons conclure que plus un oursin est gros, plus ses dimensions et ses masses sont importante. **C'est trivial et d'un intérêt très limité**, avouons-le.

\BeginKnitrBlock{warning}<div class="warning">Puisque l'ACP optimise la variance sur le premier axe, un effet saturant aura tendance à occulter d'autres effets intéressants. Nous pouvons nous en débarrasser en identifiant une des variables représentant le mieux cet effet, et en calculant les ratios entre toutes les autres variables et celle-là. Ainsi, nous passons de quantification de la taille sur toutes les variables à des ratios qui quantifient beaucoup mieux des effets de forme plus subtils.</div>\EndKnitrBlock{warning}

Notez aussi les valeurs relativement faibles, mais homogènes de toutes les variables sur l'axe PC1 dans le tableau des vecteurs propres, avec des valeurs comprises entre 0,26 et 0,30. Le graphique des variables est également très moche dans le premier plan de l'ACP, même si un effet différent relatif aux gonades apparait tout de même sur l'axe PC2, il ne compte que pour 4,2% de la variance totale\ :


```r
chart$loadings(urchin2_pca)
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-32-1.png" width="672" style="display: block; margin: auto;" />

Recommençons tout de suite l'analyse en éliminant l'effet saturant. Nous pourrons considérer comme référence de la taille, par exemple, la masse immergée (`buoyant weight`) connue comme étant une mesure pouvant être mesurée très précisément. **Elle fait partie des variables les mieux corrélées sur l'axe PC1, représentant ainsi très bien cet effet saturant que nous voulons éliminer.** Voici notre calcul\ :


```r
urchin2 %>.%
  select(., -origin, -maturity, -buoyant_weight) %>.% # Elimination des variables inutiles
  (. / urchin2$buoyant_weight) %>.% # Division par buoyant_weight
  log1p(.) -> urchin3 # Transformation log(x + 1)
urchin3
```

```
#     diameter1 diameter2   height   weight solid_parts integuments
# 1    3.380877  3.386644 2.726760 1.683990    1.497119    1.388714
# 2    2.953357  2.958109 2.240741 1.587131    1.468302    1.345385
# 3    3.485925  3.451231 2.675524 1.733496    1.582091    1.478861
# 4    3.247200  3.271795 2.643585 1.749467    1.600897    1.441601
# 5    3.247291  3.306831 2.582396 1.744070    1.595668    1.424509
# 6    3.000850  2.973973 2.254397 1.690963    1.594759    1.406850
# 7    3.624570  3.572914 2.854510 1.773052    1.635777    1.457258
# 8    2.896813  2.905770 2.282655 1.686646    1.549377    1.412300
# 9    3.511709  3.564779 2.815702 1.742476    1.564468    1.416512
# 10   3.172056  3.167181 2.505869 1.641947    1.520279    1.364092
# 11   3.268178  3.231506 2.577097 1.666978    1.476470    1.358695
# 12   3.279373  3.274129 2.643147 1.634242    1.494883    1.332908
# 13   3.332233  3.296974 2.634369 1.601647    1.495145    1.342715
# 14   3.068462  3.068462 2.425726 1.568443    1.422907    1.308572
# 15   3.242208  3.237191 2.548099 1.706237    1.580013    1.395093
# 16   3.012068  2.976515 2.261722 1.726914    1.525542    1.371891
# 17   2.647038  2.708633 2.051586 1.542420    1.417490    1.241250
# 18   2.985775  2.942983 2.359462 1.458508    1.265194    1.172611
# 19   2.792739  2.788383 2.082470 1.453483    1.404191    1.193574
# 20   3.092519  3.064731 2.347520 1.571534    1.449396    1.176850
# 21   2.917219  2.848220 2.256780 1.622968    1.488087    1.332134
# 22   3.003502  3.053961 2.367701 1.682916    1.514403    1.331852
# 23   3.154718  3.125279 2.436211 1.641039    1.466250    1.289257
# 24   2.922469  2.909115 2.230064 1.606607    1.480274    1.261646
# 25   2.910184  2.841665 2.208924 1.699601    1.547925    1.394502
# 26   2.761460  2.746772 2.211831 1.779281    1.579558    1.401832
# 27   2.562093  2.572605 2.083389 1.696792    1.554015    1.328699
# 28   3.089122  3.080447 2.367701 1.710722    1.547230    1.413540
# 29   2.884093  2.861174 2.201479 1.710650    1.573456    1.367337
# 30   2.761672  2.789784 2.159982 1.661535    1.526802    1.296325
# 31   2.365524  2.375505 1.727927 1.605613    1.438734    1.350480
# 32   2.433933  2.400808 1.947694 1.732473    1.562977    1.365690
# 33   2.630112  2.654328 2.204454 1.787621    1.577568    1.428557
# 34   2.590803  2.576689 2.098248 1.749876    1.526922    1.388874
# 35   2.612355  2.595847 2.017088 1.762127    1.534193    1.373775
# 36   2.719421  2.741591 2.215152 1.758085    1.525963    1.387385
# 37   2.617050  2.627335 2.089286 1.688742    1.489696    1.367633
# 38   2.649809  2.663792 1.970902 1.713468    1.536116    1.386020
# 39   2.705668  2.705668 2.118018 1.697599    1.498384    1.371810
# 40   2.544158  2.526946 2.023071 1.709804    1.522775    1.349347
# 41   2.676719  2.700662 2.078831 1.765802    1.581017    1.413438
# 42   2.395575  2.378601 1.927728 1.624282    1.465228    1.339916
# 43   2.436368  2.411224 1.832774 1.662396    1.491829    1.328301
# 44   2.689601  2.653517 2.045423 1.656893    1.478666    1.337048
# 45   2.609943  2.624204 2.014164 1.735732    1.578534    1.402796
# 46   2.411547  2.411547 1.850846 1.790641    1.568265    1.363072
# 47   2.245242  2.250957 1.751380 1.693380    1.543330    1.381785
# 48   2.188728  2.176072 1.630776 1.711782    1.548368    1.343019
# 49   2.196106  2.201275 1.673864 1.730336    1.569440    1.375833
# 50   2.222542  2.194723 1.696950 1.670142    1.503509    1.322522
# 51   2.278937  2.293786 1.766558 1.661875    1.487862    1.322453
# 52   2.349861  2.364163 1.747817 1.723862    1.564136    1.354659
# 53   2.088174  2.103977 1.584710 1.655495    1.491320    1.346270
# 54   2.229197  2.234865 1.740604 1.704762    1.528822    1.359544
# 55   2.306144  2.292533 1.744202 1.753642    1.555270    1.372896
# 56   2.458580  2.497064 1.877649 1.765883    1.589515    1.405968
# 57   2.195454  2.192749 1.638954 1.685001    1.522647    1.326152
# 58   2.257040  2.204114 1.716055 1.692000    1.557292    1.365834
# 59   2.300742  2.303727 1.778921 1.674407    1.546940    1.333929
# 60   2.209733  2.206964 1.697410 1.635541    1.502351    1.367447
# 61   2.020478  1.995449 1.531881 1.687285    1.523051    1.299542
# 62   1.947060  1.929225 1.447338 1.614011    1.495333    1.303020
# 63   2.121998  2.087939 1.639600 1.733118    1.547107    1.362814
# 64   2.164833  2.152896 1.547878 1.776146    1.570922    1.358635
# 65   2.127971  2.137876 1.676486 1.741183    1.562775    1.355222
# 66   2.121115  2.148247 1.646028 1.843704    1.470841    1.323484
# 67   2.015875  2.018154 1.517737 1.748363    1.509444    1.327271
# 68   1.879188  1.901281 1.321214 1.986651    1.678803    1.447539
# 69   1.768580  1.768580 1.220753 1.921925    1.667299    1.324348
# 70   1.799822  1.802908 1.191056 1.948223    1.636655    1.340177
# 71   1.715082  1.718286 1.272118 1.882412    1.611882    1.323507
# 72   1.740742  1.764206 1.339079 1.908900    1.603726    1.370051
# 73   1.653429  1.639416 1.197053 1.878428    1.677563    1.314784
# 74   1.834261  1.832619 1.333591 1.974167    1.642153    1.390429
# 75   3.608259  3.625183 3.011387 1.812023    1.640364    1.397308
# 76   3.601205  3.578978 2.982006 1.780622    1.574021    1.386395
# 77   3.494505  3.483550 2.754371 1.795639    1.616731    1.398686
# 78   3.725105  3.731186 3.006964 1.776251    1.618905    1.413803
# 79   3.712397  3.718257 3.099160 1.847346    1.660117    1.435314
# 80   3.373863  3.379214 2.693431 1.804082    1.566558    1.395841
# 81   3.578590  3.528765 2.944011 1.825211    1.640823    1.446953
# 82   3.494505  3.505341 2.742887 1.765401    1.573046    1.407507
# 83   3.695022  3.683815 2.854824 1.835815    1.575869    1.420628
# 84   3.446289  3.431205 2.745050 1.864123    1.629901    1.463321
# 85   3.405083  3.421068 2.871374 1.788789    1.605288    1.372203
# 86   3.533808  3.533808 2.837846 1.810912    1.659611    1.447467
# 87   3.337027  3.353333 2.779407 1.740893    1.554161    1.360605
# 88   3.363323  3.378066 2.746374 1.813799    1.665819    1.452715
# 89   3.745575  3.788725 3.044522 1.802809    1.640937    1.402824
# 90   2.817112  2.812980 2.249722 1.646225    1.563070    1.337018
# 91   2.783566  2.792059 2.166338 1.651531    1.547046    1.277294
# 92   2.765419  2.742856 2.075721 1.625656    1.540787    1.293726
# 93   2.850230  2.845899 2.229514 1.603080    1.528140    1.299636
# 94   2.771332  2.782699 2.124435 1.662932    1.560929    1.325637
# 95   2.812820  2.850318 2.156110 1.743294    1.604998    1.371289
# 96   3.069433  3.065296 2.291252 1.769896    1.696082    1.411394
# 97   2.909527  2.883229 2.185118 1.761894    1.680389    1.411536
# 98   3.120679  3.107124 2.394840 1.731233    1.675697    1.381222
# 99   3.048257  3.052390 2.327629 1.765594    1.655168    1.387058
# 100  3.002424  3.027450 2.286892 1.746901    1.681036    1.388906
# 101  3.013782  3.022712 2.260308 1.711930    1.654236    1.373794
# 102  3.214530  3.238126 2.537550 1.745851    1.712101    1.423081
# 103  2.883797  2.872073 2.164649 1.756765    1.682098    1.365551
# 104  3.159759  3.150089 2.446288 1.747419    1.673643    1.341153
# 105  2.646773  2.646773 1.967792 1.704936    1.606276    1.404314
# 106  2.938253  2.919944 2.279329 1.861498    1.732059    1.499767
# 107  2.493174  2.476597 1.918157 1.699585    1.604301    1.409806
# 108  2.591035  2.610058 2.036399 1.824221    1.658185    1.434865
# 109  2.675686  2.668765 2.061955 1.739192    1.621703    1.438722
# 110  2.658884  2.649011 2.014887 1.756664    1.608128    1.406489
# 111  2.909720  2.913421 2.207993 1.807735    1.676745    1.469681
# 112  2.848487  2.876323 2.194840 1.867064    1.685907    1.441721
# 113  2.670386  2.690275 2.087078 1.825990    1.674013    1.442391
# 114  2.819240  2.791807 2.083368 1.863106    1.715439    1.501270
# 115  2.566539  2.576417 2.015678 1.740564    1.604774    1.367484
# 116  2.533305  2.499255 1.962031 1.689091    1.533533    1.349239
# 117  2.556562  2.556562 1.974667 1.666055    1.508520    1.330592
# 118  2.656010  2.663056 2.053594 1.723854    1.586421    1.386357
# 119  2.590238  2.613510 1.991080 1.662775    1.534089    1.359087
# 120  2.364657  2.369957 1.824848 1.843366    1.570301    1.411310
# 121  2.481513  2.473418 1.905194 1.906150    1.610979    1.423788
# 122  2.339788  2.364355 1.846068 1.831014    1.636824    1.421554
# 123  2.284389  2.292352 1.719180 1.736846    1.526875    1.348763
# 124  2.427869  2.414001 1.791031 1.751191    1.463017    1.369894
# 125  2.256112  2.277583 1.738235 1.780583    1.581817    1.347718
# 126  2.288071  2.309428 1.746377 1.756529    1.540593    1.379583
# 127  2.134555  2.144610 1.595243 1.712436    1.576879    1.326479
# 128  2.344268  2.362604 1.733018 1.781974    1.520515    1.367853
# 129  2.228090  2.222904 1.625517 1.696192    1.538186    1.336032
# 130  2.479983  2.477252 1.852046 1.831395    1.642490    1.394160
# 131  2.663501  2.628657 1.982905 1.811822    1.696399    1.399734
# 132  2.693759  2.682036 2.038379 1.893083    1.712979    1.458685
# 133  2.570443  2.590259 1.921361 1.842848    1.681328    1.463833
# 134  2.603141  2.618061 1.917863 1.796889    1.660013    1.393536
# 135  2.291693  2.305677 1.688436 1.814067    1.518119    1.354919
# 136  2.233436  2.197556 1.719642 1.788685    1.574943    1.356677
# 137  2.454722  2.419521 1.802157 1.857979    1.644730    1.397792
# 138  2.233686  2.244868 1.680183 1.845650    1.643936    1.390981
# 139  2.476227  2.466565 1.839097 1.976108    1.733719    1.438543
# 140  2.178356  2.157186 1.690004 1.804897    1.533922    1.398837
# 141  2.223935  2.195819 1.678356 1.799842    1.584808    1.425701
# 142  2.288900  2.272597 1.694410 1.835144    1.631791    1.419600
# 143  2.286359  2.297899 1.810931 1.886430    1.601495    1.408045
# 144  2.263796  2.275307 1.774952 1.885959    1.613031    1.434163
# 145  2.169351  2.216999 1.616707 1.796323    1.605621    1.334414
# 146  2.087019  2.047573 1.609056 1.739912    1.528471    1.328353
# 147  2.142707  2.161662 1.628772 1.797975    1.560079    1.352295
# 148  2.133212  2.135685 1.672825 1.759235    1.560392    1.374143
# 149  2.236930  2.220026 1.697922 1.776070    1.616268    1.357984
# 150  1.971509  1.967553 1.485425 1.808014    1.480921    1.327164
# 151  1.960594  1.939556 1.475776 1.869096    1.557476    1.343553
# 152  1.852069  1.852069 1.369891 1.774040    1.561205    1.323308
# 153  1.856445  1.831822 1.345000 1.770131    1.502575    1.314590
# 154  1.844963  1.854796 1.358166 1.766674    1.524845    1.308725
# 155  1.721121  1.714182 1.294223 1.772503    1.555078    1.270905
# 156  2.098525  2.118494 1.551143 1.938780    1.752857    1.502929
# 157  2.029484  2.016938 1.467115 1.749005    1.551682    1.354123
# 158  1.953712  1.928662 1.388112 1.742621    1.569044    1.338666
# 159  2.283980  2.262312 1.716611 1.936324    1.674138    1.429091
# 160  2.140247  2.121518 1.629497 1.895511    1.662280    1.398207
# 161  2.128318  2.091184 1.622581 1.914637    1.647128    1.379748
# 162  2.146134  2.118276 1.527804 1.878391    1.636405    1.396450
# 163  2.071366  2.027272 1.506847 1.891439    1.553563    1.359647
# 164  1.851762  1.853542 1.385741 1.844068    1.481550    1.342201
# 165  1.857206  1.843720 1.320638 1.872822    1.527364    1.361043
# 166  1.956337  1.961549 1.453231 1.978386    1.611663    1.392980
# 167  1.918802  1.899699 1.412358 1.960743    1.657897    1.377209
# 168  1.917694  1.914081 1.421167 1.957471    1.605720    1.410443
# 169  1.892677  1.879809 1.357183 1.862641    1.656429    1.422114
# 170  1.839930  1.829087 1.363109 1.838311    1.635203    1.381470
# 171  1.749717  1.764135 1.148574 1.804046    1.532695    1.329805
# 172  1.966131  1.971606 1.457146 1.944677    1.719620    1.435607
# 173  1.980796  1.973480 1.355800 1.903923    1.637336    1.403991
# 174  2.231408  2.231408 1.623966 1.988697    1.756396    1.477638
# 175  1.979366  2.010454 1.441405 1.936676    1.671913    1.387494
# 176  1.979411  1.993029 1.449509 1.875036    1.648588    1.387038
# 177  2.065546  2.050175 1.439636 1.866876    1.676605    1.389182
# 178  1.894631  1.875951 1.366130 1.908831    1.652756    1.504012
# 179  1.546227  1.546227 1.043149 1.907213    1.659081    1.315437
# 180  1.599104  1.599104 1.077791 1.953013    1.652625    1.341666
# 181  1.519627  1.530673 1.019266 1.885362    1.620834    1.298421
# 182  1.605434  1.615452 1.011128 1.923104    1.707864    1.340049
# 183  1.605266  1.617424 1.055895 1.917758    1.644316    1.339704
# 184  1.808600  1.802610 1.290044 2.005923    1.669457    1.357437
# 185  1.582432  1.566414 1.093581 1.909814    1.661826    1.308588
# 186  1.650581  1.653354 1.077943 1.864956    1.586954    1.309596
# 187  1.585834  1.588650 1.026114 1.823636    1.592859    1.293079
# 188  1.611490  1.598449 1.121694 1.857787    1.531656    1.273141
# 189  1.624517  1.628612 1.116891 1.926994    1.579469    1.309007
# 190  1.640385  1.636214 1.111602 1.901652    1.591190    1.279196
# 191  1.648208  1.646813 1.091657 1.875657    1.569297    1.273194
# 192  1.502393  1.499699 1.044821 1.813273    1.615337    1.295261
# 193  2.010474  2.008092 1.463778 1.650688    1.494539    1.284647
# 194  2.125830  2.120966 1.581600 1.698645    1.447584    1.304989
# 195  2.040036  2.015642 1.533145 1.664898    1.491064    1.301921
# 196  1.977892  1.946957 1.477297 1.667906    1.488008    1.285591
# 197  2.060489  2.041948 1.503412 1.662090    1.498993    1.279743
# 198  2.050367  2.052745 1.469803 1.651730    1.516351    1.304068
# 199  2.027841  2.049085 1.617799 1.652422    1.480445    1.296696
# 200  2.024645  1.988626 1.447151 1.576270    1.461337    1.260604
# 201  2.003981  2.017802 1.467488 1.685770    1.474208    1.274730
# 202  2.072665  2.070150 1.545974 1.658117    1.473207    1.270553
# 203  2.060864  2.065768 1.519765 1.674396    1.565978    1.359334
# 204  2.185066  2.187736 1.579090 1.640847    1.558795    1.316048
# 205  2.142823  2.162248 1.597540 1.547485    1.415092    1.304152
# 206  2.183501  2.178221 1.575566 1.653680    1.565847    1.328293
# 207  2.277763  2.268499 1.576447 1.557679    1.460253    1.307091
# 208  2.120913  2.109418 1.576122 1.564726    1.484878    1.293304
# 209  2.116269  2.121820 1.556557 1.587256    1.455132    1.300341
# 210  1.996458  2.015302 1.504435 1.581815    1.511522    1.295609
# 211  2.222294  2.154894 1.603395 1.634847    1.583347    1.308202
# 212  2.339041  2.344620 1.755646 1.696796    1.557832    1.378276
# 213  2.228925  2.211503 1.678753 1.630710    1.559094    1.296946
# 214  2.378728  2.420380 1.818615 1.559372    1.525843    1.310458
# 215  2.604358  2.611739 2.009774 1.636556    1.533970    1.366157
# 216  2.392452  2.412224 1.967603 1.658769    1.554381    1.348334
# 217  2.221159  2.180576 1.628010 1.471891    1.411938    1.264629
# 218  2.182163  2.188594 1.644002 1.525019    1.422878    1.260030
# 219  2.595617  2.609540 2.037617 1.648800    1.599730    1.371549
# 220  2.392439  2.382441 1.834227 1.605367    1.503452    1.322687
# 221  2.582659  2.599696 1.959165 1.708383    1.563241    1.372870
# 222  2.637729  2.613454 2.054389 1.715645    1.658930    1.382359
# 223  2.775667  2.762309 2.199620 1.551754    1.462124    1.339779
# 224  2.756138  2.756138 2.168390 1.627412    1.408895    1.321206
# 225  2.903047  2.931775 2.385747 1.693821    1.564417    1.347112
# 226  3.136325  3.085418 2.528036 1.659829    1.581806    1.404715
# 227  3.119920  3.115223 2.496554 1.665151    1.480519    1.362811
# 228  2.792434  2.792434 2.188298 1.619566    1.467136    1.359784
# 229  3.069552  3.127986 2.403948 1.702053    1.578423    1.373782
# 230  2.713331  2.733085 2.022655 1.541218    1.413596    1.337594
# 231  2.966280  2.962036 2.185387 1.667826    1.531767    1.365879
# 232  3.153491  3.130323 2.488894 1.685978    1.620277    1.357210
# 233  2.891468  2.886431 2.385807 1.532616    1.500833    1.353207
# 234  3.393412  3.404587 2.816707 1.676826    1.550324    1.397231
# 235  3.116710  3.111614 2.396329 1.580096    1.454350    1.331916
# 236  3.135277  3.130093 2.518795 1.531737    1.440302    1.342768
# 237  3.088249  3.125926 2.500726 1.525638    1.412131    1.263875
# 238  3.031841  3.020966 2.396059 1.491212    1.356302    1.259309
# 239  3.611291  3.637599 3.031755 1.603699    1.526357    1.305275
# 240  3.182663  3.177264 2.509255 1.580244    1.490804    1.308494
# 241  3.459072  3.452844 2.772527 1.562239    1.453041    1.319952
# 242  1.660618  1.646662 1.191907 1.765069    1.613966    1.351456
# 243  1.646236  1.663232 1.039195 1.802720    1.512727    1.305487
# 244  1.647636  1.614739 1.021421 1.737760    1.543117    1.307285
# 245  1.661071  1.696796 1.096347 1.726499    1.469612    1.324427
# 246  1.612910  1.624079 1.064714 1.740485    1.504730    1.290181
# 247  1.543843  1.523250 1.084042 1.601273    1.446889    1.243356
# 248  1.646901  1.648558 1.158701 1.724390    1.530093    1.341155
# 249  1.595139  1.580899 1.150738 1.618374    1.446185    1.255081
# 250  1.776262  1.749839 1.242675 1.668597    1.478394    1.299467
# 251  1.588256  1.586663 1.205109 1.748625    1.497703    1.320512
# 252  1.828092  1.836619 1.434404 1.629045    1.526805    1.282508
# 253  1.859635  1.861825 1.420740 1.660870    1.553074    1.310840
# 254  1.800213  1.808207 1.352074 1.680562    1.519968    1.334967
# 255  1.809594  1.811683 1.334109 1.674601    1.484646    1.296060
# 256  1.705956  1.728231 1.296490 1.654934    1.570692    1.293943
# 257  1.689149  1.689149 1.223595 1.957958    1.541728    1.393475
# 258  1.795574  1.753062 1.278259 1.903736    1.595180    1.398939
# 259  1.734404  1.714602 1.196084 1.851143    1.585200    1.349338
# 260  1.915467  1.920300 1.359708 1.987428    1.674477    1.419539
# 261  1.965932  1.984169 1.406766 1.993001    1.577711    1.421293
# 262  2.003606  1.986421 1.455141 1.895345    1.589519    1.411274
# 263  1.874984  1.864141 1.353696 1.868673    1.646802    1.417323
# 264  1.908640  1.917622 1.378592 1.908820    1.666677    1.467983
# 265  1.890713  1.896176 1.377060 1.895995    1.652243    1.387684
# 266  1.989355  1.991254 1.498533 1.958265    1.674133    1.448863
# 267  2.453650  2.451537 1.841039 2.104544    1.756126    1.503212
# 268  2.040859  2.046814 1.453870 1.877022    1.597400    1.394487
# 269  2.126471  2.118320 1.501499 1.887354    1.545264    1.397303
# 270  2.315544  2.321961 1.688811 2.000892    1.743513    1.495344
# 271  2.183255  2.185569 1.594790 1.824228    1.578381    1.390388
# 272  2.723641  2.782175 2.168970 1.748703    1.544410    1.396815
# 273  2.700369  2.721427 2.095241 1.699604    1.371361    1.306323
# 274  2.828062  2.810744 2.236750 1.814771    1.539993    1.423950
# 275  2.760587  2.750125 2.122227 1.803031    1.577924    1.414736
# 276  2.649545  2.639127 2.090258 1.718253    1.529705    1.355523
# 277  2.966035  2.969635 2.139118 1.852760    1.622958    1.488601
# 278  2.983600  2.994250 2.269461 1.851083    1.672035    1.440510
# 279  3.067719  3.094382 2.347389 1.883332    1.630585    1.466784
# 280  2.814232  2.810744 2.192360 1.893854    1.665055    1.430917
# 281  2.941162  2.970110 2.204429 1.836902    1.604811    1.409505
# 282  2.275285  2.259050 1.688271 1.835023    1.674941    1.435916
# 283  2.535313  2.537874 1.794615 1.922215    1.655144    1.587340
# 284  2.212830  2.224450 1.595719 1.810512    1.622228    1.415300
# 285  2.492121  2.479155 1.854379 1.860697    1.643483    1.423510
# 286  2.264815  2.249809 1.787087 1.824407    1.631023    1.424967
# 287  2.147163  2.131548 1.578276 1.804406    1.567663    1.372356
# 288  2.279071  2.257174 1.757021 1.831261    1.628658    1.429823
# 289  2.162368  2.171450 1.682892 1.881991    1.693154    1.458314
# 290  2.193368  2.186699 1.603972 1.910212    1.696612    1.428453
# 291  2.458662  2.449859 1.971476 1.784275    1.519986    1.396944
# 292  2.433209  2.435921 1.901020 1.830150    1.537032    1.398854
# 293  2.399421  2.394058 1.812581 1.778836    1.555124    1.388288
# 294  2.425652  2.447642 1.891043 1.846918    1.555453    1.415374
# 295  2.357342  2.385963 1.712796 1.803585    1.634373    1.445797
# 296  2.569766  2.561044 1.963739 1.884821    1.677551    1.480633
# 297  2.587196  2.581821 1.858765 1.887234    1.711225    1.506205
# 298  2.704325  2.695195 2.029453 1.930739    1.576747    1.443868
# 299  2.705862  2.717569 2.005910 1.934032    1.756788    1.547795
# 300  2.532709  2.553892 1.925875 1.900518    1.537175    1.439099
# 301  3.158297  3.138316 2.427828 1.844571    1.600761    1.483309
# 302  3.157032  3.161222 2.505556 1.769052    1.598882    1.469142
# 303  3.252162  3.220212 2.622593 1.800142    1.615180    1.471415
# 304  3.245438  3.249898 2.620175 1.796149    1.598917    1.480884
# 305  3.039890  3.055506 2.221848 1.804088    1.597258    1.473335
# 306  3.040826  3.075485 2.315615 1.711617    1.558970    1.403984
# 307  2.920208  2.888164 2.126275 1.718536    1.547371    1.409115
# 308  2.966115  2.958079 2.231578 1.775504    1.558819    1.433588
# 309  3.163236  3.184964 2.385924 1.779040    1.617632    1.432529
# 310  3.001217  3.036286 2.317598 1.724822    1.527506    1.392154
# 311  3.486070  3.491756 2.693293 1.722700    1.598955    1.431353
# 312  3.488407  3.476514 2.848036 1.737420    1.587053    1.405270
# 313  4.149472  4.139378 3.239246 1.853239    1.780822    1.585988
# 314  3.552257  3.525210 2.718179 1.780837    1.660973    1.491795
# 315  3.415511  3.444056 2.771376 1.664436    1.616603    1.445238
# 316  3.438927  3.438927 2.736776 1.656277    1.606216    1.456164
# 317  3.570644  3.564842 2.882483 1.716343    1.682523    1.470440
# 318  3.599158  3.593520 2.916436 1.821467    1.759728    1.512086
# 319  3.689321  3.666112 2.951349 1.741356    1.712790    1.454438
#     dry_integuments digestive_tract dry_digestive_tract      gonads
# 1         1.0304814      0.10391414          0.02667720 0.009140213
# 2         1.0226296      0.18061570          0.04368131 0.040178983
# 3         1.0511653      0.16838684          0.03608357 0.000000000
# 4         1.0497965      0.20619750          0.04764294 0.023167059
# 5         1.0487367      0.31540084          0.06613980 0.028901124
# 6         1.0340836      0.34644959          0.05538973 0.016565715
# 7         1.0534036      0.15247071          0.03504791 0.000000000
# 8         1.0300792      0.25590445          0.07661986 0.046812654
# 9         1.0707127      0.19799430          0.03640025 0.000000000
# 10        1.0454793      0.16409560          0.04103081 0.049568871
# 11        1.0428674      0.13288947          0.03765026 0.014352067
# 12        1.0366295      0.28102711          0.06480762 0.016121603
# 13        1.0042986      0.25179936          0.04272802 0.008692907
# 14        0.9986704      0.14874071          0.03620753 0.006701689
# 15        1.0179148      0.22508838          0.04776505 0.012346622
# 16        1.0121387      0.28353146          0.05674193 0.059942175
# 17        0.9259048      0.28794820          0.06092123 0.026159379
# 18        0.8796841      0.03067670          0.01754110 0.009265696
# 19        0.9151264      0.23821320          0.04032805 0.044680017
# 20        0.8707698      0.39788893          0.06033452 0.000000000
# 21        0.9971166      0.29353976          0.05719694 0.031399516
# 22        0.9887988      0.37056569          0.07065421 0.025836877
# 23        0.9934323      0.25287514          0.04000089 0.003056203
# 24        0.9875865      0.26543193          0.04969768 0.005572871
# 25        1.0065033      0.37389943          0.05211706 0.004313838
# 26        1.0107301      0.38910562          0.07187256 0.106219531
# 27        0.9977788      0.34732587          0.05386293 0.165739618
# 28        1.0166316      0.17566121          0.02565203 0.001894837
# 29        1.0040314      0.30954680          0.05006314 0.080001138
# 30        0.9915289      0.42164245          0.06977041 0.008984210
# 31        1.0018815      0.16987678          0.03154393 0.120611802
# 32        0.9991331      0.46231071          0.08106852 0.156054914
# 33        1.0165620      0.41221197          0.07486424 0.083840157
# 34        1.0092567      0.33120970          0.07212606 0.093638739
# 35        1.0033095      0.33129419          0.06897799 0.125661857
# 36        1.0018291      0.37396689          0.08160689 0.016862858
# 37        1.0015843      0.30303159          0.05508759 0.090195849
# 38        1.0074349      0.29173383          0.04385813 0.112933640
# 39        1.0079571      0.28704461          0.05393305 0.021923491
# 40        1.0009137      0.34865171          0.07127015 0.114818967
# 41        1.0047286      0.45342977          0.09131743 0.010014605
# 42        1.0039041      0.32779008          0.05393305 0.036276603
# 43        0.9970188      0.27712412          0.05521690 0.181232877
# 44        0.9976165      0.31139279          0.06074434 0.030833335
# 45        1.0162379      0.38109717          0.06612730 0.070687481
# 46        0.9902349      0.40542446          0.07627435 0.030028680
# 47        1.0010039      0.32906135          0.05013581 0.114907545
# 48        0.9888082      0.36617645          0.06936842 0.193301336
# 49        0.9925674      0.36644098          0.07189423 0.196090308
# 50        0.9911278      0.38597849          0.06376897 0.073228431
# 51        0.9858936      0.30705661          0.05955491 0.159056100
# 52        0.9900742      0.42630874          0.08396619 0.143269973
# 53        0.9912859      0.26889393          0.05014105 0.165487203
# 54        1.0022857      0.31653605          0.05896725 0.093228914
# 55        0.9977459      0.34146792          0.05284687 0.141456976
# 56        1.0011469      0.48728354          0.08126572 0.007030153
# 57        0.9882789      0.34539748          0.05890013 0.125285902
# 58        0.9963434      0.35971633          0.06279527 0.142322835
# 59        0.9902401      0.38422439          0.06637753 0.133855841
# 60        0.9953996      0.22672512          0.04194211 0.147331618
# 61        0.9735576      0.42002428          0.06400452 0.250545823
# 62        0.9765619      0.34163519          0.05277872 0.228764146
# 63        0.9981126      0.37791415          0.05628948 0.102395038
# 64        0.9916398      0.43695400          0.06407444 0.209992670
# 65        0.9913273      0.44873499          0.07253642 0.165277805
# 66        0.9866307      0.30843048          0.04850521 0.184507913
# 67        0.9902460      0.44685958          0.06464982 0.117790646
# 68        1.0055519      0.37251320          0.05092387 0.470498507
# 69        0.9727230      0.43043224          0.05300801 0.652279051
# 70        0.9693045      0.32109945          0.04836224 0.600335257
# 71        0.9653713      0.40685536          0.06390737 0.517467095
# 72        0.9825536      0.46327839          0.06040483 0.331784187
# 73        0.9630466      0.51890015          0.06792006 0.588536752
# 74        0.9846172      0.50782510          0.07041042 0.319212817
# 75        1.0012214      0.23404341          0.05454226 0.284738420
# 76        0.9900995      0.20485408          0.04988271 0.280507788
# 77        1.0014705      0.35102572          0.07747255 0.214422461
# 78        1.0099668      0.29681053          0.06843618 0.231958619
# 79        1.0228058      0.33337306          0.07709474 0.267810851
# 80        1.0130948      0.21594464          0.05795858 0.087807197
# 81        1.0223939      0.30736520          0.07399742 0.104257179
# 82        1.0210483      0.17770621          0.04986601 0.122119194
# 83        1.0194600      0.19719875          0.05497567 0.053063424
# 84        1.0411314      0.23395413          0.06473427 0.151574369
# 85        0.9986650      0.23721426          0.05965359 0.337650079
# 86        1.0130782      0.24822380          0.06370454 0.162701628
# 87        1.0031274      0.21129441          0.05445665 0.210545253
# 88        1.0148327      0.27222648          0.06384115 0.132529890
# 89        1.0248461      0.29007920          0.06078147 0.124457011
# 90        0.9233299      0.35997386          0.06709630 0.212825603
# 91        0.9150759      0.35054416          0.05567798 0.303084369
# 92        0.9235066      0.32122619          0.04549650 0.337079632
# 93        0.9051233      0.26664933          0.05113261 0.311917070
# 94        0.9251902      0.34147042          0.05173772 0.279105022
# 95        0.9930066      0.33602662          0.04917802 0.294895681
# 96        0.9857054      0.43879376          0.07083973 0.300157257
# 97        0.9859610      0.41433780          0.06058060 0.341145777
# 98        0.9736209      0.41319811          0.07614712 0.239984569
# 99        0.9745349      0.25403636          0.04062992 0.359891067
# 100       0.9935572      0.31376549          0.06270846 0.394720744
# 101       0.9860676      0.34036621          0.05862387 0.311420720
# 102       0.9975524      0.38581937          0.07265508 0.266752577
# 103       0.9965963      0.34119518          0.06508547 0.384118981
# 104       0.9851430      0.32984119          0.06348200 0.291106297
# 105       0.9970242      0.27009743          0.03689775 0.372282634
# 106       1.0030623      0.41467233          0.06627135 0.414672334
# 107       0.9927919      0.23943282          0.03901259 0.411599448
# 108       1.0042314      0.49249784          0.07009422 0.287543251
# 109       1.0144697      0.34407017          0.04886886 0.227471329
# 110       0.9853491      0.30063269          0.05459623 0.393153641
# 111       1.0111554      0.48613573          0.07849216 0.246068733
# 112       0.9974769      0.48693572          0.08182610 0.293466516
# 113       0.9923701      0.39450698          0.06551744 0.432879768
# 114       1.0061993      0.42724676          0.06050787 0.330215106
# 115       0.9923418      0.41105051          0.06666358 0.318336220
# 116       0.9882593      0.31557094          0.05735282 0.281157057
# 117       0.9907433      0.27674353          0.04630457 0.240145969
# 118       0.9993763      0.34556386          0.04911372 0.260057593
# 119       0.9931421      0.29329460          0.04815356 0.309879828
# 120       0.9991774      0.45893102          0.06033855 0.000000000
# 121       1.0027824      0.54413886          0.07430308 0.047082287
# 122       0.9999167      0.55388346          0.07192370 0.156340409
# 123       0.9907576      0.42780807          0.06091708 0.065829888
# 124       0.9857348      0.24520391          0.03074517 0.027712839
# 125       0.9822041      0.44594003          0.06273222 0.278478905
# 126       0.9910489      0.34515735          0.04929406 0.190699041
# 127       0.9764819      0.34210048          0.04998438 0.438891967
# 128       0.9903609      0.37103502          0.05103971 0.024501688
# 129       0.9832869      0.32157402          0.04230352 0.309315962
# 130       0.9810327      0.44264649          0.06918823 0.382296353
# 131       0.9714201      0.42938589          0.06728285 0.530671598
# 132       0.9971219      0.42804328          0.06262480 0.343097621
# 133       1.0044438      0.44984916          0.06860162 0.281524573
# 134       0.9740048      0.43968422          0.07807929 0.413177159
# 135       0.9598025      0.36056300          0.05201969 0.195490718
# 136       0.9565836      0.34204467          0.05751404 0.415587359
# 137       0.9771549      0.43973479          0.06406378 0.398112558
# 138       0.9777411      0.30518827          0.03911171 0.566474905
# 139       0.9767814      0.44503482          0.07172323 0.605537641
# 140       0.9993118      0.30879458          0.05206590 0.162782010
# 141       0.9940361      0.33296361          0.04387218 0.198828452
# 142       0.9855889      0.28377219          0.05108249 0.406338222
# 143       0.9946795      0.39488008          0.05552911 0.162006042
# 144       1.0053328      0.42241130          0.05207476 0.210150242
# 145       0.9692429      0.34514211          0.05807990 0.482027142
# 146       0.9809897      0.33502496          0.04833069 0.267180268
# 147       0.9872683      0.39416880          0.05358106 0.268369509
# 148       0.9938485      0.31786253          0.04695986 0.223824800
# 149       0.9781206      0.29997917          0.04387342 0.499613268
# 150       0.9900641      0.33431083          0.05112811 0.169161078
# 151       0.9814760      0.41712446          0.07133564 0.298768768
# 152       0.9745701      0.37756359          0.05507381 0.359983386
# 153       0.9791513      0.29720499          0.04678997 0.288268699
# 154       0.9730055      0.34463254          0.05238266 0.333046681
# 155       0.9571873      0.37315798          0.05819523 0.510369340
# 156       1.1266778      0.34090758          0.04659178 0.505522461
# 157       0.9832380      0.23310077          0.03572585 0.397866786
# 158       0.9674320      0.35646280          0.05230908 0.355522386
# 159       0.9948145      0.41155508          0.07691584 0.346764651
# 160       0.9782517      0.33994655          0.05622808 0.484867203
# 161       0.9740221      0.38264729          0.06104686 0.466095735
# 162       0.9813640      0.27037821          0.04431385 0.440699351
# 163       0.9786968      0.28548809          0.04368206 0.313071783
# 164       0.9689423      0.36268083          0.05307165 0.068034787
# 165       0.9908809      0.36967662          0.04508905 0.138035462
# 166       0.9923580      0.47381224          0.07016333 0.214298776
# 167       0.9776821      0.48888719          0.07667540 0.455265416
# 168       0.9954649      0.49210413          0.06302009 0.156465889
# 169       0.9849954      0.38798155          0.05541323 0.429070313
# 170       0.9753897      0.44479086          0.06675197 0.422048765
# 171       0.9803534      0.40322846          0.05508199 0.253638703
# 172       0.9905033      0.52942045          0.07549805 0.405857137
# 173       0.9923384      0.41816300          0.05524989 0.375527008
# 174       0.9951672      0.56374390          0.07316594 0.461786595
# 175       0.9804847      0.52215657          0.07882060 0.411804125
# 176       0.9874534      0.49349209          0.06318180 0.357397779
# 177       0.9740654      0.44442427          0.06971767 0.471072976
# 178       0.9822338      0.44143709          0.05337272 0.413714071
# 179       0.9748747      0.41143285          0.05268460 0.611007466
# 180       0.9688130      0.43718935          0.05958488 0.562016286
# 181       0.9648915      0.31524111          0.03872487 0.641092110
# 182       0.9646895      0.50075421          0.05902064 0.609017659
# 183       0.9783349      0.46507930          0.05015269 0.461673139
# 184       0.9783036      0.46625103          0.06267881 0.486057801
# 185       0.9545969      0.37985995          0.04412443 0.682886428
# 186       0.9572149      0.35745154          0.04456435 0.500524050
# 187       0.9674208      0.41724975          0.05231353 0.466809737
# 188       0.9686561      0.39423702          0.04650559 0.328415202
# 189       0.9655454      0.46255923          0.06192470 0.326473051
# 190       0.9292616      0.55256152          0.06106268 0.376036378
# 191       0.9580334      0.41476437          0.05157092 0.379678863
# 192       0.9668308      0.43542256          0.05295324 0.530448033
# 193       0.9800850      0.12251144          0.02612717 0.448961938
# 194       0.9870613      0.17523858          0.03132104 0.169163068
# 195       0.9841368      0.20170471          0.03559354 0.324223130
# 196       0.9738363      0.19992558          0.04008889 0.365394621
# 197       0.9756787      0.26222837          0.04447527 0.379029383
# 198       0.9808755      0.16541912          0.02764094 0.430589579
# 199       0.9838720      0.23242768          0.04072121 0.292004702
# 200       0.9760957      0.21682529          0.04120561 0.330673086
# 201       0.9830530      0.29873276          0.05275403 0.264256660
# 202       0.9800680      0.20217923          0.04446478 0.324779048
# 203       0.9798153      0.33305243          0.04827159 0.338571359
# 204       0.9689939      0.23908353          0.04791600 0.525323556
# 205       0.9881319      0.20009949          0.03989200 0.129483828
# 206       0.9834409      0.33881699          0.05520210 0.409702551
# 207       0.9876829      0.18318068          0.03334773 0.258484417
# 208       0.9786079      0.22470063          0.03025094 0.367462648
# 209       0.9804393      0.26143799          0.04468175 0.205660776
# 210       0.9778610      0.33664659          0.04715283 0.311116242
# 211       0.9765646      0.20231853          0.04103466 0.548371850
# 212       1.0018765      0.29096244          0.04454384 0.288486866
# 213       0.9765582      0.32246905          0.05400388 0.504726032
# 214       0.9868826      0.17389738          0.02798388 0.428826594
# 215       1.0091024      0.29897600          0.05078224 0.072958816
# 216       0.9947503      0.41332327          0.06873068 0.216274690
# 217       0.9767540      0.20910851          0.03665238 0.236128263
# 218       0.9832753      0.24613047          0.04241898 0.182402592
# 219       0.9901248      0.31883811          0.05738550 0.403364102
# 220       0.9897455      0.29628814          0.04995556 0.235911323
# 221       1.0035613      0.34593125          0.04989174 0.042597508
# 222       0.9987507      0.28229058          0.05025808 0.451889315
# 223       1.0024665      0.31557784          0.05328605 0.026202984
# 224       0.9973975      0.12863245          0.02929440 0.068543108
# 225       0.9876843      0.35399203          0.06742808 0.210017938
# 226       1.0210362      0.24142168          0.03690907 0.227868557
# 227       1.0271533      0.12500710          0.03051120 0.031231323
# 228       1.0102544      0.20928531          0.04101933 0.093396213
# 229       1.0239291      0.46558343          0.09083310 0.065391591
# 230       1.0075450      0.15337679          0.02969899 0.040041907
# 231       1.0091534      0.30609550          0.05456481 0.102824617
# 232       1.0110059      0.34416540          0.05397125 0.099756513
# 233       1.0039839      0.38454652          0.05315197 0.069837466
# 234       1.0258971      0.28755663          0.05830993 0.045926855
# 235       0.9968364      0.13963721          0.03362230 0.186424812
# 236       1.0150859      0.19832253          0.03791205 0.045097706
# 237       0.9556966      0.24008704          0.04240918 0.016706833
# 238       0.9611539      0.12436470          0.02808366 0.064692692
# 239       1.0030883      0.08076348          0.05564985 0.164318001
# 240       0.9741529      0.36780468          0.06411646 0.011360824
# 241       0.9848914      0.25869889          0.05845135 0.024177430
# 242       0.9775825      0.30562985          0.03695846 0.566050240
# 243       0.9772660      0.39545048          0.05863307 0.281614992
# 244       0.9802287      0.33599181          0.04752379 0.413660699
# 245       0.9834913      0.27400643          0.03156209 0.213336721
# 246       0.9771493      0.22473636          0.03399413 0.416747597
# 247       0.9695914      0.22673217          0.03773503 0.398692202
# 248       0.9821998      0.30701762          0.04906673 0.305490797
# 249       0.9785016      0.21121948          0.03426215 0.359674266
# 250       0.9805748      0.32722724          0.05201434 0.225532939
# 251       0.9815934      0.19419222          0.02044963 0.351674264
# 252       0.9732490      0.18301625          0.03429401 0.558776185
# 253       0.9773468      0.19009478          0.03237865 0.547470435
# 254       0.9877789      0.29046774          0.04243262 0.320241373
# 255       0.9802784      0.23547699          0.03959686 0.325220899
# 256       0.9794976      0.20284929          0.03284117 0.639722329
# 257       0.9963443      0.32342047          0.03867240 0.203344251
# 258       0.9935949      0.45852029          0.04990987 0.227956196
# 259       0.9755388      0.27648554          0.03005777 0.491978380
# 260       0.9825412      0.36309905          0.05002700 0.515274535
# 261       0.9995320      0.32609329          0.04109534 0.151087157
# 262       0.9953673      0.37776061          0.04804717 0.176652549
# 263       0.9965651      0.50010198          0.05330381 0.263945385
# 264       1.0017977      0.51738522          0.06436433 0.169124465
# 265       0.9849029      0.53772473          0.06469725 0.313931051
# 266       0.9983971      0.46211010          0.05173801 0.283911420
# 267       1.0001144      0.69073825          0.08975735 0.176798653
# 268       0.9948018      0.55117978          0.06176076 0.045873352
# 269       0.9938680      0.39158575          0.04614016 0.015332958
# 270       0.9996723      0.64737637          0.08370156 0.176215999
# 271       0.9979789      0.47156510          0.05584045 0.067977977
# 272       0.9992778      0.21598869          0.04489840 0.243361413
# 273       0.9873742      0.14702815          0.02089831 0.000000000
# 274       1.0002757      0.28595284          0.04541333 0.028624085
# 275       1.0034703      0.35998233          0.06374290 0.133311974
# 276       0.9881339      0.23035114          0.04304825 0.275882526
# 277       1.0269251      0.22479462          0.04115243 0.190665494
# 278       1.0034901      0.41062338          0.06153684 0.106910794
# 279       0.9996810      0.26771723          0.05645295 0.202132860
# 280       1.0024093      0.32865970          0.05094769 0.208883581
# 281       1.0047449      0.26690525          0.03417267 0.172036863
# 282       0.9937574      0.46576925          0.07335615 0.389897958
# 283       1.0131945      0.47733250          0.06879001 0.053566232
# 284       0.9999180      0.42662981          0.05811125 0.263364814
# 285       0.9934017      0.32429716          0.04565613 0.366190412
# 286       1.0027872      0.49216037          0.06702443 0.055784318
# 287       0.9853902      0.32271328          0.04809152 0.303315102
# 288       1.0331608      0.41991026          0.05265217 0.261226306
# 289       1.0408893      0.41438563          0.05963994 0.297820026
# 290       0.9955805      0.46619995          0.05600676 0.355623600
# 291       0.9935123      0.32388444          0.04995711 0.106786996
# 292       0.9984246      0.37273854          0.06003844 0.027464194
# 293       0.9997438      0.39616134          0.06280065 0.167765432
# 294       0.9992234      0.38412744          0.05800497 0.072726268
# 295       1.0037741      0.36692666          0.05161665 0.263213172
# 296       1.0049960      0.42734536          0.05517119 0.261389343
# 297       1.0101033      0.51059756          0.06883609 0.178954771
# 298       1.0040986      0.37032636          0.05290373 0.004517128
# 299       1.0240523      0.43874374          0.05985715 0.282203575
# 300       1.0048134      0.22793550          0.02659112 0.033129873
# 301       1.0117109      0.31852990          0.04362971 0.114415384
# 302       0.9995879      0.38031514          0.07158827 0.066337947
# 303       1.0127507      0.26833911          0.05401222 0.100017111
# 304       1.0302258      0.30393102          0.05655256 0.053728930
# 305       1.0152354      0.29315682          0.05923519 0.106100434
# 306       0.9989947      0.24903673          0.04781464 0.068180884
# 307       0.9998823      0.28043679          0.04900249 0.023052652
# 308       1.0081831      0.25192726          0.04794244 0.013890184
# 309       0.9984632      0.32164076          0.06368760 0.020869429
# 310       1.0015840      0.14288550          0.03291716 0.006617211
# 311       0.9981367      0.28193936          0.04390105 0.058760688
# 312       0.9942832      0.29639977          0.04409836 0.140401124
# 313       1.0461062      0.34656350          0.04573826 0.042387149
# 314       1.0252203      0.39861363          0.05559502 0.008898835
# 315       1.0287882      0.25676308          0.03956636 0.037191868
# 316       1.0458302      0.16428090          0.02650953 0.045198369
# 317       1.0329035      0.30034823          0.04187407 0.031170303
# 318       1.0401953      0.25995752          0.03375763 0.047575454
# 319       1.0664267      0.28934468          0.04751062 0.059260098
#       dry_gonads
# 1   0.0015291822
# 2   0.0078217770
# 3   0.0000000000
# 4   0.0021748873
# 5   0.0049842705
# 6   0.0031049038
# 7   0.0000000000
# 8   0.0103787520
# 9   0.0000000000
# 10  0.0115354838
# 11  0.0037756248
# 12  0.0027743117
# 13  0.0029060401
# 14  0.0001120637
# 15  0.0010240656
# 16  0.0084846137
# 17  0.0011200718
# 18  0.0045977092
# 19  0.0139507653
# 20  0.0000000000
# 21  0.0053788358
# 22  0.0018948371
# 23  0.0000000000
# 24  0.0017243507
# 25  0.0000000000
# 26  0.0215675616
# 27  0.0477002529
# 28  0.0051078431
# 29  0.0175431349
# 30  0.0070658053
# 31  0.0384205947
# 32  0.0429747171
# 33  0.0144706365
# 34  0.0185250516
# 35  0.0350836245
# 36  0.0000000000
# 37  0.0187017466
# 38  0.0294514268
# 39  0.0000000000
# 40  0.0299641473
# 41  0.0000000000
# 42  0.0110218244
# 43  0.0652390279
# 44  0.0052053629
# 45  0.0145428905
# 46  0.0151270507
# 47  0.0240621153
# 48  0.0588417827
# 49  0.0675491376
# 50  0.0175724976
# 51  0.0512611508
# 52  0.0326845898
# 53  0.0602779850
# 54  0.0234899620
# 55  0.0398946639
# 56  0.0000000000
# 57  0.0357546964
# 58  0.0356498605
# 59  0.0322959108
# 60  0.0515603420
# 61  0.0719217127
# 62  0.0745724044
# 63  0.0158472504
# 64  0.0562860827
# 65  0.0409457017
# 66  0.0539535267
# 67  0.0270304990
# 68  0.1767211503
# 69  0.2465255507
# 70  0.2376948448
# 71  0.2049067078
# 72  0.1237671446
# 73  0.2395029332
# 74  0.1153496038
# 75  0.0942980823
# 76  0.0933678267
# 77  0.0696385188
# 78  0.0736199777
# 79  0.0822064301
# 80  0.0271915409
# 81  0.0325464068
# 82  0.0369825558
# 83  0.0171210408
# 84  0.0506901368
# 85  0.1212637963
# 86  0.0495468406
# 87  0.0779769592
# 88  0.0426906893
# 89  0.0494248833
# 90  0.0616629087
# 91  0.0881604590
# 92  0.1076439248
# 93  0.0961183265
# 94  0.0932487347
# 95  0.0954721737
# 96  0.0937577360
# 97  0.1057754238
# 98  0.0751099253
# 99  0.1159286572
# 100 0.1162038018
# 101 0.0919165114
# 102 0.0808013186
# 103 0.1203579525
# 104 0.0883900233
# 105 0.1277356407
# 106 0.1223783004
# 107 0.1442353579
# 108 0.0858967025
# 109 0.0677667650
# 110 0.1353709435
# 111 0.0784921573
# 112 0.0874224016
# 113 0.1439001630
# 114 0.1074281243
# 115 0.1061787934
# 116 0.0925580799
# 117 0.0866163543
# 118 0.0913437857
# 119 0.1062660164
# 120 0.0000000000
# 121 0.0064071980
# 122 0.0393696281
# 123 0.0155806839
# 124 0.0062251197
# 125 0.0962161528
# 126 0.0518223816
# 127 0.1649985010
# 128 0.0054969356
# 129 0.1090602122
# 130 0.1508424272
# 131 0.1997323576
# 132 0.1291595060
# 133 0.0995244811
# 134 0.1539882499
# 135 0.0586096414
# 136 0.1503597058
# 137 0.1443887034
# 138 0.1606225171
# 139 0.2199008689
# 140 0.0442300925
# 141 0.0580786957
# 142 0.1460284593
# 143 0.0446672290
# 144 0.0667660930
# 145 0.1793877855
# 146 0.0927626969
# 147 0.0933567661
# 148 0.0783737942
# 149 0.1763499127
# 150 0.0443697748
# 151 0.1021516838
# 152 0.1339010784
# 153 0.1055952074
# 154 0.1233856095
# 155 0.2280258802
# 156 0.1937573544
# 157 0.1630659365
# 158 0.1224174484
# 159 0.1318757696
# 160 0.2052649329
# 161 0.1801130791
# 162 0.1862602677
# 163 0.1191232958
# 164 0.0191199661
# 165 0.0409794371
# 166 0.0655589553
# 167 0.1728446731
# 168 0.0490951956
# 169 0.1567640611
# 170 0.1654563201
# 171 0.0887720783
# 172 0.1331882413
# 173 0.1275723002
# 174 0.1667350163
# 175 0.1587858844
# 176 0.1251131004
# 177 0.1748298294
# 178 0.1665737326
# 179 0.2366483616
# 180 0.2162155844
# 181 0.2194307702
# 182 0.2287100836
# 183 0.1376096800
# 184 0.1416262793
# 185 0.2793207548
# 186 0.2064514298
# 187 0.1734055080
# 188 0.1239795571
# 189 0.1183884991
# 190 0.1249450968
# 191 0.1435995086
# 192 0.2037409383
# 193 0.1520291901
# 194 0.0453688510
# 195 0.1126254583
# 196 0.1188780133
# 197 0.1163800442
# 198 0.1136305753
# 199 0.0740587733
# 200 0.1096703831
# 201 0.0811103482
# 202 0.1052746853
# 203 0.1146982975
# 204 0.1874109480
# 205 0.0366277755
# 206 0.1056246167
# 207 0.0678348841
# 208 0.1071982789
# 209 0.0574792183
# 210 0.1013185021
# 211 0.1520828324
# 212 0.0799824476
# 213 0.1527629774
# 214 0.1277201547
# 215 0.0170842455
# 216 0.0553159239
# 217 0.0618081910
# 218 0.0443361057
# 219 0.1013066374
# 220 0.0576218572
# 221 0.0104495104
# 222 0.1095956429
# 223 0.0049568151
# 224 0.0195429844
# 225 0.0429358877
# 226 0.0701831336
# 227 0.0075049240
# 228 0.0266477464
# 229 0.0160177441
# 230 0.0103553435
# 231 0.0303703951
# 232 0.0185967848
# 233 0.0166984758
# 234 0.0089910696
# 235 0.0616298199
# 236 0.0115840539
# 237 0.0031238761
# 238 0.0198729826
# 239 0.0022175690
# 240 0.0018160597
# 241 0.0078632212
# 242 0.1731426058
# 243 0.0890780225
# 244 0.1251207620
# 245 0.0646826704
# 246 0.1299938249
# 247 0.1331057803
# 248 0.0922034740
# 249 0.0981043404
# 250 0.0716369075
# 251 0.1063377341
# 252 0.1510610486
# 253 0.1950311713
# 254 0.0879974013
# 255 0.1086879239
# 256 0.1864115420
# 257 0.0546563705
# 258 0.0617370546
# 259 0.2021502096
# 260 0.1865892221
# 261 0.0387417364
# 262 0.0374475495
# 263 0.0705677506
# 264 0.0497435508
# 265 0.1119281612
# 266 0.0834499969
# 267 0.0403711306
# 268 0.0066444412
# 269 0.0019210342
# 270 0.0466553820
# 271 0.0148175048
# 272 0.0719695134
# 273 0.0000000000
# 274 0.0057907394
# 275 0.0323792586
# 276 0.0888018878
# 277 0.0478494774
# 278 0.0346535021
# 279 0.0485807496
# 280 0.0673693152
# 281 0.0541294273
# 282 0.1081102575
# 283 0.0096637783
# 284 0.0621356429
# 285 0.1338816894
# 286 0.0071456097
# 287 0.0853140264
# 288 0.0659434642
# 289 0.0781439053
# 290 0.1108220062
# 291 0.0169311914
# 292 0.0030890918
# 293 0.0432431251
# 294 0.0124943059
# 295 0.0777496310
# 296 0.0587432584
# 297 0.0418675533
# 298 0.0045171276
# 299 0.0722304869
# 300 0.0067143637
# 301 0.0275474218
# 302 0.0171228683
# 303 0.0251087997
# 304 0.0109565123
# 305 0.0233296498
# 306 0.0167379405
# 307 0.0056859754
# 308 0.0022509421
# 309 0.0034223523
# 310 0.0000000000
# 311 0.0471038980
# 312 0.1025314402
# 313 0.0142249909
# 314 0.0000000000
# 315 0.0059743632
# 316 0.0092770106
# 317 0.0041034116
# 318 0.0065574005
# 319 0.0116373247
```

Refaisons notre ACP sur `urchin3` ainsi calculé\ :


```r
urchin3_pca <- pca(urchin3, scale = TRUE)
summary(urchin3_pca)
```

```
# Importance of components (eigenvalues):
#                          PC1   PC2   PC3    PC4    PC5    PC6     PC7
# Variance               4.687 3.353 1.273 0.9666 0.3668 0.1724 0.10547
# Proportion of Variance 0.426 0.305 0.116 0.0879 0.0333 0.0157 0.00959
# Cumulative Proportion  0.426 0.731 0.847 0.9345 0.9678 0.9835 0.99308
#                            PC8     PC9    PC10    PC11
# Variance               0.04761 0.01943 0.00834 0.00068
# Proportion of Variance 0.00433 0.00177 0.00076 0.00006
# Cumulative Proportion  0.99741 0.99918 0.99994 1.00000
# 
# Loadings (eigenvectors, rotation matrix):
#                     PC1    PC2    PC3    PC4    PC5    PC6    PC7   
# diameter1           -0.425         0.145 -0.297               -0.101
# diameter2           -0.425 -0.101  0.143 -0.296               -0.103
# height              -0.427         0.131 -0.300                     
# weight               0.189 -0.428         0.216  0.497 -0.672 -0.145
# solid_parts                -0.495  0.254 -0.117         0.335 -0.245
# integuments         -0.142 -0.463  0.152  0.283  0.165  0.345  0.667
# dry_integuments     -0.259 -0.242  0.173  0.533 -0.669 -0.161 -0.277
# digestive_tract      0.214 -0.370 -0.448 -0.102         0.388 -0.462
# dry_digestive_tract        -0.360 -0.485 -0.401 -0.429 -0.338  0.382
# gonads               0.374         0.438 -0.257 -0.223              
# dry_gonads           0.371         0.440 -0.269 -0.162 -0.122       
#                     PC8    PC9    PC10   PC11  
# diameter1            0.101         0.400  0.714
# diameter2                          0.425 -0.700
# height               0.141  0.197 -0.790       
# weight                      0.109              
# solid_parts         -0.662 -0.223              
# integuments          0.266                     
# dry_integuments                                
# digestive_tract      0.468  0.148              
# dry_digestive_tract -0.155                     
# gonads                      0.725  0.130       
# dry_gonads           0.446 -0.589
```


```r
chart$scree(urchin3_pca)
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-35-1.png" width="672" style="display: block; margin: auto;" />

Maintenant que l'effet saturant est éliminé, la répartition des variances sur les axes se fait mieux. L'axe PC1 contraste les diamètres avec les gonades, l'axe PC2 représente les masses somatiques (dans l'ordre inverse), et l'axe PC3 contraste de manière intéressante les masses du tube digestif avec celles des gonades (le tout en ratios sur la masse immergée, ne l'oublions pas). Les deux premiers axes reprennent 73% de la variance, mais il semble qu'un effet intéressant se marque également sur PC3 avec 85% de la variance totale sur les trois premiers axes.

Tout ceci est également visible sur les graphiques dans l'espace des variables (plans PC1 - PC2 et PC2 - PC3 représentés ici).


```r
chart$loadings(urchin3_pca, choices = c(1, 2))
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-36-1.png" width="672" style="display: block; margin: auto;" />


```r
chart$loadings(urchin3_pca, choices = c(2, 3))
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-37-1.png" width="672" style="display: block; margin: auto;" />

Enfin, dans l'espace des individus, avec l'origine reprise en couleur, nous observons ceci dans le prmeier plan de l'ACP\ :


```r
chart$scores(urchin3_pca, choices = c(1, 2),
  col = urchin2$origin, labels = urchin2$maturity, aspect.ratio = 3/5) +
  theme(legend.position = "right") +
  stat_ellipse()
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-38-1.png" width="672" style="display: block; margin: auto;" />

Et pour le plan PC2 - PC3\ :


```r
chart$scores(urchin3_pca, choices = c(2, 3),
  col = urchin2$origin, labels = urchin2$maturity, aspect.ratio = 3/5) +
  theme(legend.position = "right") +
  stat_ellipse()
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-39-1.png" width="672" style="display: block; margin: auto;" />

Vous devriez pouvoir interpréter ces résultats par vous-même maintenant.


### Visualisation de données quantitatives

#### Deux dimensions

**Le nuage de points** est le graphe idéal pour visualiser la distribution des données bivariées pour deux vaeriavbles quantitatives. Il permet de visualiser également une **association** entre deux variables. Il permet aussi de visualiser comment deux ou plusieurs groupes peuvent être séparés en fonction de ces deux variables.


```r
chart(data = pima, glucose ~ insulin %col=% diabetes) +
  geom_point()
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-40-1.png" width="672" style="display: block; margin: auto;" />


#### Trois dimensions

**Le nuage de points en pseudo-3D** est l’équivalent pour visualiser trois variables quantitatives simultanément. Il est nécessaire de rendre l’effet de la **troisième dimension** (perspective, variation de taille des objets, ...).  La possibilité de **faire tourner l’objet 3D virtuel** est indispensable pour concrétiser l’effet 3D et pour le visionner sous différents angles

Le package `rgl` permet de réaliser ce genre de graphique 3D interactif (que vous pouvez faire tourner dans l'orientation que vous voulez à la souris)\ :




```r
rgl::plot3d(pima$insulin, pima$glucose, pima$triceps,
  col = as.integer(pima$diabetes))
```

<!--html_preserve--><div id="rgl66847" style="width:800px;height:800px;" class="rglWebGL html-widget"></div>
<script type="application/json" data-for="rgl66847">{"x":{"material":{"color":"#000000","alpha":1,"lit":true,"ambient":"#000000","specular":"#FFFFFF","emission":"#000000","shininess":50,"smooth":true,"front":"filled","back":"filled","size":3,"lwd":1,"fog":false,"point_antialias":false,"line_antialias":false,"texture":null,"textype":"rgb","texmipmap":false,"texminfilter":"linear","texmagfilter":"linear","texenvmap":false,"depth_mask":true,"depth_test":"less","isTransparent":false,"polygon_offset":[0,0]},"rootSubscene":1,"objects":{"7":{"id":7,"type":"points","material":{"lit":false},"vertices":[[94,89,23],[168,137,35],[88,78,32],[543,197,45],[846,189,23],[175,166,19],[230,118,47],[83,103,38],[96,115,30],[235,126,41],[146,143,33],[115,125,26],[140,97,15],[110,145,19],[245,158,36],[54,88,11],[192,103,33],[207,111,47],[70,180,25],[240,171,24],[82,103,11],[36,101,15],[23,88,21],[300,176,34],[342,150,42],[304,187,39],[110,100,60],[142,105,41],[128,141,34],[38,95,13],[100,146,27],[90,100,20],[140,139,35],[270,129,20],[71,83,26],[125,110,29],[71,100,25],[110,136,32],[176,123,15],[48,81,40],[64,142,18],[228,144,27],[76,71,18],[64,93,30],[220,122,51],[40,81,18],[152,126,29],[140,144,28],[18,83,31],[36,95,25],[135,171,33],[495,155,26],[37,89,34],[175,160,32],[51,99,15],[100,162,56],[100,107,30],[99,88,42],[135,120,30],[94,118,36],[145,117,24],[168,173,14],[225,170,37],[49,96,13],[140,125,20],[50,100,26],[92,93,25],[325,105,29],[63,108,26],[284,154,31],[119,106,35],[204,136,50],[155,156,28],[485,153,42],[94,99,15],[135,109,21],[53,88,19],[114,163,41],[105,102,40],[285,114,34],[156,104,18],[78,111,12],[130,134,23],[48,79,42],[55,75,24],[130,179,42],[130,129,46],[92,119,18],[495,181,36],[58,128,41],[114,109,39],[160,139,35],[94,123,44],[210,158,41],[48,107,13],[99,109,44],[318,148,27],[44,99,16],[190,103,32],[280,196,29],[87,96,27],[130,140,26],[175,112,32],[271,151,40],[129,109,41],[120,125,30],[478,177,29],[190,142,33],[56,100,15],[32,87,27],[744,197,39],[53,117,31],[370,134,37],[37,79,25],[45,74,28],[192,181,21],[88,91,32],[176,119,22],[194,146,35],[680,165,33],[402,124,33],[55,90,14],[258,92,7],[375,193,16],[150,155,28],[130,191,15],[67,96,18],[56,108,32],[45,71,50],[57,100,52],[116,104,23],[278,108,10],[122,129,28],[155,133,15],[135,136,26],[545,155,44],[220,119,39],[49,96,17],[75,108,43],[40,78,29],[74,107,30],[182,128,37],[194,128,45],[120,151,31],[360,146,38],[215,126,29],[184,100,25],[135,144,33],[42,77,41],[105,120,37],[132,161,23],[148,137,14],[180,128,19],[205,124,28],[148,106,37],[96,155,17],[85,113,10],[94,112,22],[64,99,11],[140,115,39],[231,129,12],[29,152,33],[168,157,21],[156,122,32],[120,102,36],[68,105,32],[52,87,16],[58,95,18],[255,165,43],[171,152,34],[105,130,13],[73,95,21],[108,126,36],[83,139,19],[74,99,19],[43,90,12],[167,125,40],[54,88,40],[249,196,36],[325,189,33],[293,147,25],[83,99,28],[66,81,16],[140,133,28],[465,173,48],[66,84,22],[94,105,40],[158,122,43],[325,140,43],[84,98,15],[75,87,37],[72,93,39],[82,107,30],[182,109,8],[59,90,18],[110,125,24],[50,119,13],[285,144,26],[81,100,23],[196,100,29],[415,131,14],[87,116,12],[275,127,24],[115,96,34],[88,136,41],[165,123,32],[579,172,49],[176,112,30],[310,143,23],[61,143,22],[167,138,35],[474,173,33],[115,129,29],[170,119,41],[76,94,18],[78,102,46],[210,151,32],[277,184,39],[180,181,30],[145,135,46],[180,95,25],[85,89,16],[60,80,11],[50,83,23],[120,117,27],[14,180,63],[70,100,12],[92,95,45],[64,104,37],[63,120,18],[95,82,13],[210,91,32],[105,100,28],[71,86,28],[237,148,48],[60,134,33],[56,120,22],[49,74,40],[105,124,13],[36,74,10],[100,97,36],[140,154,41],[191,105,45],[110,114,17],[75,126,38],[328,158,30],[49,85,22],[125,84,31],[250,135,42],[480,139,41],[265,173,32],[66,83,28],[122,125,18],[76,81,15],[145,195,33],[193,154,32],[71,117,19],[79,94,25],[90,180,26],[170,130,23],[76,84,23],[210,139,17],[86,99,19],[105,163,18],[165,145,34],[326,129,7],[66,68,32],[130,124,33],[82,97,19],[105,116,15],[188,117,31],[106,122,18],[65,86,52],[56,77,30],[210,127,37],[155,129,49],[215,100,40],[190,128,25],[56,84,23],[76,88,29],[225,186,35],[207,187,27],[166,131,21],[67,164,43],[106,84,30],[44,88,24],[115,84,23],[215,124,33],[274,198,32],[77,87,34],[54,99,19],[88,95,14],[18,99,30],[126,92,32],[126,154,29],[165,121,30],[44,111,31],[120,98,17],[330,143,30],[63,119,47],[130,108,20],[600,124,24],[156,176,27],[140,112,50],[115,82,22],[230,123,45],[185,188,14],[25,89,19],[120,109,18],[126,150,29],[293,181,42],[41,92,25],[272,152,39],[182,111,13],[158,106,21],[194,174,22],[321,168,42],[144,138,26],[15,68,13],[160,112,42],[115,94,27],[54,90,47],[90,102,40],[183,128,17],[66,94,18],[91,97,32],[46,100,12],[105,102,17],[152,103,30],[440,157,35],[144,167,17],[159,179,36],[130,136,35],[100,91,25],[106,117,23],[77,123,40],[135,106,28],[540,155,27],[90,101,35],[200,120,48],[70,80,31],[231,167,46],[130,145,46],[132,112,45],[190,98,33],[100,154,30],[168,165,26],[49,68,23],[240,123,35],[265,101,17],[45,56,28],[105,95,39],[205,129,26],[180,140,26],[180,144,46],[95,121,32],[125,129,49],[480,142,24],[125,169,19],[155,127,11],[200,122,27],[100,110,20],[335,127,21],[160,93,32],[387,158,13],[22,126,27],[291,134,20],[392,187,33],[185,173,39],[178,108,46],[200,114,36],[127,149,29],[105,117,30],[180,116,29],[79,130,23],[120,174,37],[165,106,27],[120,126,27],[160,99,17],[150,120,37],[94,102,20],[116,109,18],[140,153,37],[105,100,33],[57,81,41],[200,187,22],[74,121,39],[510,181,44],[110,128,39],[16,88,26],[180,101,48],[112,121,23]],"colors":[[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1]],"centers":[[94,89,23],[168,137,35],[88,78,32],[543,197,45],[846,189,23],[175,166,19],[230,118,47],[83,103,38],[96,115,30],[235,126,41],[146,143,33],[115,125,26],[140,97,15],[110,145,19],[245,158,36],[54,88,11],[192,103,33],[207,111,47],[70,180,25],[240,171,24],[82,103,11],[36,101,15],[23,88,21],[300,176,34],[342,150,42],[304,187,39],[110,100,60],[142,105,41],[128,141,34],[38,95,13],[100,146,27],[90,100,20],[140,139,35],[270,129,20],[71,83,26],[125,110,29],[71,100,25],[110,136,32],[176,123,15],[48,81,40],[64,142,18],[228,144,27],[76,71,18],[64,93,30],[220,122,51],[40,81,18],[152,126,29],[140,144,28],[18,83,31],[36,95,25],[135,171,33],[495,155,26],[37,89,34],[175,160,32],[51,99,15],[100,162,56],[100,107,30],[99,88,42],[135,120,30],[94,118,36],[145,117,24],[168,173,14],[225,170,37],[49,96,13],[140,125,20],[50,100,26],[92,93,25],[325,105,29],[63,108,26],[284,154,31],[119,106,35],[204,136,50],[155,156,28],[485,153,42],[94,99,15],[135,109,21],[53,88,19],[114,163,41],[105,102,40],[285,114,34],[156,104,18],[78,111,12],[130,134,23],[48,79,42],[55,75,24],[130,179,42],[130,129,46],[92,119,18],[495,181,36],[58,128,41],[114,109,39],[160,139,35],[94,123,44],[210,158,41],[48,107,13],[99,109,44],[318,148,27],[44,99,16],[190,103,32],[280,196,29],[87,96,27],[130,140,26],[175,112,32],[271,151,40],[129,109,41],[120,125,30],[478,177,29],[190,142,33],[56,100,15],[32,87,27],[744,197,39],[53,117,31],[370,134,37],[37,79,25],[45,74,28],[192,181,21],[88,91,32],[176,119,22],[194,146,35],[680,165,33],[402,124,33],[55,90,14],[258,92,7],[375,193,16],[150,155,28],[130,191,15],[67,96,18],[56,108,32],[45,71,50],[57,100,52],[116,104,23],[278,108,10],[122,129,28],[155,133,15],[135,136,26],[545,155,44],[220,119,39],[49,96,17],[75,108,43],[40,78,29],[74,107,30],[182,128,37],[194,128,45],[120,151,31],[360,146,38],[215,126,29],[184,100,25],[135,144,33],[42,77,41],[105,120,37],[132,161,23],[148,137,14],[180,128,19],[205,124,28],[148,106,37],[96,155,17],[85,113,10],[94,112,22],[64,99,11],[140,115,39],[231,129,12],[29,152,33],[168,157,21],[156,122,32],[120,102,36],[68,105,32],[52,87,16],[58,95,18],[255,165,43],[171,152,34],[105,130,13],[73,95,21],[108,126,36],[83,139,19],[74,99,19],[43,90,12],[167,125,40],[54,88,40],[249,196,36],[325,189,33],[293,147,25],[83,99,28],[66,81,16],[140,133,28],[465,173,48],[66,84,22],[94,105,40],[158,122,43],[325,140,43],[84,98,15],[75,87,37],[72,93,39],[82,107,30],[182,109,8],[59,90,18],[110,125,24],[50,119,13],[285,144,26],[81,100,23],[196,100,29],[415,131,14],[87,116,12],[275,127,24],[115,96,34],[88,136,41],[165,123,32],[579,172,49],[176,112,30],[310,143,23],[61,143,22],[167,138,35],[474,173,33],[115,129,29],[170,119,41],[76,94,18],[78,102,46],[210,151,32],[277,184,39],[180,181,30],[145,135,46],[180,95,25],[85,89,16],[60,80,11],[50,83,23],[120,117,27],[14,180,63],[70,100,12],[92,95,45],[64,104,37],[63,120,18],[95,82,13],[210,91,32],[105,100,28],[71,86,28],[237,148,48],[60,134,33],[56,120,22],[49,74,40],[105,124,13],[36,74,10],[100,97,36],[140,154,41],[191,105,45],[110,114,17],[75,126,38],[328,158,30],[49,85,22],[125,84,31],[250,135,42],[480,139,41],[265,173,32],[66,83,28],[122,125,18],[76,81,15],[145,195,33],[193,154,32],[71,117,19],[79,94,25],[90,180,26],[170,130,23],[76,84,23],[210,139,17],[86,99,19],[105,163,18],[165,145,34],[326,129,7],[66,68,32],[130,124,33],[82,97,19],[105,116,15],[188,117,31],[106,122,18],[65,86,52],[56,77,30],[210,127,37],[155,129,49],[215,100,40],[190,128,25],[56,84,23],[76,88,29],[225,186,35],[207,187,27],[166,131,21],[67,164,43],[106,84,30],[44,88,24],[115,84,23],[215,124,33],[274,198,32],[77,87,34],[54,99,19],[88,95,14],[18,99,30],[126,92,32],[126,154,29],[165,121,30],[44,111,31],[120,98,17],[330,143,30],[63,119,47],[130,108,20],[600,124,24],[156,176,27],[140,112,50],[115,82,22],[230,123,45],[185,188,14],[25,89,19],[120,109,18],[126,150,29],[293,181,42],[41,92,25],[272,152,39],[182,111,13],[158,106,21],[194,174,22],[321,168,42],[144,138,26],[15,68,13],[160,112,42],[115,94,27],[54,90,47],[90,102,40],[183,128,17],[66,94,18],[91,97,32],[46,100,12],[105,102,17],[152,103,30],[440,157,35],[144,167,17],[159,179,36],[130,136,35],[100,91,25],[106,117,23],[77,123,40],[135,106,28],[540,155,27],[90,101,35],[200,120,48],[70,80,31],[231,167,46],[130,145,46],[132,112,45],[190,98,33],[100,154,30],[168,165,26],[49,68,23],[240,123,35],[265,101,17],[45,56,28],[105,95,39],[205,129,26],[180,140,26],[180,144,46],[95,121,32],[125,129,49],[480,142,24],[125,169,19],[155,127,11],[200,122,27],[100,110,20],[335,127,21],[160,93,32],[387,158,13],[22,126,27],[291,134,20],[392,187,33],[185,173,39],[178,108,46],[200,114,36],[127,149,29],[105,117,30],[180,116,29],[79,130,23],[120,174,37],[165,106,27],[120,126,27],[160,99,17],[150,120,37],[94,102,20],[116,109,18],[140,153,37],[105,100,33],[57,81,41],[200,187,22],[74,121,39],[510,181,44],[110,128,39],[16,88,26],[180,101,48],[112,121,23]],"ignoreExtent":false,"flags":4096},"9":{"id":9,"type":"text","material":{"lit":false},"vertices":[[430,31.9309997558594,-2.49200010299683]],"colors":[[0,0,0,1]],"texts":[["pima$insulin"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[430,31.9309997558594,-2.49200010299683]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"10":{"id":10,"type":"text","material":{"lit":false},"vertices":[[-127.024002075195,127,-2.49200010299683]],"colors":[[0,0,0,1]],"texts":[["pima$glucose"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-127.024002075195,127,-2.49200010299683]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"11":{"id":11,"type":"text","material":{"lit":false},"vertices":[[-127.024002075195,31.9309997558594,35]],"colors":[[0,0,0,1]],"texts":[["pima$triceps"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-127.024002075195,31.9309997558594,35]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"5":{"id":5,"type":"light","vertices":[[0,0,1]],"colors":[[1,1,1,1],[1,1,1,1],[1,1,1,1]],"viewpoint":true,"finite":false},"4":{"id":4,"type":"background","material":{"fog":true},"colors":[[0.298039227724075,0.298039227724075,0.298039227724075,1]],"centers":[[0,0,0]],"sphere":false,"fogtype":"none","flags":0},"6":{"id":6,"type":"background","material":{"lit":false,"back":"lines"},"colors":[[1,1,1,1]],"centers":[[0,0,0]],"sphere":false,"fogtype":"none","flags":0},"8":{"id":8,"type":"bboxdeco","material":{"front":"lines","back":"lines"},"vertices":[[200,"NA","NA"],[400,"NA","NA"],[600,"NA","NA"],[800,"NA","NA"],["NA",100,"NA"],["NA",150,"NA"],["NA","NA",10],["NA","NA",20],["NA","NA",30],["NA","NA",40],["NA","NA",50],["NA","NA",60]],"colors":[[0,0,0,1]],"draw_front":true,"newIds":[19,20,21,22,23,24,25]},"1":{"id":1,"type":"subscene","par3d":{"antialias":0,"FOV":30,"ignoreExtent":false,"listeners":1,"mouseMode":{"left":"trackball","right":"zoom","middle":"fov","wheel":"pull"},"observer":[0,0,2069.89453125],"modelMatrix":[[0.586986541748047,0,0,-252.404205322266],[0,1.1762912273407,8.19500637054443,-436.214233398438],[0,-3.23183345794678,2.98273849487305,-1763.84753417969],[0,0,0,1]],"projMatrix":[[3.73205089569092,0,0,0],[0,3.73205089569092,0,0],[0,0,-3.86370348930359,-7461.73046875],[0,0,-1,0]],"skipRedraw":false,"userMatrix":[[1,0,0,0],[0,0.342020143325668,0.939692620785909,0],[0,-0.939692620785909,0.342020143325668,0],[0,0,0,1]],"userProjection":[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],"scale":[0.586986541748047,3.43924522399902,8.72094345092773],"viewport":{"x":0,"y":0,"width":1,"height":1},"zoom":1,"bbox":[14,846,56,198,7,63],"windowRect":[463,385,719,641],"family":"sans","font":1,"cex":1,"useFreeType":true,"fontname":"/usr/local/lib/R/site-library/rgl/fonts/FreeSans.ttf","maxClipPlanes":6,"glVersion":2.09999990463257,"activeSubscene":0},"embeddings":{"viewport":"replace","projection":"replace","model":"replace","mouse":"replace"},"objects":[6,8,7,9,10,11,5,19,20,21,22,23,24,25],"subscenes":[],"flags":6736},"19":{"id":19,"type":"lines","material":{"lit":false},"vertices":[[200,53.8699989318848,6.15999984741211],[800,53.8699989318848,6.15999984741211],[200,53.8699989318848,6.15999984741211],[200,50.2135009765625,4.71799993515015],[400,53.8699989318848,6.15999984741211],[400,50.2135009765625,4.71799993515015],[600,53.8699989318848,6.15999984741211],[600,50.2135009765625,4.71799993515015],[800,53.8699989318848,6.15999984741211],[800,50.2135009765625,4.71799993515015]],"colors":[[0,0,0,1]],"centers":[[500,53.8699989318848,6.15999984741211],[200,52.041748046875,5.43900012969971],[400,52.041748046875,5.43900012969971],[600,52.041748046875,5.43900012969971],[800,52.041748046875,5.43900012969971]],"ignoreExtent":true,"origId":8,"flags":64},"20":{"id":20,"type":"text","material":{"lit":false},"vertices":[[200,42.9005012512207,1.83399999141693],[400,42.9005012512207,1.83399999141693],[600,42.9005012512207,1.83399999141693],[800,42.9005012512207,1.83399999141693]],"colors":[[0,0,0,1]],"texts":[["200"],["400"],["600"],["800"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[200,42.9005012512207,1.83399999141693],[400,42.9005012512207,1.83399999141693],[600,42.9005012512207,1.83399999141693],[800,42.9005012512207,1.83399999141693]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":8,"flags":2064},"21":{"id":21,"type":"lines","material":{"lit":false},"vertices":[[1.51999998092651,100,6.15999984741211],[1.51999998092651,150,6.15999984741211],[1.51999998092651,100,6.15999984741211],[-19.9039993286133,100,4.71799993515015],[1.51999998092651,150,6.15999984741211],[-19.9039993286133,150,4.71799993515015]],"colors":[[0,0,0,1]],"centers":[[1.51999998092651,125,6.15999984741211],[-9.1919994354248,100,5.43900012969971],[-9.1919994354248,150,5.43900012969971]],"ignoreExtent":true,"origId":8,"flags":64},"22":{"id":22,"type":"text","material":{"lit":false},"vertices":[[-62.7519989013672,100,1.83399999141693],[-62.7519989013672,150,1.83399999141693]],"colors":[[0,0,0,1]],"texts":[["100"],["150"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-62.7519989013672,100,1.83399999141693],[-62.7519989013672,150,1.83399999141693]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":8,"flags":2064},"23":{"id":23,"type":"lines","material":{"lit":false},"vertices":[[1.51999998092651,53.8699989318848,10],[1.51999998092651,53.8699989318848,60],[1.51999998092651,53.8699989318848,10],[-19.9039993286133,50.2135009765625,10],[1.51999998092651,53.8699989318848,20],[-19.9039993286133,50.2135009765625,20],[1.51999998092651,53.8699989318848,30],[-19.9039993286133,50.2135009765625,30],[1.51999998092651,53.8699989318848,40],[-19.9039993286133,50.2135009765625,40],[1.51999998092651,53.8699989318848,50],[-19.9039993286133,50.2135009765625,50],[1.51999998092651,53.8699989318848,60],[-19.9039993286133,50.2135009765625,60]],"colors":[[0,0,0,1]],"centers":[[1.51999998092651,53.8699989318848,35],[-9.1919994354248,52.041748046875,10],[-9.1919994354248,52.041748046875,20],[-9.1919994354248,52.041748046875,30],[-9.1919994354248,52.041748046875,40],[-9.1919994354248,52.041748046875,50],[-9.1919994354248,52.041748046875,60]],"ignoreExtent":true,"origId":8,"flags":64},"24":{"id":24,"type":"text","material":{"lit":false},"vertices":[[-62.7519989013672,42.9005012512207,10],[-62.7519989013672,42.9005012512207,20],[-62.7519989013672,42.9005012512207,30],[-62.7519989013672,42.9005012512207,40],[-62.7519989013672,42.9005012512207,50],[-62.7519989013672,42.9005012512207,60]],"colors":[[0,0,0,1]],"texts":[["10"],["20"],["30"],["40"],["50"],["60"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-62.7519989013672,42.9005012512207,10],[-62.7519989013672,42.9005012512207,20],[-62.7519989013672,42.9005012512207,30],[-62.7519989013672,42.9005012512207,40],[-62.7519989013672,42.9005012512207,50],[-62.7519989013672,42.9005012512207,60]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":8,"flags":2064},"25":{"id":25,"type":"lines","material":{"lit":false},"vertices":[[1.51999998092651,53.8699989318848,6.15999984741211],[1.51999998092651,200.130004882812,6.15999984741211],[1.51999998092651,53.8699989318848,63.8400001525879],[1.51999998092651,200.130004882812,63.8400001525879],[1.51999998092651,53.8699989318848,6.15999984741211],[1.51999998092651,53.8699989318848,63.8400001525879],[1.51999998092651,200.130004882812,6.15999984741211],[1.51999998092651,200.130004882812,63.8400001525879],[1.51999998092651,53.8699989318848,6.15999984741211],[858.47998046875,53.8699989318848,6.15999984741211],[1.51999998092651,53.8699989318848,63.8400001525879],[858.47998046875,53.8699989318848,63.8400001525879],[1.51999998092651,200.130004882812,6.15999984741211],[858.47998046875,200.130004882812,6.15999984741211],[1.51999998092651,200.130004882812,63.8400001525879],[858.47998046875,200.130004882812,63.8400001525879],[858.47998046875,53.8699989318848,6.15999984741211],[858.47998046875,200.130004882812,6.15999984741211],[858.47998046875,53.8699989318848,63.8400001525879],[858.47998046875,200.130004882812,63.8400001525879],[858.47998046875,53.8699989318848,6.15999984741211],[858.47998046875,53.8699989318848,63.8400001525879],[858.47998046875,200.130004882812,6.15999984741211],[858.47998046875,200.130004882812,63.8400001525879]],"colors":[[0,0,0,1]],"centers":[[1.51999998092651,127,6.15999984741211],[1.51999998092651,127,63.8400001525879],[1.51999998092651,53.8699989318848,35],[1.51999998092651,200.130004882812,35],[430,53.8699989318848,6.15999984741211],[430,53.8699989318848,63.8400001525879],[430,200.130004882812,6.15999984741211],[430,200.130004882812,63.8400001525879],[858.47998046875,127,6.15999984741211],[858.47998046875,127,63.8400001525879],[858.47998046875,53.8699989318848,35],[858.47998046875,200.130004882812,35]],"ignoreExtent":true,"origId":8,"flags":64},"32":{"id":32,"type":"linestrip","material":{"alpha":0.498039215803146,"lit":false,"depth_test":"always","isTransparent":true},"vertices":[[0,0,-0.999000012874603],[1,0,-0.999000012874603],[1,1,-0.999000012874603],[0,1,-0.999000012874603],[0,0,-0.999000012874603]],"colors":[[0,0,0,0.498039215803146]],"centers":[[0,0,-0.999000012874603],[1,0,-0.999000012874603],[1,1,-0.999000012874603],[0,1,-0.999000012874603],[0,0,-0.999000012874603]],"ignoreExtent":false,"flags":32864}},"crosstalk":{"key":[null],"group":"","id":0,"options":[null]},"width":800,"height":800,"brushId":32,"sphereVerts":{"vb":[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0746578340503426,0.146446609406726,0.212607523691814,0.270598050073099,0.318189645143208,0.353553390593274,0.375330277517865,0.38268343236509,0.375330277517865,0.353553390593274,0.318189645143208,0.270598050073099,0.212607523691814,0.146446609406726,0.0746578340503426,0,0,0.137949689641472,0.270598050073099,0.392847479193551,0.5,0.587937801209679,0.653281482438188,0.693519922661074,0.707106781186547,0.693519922661074,0.653281482438188,0.587937801209679,0.5,0.392847479193551,0.270598050073099,0.137949689641472,0,0,0.180239955501737,0.353553390593274,0.513279967159337,0.653281482438188,0.768177756711416,0.853553390593274,0.906127446352888,0.923879532511287,0.906127446352888,0.853553390593274,0.768177756711416,0.653281482438188,0.513279967159337,0.353553390593274,0.180239955501737,0,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186548,0.831469612302545,0.923879532511287,0.98078528040323,1,0.98078528040323,0.923879532511287,0.831469612302545,0.707106781186548,0.555570233019602,0.38268343236509,0.195090322016128,0,0,0.180239955501737,0.353553390593274,0.513279967159337,0.653281482438188,0.768177756711416,0.853553390593274,0.906127446352888,0.923879532511287,0.906127446352888,0.853553390593274,0.768177756711416,0.653281482438188,0.513279967159337,0.353553390593274,0.180239955501737,0,0,0.137949689641472,0.270598050073099,0.392847479193551,0.5,0.587937801209679,0.653281482438188,0.693519922661074,0.707106781186548,0.693519922661074,0.653281482438188,0.587937801209679,0.5,0.392847479193551,0.270598050073099,0.137949689641472,0,0,0.0746578340503427,0.146446609406726,0.212607523691814,0.270598050073099,0.318189645143209,0.353553390593274,0.375330277517865,0.38268343236509,0.375330277517865,0.353553390593274,0.318189645143209,0.270598050073099,0.212607523691814,0.146446609406726,0.0746578340503427,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0.0746578340503427,-0.146446609406726,-0.212607523691814,-0.270598050073099,-0.318189645143209,-0.353553390593274,-0.375330277517865,-0.38268343236509,-0.375330277517865,-0.353553390593274,-0.318189645143209,-0.270598050073099,-0.212607523691814,-0.146446609406726,-0.0746578340503427,-0,-0,-0.137949689641472,-0.270598050073099,-0.392847479193551,-0.5,-0.587937801209679,-0.653281482438188,-0.693519922661074,-0.707106781186548,-0.693519922661074,-0.653281482438188,-0.587937801209679,-0.5,-0.392847479193551,-0.270598050073099,-0.137949689641472,-0,-0,-0.180239955501737,-0.353553390593274,-0.513279967159337,-0.653281482438188,-0.768177756711416,-0.853553390593274,-0.906127446352888,-0.923879532511287,-0.906127446352888,-0.853553390593274,-0.768177756711416,-0.653281482438188,-0.513279967159337,-0.353553390593274,-0.180239955501737,-0,-0,-0.195090322016128,-0.38268343236509,-0.555570233019602,-0.707106781186548,-0.831469612302545,-0.923879532511287,-0.98078528040323,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186548,-0.555570233019602,-0.38268343236509,-0.195090322016128,-0,-0,-0.180239955501737,-0.353553390593274,-0.513279967159337,-0.653281482438188,-0.768177756711416,-0.853553390593274,-0.906127446352888,-0.923879532511287,-0.906127446352888,-0.853553390593274,-0.768177756711416,-0.653281482438188,-0.513279967159337,-0.353553390593274,-0.180239955501737,-0,-0,-0.137949689641472,-0.270598050073099,-0.392847479193551,-0.5,-0.587937801209679,-0.653281482438188,-0.693519922661074,-0.707106781186547,-0.693519922661074,-0.653281482438188,-0.587937801209679,-0.5,-0.392847479193551,-0.270598050073099,-0.137949689641472,-0,-0,-0.0746578340503426,-0.146446609406726,-0.212607523691814,-0.270598050073099,-0.318189645143208,-0.353553390593274,-0.375330277517865,-0.38268343236509,-0.375330277517865,-0.353553390593274,-0.318189645143208,-0.270598050073099,-0.212607523691814,-0.146446609406726,-0.0746578340503426,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1],[0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186548,0.831469612302545,0.923879532511287,0.98078528040323,1,0.98078528040323,0.923879532511287,0.831469612302545,0.707106781186548,0.555570233019602,0.38268343236509,0.195090322016128,0,0,0.180239955501737,0.353553390593274,0.513279967159337,0.653281482438188,0.768177756711416,0.853553390593274,0.906127446352888,0.923879532511287,0.906127446352888,0.853553390593274,0.768177756711416,0.653281482438188,0.513279967159337,0.353553390593274,0.180239955501737,0,0,0.137949689641472,0.270598050073099,0.392847479193551,0.5,0.587937801209679,0.653281482438188,0.693519922661074,0.707106781186548,0.693519922661074,0.653281482438188,0.587937801209679,0.5,0.392847479193551,0.270598050073099,0.137949689641472,0,0,0.0746578340503426,0.146446609406726,0.212607523691814,0.270598050073099,0.318189645143209,0.353553390593274,0.375330277517865,0.38268343236509,0.375330277517865,0.353553390593274,0.318189645143209,0.270598050073099,0.212607523691814,0.146446609406726,0.0746578340503426,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0.0746578340503426,-0.146446609406726,-0.212607523691814,-0.270598050073099,-0.318189645143208,-0.353553390593274,-0.375330277517865,-0.38268343236509,-0.375330277517865,-0.353553390593274,-0.318189645143208,-0.270598050073099,-0.212607523691814,-0.146446609406726,-0.0746578340503426,-0,-0,-0.137949689641472,-0.270598050073099,-0.392847479193551,-0.5,-0.587937801209679,-0.653281482438188,-0.693519922661074,-0.707106781186547,-0.693519922661074,-0.653281482438188,-0.587937801209679,-0.5,-0.392847479193551,-0.270598050073099,-0.137949689641472,-0,-0,-0.180239955501737,-0.353553390593274,-0.513279967159337,-0.653281482438188,-0.768177756711416,-0.853553390593274,-0.906127446352888,-0.923879532511287,-0.906127446352888,-0.853553390593274,-0.768177756711416,-0.653281482438188,-0.513279967159337,-0.353553390593274,-0.180239955501737,-0,-0,-0.195090322016128,-0.38268343236509,-0.555570233019602,-0.707106781186548,-0.831469612302545,-0.923879532511287,-0.98078528040323,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186548,-0.555570233019602,-0.38268343236509,-0.195090322016128,-0,-0,-0.180239955501737,-0.353553390593274,-0.513279967159337,-0.653281482438188,-0.768177756711416,-0.853553390593274,-0.906127446352888,-0.923879532511287,-0.906127446352888,-0.853553390593274,-0.768177756711416,-0.653281482438188,-0.513279967159337,-0.353553390593274,-0.180239955501737,-0,-0,-0.137949689641472,-0.270598050073099,-0.392847479193551,-0.5,-0.58793780120968,-0.653281482438188,-0.693519922661074,-0.707106781186548,-0.693519922661074,-0.653281482438188,-0.58793780120968,-0.5,-0.392847479193551,-0.270598050073099,-0.137949689641472,-0,-0,-0.0746578340503427,-0.146446609406726,-0.212607523691814,-0.270598050073099,-0.318189645143209,-0.353553390593274,-0.375330277517866,-0.38268343236509,-0.375330277517866,-0.353553390593274,-0.318189645143209,-0.270598050073099,-0.212607523691814,-0.146446609406726,-0.0746578340503427,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0746578340503427,0.146446609406726,0.212607523691814,0.270598050073099,0.318189645143209,0.353553390593274,0.375330277517865,0.38268343236509,0.375330277517865,0.353553390593274,0.318189645143209,0.270598050073099,0.212607523691814,0.146446609406726,0.0746578340503427,0,0,0.137949689641472,0.270598050073098,0.392847479193551,0.5,0.587937801209679,0.653281482438188,0.693519922661074,0.707106781186547,0.693519922661074,0.653281482438188,0.587937801209679,0.5,0.392847479193551,0.270598050073098,0.137949689641472,0,0,0.180239955501737,0.353553390593274,0.513279967159337,0.653281482438188,0.768177756711416,0.853553390593274,0.906127446352888,0.923879532511287,0.906127446352888,0.853553390593274,0.768177756711416,0.653281482438188,0.513279967159337,0.353553390593274,0.180239955501737,0,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186548,0.831469612302545,0.923879532511287,0.98078528040323,1,0.98078528040323,0.923879532511287,0.831469612302545,0.707106781186548,0.555570233019602,0.38268343236509,0.195090322016128,0]],"it":[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270],[17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288],[18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271]],"material":[],"normals":null,"texcoords":[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1]]},"context":{"shiny":false,"rmarkdown":null},"players":[],"webGLoptions":{"preserveDrawingBuffer":true}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->


#### Plus de trois dimensions

Déjà à trois dimensions la visualisation devient délicate, mais au delà, cela devient pratiquement mission impossible. La **matrice de nuages de points** peut rendre service ici, mais dans certaines limites (tous les angles de vue ne sont pas accessibles).


```r
GGally::ggscatmat(pima, 2:6, color = "diabetes")
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-44-1.png" width="672" style="display: block; margin: auto;" />

Nous voyons qu'ici nous atteignons les limites des possibilités. C'est pour cela que, pour des données multivariées comportant beaucoup de variables quantitatives, les techniques de réduction des dimensions comme l'ACP sont indispensables.


### ACP : mécanisme

Nous allons partir d'un exemple presque trivial pour illuster le principe de l'ACP. Comment réduire un tableau bivarié en une représentation des individus en une seule dimension (classement sur une droite) avec perte minimale d'information\ ? Par exemple, en partant de ces données fictives\ :

<img src="07-acp-afc_files/figure-html/unnamed-chunk-45-1.png" width="864" style="display: block; margin: auto;" />

Voic une représentation graphique 2D de ces données\ :

<img src="07-acp-afc_files/figure-html/unnamed-chunk-46-1.png" width="864" style="display: block; margin: auto;" />

Si nous réduisons à une seule dimension en laissant tomber une des deux variables, voici ce que cela donne (ici on ne garde que `Var1`, donc, on projette les points sur l'axe des abscisses).


<img src="07-acp-afc_files/figure-html/unnamed-chunk-47-1.png" width="432" style="display: block; margin: auto;" />

Au final, nous avons ordonné nos individus en une dimension comme suit\ :

<img src="07-acp-afc_files/figure-html/unnamed-chunk-48-1.png" width="432" style="display: block; margin: auto;" />

C'est une mauvaise solution car il y a trop de perte d'information. Regardez l'écart entre 7 et 9 sur le graphqie en deux dimensions et dans celui à une dimension\ : les points sont trop près. Comparez sur les deux graphiques les distances 7 - 9 avec  9 - 8 et 1 - 2 *versus* 1 - 3. Tout cela est très mal représenté en une seule dimension.

Une autre solution serait de projeter le long de la droite de "tendance générale", c'est-à-dire le long de l'axe de plus grand allongement du nuage de points.


<img src="07-acp-afc_files/figure-html/unnamed-chunk-49-1.png" width="432" style="display: block; margin: auto;" />

Cela donne ceci en une seule dimension\ :

<img src="07-acp-afc_files/figure-html/unnamed-chunk-50-1.png" width="432" style="display: block; margin: auto;" />

C'est une bien meilleure solution car la perte d'information est ici minimale. Regardez à nouveau la distance entre 7 et 9 sur le graphique initial à deux dimensions et sur le nouveau graphique réduit à une dimension\ : c'est mieux qu'avant. Comparez aussi les distances respectives entre les paires 7 - 9 et 9 - 8, ainsi que 1 - 2 par rapport à 1 - 3. Tout cela est bien pieux représenté à présent.

**L’ACP effectue précisément la projection que nous venons d’imaginer.**

- La droite de projection est appelée **composante principale 1**. La composante principale 1 présente la plus grande variabilité possible sur un seul axe.

- Ensuite on calcule la composante 2 comme étant orthogonale (i.e., perpendiculaire) à PC1 et présentant la plus grande variabilité non encore capturée par la composante 1.

- Le mécanisme revient à projeter les points sur des axes orientés différemment dans l'e plan'espace à N dimensions (pour N variables intiales). En effet, mathématiquement ce mécanisme se généralise facilement à trois, puis à N dimensions.


### Calcul matriciel ACP

La rotation optimale des axes vers les PC1 à PCN se résoud par un calcul matriciel. Nous allons maintenant le détailler. Mais auparavant, nous devons nous rafraîchir l'esprit concernant quelques notions.

- Multiplication matricielle\ : $\begin{pmatrix}
2 & 3\\
2 & 1
\end{pmatrix}
\times
\begin{pmatrix}
1\\
3
\end{pmatrix}
=
\begin{pmatrix}
11\\
5
\end{pmatrix}$ 

- Vecteurs propres et valeurs propres (il en existe autant qu’il y a de colonnes dans la matrice de départ) :

$$
\begin{pmatrix}
2 & 3\\
2 & 1
\end{pmatrix}
\times
\begin{pmatrix}
3\\
2
\end{pmatrix}
=
\begin{pmatrix}
12\\
8
\end{pmatrix}
= 4 \times
\begin{pmatrix}
3\\
2
\end{pmatrix}
$$

- La constante (4) est une **valeur propre** et la matrice multipliée (à droite) est la matrice des **vecteurs propres**.

- La rotation d'un système d'axes à deux dimensions d'un angle $\alpha$ peut se représenter sous forme d'un calcul matriciel\ :


$$
\begin{pmatrix}
\cos \alpha & \sin \alpha\\
-\sin \alpha & \cos \alpha
\end{pmatrix}
\times
\begin{pmatrix}
x\\
y
\end{pmatrix}
=
\begin{pmatrix}
x'\\
y'
\end{pmatrix}
$$

Dans le cas particulier de l’ACP, la matrice de transformation qui effectue la rotation voulue pour obtenir les axes principaux est **la matrice rassemblant tous les vecteurs propres calculés après diagonalisation de la matrice de corrélation ou de variance/covariance** (réduction ou non, respectivement). Le schéma suivant visualise la rotation depuis les axes initiaux X et Y (variables de départ) en bleu royal vers les PC1, PC2 en rouge. Un individu p est représenté par les coordonnées {x, y} dans le système d'axes initial XY. Les nouvelles coordonnées {x', y'} sont recalculées par projection sur les nouveaux axes PC1-PC2. Les flèches bleues sont représentées dans l'espace des variables, tandis que les points reprojettés sur PC1-PC2 sont représentés dans l'espace des individus selon les coordonnées primes en rouge.

<img src="07-acp-afc_files/figure-html/unnamed-chunk-51-1.png" width="432" style="display: block; margin: auto;" />


#### Résolution numérique simple

Effectuons une ACP sur matrice var/covar sans réduction des données (mais calcul très similaire lorsque les données sont réduites) sur un exemple numérique simple.

- Étape 1\ :  centrage des données

$$
\mathop{\begin{pmatrix}
2 & 1 \\
3 & 4 \\
5 & 0 \\
7 & 6 \\
9 & 2
\end{pmatrix}}_{\text{Tableau brut}}
\xrightarrow{\phantom{---}\text{centrage}\phantom{---}}
\mathop{\begin{pmatrix}
-3.2 & -1.8 \\
-2.2 & \phantom{-}1.4 \\
-0.2 & -2.6 \\
\phantom{-}1.8 & \phantom{-}3.4 \\
\phantom{-}3.8 & -0.6
\end{pmatrix}}_{\text{Tableau centré (X)}}
$$

- Étape 2\ :  calcul de la matrice de variance/covariance

$$
\mathop{\begin{pmatrix}
-3.2 & -1.8 \\
-2.2 & \phantom{-}1.4 \\
-0.2 & -2.6 \\
\phantom{-}1.8 & \phantom{-}3.4 \\
\phantom{-}3.8 & -0.6
\end{pmatrix}}_{\text{Tableau centré (X)}}
\xrightarrow{\phantom{---}\text{var/covar}\phantom{---}}
\mathop{\begin{pmatrix}
8.2 & 1.6 \\
1.6 & 5.8
\end{pmatrix}}_{\text{Matrice carrée (A)}}
$$

- Étape 3\ : diagonalisation de la matrice var/covar

$$
\mathop{\begin{pmatrix}
8.2 & 1.6 \\
1.6 & 5.8
\end{pmatrix}}_{\text{Matrice carrée (A)}}
\xrightarrow{\phantom{---}\text{diagonalisation}\phantom{---}}
\mathop{\begin{pmatrix}
9 & 0 \\
0 & 5
\end{pmatrix}}_{\text{Matrice diagonalisée (B)}}
$$

- La **trace** des deux matrices A et B (somme des éléments sur la diagonale) est égale à : 8.2 + 5.8 = 14 = 9 + 5.
- 8.2 est la **part de variance** exprimée sur le premier axe initial (X)
- 5.8 est la **part de variance** exprimée sur le second axe initial (Y)
- 14 est la **variance totale** du jeu de données
- La matrice diagonale B est la solution exprimant **la plus grande part de variance possible sur le premier axe de l’ACP**\ : 9, soit 64,3% de la variance totale.
- Les éléments sur la diagonale sont les valeurs propres $\lambda_i$\ ! Vous vous rappelez les fameuses "eigenvalues" dans la sortie de `summary(pima_pca)`.


- Étape 4\ : calcul de la matrice de rotation des axes (en utilisant la propriété des valeurs propres $\text{A}.\text{U} = \text{B}.\text{U}$).

$$
\mathop{\begin{pmatrix}
8.2 & 1.6 \\
1.6 & 5.8
\end{pmatrix}}_{\text{Matrice A}}
\times
\text{U}
=
\mathop{\begin{pmatrix}
9 & 0 \\
0 & 5
\end{pmatrix}}_{\text{Matrice B}}
\times
\text{U}
\rightarrow
\text{U}
=
\mathop{\begin{pmatrix}
\phantom{-}0.894 & -0.447 \\
\phantom{-}0.447 & \phantom{-}0.894
\end{pmatrix}}_{\text{Matrice des vecteur propres (U)}}
$$

- La **matrice des vecteurs propres (U)** ("eigenvectors" en anglais) effectue la transformation (**rotation des axes**) pour obtenir les **composantes principales**.
- L'angle de rotation se déduit en considérant que cette matrice contient des sin et cos d'angles de rotation des axes\ :

$$
\begin{pmatrix}
\phantom{-}0.894 & -0.447 \\
\phantom{-}0.447 & \phantom{-}0.894
\end{pmatrix}
=
\begin{pmatrix}
\phantom{-}\cos(-26.6°) & \phantom{-}\sin(-26.6°) \\
-\sin(-26.6°) & \phantom{-}\cos(-26.6°)
\end{pmatrix}
$$

- Étape 5\ : représentation dans l'espace des variables. C'est une représentation dans un cercle de la matrice des vecteurs propres U sous forme de vecteurs.

- Étape 6\ : représentation dans l'espace des individus. On recalcule les coordonnées des individus dans le système d'axe après rotation.

$$
\mathop{\begin{pmatrix}
-3.2 & -1.8 \\
-2.2 & \phantom{-}1.4 \\
-0.2 & -2.6 \\
\phantom{-}1.8 & \phantom{-}3.4 \\
\phantom{-}3.8 & -0.6
\end{pmatrix}}_{\text{Tableau centré (X)}}
\times
\mathop{\begin{pmatrix}
\phantom{-}0.894 & -0.447 \\
\phantom{-}0.447 & \phantom{-}0.894
\end{pmatrix}}_{\text{Matrice des vecteur propres (U)}}
\xrightarrow{\phantom{---}\text{X}.\text{U} = \text{X'}\phantom{---}}
\mathop{\begin{pmatrix}
-3.58 & \phantom{-}0.00 \\
-1.34 & \phantom{-}2.24 \\
-1.34 & -2.24 \\
\phantom{-}3.13 & \phantom{-}2.24 \\
\phantom{-}3.13 & -2.24
\end{pmatrix}}_{\text{Tableau avec rotation (X')}}
$$

- Ensuite, on représente ces individus à l'aide d'un graphique en nuage de points.


Tout ces calculs se généralisent facilement à trois, puis à N dimensions.


##### Pour aller plus loin {-}

- N'hésitez pas à **combiner** plusieurs techniques. Par exemple, vous pouvez représenter les groupes créés par classification ascendante hiérarchiques sur un graphique de l'ACP dans l'espace des individus en faisant varier les couleurs ou les labels des individus en fonction des groupes de la CAH.

- Une [autre explication de l'ACP](http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/) en utilsant quelques autres fonctions de R pour visualiser le résultat

- Une [explication détaillée de la PCA en anglais](http://staff.ustc.edu.cn/~zwp/teach/MVA/abdi-awPCA2010.pdf).

- Une page qui reprend [une série de vidéos](http://www.sthda.com/english/articles/21-courses/65-principal-component-analysis-course-using-factominer/) qui présentent les différentes facettes de l'A CP (en franglais).


## Analyse factorielle des correspondances

Comme l'ACP s'intéresse à des corrélations linéaires entre variables quantitatives, elle n'est absolument pas utilisable pour traiter des variables qualitatives. L'**Analyse Factorielle des Correspondances** sera utile dans ce dernier cas (AFC, ou en anglais "Correspondence Analysis" ou CA).


### AFC dans SciViews::R

L'AFC utilises la fonction `ca()` du package `ca` dans `SciViews::R`, mais au stade actuel, tout le code nécessaire (en particulier pour réaliser les graphiques avec `chart()`) n'est pas encore complètement intégré dans les packages. Ainsi, vous pouvez copier-coller le code du chunk suivant au début de vos scripts ou dans un chunk de `setup` dans vos documenbts R Markdown/Notebook.


```r
SciViews::R
library(broom)

# Au lieu de MASS::corresp(, nf = 2), nous préférons ca::ca()
ca <- ca::ca

scale_axes <- function(data, aspect.ratio = 1) {
  range_x <- range(data[, 1])
  span_x <- max(range_x) - min(range_x)
  range_y <- range(data[, 2])
  span_y <- max(range_y) - min(range_y)
  if (span_y * aspect.ratio < span_x) {
    # Adjust range_x
    span_x_2 <- span_y / aspect.ratio / 2
    range_x_mid <- sum(range_x) / 2
    range_x <- c(range_x_mid - span_x_2, range_x_mid + span_x_2)
  } else {
    # Adjust range_y
    span_y_2 <- span_x * aspect.ratio / 2
    range_y_mid <- sum(range_y) / 2
    range_y <- c(range_y_mid - span_y_2, range_y_mid + span_y_2)
  }
  list(x = range_x, y = range_y)
}

plot3d <- rgl::plot3d
plot3d.ca <- ca:::plot3d.ca

autoplot.ca <- function(object, choices = 1L:2L,
type = c("screeplot", "altscreeplot", "biplot"), col = "black", fill = "gray",
aspect.ratio = 1, repel = FALSE, ...) {
  type = match.arg(type)

  res <- switch(type,
    screeplot = object %>.% # Classical screeplot
      `[[`(., "sv") %>.%
      tibble(Dimension = 1:length(.), sv = .) %>.%
      chart(data = ., sv^2 ~ Dimension) +
      geom_col(col = col, fill = fill) +
      labs(y = "Inertia"),

    altscreeplot = object %>.% # screeplot represented by dots and lines
      `[[`(., "sv") %>.%
      tibble(Dimension = 1:length(.), sv = .) %>.%
      chart(data = ., sv^2 ~ Dimension) +
      geom_line(col = col) +
      geom_point(col = "white", fill = col, size = 2, shape = 21, stroke = 3) +
      labs(y = "Inertia"),

    biplot = {
      # We want to use the function plot.ca(), but without plotting the base plot
      # So, we place it in a specific environment where all base plot functions are
      # fake and do nothing (we just want to collect points coordinates at the end)
      env <- new.env()
      env$plot_ca <- ca:::plot.ca
      environment(env$plot_ca) <- env
      env$plot <- function(...) NULL
      env$box <- function(...) NULL
      env$abline <- function(...) NULL
      env$axis <- function(...) NULL
      env$par <- function(...) NULL
      env$points <- function(...) NULL
      env$lines <- function(...) NULL
      env$.arrows <- function(...) NULL
      env$text <- function(...) NULL
      env$strwidth <- function(...) NULL
      env$strheight <- function(...) NULL

      contribs <- paste0("Dimension ", 1:length(object$sv), " (",
        round(object$sv^2 / sum(object$sv^2) * 100, 1), "%)")[choices]

      res <- env$plot_ca(object, dim = choices, ...)

      rows <- as.data.frame(res$rows)
      rows$Type <- "rows"
      rows$Labels <- object$rownames
      cols <- as.data.frame(res$cols)
      cols$Type <- "cols"
      cols$Labels <- object$colnames
      res <- bind_rows(rows, cols)
      names(res) <- c("x", "y", "type", "labels")

      lims <- scale_axes(res, aspect.ratio = aspect.ratio)
      nudge <- (lims$x[2] - lims$x[1]) / 100

      res <- chart(data = res, y ~ x %col=% type %label=% labels) +
        geom_hline(yintercept = 0, col = "gray") +
        geom_vline(xintercept = 0, col = "gray") +
        coord_fixed(ratio = 1, xlim = lims$x, ylim = lims$y, expand = TRUE) +
        theme(legend.position = "none") +
        labs(x = contribs[1], y = contribs[2])

      if (isTRUE(repel)) {
        res <- res + geom_point() + ggrepel::geom_text_repel()
      } else {# Use text
        res <- res + geom_point() +
          geom_text(hjust = 0, vjust = 0, nudge_x = nudge, nudge_y = nudge)
      }
      res
    }
  )
  res
}
chart.ca <- function(data, choices = 1L:2L, ...,
type = c("screeplot", "altscreeplot", "boiplot"), env = parent.frame())
  autoplot.ca(data, choices = choices, ..., type = type, env = env)
class(chart.ca) <- c("function", "subsettable_type")
```

