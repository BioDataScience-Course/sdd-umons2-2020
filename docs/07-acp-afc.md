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

- Pour **visualiser** cette structure, les données sont simplifiées (réduites) de **N variables à n (n < N et n = 2 ou 3)**. La représentation sous forme d'un nuage de points s'appelle une **carte**.

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
  if (span_y * aspect.ratio > span_x) {
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

Les indiens Pimas sont des amérindiens originaires du nord du Mexique qui sont connus pour compter le plus haut pourcentage d'obèses et de diabétiques de toutes les éthnies. Ils ont fait l'objet de plusieurs études scientifiques d'autant plus que les Pimas en Arizona développent principalement cette obésité et ce diabète, alors que les Pimas mexicains les ont plus rarement. Il est supposé que leur mode de vie différent aux Etats_Units pourrait en être la raison. Voici un jeu de données qui permet d'explorer un peu ceci (nous éliminons les lignes qui contiennent des données manquantes dans le tableau initial à l'aide de `drop_na()`\ :


```r
read("PimaIndiansDiabetes2", package = "mlbench") %>.%
  drop_na(.) -> pima
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

Avant de nous lancer dans une ACP, nous devons décrire les données, repérer les variables quantitatives d'intérêt, et synthétiser les corrélations linéaires (coefficients de corrélation de Pearson) entre ces variables.


```r
skimr::skim(pima)
```

```
# Skim summary statistics
#  n obs: 392 
#  n variables: 9 
# 
# ── Variable type:factor ──────────────────────────────────────────────────────────────────────────────
#  variable missing complete   n n_unique                top_counts ordered
#  diabetes       0      392 392        2 neg: 262, pos: 130, NA: 0   FALSE
# 
# ── Variable type:numeric ─────────────────────────────────────────────────────────────────────────────
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

<img src="07-acp-afc_files/figure-html/unnamed-chunk-4-1.png" width="672" style="display: block; margin: auto;" />

Quelques corrélations positives d'intensités moyennes se dégagent ici, notamment entre `mass` et `triceps` (épaisseur du pli cutané au niveau du triceps), ainsi qu'entre `glucose` (taux de glucose dans le sang), `insulin` (taux d'insuline dans le sang) et `age`. Par contre, la pression artérielle (`pressure`) et le `pedigree` (variable qui quantifie la susceptibilité au diabète en fonction de la parenté) semblent peu corrélés avec les autres variables.

L'ACP est en fait équivalente à une Analyse en Coordonnées Principales sur une matrice de distances euclidiennes (MDS métrique), mais en plus efficace en terme de calculs. Nous pouvons donc nous lancer dans l'analyse et en comprendre les résultats en gardant ceci à l'esprit.

Nous utiliserons la fonction `pca()` qui prend un argument `data =` et une formule du type `~ var1 + var2 + .... + varn`, ou plus simplement, directement un tableau contenant uniquement les variables à analyser comme argument unique. Comme les différentes variables sont mesurées dans des unités différentes, nous devons les standardiser (écart type ramené à un pour toutes). Ceci est réalisé par la fonction `pca()` en lui indiquant `scale = TRUE`. Donc\ :


```r
pima_pca <- pca(data = pima, ~ glucose + pressure + triceps + insulin + mass +
  pedigree + age, scale = TRUE)
```

Ou alors, nous selectionnons les variables d'intérêt avec `select()` et appliquons `pca()` directement sur ce tableau, ce qui donnera le même résultat.


```r
pima %>.%
  select(., glucose:age) %>.%
  pca(., scale = TRUE) -> pima_pca
```

Le nuage de points dans l'espace initial à sept dimensions a été centré (origine ramenée au centre de gravité du nuage de points = moyenne des variables). Ensuite une rotation des axes a été réalisée pour orienter son plus grand axe selon un premier **axe principal 1** ou **PC1** . Ensuite **PC2** est construit orthogonal au premier et dans la seconde direction de plus grande variabilité du nuage de points, et ainsi de suite pour les autres axes. Ainsi les axes PC1, PC2, PC3, ... représentent une **part de variance** de plus en plus faible par rapport à la variance totale du jeu de données. Ceci est présenté dans le résumé\ :


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

Le premier tableau `Importance of components (eigenvalues):` montre la part de variance présentée sur chacun des sept axes de l'ACP (PC1, PC2, ..., PC7). Le fait qu'il s'agit de *valeurs propres* (*eigenvalues* en anglais) apparaitra plus clair lorsque vous aurez lu les explications détaillées plus bas. Ces parts de variance s'additionnent pour donner la variance totale du nuage de points dans les sept dimensions (propriété d'additivité des variances). Pour facilité la lecture, la `Proportion de Variance` en %  est reprise également, ainsi que les proportions cumulées. Ainsi, les deux premiers axes de l'ACP capturent ici 53% de la variance totale. Et il faudrait considérer les cinq premiers axes pour capturer 90% de la variance totale.

Le second tableau `Loadings (eigenvectors, rotation matrix):` est la matrice de transformation des coordonnées initiales sur les lignes en coordonnées PC1 à PC7 en colonnes. Nous pouvons y lire l'**importante** des variables initiales sur les axes de l'ACP. Par exemple, l'axe PC3 contraste essentiellement `pressure` et `pedigree`.

Le **graphique des éboulis** sert à visualiser la "chute" de la variance d'un axe principal à l'autre, et aide à choisir le nombre d'axes à conserver (espace à dimensions réduites avec perte minimale d'information). Deux variantes en diagramme en barres versticales `chart$screeplot()` ou `chart$scree()` ou sous forme d'une ligne brisée `chart$altscree()` sont disponibles\ :


```r
chart$scree(pima_pca, fill = "cornsilk")
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-8-1.png" width="672" style="display: block; margin: auto;" />


```r
chart$altscree(pima_pca)
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-9-1.png" width="672" style="display: block; margin: auto;" />

La diminution est importante entre le premier et le second axe, mais plus progressive ensuite. Ceci traduit une structure plus complexe dans les données qui ne se réduit pas facilement à un très petit nombre d'axes. Nous pouvons visualiser le **premier plan principal** constitué par PC1 et PC2, tout en gardant à l'esprit que seulement 53% de la variance totale y est capturée. Donc, nous pouvons nous attendre à des déformations non négligeables des données dans ce plan. Nous verrons qu'il est porteur, toutefois, d'information utile.

Deux types de représentations peuvent être réalisées à partir d'ici\ : la représentation dans **l'espace des variables**, et la représentation complémentaire dans **l'espace des individus**. Ces deux représentations sont complémentaires et s'analysent conjointement. L'espace des variables représente les axes initiaux projettés comme des ombres dans le plan choisi de l'ACP. Il se réalise à l'aide de `chart$loadings()`. Par exemple pour PC1 et PC2 (`choices = c(1, 2)`)\ :


```r
chart$loadings(pima_pca, choices = c(1, 2))
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-10-1.png" width="672" style="display: block; margin: auto;" />

Ce graphique s'interpète selon les critères suivants\ :

- Plus la norme (longueur) du vecteur qui représente une variable est grande et se rapporche de un (matérialisé par le cer cle gris), plus la variable est bien représentée dans le plan choisi. On évitera d'interpréter ici les variables qui ont des normes petites, comme `pedigree` ou `pressure`.

- Des vecteurs qui pointent dans la même direction représentent des variables **directement corrélés** entre elles. C'est le cas de `glucose`, `insulin` et `age`d'une part, et par ailleurs aussi de `mass` et `triceps`.

- Des vecteurs qui pointent en directions opposées représentent des variables **inversément proportionnelles**. Il n'y en a pas ici.

- Des vecteurs orthogonaux représentent des variables **non corrélées** entre elles. ainsi le groupoe `glucose`/`insulin`/`age` n'est pas corrélé avec le groupe `mass`/`triceps`

- Les PCs sont orientés en fonction des variables initiales, ou à défaut, les zones du graphique sont orientés. Ici, les gros sont dans le haut à droite du graphique, alors que ceux qui sont agés, et ont beaucoup de sucre et d'insuline dans le sang sont en bas à droite. A l'opposé, on trouve les plus maigres en bas à gauche et les jeunes ayant moins de glucose et d'insuline dans le sang en haut à gauche du graphique.

Cela donne déjà une vision synthétique des différentes corrélations entre la variables. Naturellement, on peut très bien choisir d'autres axes, pour peu qu'ils représentent une part de variance relativement importante. Par exemple, ici, nous pouvons représente le plan constitué par PC1 et PC3\ :


```r
chart$loadings(pima_pca, choices = c(1, 3))
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-11-1.png" width="672" style="display: block; margin: auto;" />

Nous voyons que `pedigree` et `pressure` (inversément proportionnels) sont bien mieux représentés le long de PC3. Ici l'axe PC3 est plus facile à orienter\ : en haut les pédigrées élevés et les pressions qartérielles basses, et en bas le contraire.

La seconde représentation se fait dans **l'espace des individus**. Ici, nous allons projeter les points relatifs à chaque individu dans le plan de l'ACP choisi. Cela se réalise à l'aide de `chart$scores()` (l'aspect ratio est le rapport hauteur/largeur peut s'adapter)\ :


```r
chart$scores(pima_pca, choices = c(1, 2), aspect.ratio = 3/5)
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-12-1.png" width="672" style="display: block; margin: auto;" />

Ce graphique est peu lisible tel quel. Généralement, nous représentons d'autres informations utiles sous forme de labels et ou de couleurs différentes. Nous pouvons ainsi contraster les individus qui ont le diabète de ceux qui ne l'ont pas sur ce graphique et aussi ajouter des ellipses de confiance à 95% autour des deux groupes pour aider à la cerner à l'aide de `stat_ellipse()`\ :


```r
chart$scores(pima_pca, choices = c(1, 2),
  labels = pima$diabetes) +
  stat_ellipse()
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-13-1.png" width="672" style="display: block; margin: auto;" />

Ce graphique est nettement plus intéressant. Il s'interprète comme suit\ :

- Nous savons que les individus plus âgés et ayant plus de glucose et d'insuline dans le sang sont dans le bas à droite du graphique. Or le groupe des diabétique, s'il ne se détache pas complètement tend à s'étaler plus dans cette région.

- A l'inverse, le groupe des non diabétiques s'étale vers la gauche, c'est-à-dire dans une région reprenant les individus les plus jeunes et les plus maigres.

Le graphique entre PC1 et PC3 donne ceci\ :


```r
chart$scores(pima_pca, choices = c(1, 3),
  labels = pima$diabetes) +
  stat_ellipse()
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-14-1.png" width="672" style="display: block; margin: auto;" />

Ici, la séparation se fait essentiellement sur l'axe horizontal (PC1). Donc, les différentes de pédigrée (élevé dans le haut du graphique) et de pression artérielle (élevée dans le bas du graphique) semblent être moins liés au diabète.

Etant donné que les deux graphiques (variables et individus) s'interprètent conjointement, nous pourrions être tentés de les superposer, cela s'appelle un **biplot**. Mais se pose alors un problème\ : celui de mettre à l'échelle les deux représentations pour qu'elles soient cohérentes entre elles. Ceci n'est pas facile et différentes représentations coexistent. L'argument `scale =` de la fonction `chart$biplot()` permet d'utiliser différentes mises à l'échelle. Enfin, ce type de graphique tend à être souvent bien trop encombré. Il est donc plus difficile à lire que les deux graphiques des variables et individus séparés. Voici ce que cela donne pour notre jeu de données exemple\ :


```r
chart$biplot(pima_pca)
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-15-1.png" width="672" style="display: block; margin: auto;" />

Bien moins lisible\ !


### Visualisation de données quantitatives

#### Deux dimensions

**Le nuage de points** est le graphe idéal pour visualiser la distribution des données bivariées pour deux vaeriavbles quantitatives. Il permet de visualiser également une **association** entre deux variables. Il permet aussi de visualiser comment deux ou plusieurs groupes peuvent être séparés en fonction de ces deux variables.


```r
chart(data = pima, glucose ~ insulin %col=% diabetes) +
  geom_point()
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-16-1.png" width="672" style="display: block; margin: auto;" />


#### Trois dimensions

**Le nuage de points en pseudo-3D** est l’équivalent pour visualiser trois variables quantitatives simultanément. Il est nécessaire de rendre l’effet de la **troisième dimension** (perspective, variation de taille des objets, ...).  La possibilité de **faire tourner l’objet 3D virtuel** est indispensable pour concrétiser l’effet 3D et pour le visionner sous différents angles

Le package `rgl` permet de réaliser ce genre de graphique 3D interactif (que vous pouvez faire tourner dans l'orientation que vous voulez à la souris)\ :




```r
rgl::plot3d(pima$insulin, pima$glucose, pima$triceps,
  col = as.integer(pima$diabetes))
```

<!--html_preserve--><div id="rgl11470" style="width:800px;height:800px;" class="rglWebGL html-widget"></div>
<script type="application/json" data-for="rgl11470">{"x":{"material":{"color":"#000000","alpha":1,"lit":true,"ambient":"#000000","specular":"#FFFFFF","emission":"#000000","shininess":50,"smooth":true,"front":"filled","back":"filled","size":3,"lwd":1,"fog":false,"point_antialias":false,"line_antialias":false,"texture":null,"textype":"rgb","texmipmap":false,"texminfilter":"linear","texmagfilter":"linear","texenvmap":false,"depth_mask":true,"depth_test":"less","isTransparent":false,"polygon_offset":[0,0]},"rootSubscene":1,"objects":{"7":{"id":7,"type":"points","material":{"lit":false},"vertices":[[94,89,23],[168,137,35],[88,78,32],[543,197,45],[846,189,23],[175,166,19],[230,118,47],[83,103,38],[96,115,30],[235,126,41],[146,143,33],[115,125,26],[140,97,15],[110,145,19],[245,158,36],[54,88,11],[192,103,33],[207,111,47],[70,180,25],[240,171,24],[82,103,11],[36,101,15],[23,88,21],[300,176,34],[342,150,42],[304,187,39],[110,100,60],[142,105,41],[128,141,34],[38,95,13],[100,146,27],[90,100,20],[140,139,35],[270,129,20],[71,83,26],[125,110,29],[71,100,25],[110,136,32],[176,123,15],[48,81,40],[64,142,18],[228,144,27],[76,71,18],[64,93,30],[220,122,51],[40,81,18],[152,126,29],[140,144,28],[18,83,31],[36,95,25],[135,171,33],[495,155,26],[37,89,34],[175,160,32],[51,99,15],[100,162,56],[100,107,30],[99,88,42],[135,120,30],[94,118,36],[145,117,24],[168,173,14],[225,170,37],[49,96,13],[140,125,20],[50,100,26],[92,93,25],[325,105,29],[63,108,26],[284,154,31],[119,106,35],[204,136,50],[155,156,28],[485,153,42],[94,99,15],[135,109,21],[53,88,19],[114,163,41],[105,102,40],[285,114,34],[156,104,18],[78,111,12],[130,134,23],[48,79,42],[55,75,24],[130,179,42],[130,129,46],[92,119,18],[495,181,36],[58,128,41],[114,109,39],[160,139,35],[94,123,44],[210,158,41],[48,107,13],[99,109,44],[318,148,27],[44,99,16],[190,103,32],[280,196,29],[87,96,27],[130,140,26],[175,112,32],[271,151,40],[129,109,41],[120,125,30],[478,177,29],[190,142,33],[56,100,15],[32,87,27],[744,197,39],[53,117,31],[370,134,37],[37,79,25],[45,74,28],[192,181,21],[88,91,32],[176,119,22],[194,146,35],[680,165,33],[402,124,33],[55,90,14],[258,92,7],[375,193,16],[150,155,28],[130,191,15],[67,96,18],[56,108,32],[45,71,50],[57,100,52],[116,104,23],[278,108,10],[122,129,28],[155,133,15],[135,136,26],[545,155,44],[220,119,39],[49,96,17],[75,108,43],[40,78,29],[74,107,30],[182,128,37],[194,128,45],[120,151,31],[360,146,38],[215,126,29],[184,100,25],[135,144,33],[42,77,41],[105,120,37],[132,161,23],[148,137,14],[180,128,19],[205,124,28],[148,106,37],[96,155,17],[85,113,10],[94,112,22],[64,99,11],[140,115,39],[231,129,12],[29,152,33],[168,157,21],[156,122,32],[120,102,36],[68,105,32],[52,87,16],[58,95,18],[255,165,43],[171,152,34],[105,130,13],[73,95,21],[108,126,36],[83,139,19],[74,99,19],[43,90,12],[167,125,40],[54,88,40],[249,196,36],[325,189,33],[293,147,25],[83,99,28],[66,81,16],[140,133,28],[465,173,48],[66,84,22],[94,105,40],[158,122,43],[325,140,43],[84,98,15],[75,87,37],[72,93,39],[82,107,30],[182,109,8],[59,90,18],[110,125,24],[50,119,13],[285,144,26],[81,100,23],[196,100,29],[415,131,14],[87,116,12],[275,127,24],[115,96,34],[88,136,41],[165,123,32],[579,172,49],[176,112,30],[310,143,23],[61,143,22],[167,138,35],[474,173,33],[115,129,29],[170,119,41],[76,94,18],[78,102,46],[210,151,32],[277,184,39],[180,181,30],[145,135,46],[180,95,25],[85,89,16],[60,80,11],[50,83,23],[120,117,27],[14,180,63],[70,100,12],[92,95,45],[64,104,37],[63,120,18],[95,82,13],[210,91,32],[105,100,28],[71,86,28],[237,148,48],[60,134,33],[56,120,22],[49,74,40],[105,124,13],[36,74,10],[100,97,36],[140,154,41],[191,105,45],[110,114,17],[75,126,38],[328,158,30],[49,85,22],[125,84,31],[250,135,42],[480,139,41],[265,173,32],[66,83,28],[122,125,18],[76,81,15],[145,195,33],[193,154,32],[71,117,19],[79,94,25],[90,180,26],[170,130,23],[76,84,23],[210,139,17],[86,99,19],[105,163,18],[165,145,34],[326,129,7],[66,68,32],[130,124,33],[82,97,19],[105,116,15],[188,117,31],[106,122,18],[65,86,52],[56,77,30],[210,127,37],[155,129,49],[215,100,40],[190,128,25],[56,84,23],[76,88,29],[225,186,35],[207,187,27],[166,131,21],[67,164,43],[106,84,30],[44,88,24],[115,84,23],[215,124,33],[274,198,32],[77,87,34],[54,99,19],[88,95,14],[18,99,30],[126,92,32],[126,154,29],[165,121,30],[44,111,31],[120,98,17],[330,143,30],[63,119,47],[130,108,20],[600,124,24],[156,176,27],[140,112,50],[115,82,22],[230,123,45],[185,188,14],[25,89,19],[120,109,18],[126,150,29],[293,181,42],[41,92,25],[272,152,39],[182,111,13],[158,106,21],[194,174,22],[321,168,42],[144,138,26],[15,68,13],[160,112,42],[115,94,27],[54,90,47],[90,102,40],[183,128,17],[66,94,18],[91,97,32],[46,100,12],[105,102,17],[152,103,30],[440,157,35],[144,167,17],[159,179,36],[130,136,35],[100,91,25],[106,117,23],[77,123,40],[135,106,28],[540,155,27],[90,101,35],[200,120,48],[70,80,31],[231,167,46],[130,145,46],[132,112,45],[190,98,33],[100,154,30],[168,165,26],[49,68,23],[240,123,35],[265,101,17],[45,56,28],[105,95,39],[205,129,26],[180,140,26],[180,144,46],[95,121,32],[125,129,49],[480,142,24],[125,169,19],[155,127,11],[200,122,27],[100,110,20],[335,127,21],[160,93,32],[387,158,13],[22,126,27],[291,134,20],[392,187,33],[185,173,39],[178,108,46],[200,114,36],[127,149,29],[105,117,30],[180,116,29],[79,130,23],[120,174,37],[165,106,27],[120,126,27],[160,99,17],[150,120,37],[94,102,20],[116,109,18],[140,153,37],[105,100,33],[57,81,41],[200,187,22],[74,121,39],[510,181,44],[110,128,39],[16,88,26],[180,101,48],[112,121,23]],"colors":[[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,0,0,1],[0,0,0,1],[1,0,0,1],[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1]],"centers":[[94,89,23],[168,137,35],[88,78,32],[543,197,45],[846,189,23],[175,166,19],[230,118,47],[83,103,38],[96,115,30],[235,126,41],[146,143,33],[115,125,26],[140,97,15],[110,145,19],[245,158,36],[54,88,11],[192,103,33],[207,111,47],[70,180,25],[240,171,24],[82,103,11],[36,101,15],[23,88,21],[300,176,34],[342,150,42],[304,187,39],[110,100,60],[142,105,41],[128,141,34],[38,95,13],[100,146,27],[90,100,20],[140,139,35],[270,129,20],[71,83,26],[125,110,29],[71,100,25],[110,136,32],[176,123,15],[48,81,40],[64,142,18],[228,144,27],[76,71,18],[64,93,30],[220,122,51],[40,81,18],[152,126,29],[140,144,28],[18,83,31],[36,95,25],[135,171,33],[495,155,26],[37,89,34],[175,160,32],[51,99,15],[100,162,56],[100,107,30],[99,88,42],[135,120,30],[94,118,36],[145,117,24],[168,173,14],[225,170,37],[49,96,13],[140,125,20],[50,100,26],[92,93,25],[325,105,29],[63,108,26],[284,154,31],[119,106,35],[204,136,50],[155,156,28],[485,153,42],[94,99,15],[135,109,21],[53,88,19],[114,163,41],[105,102,40],[285,114,34],[156,104,18],[78,111,12],[130,134,23],[48,79,42],[55,75,24],[130,179,42],[130,129,46],[92,119,18],[495,181,36],[58,128,41],[114,109,39],[160,139,35],[94,123,44],[210,158,41],[48,107,13],[99,109,44],[318,148,27],[44,99,16],[190,103,32],[280,196,29],[87,96,27],[130,140,26],[175,112,32],[271,151,40],[129,109,41],[120,125,30],[478,177,29],[190,142,33],[56,100,15],[32,87,27],[744,197,39],[53,117,31],[370,134,37],[37,79,25],[45,74,28],[192,181,21],[88,91,32],[176,119,22],[194,146,35],[680,165,33],[402,124,33],[55,90,14],[258,92,7],[375,193,16],[150,155,28],[130,191,15],[67,96,18],[56,108,32],[45,71,50],[57,100,52],[116,104,23],[278,108,10],[122,129,28],[155,133,15],[135,136,26],[545,155,44],[220,119,39],[49,96,17],[75,108,43],[40,78,29],[74,107,30],[182,128,37],[194,128,45],[120,151,31],[360,146,38],[215,126,29],[184,100,25],[135,144,33],[42,77,41],[105,120,37],[132,161,23],[148,137,14],[180,128,19],[205,124,28],[148,106,37],[96,155,17],[85,113,10],[94,112,22],[64,99,11],[140,115,39],[231,129,12],[29,152,33],[168,157,21],[156,122,32],[120,102,36],[68,105,32],[52,87,16],[58,95,18],[255,165,43],[171,152,34],[105,130,13],[73,95,21],[108,126,36],[83,139,19],[74,99,19],[43,90,12],[167,125,40],[54,88,40],[249,196,36],[325,189,33],[293,147,25],[83,99,28],[66,81,16],[140,133,28],[465,173,48],[66,84,22],[94,105,40],[158,122,43],[325,140,43],[84,98,15],[75,87,37],[72,93,39],[82,107,30],[182,109,8],[59,90,18],[110,125,24],[50,119,13],[285,144,26],[81,100,23],[196,100,29],[415,131,14],[87,116,12],[275,127,24],[115,96,34],[88,136,41],[165,123,32],[579,172,49],[176,112,30],[310,143,23],[61,143,22],[167,138,35],[474,173,33],[115,129,29],[170,119,41],[76,94,18],[78,102,46],[210,151,32],[277,184,39],[180,181,30],[145,135,46],[180,95,25],[85,89,16],[60,80,11],[50,83,23],[120,117,27],[14,180,63],[70,100,12],[92,95,45],[64,104,37],[63,120,18],[95,82,13],[210,91,32],[105,100,28],[71,86,28],[237,148,48],[60,134,33],[56,120,22],[49,74,40],[105,124,13],[36,74,10],[100,97,36],[140,154,41],[191,105,45],[110,114,17],[75,126,38],[328,158,30],[49,85,22],[125,84,31],[250,135,42],[480,139,41],[265,173,32],[66,83,28],[122,125,18],[76,81,15],[145,195,33],[193,154,32],[71,117,19],[79,94,25],[90,180,26],[170,130,23],[76,84,23],[210,139,17],[86,99,19],[105,163,18],[165,145,34],[326,129,7],[66,68,32],[130,124,33],[82,97,19],[105,116,15],[188,117,31],[106,122,18],[65,86,52],[56,77,30],[210,127,37],[155,129,49],[215,100,40],[190,128,25],[56,84,23],[76,88,29],[225,186,35],[207,187,27],[166,131,21],[67,164,43],[106,84,30],[44,88,24],[115,84,23],[215,124,33],[274,198,32],[77,87,34],[54,99,19],[88,95,14],[18,99,30],[126,92,32],[126,154,29],[165,121,30],[44,111,31],[120,98,17],[330,143,30],[63,119,47],[130,108,20],[600,124,24],[156,176,27],[140,112,50],[115,82,22],[230,123,45],[185,188,14],[25,89,19],[120,109,18],[126,150,29],[293,181,42],[41,92,25],[272,152,39],[182,111,13],[158,106,21],[194,174,22],[321,168,42],[144,138,26],[15,68,13],[160,112,42],[115,94,27],[54,90,47],[90,102,40],[183,128,17],[66,94,18],[91,97,32],[46,100,12],[105,102,17],[152,103,30],[440,157,35],[144,167,17],[159,179,36],[130,136,35],[100,91,25],[106,117,23],[77,123,40],[135,106,28],[540,155,27],[90,101,35],[200,120,48],[70,80,31],[231,167,46],[130,145,46],[132,112,45],[190,98,33],[100,154,30],[168,165,26],[49,68,23],[240,123,35],[265,101,17],[45,56,28],[105,95,39],[205,129,26],[180,140,26],[180,144,46],[95,121,32],[125,129,49],[480,142,24],[125,169,19],[155,127,11],[200,122,27],[100,110,20],[335,127,21],[160,93,32],[387,158,13],[22,126,27],[291,134,20],[392,187,33],[185,173,39],[178,108,46],[200,114,36],[127,149,29],[105,117,30],[180,116,29],[79,130,23],[120,174,37],[165,106,27],[120,126,27],[160,99,17],[150,120,37],[94,102,20],[116,109,18],[140,153,37],[105,100,33],[57,81,41],[200,187,22],[74,121,39],[510,181,44],[110,128,39],[16,88,26],[180,101,48],[112,121,23]],"ignoreExtent":false,"flags":4096},"9":{"id":9,"type":"text","material":{"lit":false},"vertices":[[430,31.9309997558594,-2.49200010299683]],"colors":[[0,0,0,1]],"texts":[["pima$insulin"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[430,31.9309997558594,-2.49200010299683]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"10":{"id":10,"type":"text","material":{"lit":false},"vertices":[[-127.024002075195,127,-2.49200010299683]],"colors":[[0,0,0,1]],"texts":[["pima$glucose"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-127.024002075195,127,-2.49200010299683]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"11":{"id":11,"type":"text","material":{"lit":false},"vertices":[[-127.024002075195,31.9309997558594,35]],"colors":[[0,0,0,1]],"texts":[["pima$triceps"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-127.024002075195,31.9309997558594,35]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"5":{"id":5,"type":"light","vertices":[[0,0,1]],"colors":[[1,1,1,1],[1,1,1,1],[1,1,1,1]],"viewpoint":true,"finite":false},"4":{"id":4,"type":"background","material":{"fog":true},"colors":[[0.298039227724075,0.298039227724075,0.298039227724075,1]],"centers":[[0,0,0]],"sphere":false,"fogtype":"none","flags":0},"6":{"id":6,"type":"background","material":{"lit":false,"back":"lines"},"colors":[[1,1,1,1]],"centers":[[0,0,0]],"sphere":false,"fogtype":"none","flags":0},"8":{"id":8,"type":"bboxdeco","material":{"front":"lines","back":"lines"},"vertices":[[200,"NA","NA"],[400,"NA","NA"],[600,"NA","NA"],[800,"NA","NA"],["NA",100,"NA"],["NA",150,"NA"],["NA","NA",10],["NA","NA",20],["NA","NA",30],["NA","NA",40],["NA","NA",50],["NA","NA",60]],"colors":[[0,0,0,1]],"draw_front":true,"newIds":[19,20,21,22,23,24,25]},"1":{"id":1,"type":"subscene","par3d":{"antialias":0,"FOV":30,"ignoreExtent":false,"listeners":1,"mouseMode":{"left":"trackball","right":"zoom","middle":"fov","wheel":"pull"},"observer":[0,0,2069.89453125],"modelMatrix":[[0.586986541748047,0,0,-252.404205322266],[0,1.1762912273407,8.19500637054443,-436.214233398438],[0,-3.23183345794678,2.98273849487305,-1763.84753417969],[0,0,0,1]],"projMatrix":[[3.73205089569092,0,0,0],[0,3.73205089569092,0,0],[0,0,-3.86370348930359,-7461.73046875],[0,0,-1,0]],"skipRedraw":false,"userMatrix":[[1,0,0,0],[0,0.342020143325668,0.939692620785909,0],[0,-0.939692620785909,0.342020143325668,0],[0,0,0,1]],"userProjection":[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],"scale":[0.586986541748047,3.43924522399902,8.72094345092773],"viewport":{"x":0,"y":0,"width":1,"height":1},"zoom":1,"bbox":[14,846,56,198,7,63],"windowRect":[463,385,719,641],"family":"sans","font":1,"cex":1,"useFreeType":true,"fontname":"/usr/local/lib/R/site-library/rgl/fonts/FreeSans.ttf","maxClipPlanes":6,"glVersion":2.09999990463257,"activeSubscene":0},"embeddings":{"viewport":"replace","projection":"replace","model":"replace","mouse":"replace"},"objects":[6,8,7,9,10,11,5,19,20,21,22,23,24,25],"subscenes":[],"flags":6736},"19":{"id":19,"type":"lines","material":{"lit":false},"vertices":[[200,53.8699989318848,6.15999984741211],[800,53.8699989318848,6.15999984741211],[200,53.8699989318848,6.15999984741211],[200,50.2135009765625,4.71799993515015],[400,53.8699989318848,6.15999984741211],[400,50.2135009765625,4.71799993515015],[600,53.8699989318848,6.15999984741211],[600,50.2135009765625,4.71799993515015],[800,53.8699989318848,6.15999984741211],[800,50.2135009765625,4.71799993515015]],"colors":[[0,0,0,1]],"centers":[[500,53.8699989318848,6.15999984741211],[200,52.041748046875,5.43900012969971],[400,52.041748046875,5.43900012969971],[600,52.041748046875,5.43900012969971],[800,52.041748046875,5.43900012969971]],"ignoreExtent":true,"origId":8,"flags":64},"20":{"id":20,"type":"text","material":{"lit":false},"vertices":[[200,42.9005012512207,1.83399999141693],[400,42.9005012512207,1.83399999141693],[600,42.9005012512207,1.83399999141693],[800,42.9005012512207,1.83399999141693]],"colors":[[0,0,0,1]],"texts":[["200"],["400"],["600"],["800"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[200,42.9005012512207,1.83399999141693],[400,42.9005012512207,1.83399999141693],[600,42.9005012512207,1.83399999141693],[800,42.9005012512207,1.83399999141693]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":8,"flags":2064},"21":{"id":21,"type":"lines","material":{"lit":false},"vertices":[[1.51999998092651,100,6.15999984741211],[1.51999998092651,150,6.15999984741211],[1.51999998092651,100,6.15999984741211],[-19.9039993286133,100,4.71799993515015],[1.51999998092651,150,6.15999984741211],[-19.9039993286133,150,4.71799993515015]],"colors":[[0,0,0,1]],"centers":[[1.51999998092651,125,6.15999984741211],[-9.1919994354248,100,5.43900012969971],[-9.1919994354248,150,5.43900012969971]],"ignoreExtent":true,"origId":8,"flags":64},"22":{"id":22,"type":"text","material":{"lit":false},"vertices":[[-62.7519989013672,100,1.83399999141693],[-62.7519989013672,150,1.83399999141693]],"colors":[[0,0,0,1]],"texts":[["100"],["150"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-62.7519989013672,100,1.83399999141693],[-62.7519989013672,150,1.83399999141693]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":8,"flags":2064},"23":{"id":23,"type":"lines","material":{"lit":false},"vertices":[[1.51999998092651,53.8699989318848,10],[1.51999998092651,53.8699989318848,60],[1.51999998092651,53.8699989318848,10],[-19.9039993286133,50.2135009765625,10],[1.51999998092651,53.8699989318848,20],[-19.9039993286133,50.2135009765625,20],[1.51999998092651,53.8699989318848,30],[-19.9039993286133,50.2135009765625,30],[1.51999998092651,53.8699989318848,40],[-19.9039993286133,50.2135009765625,40],[1.51999998092651,53.8699989318848,50],[-19.9039993286133,50.2135009765625,50],[1.51999998092651,53.8699989318848,60],[-19.9039993286133,50.2135009765625,60]],"colors":[[0,0,0,1]],"centers":[[1.51999998092651,53.8699989318848,35],[-9.1919994354248,52.041748046875,10],[-9.1919994354248,52.041748046875,20],[-9.1919994354248,52.041748046875,30],[-9.1919994354248,52.041748046875,40],[-9.1919994354248,52.041748046875,50],[-9.1919994354248,52.041748046875,60]],"ignoreExtent":true,"origId":8,"flags":64},"24":{"id":24,"type":"text","material":{"lit":false},"vertices":[[-62.7519989013672,42.9005012512207,10],[-62.7519989013672,42.9005012512207,20],[-62.7519989013672,42.9005012512207,30],[-62.7519989013672,42.9005012512207,40],[-62.7519989013672,42.9005012512207,50],[-62.7519989013672,42.9005012512207,60]],"colors":[[0,0,0,1]],"texts":[["10"],["20"],["30"],["40"],["50"],["60"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-62.7519989013672,42.9005012512207,10],[-62.7519989013672,42.9005012512207,20],[-62.7519989013672,42.9005012512207,30],[-62.7519989013672,42.9005012512207,40],[-62.7519989013672,42.9005012512207,50],[-62.7519989013672,42.9005012512207,60]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":8,"flags":2064},"25":{"id":25,"type":"lines","material":{"lit":false},"vertices":[[1.51999998092651,53.8699989318848,6.15999984741211],[1.51999998092651,200.130004882812,6.15999984741211],[1.51999998092651,53.8699989318848,63.8400001525879],[1.51999998092651,200.130004882812,63.8400001525879],[1.51999998092651,53.8699989318848,6.15999984741211],[1.51999998092651,53.8699989318848,63.8400001525879],[1.51999998092651,200.130004882812,6.15999984741211],[1.51999998092651,200.130004882812,63.8400001525879],[1.51999998092651,53.8699989318848,6.15999984741211],[858.47998046875,53.8699989318848,6.15999984741211],[1.51999998092651,53.8699989318848,63.8400001525879],[858.47998046875,53.8699989318848,63.8400001525879],[1.51999998092651,200.130004882812,6.15999984741211],[858.47998046875,200.130004882812,6.15999984741211],[1.51999998092651,200.130004882812,63.8400001525879],[858.47998046875,200.130004882812,63.8400001525879],[858.47998046875,53.8699989318848,6.15999984741211],[858.47998046875,200.130004882812,6.15999984741211],[858.47998046875,53.8699989318848,63.8400001525879],[858.47998046875,200.130004882812,63.8400001525879],[858.47998046875,53.8699989318848,6.15999984741211],[858.47998046875,53.8699989318848,63.8400001525879],[858.47998046875,200.130004882812,6.15999984741211],[858.47998046875,200.130004882812,63.8400001525879]],"colors":[[0,0,0,1]],"centers":[[1.51999998092651,127,6.15999984741211],[1.51999998092651,127,63.8400001525879],[1.51999998092651,53.8699989318848,35],[1.51999998092651,200.130004882812,35],[430,53.8699989318848,6.15999984741211],[430,53.8699989318848,63.8400001525879],[430,200.130004882812,6.15999984741211],[430,200.130004882812,63.8400001525879],[858.47998046875,127,6.15999984741211],[858.47998046875,127,63.8400001525879],[858.47998046875,53.8699989318848,35],[858.47998046875,200.130004882812,35]],"ignoreExtent":true,"origId":8,"flags":64},"32":{"id":32,"type":"linestrip","material":{"alpha":0.498039215803146,"lit":false,"depth_test":"always","isTransparent":true},"vertices":[[0,0,-0.999000012874603],[1,0,-0.999000012874603],[1,1,-0.999000012874603],[0,1,-0.999000012874603],[0,0,-0.999000012874603]],"colors":[[0,0,0,0.498039215803146]],"centers":[[0,0,-0.999000012874603],[1,0,-0.999000012874603],[1,1,-0.999000012874603],[0,1,-0.999000012874603],[0,0,-0.999000012874603]],"ignoreExtent":false,"flags":32864}},"crosstalk":{"key":[null],"group":"","id":0,"options":[null]},"width":800,"height":800,"brushId":32,"sphereVerts":{"vb":[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0746578340503426,0.146446609406726,0.212607523691814,0.270598050073099,0.318189645143208,0.353553390593274,0.375330277517865,0.38268343236509,0.375330277517865,0.353553390593274,0.318189645143208,0.270598050073099,0.212607523691814,0.146446609406726,0.0746578340503426,0,0,0.137949689641472,0.270598050073099,0.392847479193551,0.5,0.587937801209679,0.653281482438188,0.693519922661074,0.707106781186547,0.693519922661074,0.653281482438188,0.587937801209679,0.5,0.392847479193551,0.270598050073099,0.137949689641472,0,0,0.180239955501737,0.353553390593274,0.513279967159337,0.653281482438188,0.768177756711416,0.853553390593274,0.906127446352888,0.923879532511287,0.906127446352888,0.853553390593274,0.768177756711416,0.653281482438188,0.513279967159337,0.353553390593274,0.180239955501737,0,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186548,0.831469612302545,0.923879532511287,0.98078528040323,1,0.98078528040323,0.923879532511287,0.831469612302545,0.707106781186548,0.555570233019602,0.38268343236509,0.195090322016128,0,0,0.180239955501737,0.353553390593274,0.513279967159337,0.653281482438188,0.768177756711416,0.853553390593274,0.906127446352888,0.923879532511287,0.906127446352888,0.853553390593274,0.768177756711416,0.653281482438188,0.513279967159337,0.353553390593274,0.180239955501737,0,0,0.137949689641472,0.270598050073099,0.392847479193551,0.5,0.587937801209679,0.653281482438188,0.693519922661074,0.707106781186548,0.693519922661074,0.653281482438188,0.587937801209679,0.5,0.392847479193551,0.270598050073099,0.137949689641472,0,0,0.0746578340503427,0.146446609406726,0.212607523691814,0.270598050073099,0.318189645143209,0.353553390593274,0.375330277517865,0.38268343236509,0.375330277517865,0.353553390593274,0.318189645143209,0.270598050073099,0.212607523691814,0.146446609406726,0.0746578340503427,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0.0746578340503427,-0.146446609406726,-0.212607523691814,-0.270598050073099,-0.318189645143209,-0.353553390593274,-0.375330277517865,-0.38268343236509,-0.375330277517865,-0.353553390593274,-0.318189645143209,-0.270598050073099,-0.212607523691814,-0.146446609406726,-0.0746578340503427,-0,-0,-0.137949689641472,-0.270598050073099,-0.392847479193551,-0.5,-0.587937801209679,-0.653281482438188,-0.693519922661074,-0.707106781186548,-0.693519922661074,-0.653281482438188,-0.587937801209679,-0.5,-0.392847479193551,-0.270598050073099,-0.137949689641472,-0,-0,-0.180239955501737,-0.353553390593274,-0.513279967159337,-0.653281482438188,-0.768177756711416,-0.853553390593274,-0.906127446352888,-0.923879532511287,-0.906127446352888,-0.853553390593274,-0.768177756711416,-0.653281482438188,-0.513279967159337,-0.353553390593274,-0.180239955501737,-0,-0,-0.195090322016128,-0.38268343236509,-0.555570233019602,-0.707106781186548,-0.831469612302545,-0.923879532511287,-0.98078528040323,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186548,-0.555570233019602,-0.38268343236509,-0.195090322016128,-0,-0,-0.180239955501737,-0.353553390593274,-0.513279967159337,-0.653281482438188,-0.768177756711416,-0.853553390593274,-0.906127446352888,-0.923879532511287,-0.906127446352888,-0.853553390593274,-0.768177756711416,-0.653281482438188,-0.513279967159337,-0.353553390593274,-0.180239955501737,-0,-0,-0.137949689641472,-0.270598050073099,-0.392847479193551,-0.5,-0.587937801209679,-0.653281482438188,-0.693519922661074,-0.707106781186547,-0.693519922661074,-0.653281482438188,-0.587937801209679,-0.5,-0.392847479193551,-0.270598050073099,-0.137949689641472,-0,-0,-0.0746578340503426,-0.146446609406726,-0.212607523691814,-0.270598050073099,-0.318189645143208,-0.353553390593274,-0.375330277517865,-0.38268343236509,-0.375330277517865,-0.353553390593274,-0.318189645143208,-0.270598050073099,-0.212607523691814,-0.146446609406726,-0.0746578340503426,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186547,-0.555570233019602,-0.38268343236509,-0.195090322016128,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186547,0.831469612302545,0.923879532511287,0.98078528040323,1],[0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186548,0.831469612302545,0.923879532511287,0.98078528040323,1,0.98078528040323,0.923879532511287,0.831469612302545,0.707106781186548,0.555570233019602,0.38268343236509,0.195090322016128,0,0,0.180239955501737,0.353553390593274,0.513279967159337,0.653281482438188,0.768177756711416,0.853553390593274,0.906127446352888,0.923879532511287,0.906127446352888,0.853553390593274,0.768177756711416,0.653281482438188,0.513279967159337,0.353553390593274,0.180239955501737,0,0,0.137949689641472,0.270598050073099,0.392847479193551,0.5,0.587937801209679,0.653281482438188,0.693519922661074,0.707106781186548,0.693519922661074,0.653281482438188,0.587937801209679,0.5,0.392847479193551,0.270598050073099,0.137949689641472,0,0,0.0746578340503426,0.146446609406726,0.212607523691814,0.270598050073099,0.318189645143209,0.353553390593274,0.375330277517865,0.38268343236509,0.375330277517865,0.353553390593274,0.318189645143209,0.270598050073099,0.212607523691814,0.146446609406726,0.0746578340503426,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0.0746578340503426,-0.146446609406726,-0.212607523691814,-0.270598050073099,-0.318189645143208,-0.353553390593274,-0.375330277517865,-0.38268343236509,-0.375330277517865,-0.353553390593274,-0.318189645143208,-0.270598050073099,-0.212607523691814,-0.146446609406726,-0.0746578340503426,-0,-0,-0.137949689641472,-0.270598050073099,-0.392847479193551,-0.5,-0.587937801209679,-0.653281482438188,-0.693519922661074,-0.707106781186547,-0.693519922661074,-0.653281482438188,-0.587937801209679,-0.5,-0.392847479193551,-0.270598050073099,-0.137949689641472,-0,-0,-0.180239955501737,-0.353553390593274,-0.513279967159337,-0.653281482438188,-0.768177756711416,-0.853553390593274,-0.906127446352888,-0.923879532511287,-0.906127446352888,-0.853553390593274,-0.768177756711416,-0.653281482438188,-0.513279967159337,-0.353553390593274,-0.180239955501737,-0,-0,-0.195090322016128,-0.38268343236509,-0.555570233019602,-0.707106781186548,-0.831469612302545,-0.923879532511287,-0.98078528040323,-1,-0.98078528040323,-0.923879532511287,-0.831469612302545,-0.707106781186548,-0.555570233019602,-0.38268343236509,-0.195090322016128,-0,-0,-0.180239955501737,-0.353553390593274,-0.513279967159337,-0.653281482438188,-0.768177756711416,-0.853553390593274,-0.906127446352888,-0.923879532511287,-0.906127446352888,-0.853553390593274,-0.768177756711416,-0.653281482438188,-0.513279967159337,-0.353553390593274,-0.180239955501737,-0,-0,-0.137949689641472,-0.270598050073099,-0.392847479193551,-0.5,-0.58793780120968,-0.653281482438188,-0.693519922661074,-0.707106781186548,-0.693519922661074,-0.653281482438188,-0.58793780120968,-0.5,-0.392847479193551,-0.270598050073099,-0.137949689641472,-0,-0,-0.0746578340503427,-0.146446609406726,-0.212607523691814,-0.270598050073099,-0.318189645143209,-0.353553390593274,-0.375330277517866,-0.38268343236509,-0.375330277517866,-0.353553390593274,-0.318189645143209,-0.270598050073099,-0.212607523691814,-0.146446609406726,-0.0746578340503427,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0746578340503427,0.146446609406726,0.212607523691814,0.270598050073099,0.318189645143209,0.353553390593274,0.375330277517865,0.38268343236509,0.375330277517865,0.353553390593274,0.318189645143209,0.270598050073099,0.212607523691814,0.146446609406726,0.0746578340503427,0,0,0.137949689641472,0.270598050073098,0.392847479193551,0.5,0.587937801209679,0.653281482438188,0.693519922661074,0.707106781186547,0.693519922661074,0.653281482438188,0.587937801209679,0.5,0.392847479193551,0.270598050073098,0.137949689641472,0,0,0.180239955501737,0.353553390593274,0.513279967159337,0.653281482438188,0.768177756711416,0.853553390593274,0.906127446352888,0.923879532511287,0.906127446352888,0.853553390593274,0.768177756711416,0.653281482438188,0.513279967159337,0.353553390593274,0.180239955501737,0,0,0.195090322016128,0.38268343236509,0.555570233019602,0.707106781186548,0.831469612302545,0.923879532511287,0.98078528040323,1,0.98078528040323,0.923879532511287,0.831469612302545,0.707106781186548,0.555570233019602,0.38268343236509,0.195090322016128,0]],"it":[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270],[17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288],[18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271]],"material":[],"normals":null,"texcoords":[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1]]},"context":{"shiny":false,"rmarkdown":null},"players":[],"webGLoptions":{"preserveDrawingBuffer":true}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->


#### Plus de trois dimensions

Déjà à trois dimensions la visualisation devient délicate, mais au delà, cela devient pratiquement mission impossible. La **matrice de nuages de points** peut rendre service ici, mais dans certaines limites (tous les angles de vue ne sont pas accessibles).


```r
GGally::ggscatmat(pima, 2:6, color = "diabetes")
```

<img src="07-acp-afc_files/figure-html/unnamed-chunk-20-1.png" width="672" style="display: block; margin: auto;" />

Nous voyons qu'ici nous atteignons les limites des possibilités. C'est pour cela que, pour des données multivariées comportant beaucoup de variables quantitatives, les techniques de réduction des dimensions comme l'ACP sont indispensables.


### ACP : mécanisme

Nous allons partir d'un exemple presque trivial pour illuster le principe de l'ACP. Comment réduire un tableau bivarié en une représentation des individus en une seule dimension (classement sur une droite) avec perte minimale d'information\ ? Par exemple, en partant de ces données fictives\ :

<img src="07-acp-afc_files/figure-html/unnamed-chunk-21-1.png" width="864" style="display: block; margin: auto;" />

Voic une représentation graphique 2D de ces données\ :

<img src="07-acp-afc_files/figure-html/unnamed-chunk-22-1.png" width="864" style="display: block; margin: auto;" />

Si nous réduisons à une seule dimension en laissant tomber une des deux variables, voici ce que cela donne (ici on ne garde que `Var1`, donc, on projette les points sur l'axe des abscisses).


<img src="07-acp-afc_files/figure-html/unnamed-chunk-23-1.png" width="432" style="display: block; margin: auto;" />

Au final, nous avons ordonné nos individus en une dimension comme suit\ :

<img src="07-acp-afc_files/figure-html/unnamed-chunk-24-1.png" width="432" style="display: block; margin: auto;" />

C'est une mauvaise solution car il y a trop de perte d'information. Regardez l'écart entre 7 et 9 sur le graphqie en deux dimensions et dans celui à une dimension\ : les points sont trop près. Comparez sur les deux graphiques les distances 7 - 9 avec  9 - 8 et 1 - 2 *versus* 1 - 3. Tout cela est très mal représenté en une seule dimension.

Une autre solution serait de projeter le long de la droite de "tendance générale", c'est-à-dire le long de l'axe de plus grand allongement du nuage de points.


<img src="07-acp-afc_files/figure-html/unnamed-chunk-25-1.png" width="432" style="display: block; margin: auto;" />

Cela donne ceci en une seule dimension\ :

<img src="07-acp-afc_files/figure-html/unnamed-chunk-26-1.png" width="432" style="display: block; margin: auto;" />

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

<img src="07-acp-afc_files/figure-html/unnamed-chunk-27-1.png" width="432" style="display: block; margin: auto;" />


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
  if (span_y * aspect.ratio > span_x) {
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

