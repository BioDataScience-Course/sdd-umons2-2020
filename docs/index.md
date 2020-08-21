--- 
title: "Science des données biologiques 2"
author: "Philippe Grosjean & Guyliann Engels (avec des contributions de Raphaël Conotte)"
date: "2020-08-12"
site: bookdown::bookdown_site
output:
  bookdown::gitbook:
    includes:
      after_body: disqus.html
documentclass: book
bibliography: [book.bib, packages.bib]
biblio-style: apalike
link-citations: yes
github-repo: biodatascience-course/sdd-umons2
url: 'https\://biodatascience-course.sciviews.org/sdd-umons2/'
description: "Tutoriel interactif pour la science des données avec R & SciViews-R, partie 2."
cover-image: "images/front-cover.png"
---

# Préambule {-}



Cet ouvrage interactif est le second d'une série de trois ouvrages traitant de la science des données biologiques. L'écriture de cette suite de livres a débuté au cours de l'année académique 2018-2019. 

Pour l'année académique 2019-2020, cet ouvrage interactif sera le support du cours suivant :

- [Science des données II : Analyse et modélisation, UMONS](http://applications.umons.ac.be/web/fr/pde/2019-2020/ue/US-B3-SCBIOL-006-M.htm) dont le responsable est Grosjean Philippe

Cet ouvrage est conçu pour être utilisé de manière interactive en ligne. En effet, nous y ajoutons des vidéos, des démonstrations interactives, et des exercices sous forme de questionnaires interactifs également. **Ces différents éléments ne sont, bien évidemment, utilisables qu'en ligne.**

![](images/front-cover.jpg)

----

_Le matériel dans cet ouvrage est distribué sous licence [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/deed.fr)._

----


## Vue générale des cours {-}

Le cours de  **Science des données II: analyse et modélisation**  est dispensé aux biologistes de troisième Bachelier en Faculté des Sciences de l'Université de Mons à partir de l'année académique 2019-2020.

La matière est divisée en huit modules de 6h chacun en présentiel. Il nécessitera environ un tiers de ce temps (voir plus, en fonction de votre rythme et de votre technique d'apprentissage) en travail à domicile. **Cette matière fait suite au premier cours dont le contenu est considéré comme assimilé (voir https://biodatascience-course.sciviews.org/sdd-umons/).**

<!-- A faire: un diagramme qui montre la relation entre ces différents modules, et avec les modules du cours 1 -->

La première moitié du cours est consacrée à la **modélisation**, un domaine particulièrement important de la science des données qui étend les concepts déjà vu au cours 1 d'analyse de variance et de corrélation entre deux variables. Ces quatre modules formeront aussi un socle sur lequel nous pourrons élaborer les techniques d'apprentissage machine (classification supervisée), et puis ensuite l'apprentissage profond à la base de l'intelligence artificielle qui seront abordées plus tard dans le cours 3. Cette partie est dense, mais *ultra* importante\ !

La seconde moitié s'intéressera à l'**exploration des données**, encore appelée **analyse des données** qui vise à découvrir des caractéristiques intéressantes dans des très gros jeux de données. Ces techniques sont d'autant plus utiles que les données volumineuses deviennent de plus en plus courantes en biologie.


## Matériel pédagogique {-}

Le matériel pédagogique, rassemblé dans ce syllabus interactif est aussi varié que possible. Vous pourrez ainsi piocher dans l'offre en fonction de vos envies et de votre profil d'apprenant pour optimiser votre travail. Vous trouverez:

- le présent ouvrage en ligne,

- des tutoriaux interactifs (réalisés avec un logiciel appelé `learnr`). Vous pourrez exécuter ces tutoriaux directement sur votre ordinateur, et vous aurez alors accès à des pages Web réactives contenant des explications, des exercices et des quizzs en ligne,

- des slides de présentations,

- des dépôts Github Classroom dans la section `BioDataScience-Course` pour réaliser et documenter vos travaux personnels.

- des renvois vers des documents externes en ligne, types vidéos youtube ou vimeo, des ouvrages en ligne en anglais ou en français, des blogs, des tutoriaux, des parties gratuites de cours Datacamp ou équivalents, des questions sur des sites comme "Stackoverflow" ou issues des "mailing lists" R, ...

<div class="info">
<p>Tout ce matériel est accessible à partir du <a href="http://biodatascience-course.sciviews.org">site Web du cours</a>, du présent syllabus interactif (et de Moodle pour les étudiants de l’UMONS). Ces derniers ont aussi accès au dossier <code>SDD</code> sur <code>StudentTemp</code> en Intranet à l’UMONS. Les aspects pratiques seront à réaliser en utilisant la <strong>‘SciViews Box’</strong>, une machine virtuelle préconfigurée. Nous installerons ensemble la nouvelle version de cette SciViews Box au premier cours. Il est donc très important que vous soyez présent à ce cours, et vous pouvez venir aussi si vous le souhaitez avec votre propre ordinateur portable comme pour le cours 1. Enfin, vous pourrez poser vos questions par mail à l’adresse <code>sdd@sciviews.org</code>.</p>
</div>


##### System information {-}


```r
sessioninfo::session_info()
```

```
# ─ Session info ──────────────────────────────────────────────────────────
#  setting  value                       
#  version  R version 3.5.3 (2019-03-11)
#  os       Ubuntu 18.04.2 LTS          
#  system   x86_64, linux-gnu           
#  ui       X11                         
#  language (EN)                        
#  collate  en_US.UTF-8                 
#  ctype    en_US.UTF-8                 
#  tz       Europe/Brussels             
#  date     2020-08-12                  
# 
# ─ Packages ──────────────────────────────────────────────────────────────
#  package     * version date       lib source        
#  assertthat    0.2.1   2019-03-21 [2] CRAN (R 3.5.3)
#  bookdown      0.9     2018-12-21 [2] CRAN (R 3.5.3)
#  cli           1.1.0   2019-03-19 [2] CRAN (R 3.5.3)
#  colorspace    1.4-1   2019-03-18 [2] CRAN (R 3.5.3)
#  crayon        1.3.4   2017-09-16 [2] CRAN (R 3.5.3)
#  digest        0.6.18  2018-10-10 [2] CRAN (R 3.5.3)
#  dplyr         0.8.0.1 2019-02-15 [2] CRAN (R 3.5.3)
#  evaluate      0.13    2019-02-12 [2] CRAN (R 3.5.3)
#  farver        1.1.0   2018-11-20 [2] CRAN (R 3.5.3)
#  gganimate     1.0.3   2019-04-02 [2] CRAN (R 3.5.3)
#  ggplot2       3.1.1   2019-04-07 [2] CRAN (R 3.5.3)
#  glue          1.3.1   2019-03-12 [2] CRAN (R 3.5.3)
#  gtable        0.3.0   2019-03-25 [2] CRAN (R 3.5.3)
#  hms           0.4.2   2018-03-10 [2] CRAN (R 3.5.3)
#  htmltools     0.3.6   2017-04-28 [2] CRAN (R 3.5.3)
#  inline        0.3.15  2018-05-18 [2] CRAN (R 3.5.3)
#  knitr         1.22    2019-03-08 [2] CRAN (R 3.5.3)
#  lazyeval      0.2.2   2019-03-15 [2] CRAN (R 3.5.3)
#  magick        2.0     2018-10-05 [2] CRAN (R 3.5.3)
#  magrittr      1.5     2014-11-22 [2] CRAN (R 3.5.3)
#  munsell       0.5.0   2018-06-12 [2] CRAN (R 3.5.3)
#  pillar        1.3.1   2018-12-15 [2] CRAN (R 3.5.3)
#  pkgconfig     2.0.2   2018-08-16 [2] CRAN (R 3.5.3)
#  plyr          1.8.4   2016-06-08 [2] CRAN (R 3.5.3)
#  prettyunits   1.0.2   2015-07-13 [2] CRAN (R 3.5.3)
#  progress      1.2.0   2018-06-14 [2] CRAN (R 3.5.3)
#  purrr         0.3.2   2019-03-15 [2] CRAN (R 3.5.3)
#  R6            2.4.0   2019-02-14 [2] CRAN (R 3.5.3)
#  Rcpp          1.0.1   2019-03-17 [2] CRAN (R 3.5.3)
#  rlang         0.3.4   2019-04-07 [2] CRAN (R 3.5.3)
#  rmarkdown     1.12    2019-03-14 [2] CRAN (R 3.5.3)
#  rstudioapi    0.10    2019-03-19 [2] CRAN (R 3.5.3)
#  scales        1.0.0   2018-08-09 [2] CRAN (R 3.5.3)
#  sessioninfo   1.1.1   2018-11-05 [2] CRAN (R 3.5.3)
#  stringi       1.4.3   2019-03-12 [2] CRAN (R 3.5.3)
#  stringr       1.4.0   2019-02-10 [2] CRAN (R 3.5.3)
#  tibble        2.1.1   2019-03-16 [2] CRAN (R 3.5.3)
#  tidyselect    0.2.5   2018-10-11 [2] CRAN (R 3.5.3)
#  tweenr        1.0.1   2018-12-14 [2] CRAN (R 3.5.3)
#  withr         2.1.2   2018-03-15 [2] CRAN (R 3.5.3)
#  xfun          0.6     2019-04-02 [2] CRAN (R 3.5.3)
#  yaml          2.2.0   2018-07-25 [2] CRAN (R 3.5.3)
# 
# [1] /home/sv/R/x86_64-pc-linux-gnu-library/3.5
# [2] /usr/local/lib/R/site-library
# [3] /usr/lib/R/site-library
# [4] /usr/lib/R/library
```
