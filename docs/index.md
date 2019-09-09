--- 
title: "Science des données biologiques 2"
author: "Philippe Grosjean & Guyliann Engels"
date: "2019-09-09"
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

![](images/front-cover.png)

----

_Le matériel dans cet ouvrage est distribué sous licence [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/deed.fr)._

----


## Vue générale des cours {-}

Le cours de  **Science des données II: analyse et modélisation**  est dispensé aux biologistes de troisième Bachelier en Faculté des Sciences de l'Université de Mons à partir de l'année académique 2019-2020.

La matière est divisée en 8 modules de sessions de 6h chacuns en présentiel. Il nécessitera environ un tiers de ce temps (voir plus, en fonction de votre rythme et de votre technique d'apprentissage) en travail à domicile. **Cette matière fait suite au premier cours dont le contenu est considéré comme assimilé (voir https://biodatascience-course.sciviews.org/sdd-umons/).**

<!-- A faire: un diagramme qui montre la relation entre ces différents modules, et avec les modules du cours 1 -->


## Matériel pédagogique {-}

Le matériel pédagogique, rassemblé dans ce syllabus interactif est aussi varié que possible. Vous pourrez ainsi piocher dans l'offre en fonction de vos envies et de votre profil d'apprenant pour optimiser votre travail. Vous trouverez:

- le présent ouvrage en ligne,

- des tutoriaux interactifs (réalisés avec un logiciel appelé `learnr`). Vous pourrez exécuter ces tutoriaux directement sur votre ordinateur, et vous aurez alors accès à des pages Web réactives contenant des explications, des exercices et des quizzs en ligne,

- des slides de présentations,

- des dépôts Github Classroom dans la section `BioDataScience-Course` pour réaliser et documenter vos travaux personnels.

- des renvois vers des documents externes en ligne, types vidéos youtube ou vimeo, des ouvrages en ligne en anglais ou en français, des blogs, des tutoriaux, des parties gratuites de cours Datacamp ou équivalents, des questions sur des sites comme "Stackoverflow" ou issues des "mailing lists" R, ...

<div class="info">
<p>Tout ce matériel est accessible à partir du <a href="http://biodatascience-course.sciviews.org">site Web du cours</a>, du présent syllabus interactif (et de Moodle pour les étudiants de l’UMONS). Ces derniers ont aussi accès au dossier <code>SDD</code> sur <code>StudentTemp</code> en Intranet à l’UMONS. Les aspects pratiques seront à réaliser en utilisant la <strong>‘SciViews Box’</strong>, une machine virtuelle préconfigurée<a href="#fn1" class="footnote-ref" id="fnref1"><sup>1</sup></a>. Enfin, vous pourrez poser vos questions par mail à l’adresse <code>sdd@sciviews.org</code>.</p>
<div class="footnotes">
<hr />
<ol>
<li id="fn1"><p>Nous installerons ensemble la nouvelle version de cette SciViews Box au premier cours. Il est donc très important que vous soyez présent à ce cours, et vous pouvez venir aussi si vous le souhaitez avec votre propre ordinateur portable comme pour le cours 1.<a href="#fnref1" class="footnote-back">↩</a></p></li>
</ol>
</div>
</div>


## Comment apprendre? {-}


```r
fortunes::fortune("brain surgery")
```

```
# 
# I wish to perform brain surgery this afternoon at 4pm and don't know where
# to start. My background is the history of great statistician sports
# legends but I am willing to learn. I know there are courses and numerous
# books on brain surgery but I don't have the time for those. Please direct
# me to the appropriate HowTos, and be on standby for solving any problem I
# may encounter while in the operating room. Some of you might ask for
# specifics of the case, but that would require my following the posting
# guide and spending even more time than I am already taking to write this
# note.
#    -- I. Ben Fooled (aka Frank Harrell)
#       R-help (April 1, 2005)
```

Version courte: **en pratiquant, en faisant des erreurs !**

Version longue: aujourd'hui --et encore plus à l'avenir-- les données sont complexes et ne se manipulent plus simplement avec un tableur comme Microsoft Excel. Vous allez apprendre à maitriser des outils professionnels, ce qui sous-entend qu'ils sont très puissants mais aussi relativement complexes. La méthode d'apprentissage que nous vous proposons a pour objectif prioritaire de vous faciliter la tâche, quelles que soient vos aptitudes au départ. Envisagez votre voyage en science des données comme l'apprentissage d'une nouvelle langue. **C'est en pratiquant, et en pratiquant encore sur le long terme que vous allez progresser.** La formation s'étale sur quatre années, et est répartie en cinq cours de difficulté croissante pour vous aider dans cet apprentissage progressif et sur la durée. N'hésitez pas à expérimenter, tester, essayer des nouvelles idées (même au delà de ce qui sera demandé dans les exercices) et **n'ayez pas peur de faire des erreurs**. Vous en ferez, ... beaucoup ... _nous vous le souhaitons!_ En fait, la meilleure manière d'apprendre, c'est justement en faisant des erreurs, et puis en mettant tout en oeuvre pour les comprendre et les corriger. Donc, si un message d'erreur, ou un "warning" apparait, ne soyez pas intimidé. Prenez une bonne respiration, lisez-le attentivement, essayez de le comprendre, et au besoin faites-vous aider: la solution est sur le Net, 'Google^[Il existe tout de même des outils plus pointus pour obtenir de l'aide sur le logiciel R comme [rseek.org](https://rseek.org), [rdocumentation.org](https://www.rdocumentation.org) ou [rdrr.io](https://rdrr.io). Rien ne sert de chercher 'R' dans Goggle.] est votre ami'!


## Evaluation {-}

L'évaluation sera basée sur une somme de petites contributions qui matérialiseront votre progression sur le long terme. Avec cette évaluation, nous souhaitons vous gratifier chaque fois que vous franchirez des étapes, plutôt que de vous sanctionner lorsque vous bloquez. Donc, pour une note finale sur 20:

<!-- - 2 points pour la progression sur base des exercices que vous réaliserez en classe inversée (donc, chez vous). -->

- 3 points pour la restitution des capsules et votre participation en présentiel. Au début de chaque séance, nous discuterons des notions que vous aurez à préparer par avance, et votre participation sera évaluée.

- 6 points pour un quizz final.

<!-- : vous aurez à répondre à cinq questions au hasard (set différent pour chaque étudiant sur base de 20 questions au total). -->

- 11 points pour l’évaluation d’un des rapports d'analyse de données (choisi au hasard en fin de cours).

- Enfin, vous pourrez éventuellement encore gagner un point bonus pour une participation remarquable, ou tout autre élément à valoriser (rapport particulièrement bien réalisé, aide des autres étudiants, etc.). Ceci étant à l'appréciation des enseignants.



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
#  tz       Europe/Madrid               
#  date     2019-09-09                  
# 
# ─ Packages ──────────────────────────────────────────────────────────────
#  package     * version date       lib source        
#  assertthat    0.2.1   2019-03-21 [2] CRAN (R 3.5.3)
#  bookdown      0.9     2018-12-21 [2] CRAN (R 3.5.3)
#  cli           1.1.0   2019-03-19 [2] CRAN (R 3.5.3)
#  crayon        1.3.4   2017-09-16 [2] CRAN (R 3.5.3)
#  digest        0.6.18  2018-10-10 [2] CRAN (R 3.5.3)
#  evaluate      0.13    2019-02-12 [2] CRAN (R 3.5.3)
#  fortunes      1.5-4   2016-12-29 [2] CRAN (R 3.5.3)
#  htmltools     0.3.6   2017-04-28 [2] CRAN (R 3.5.3)
#  inline        0.3.15  2018-05-18 [2] CRAN (R 3.5.3)
#  knitr         1.22    2019-03-08 [2] CRAN (R 3.5.3)
#  magrittr      1.5     2014-11-22 [2] CRAN (R 3.5.3)
#  Rcpp          1.0.1   2019-03-17 [2] CRAN (R 3.5.3)
#  rmarkdown     1.12    2019-03-14 [2] CRAN (R 3.5.3)
#  sessioninfo   1.1.1   2018-11-05 [2] CRAN (R 3.5.3)
#  stringi       1.4.3   2019-03-12 [2] CRAN (R 3.5.3)
#  stringr       1.4.0   2019-02-10 [2] CRAN (R 3.5.3)
#  withr         2.1.2   2018-03-15 [2] CRAN (R 3.5.3)
#  xfun          0.6     2019-04-02 [2] CRAN (R 3.5.3)
#  yaml          2.2.0   2018-07-25 [2] CRAN (R 3.5.3)
# 
# [1] /home/sv/R/x86_64-pc-linux-gnu-library/3.5
# [2] /usr/local/lib/R/site-library
# [3] /usr/lib/R/site-library
# [4] /usr/lib/R/library
```
