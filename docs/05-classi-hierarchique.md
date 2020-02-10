# (PART) Cours II: analyse {-}

# Classification hiérarchique {#hierarchique}




##### Objectifs {-}

- Comprendre la notion de distance et la matrice de distance.

- Appréhender la classification hiérarchique et le dendrogramme.

- Être capable d'effectuer un regroupement pertinent des individus d'un jeu de données multivarié à l'aide de ces techniques.


##### Prérequis {-}

Vous devez être à l'aise avec l'utilisation de R et Rstudio, en particulier pour l'importation, le remaniement et la visualisation de données multivariées. Ceci correspond au cours SDD I. Il n'est pas nécessaire d'avoir acquis toutes les notions vue dans la partie **Cours II :modélisation** pour pouvoir comprendre cette seconde partie du cours. Si vous ne vous sentez pas assez à l'aise avec R et RStudio, c'est peut-être le bon moment pour refaire le premier "learnr" du package `BioDataScience2`\ :

\BeginKnitrBlock{bdd}<div class="bdd">Démarrez la SciViews Box et RStudio. Dans la fenêtre **Console** de RStudio, entrez l'instruction suivante suivie de la touche `Entrée` pour ouvrir le tutoriel concernant les bases de R\ :

    BioDataScience2::run("01a_rappel")

N’oubliez pas d’appuyer sur la touche `ESC` pour reprendre la main dans R à la fin d’un tutoriel dans la console R.</div>\EndKnitrBlock{bdd}


## Analyse de données

L'analyse de données (on dit aussi *analyse exploratoire des données*, EAD ou **statistiques exploratoires**) mets en œuvre des méthodes statistiques multivariées visant à découvrir de l'information pertinente dans un gros jeu de données via des approches multidimensionnelles et essentiellement descriptives. Ces méthodes se regroupent en deux grandes familles\ :

- Celles visant à **réduire la dimensionnalité** (travailler avec des tableaux ayant moins de colonnes). Elles permettent ensuite de présenter les données de manière synthétique pour observer des relations entre les variables ou les individus via des **représentations graphiques**. Nous aborderons ces techniques dans les modules suivants.

- Celles cherchant à **classifier** (ou regrouper) les individus. Il s'agit ici de synthétiser le gros tableau de données dans l'autre sens, selon les lignes. L'approche via la *classification hiérarchique* sera détaillée ici.

La vidéo suivante introduit l'EAD (jusqu'à 2:11)\ :

<!--html_preserve--><iframe src="https://www.youtube.com/embed/q-IVQoh1nxA?end=131" width="770" height="433" frameborder="0" allowfullscreen=""></iframe><!--/html_preserve-->


## Distance entre individus

...
