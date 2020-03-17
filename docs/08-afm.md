# AFM {#afm}




##### Objectifs {-}

- TODO

##### Prérequis {-}

TODO

## Analyse factorielle multiple (AFM)

L'analyse factorielle multiple (AFM) se nomme *principal component analysis* (PCA) en anglais. 

## Open Data

### Gestion des données

Lors de la préparation d'une expérience, vous devez réfléchir à un plan d'expérience. Vous avez donc défini des protocoles d'expérience, le nombre de réplicas,... Vous devez cependant intégrer à votre réflexion, un plan de gestion de vos données. Dans ce plan, vous aurez à définir l'acquisition, la description ou encore le partage des données.

#### Principe FAIR

Pour assurer une gestion cohérente des données scientifiques, il faut respecter le plus possible l'acronyme en anglais **FAIR** : Findable Accessible Interoperable Reusable [@Wilkinson2016]. 

| Acronyme en anglais | Acronyme en français | Description  |
|---------------------|----------------------|-----------------------------------------------------------------------------------------------|
| Findable | Facile à trouver  | Les données ont besoin d'un code unique et persistant pour les retrouver |
| Accessible | Accessible | Les données et surtout les métadonnées avec une license sont mises à disposition. |
| Interoperable | Intéropérable | Les données et les métadonnées doivent suivre des respecter les standards internationaux |
| Reusable | Réutilisable | Les données doivent être réutilisable grâce à des métadonnées riches et des licenses claires  |

##### Facile à retrouver (*Findable*)

Vos données et vos métadonnées détaillées doivent être facile à retrouver. Vous devez donc fournir un identifiant **unique et permanent**. Il existe de nombreux identifiants comme ISBN, ISSN, DOI,...

Dans le cadre de la recherche scientifique, Le Digital Object Identifier (*DOI*) est la méthode standardisée conseillée. 

Vous avez déjà été confronté à des DOI. Par exemple, le DOI suivant <https://doi.org/10.1038/sdata.2016.18> fait référence à l'article *The FAIR Guiding Principles for scientific data management and stewardship*. Ce code est unique et persistant. Ce code va toujours renvoyer vers cet article de Scientific Data.

Imaginons que la revie Scientific Data disparaissent, le DOI sera toujours associé à cet article.

Le DOI ne couvre pas que les articles scientifiques, il est également utilisé pour les données. Par exemple, le DOI suivant <https://doi.org/10.5281/zenodo.3711592> fait référence au données intitulé. *Dataset: Number of diagnoses with coronavirus disease (COVID-19) in The Netherlands*. Nous reviendrons dans les sections suivantes sur [zenodo](https://zenodo.org). 

##### Accessible (*Accessible*)

Les données et les métadonnées que vous collectez doivent de plus en plus souvent être rendue disponible. Certaines revues scientifiques requierent la mise à disposition des données. Les recherches financées par des fonds publics (nationaux, Européen,...) requierent également la mise à disposition des données.

Les données ne doivent pas être disponible à tous. Par contre les métadonnées doivent l'être. Il est de plus important de préciser la procédure afin d'obtenir les données.

Il parait presque logique et évident de mettre à disposition ces données afin d'en faire profiter la recherche académique dans son ensemble. La recherche va progresser plus rapidement si les chercheurs collaborent. Un scientifique seul dans son laboratoire ne peut pas espérer progresser plus rapidement que 20 scientifiques qui collaborent et utilisent les données.

Il ne s'agit cependant pas de donner ces données sans aucune sécurité. En effet, il serait frustrant de travailler très dur sur un sujet précis et qu'un autre scientique "pique" le fruit de ces nombreuses heures de travail et publie un article avant vous.

Il existe une solution pour spécifier les droits d'utilisation de vos données. Vous devez associer une **license** à vos données et métadonnées. 

Vous avez très certainement déjà entendu parlé des licences [Creative Commons](https://creativecommons.org). Il est de plus en plus courant de voir apparaitre ce genre d'information sur des sites web comme [CC0](https://creativecommons.org/share-your-work/public-domain/cc0), [CC-by](https://creativecommons.org/licenses/by/4.0/), ou encore [CC-by-sa](https://creativecommons.org/licenses/by-sa/4.0/). 

Vous êtes peut être plus familié avec les logos ci-dessous 

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Licence Creative Commons" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a>

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Licence Creative Commons" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a>

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Licence Creative Commons" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a>

<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Licence Creative Commons" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a>.

**Que ce cache t'il derrière ces logos ? Nous allons détailler ensemble ce qui se cache derrière ces abréviations.** Nous pouvons résumer cela de manière simple en se posant deux questions.

- Souhaitez vous autoriser le partage des adaptations de votre Oeuvre ?
    + Oui
    + Non
    + Oui, sous condition de partage dans les mêmes conditions.

- Autorisez-vous les utilisations commerciales de votre œuvre ?
    + Oui
    + Non 

Ces deux questions proviennent de l'outil mis à disposition sur le site creative commons <https://creativecommons.org/choose/> pour définir la license la plus adpatée pour vous.

Dans le cadre de la recherche, tout n'est pas si simple. Vous devez tenir compte de votre chef de laboratoire, de la position de votre institution, de la position de l'institution qui vous finance. La bonne pratique est donc de discuter avec l'ensemble des acteurs pour décider de la bonne license à employer.

Repartons de notre jeu de données sur le COVID-19 : [Dataset: Number of diagnoses with coronavirus disease (COVID-19) in The Netherlands](https://doi.org/10.5281/zenodo.3711592). Nous pouvons voir que l'auteur a décidé d'employer la licence 

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Licence Creative Commons" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />Ce(tte) œuvre est mise à disposition selon les termes de la <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Licence Creative Commons Attribution 4.0 International</a>.

Il est donc autoriser d'adapter l'oeuvre et de l'utiliser à des fins commerciales. Il s'agit d'une license très peu contraignantes. Il suffit simplement de créditer l'auteur de l'oeuvre original.

Il existe également des licenses plus spécifiques aux bases de données, de la moins contraignante à la plus contraingante [PDDL](https://opendatacommons.org/licenses/pddl/) , [ODC-by](https://opendatacommons.org/licenses/by/) et [ODbL](https://opendatacommons.org/licenses/odbl/)

##### Intéropérable (*Interoperable*)

Les données sont associées à des métadonnées riches. Sur base des métadonnées, les données doivent être utilisable, compréhensible et combinable avec d'autres données. Ce principe est très compliqué à mettre en place et requiert une réflexion approfondie.

##### Réutilisable (*Reusable*) 

Les données doivent tant que possible être associées à des métadonnées riches avec une licence claire afin de pouvoir être réutilisé.

Vous avez certainement le sentiment que ces 4 principes se mélangent un peu. En effet, ils insistent avec des petites nuances sur des concepts particuliers. 

En Résumé : 

$$Données \ inutilisables = données \ seules$$


$$Données \ utilisables = données + contexte$$ 

Le contexte c'est :

- un code unique et persistant associé aux données et au contexte
- une description du le projet associée aux données
- des metodonnées riches sur les données
- la licence associée aux données 

Dans le cadre de vos futures recherches, un outil comme zenodo est très intéressants pour suivre au mieux le principe **FAIR**

#### DMP

Afin de respecter ces principes **FAIR**, des outils ont été développés. Il s'agit des plans de gestion des données (ou Data Management Plan : **DMP**). L'université de Mons dispose d'un [DMP](https://dmponline.be). Cet outil est partagé par l'ensemble des universités de Belgique.

Lorsque vous allez concevoir un plan d'expérience, n'oubliez pas de concevoir votre plan de gestion des données. 

Voici une checklist pour un plan de gestion des données efficaces <http://www.dcc.ac.uk/sites/default/files/documents/resource/DMP_Checklist_2013.pdf>


##### A vous de jouer ! {-}

**cet exercice va être déplacé dans l'exercice spécifique du module 8**

Suite à la lecture de la section sur la gestion des données avec le principe *FAIR*. Rédiger un document de synthèse de maximum 2 pages au format pdf qui critique votre plan de gestion des données réalisée dans le cadre des données obtenues sur la biométrie humaine. 

Dans le cadre du cours de science des données biologiques I, vous avez donc du concevoir un plan d'expérience afin d'acquérir des données. Etant donnée que vous avez été des scientifiques consciencieux, vous  avez collecté vos données en y associant des métadonnées.

Les données ont été collectée via un document googlesheet dont l’url est le suivant:

  - <https://docs.google.com/spreadsheets/d/1UfpZvx1_nd7d10vIMAfGVZ1vWyIuzeiKxPL0jfkNSQM/edit#gid=0>
  
Les métadonnées associées aux données ont été recensées dans un document googledoc dont l’url est le suivant:

  - <https://docs.google.com/document/d/1lgYD39W7vmVYyS5ea0wEl9ArE1dhuDRkIIBzZ4K6d_o/edit>
  
### Utilisation de données ouvertes (*Open Data*)

Il existe de nombreux sites qui regroupent un ensemble de données que nous allons appeler *Open Data*. Nous avons parlé précédement de Zenodo mais de nombreuses bases de données sont égalemement disponible comme comme le [Portail européen de données](https://www.europeandataportal.eu/fr), [Portail belge de données](https://data.gov.be/fr),...

Afin de connaitre la qualité des données voici une checklist très utile pour appréhender des données ouvertes 

Vous devez être capable de trouver facilement 

- But des données
- Code unique et persistant des données
- Licence des données
- Format des données
- Qualité des données

**A nouveau, vous vous rendez compte que nous revenons à notre principe FAIR expliqué plus haut.**

Prenons notre exemple sur le *Dataset: Number of diagnoses with coronavirus disease (COVID-19) in The Netherlands* et appliquans notre checklist.

- Le but des données, 

Une description des données est proposée. Le nom de l'auteur est spécifié. Il est également précisé la date de publication avec la version des données. Le 16 mars 2020, la version est `v2020.3.16`. Les données sont égalemetn associée à un dépot github qui les traitent : [J535D165/CoronaWatchNL](https://github.com/J535D165/CoronaWatchNL)

- Un code unique et persistant

Ces données ont un DOI : [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3711592.svg)](https://doi.org/10.5281/zenodo.3711592)

- La licence

Les données sont mise à disposition avec la licence  [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/legalcode)

- le format

Les données sont proposées sous le format `csv`. Ce format est un standard très employé. Il est à privilégier par rapport à format `.xlsx`

- La qualité

Ce dernier critère est le plus difficile étudier. Une première chose à vérifier sont les métadonnées associés à chaque variable monitorée. Nous pouvons voir que l'auteur peut encore améliorer les métadonnées associées à ces données. Le nom des variables reste cependant tout à fait compréhensible. 

Comme vous venez de le voir, zenodo de part sa structuration permet de réaliser très simplement cette checklist.

##### Pour en savoir plus {-}

- [Qu'est ce que l'Open data ?](https://www.europeandataportal.eu/fr/training/what-open-data) 

- [Aide sur l'interprétation et le choix des licenses](https://www.europeandataportal.eu/en/training/licensing-assistant)

- [Choisir la bonne license Open Source](https://choosealicense.com)

- Des stockages spécifiques ont été mis en place pour les données scientifiques comme [zenodo (dépot des données hebergé par le CERN)](https://zenodo.org) ou encore [Figshare](https://figshare.com)

- [Data Management Plan](https://www.unil.ch/openscience/fr/home/menuinst/open-research-data/gerer-ses-donnees-de-recherche/data-management-plan-dmp.html) 

- Le [Principe FAIR](https://ogsl.ca/fr/principes-fair) expliqué par l'Observatoire gloval du Saint-Laurent.

- L'outil institutionnel de l'université de Mons afin de réaliser un plan de gestion de données est disponible [DMPonline.be](https://dmponline.be)

- Article scientifique sur le FAIR plan : [The FAIR Guiding Principles for scientific data management and stewardship](https://www.nature.com/articles/sdata201618)

- [Guide sur les licenses open data](https://theodi.org/article/publishers-guide-to-open-data-licensing/)
