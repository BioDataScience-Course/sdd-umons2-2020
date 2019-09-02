# (APPENDIX) Appendices {-}


# Installation de la SciViews Box {#svbox}

<img src="images/sdd1_01/svBox-256.png" width="256px" style="display: block; margin: auto 0 auto auto;" />

Pour ce cours SDD 2, nous utiliserons la même SciViews Box... mais actualisée (version de l'année). Vous allez donc devoir installer la nouvelle version. La procédure n'a changé que sur des points de détails. Référez-vous à [l'appendice A1 du cours SDD 1](http://biodatascience-course.sciviews.org/sdd-umons/svbox.html). **Vous pouvez conserver l'ancienne SciViews Box en parallèle avec cette nouvelle version**, mais vérifiez si vous avez assez d'espace sur le disque dur pour contenir les deux simultanément. Comptez par sécurité 20Go par version. Si vous manquez de place, vous pouvez éliminer l'ancienne version avant d'installer la nouvelle (vos projets ne seront pas effacés).

**TODO: procédure pour éliminer une ancienne version.**


## Migration des projets

Concernant les projets réalisés dans une version précédente de la SciViews Box, ceux-ci restent disponibles, même si vous éliminez l'ancienne. Plusieurs cas de figure se présentent\ :

1. Vous conserver deux ou plusieurs version de la SciViews Box en parallèle. Dans ce cas, nous conseillons fortement de garder chaque projet accessible à partir de la version dans laquelle il a été créé. Seulement les projets que vous décidez de *migrer* explicitement (voir ci-dessous) seront à déplacer dans le dossier `shared` de la nouvelle SciViews Box. Vous aurez à faire cette manipulation, par exemple, si vous devez recommencer un cours l'année suivante afin d'être en phase (même version de la svbox) par rapport à vos nouveaux collègues.

2. Vous ne conservez que la dernière version de la SciViews Box, mais ne devez pas accéder fréquemment vos anciens projets, et dans ce cas, vous pouvez réinstaller temporairement l'ancienne version de svbox. Dans ce cas, ne migrez pas vos anciens projets. Éliminez simplement l'ancienne svbox, tout en laisant vos projets intacts dans son répertoire `shared`. Lors de la réinstallation de l'ancienne svbox, vous retrouverez alors tous vos anciens projets intactes.

3. Vous ne conservez pas d'ancienne version de la svbox et vous ne souhaitez pas devoir la réinstaller. Il est possible de *migrer* vos anciens projets en les déplaçant de l'ancien répertoire `shared` vers le nouveau. **Soyez toutefois conscients que vos documents R Markdown et scripts R ne fonctionneront pas forcément dans la nouvelle svbox et qu'une adaptation sera peut-être nécessaire\ !**


## Configuration Git et Github

A chaque nouvelle installation de la SciViews Box, vous devez la reconfigurer via la boite de dialogue `SciViews Box Configuration`. En particulier, il est très important d'indiquer correctement votre identifiant et email Git (zone encadrée en rouge dans la copie d'écran ci-dessous).

![](images/sdd2_A1/svbox_config.png)

Assurez-vous (si ce n'est déjà fait) que vous possédez un compte Github valide. Vous pouvez cliquer sur le bouton `Go to Github` par facilté dans la même boite de dialogue. **Choisissez de manière judicieuse votre login**. Vous pourriez être amenés à l'utiliser bien plus longtemps que vous ne le pensez, y compris plus tard dans votre carrière. Donc, lisez les conseils ci-dessous (inspirés et adaptés de [Happy Git and Github for the UseR - Register a Github Account](https://happygitwithr.com/github-acct.html)\ :

- Incluez votre nom réel. Les gens aiment savoir à qui ils ont affaire. Rendez aussi votre nom/login facile à deviner et à retenir. Philippe Grosjean a comme login `phgrosjean`, par exemple.

- Vous pouvez réutiliser votre login d'autres contextes, par exemple Twitter ou Slack (ou Facebook).

- Choisissez un login que vous pourrez échanger de manière confortable avec votre futur boss.

- Un login plus court est préférable.

- Soyez unique dans votre login, mais à l'aide d'aussi peu de caractères que possible. Github propose parfois des logins en auto-complétion. Examinez ce qu'il propose.

- Rendez votre login invariable dans le temps. Par exemple, n'utilisez pas un login lié à votre université (numéro de matricule, ou nom de l'université inclue dans le login). Si tout va bien votre login vous suivra dans votre carrière, ... donc, potentiellement loin de l'université où vous avez fait vos études.

- N'utilisez pas de logins qui sont aussi des mots ayant une signification particulière en programmation, par exemple, n'utilisez pas `NA`, même si c'est vos initiales\ !

Une fois votre compte Github créé, et votre login/email pour votre identification Git correctement enregistrés dans la SciViews Box, vous devez pouvoir travailler, faire des "pushs", des "pulls" et des "commits"^[Vérifiez **toujours** lors de votre premier commit que Github vous reconnait bien. Pour cela, naviguez vers le dépôt où vous avez commité avec votre explorateur web, et vérifiez l'identité prise en compte lors de votre commit.]. Cependant, RStudio vous demandera constamment vos logins et mots de passe... à la longue, c'est lassant\ ! La procédure ci-dessous vous enregistre une fois pour toutes sur votre compte Github dans RStudio.

### Compte Github dans RStudio

RStudio offre la possibilité d'enregistrer une clé publique/privée dans votre SciViews Box afin de vous enregistrer sur Github de manière permanente. L'avantage, c'est que vous ne devrez plus constamment entrer votre login et mot de passe à chaque opération sur Github\ ! Nous vous le conseillons donc vivement.

- Entrez dans Rstudio Server, et allez dans le menu `Tools -> Global Options...`. Ensuite, cliquez dans la rubrique `Git/SVN` dans la boite de dialogue.

![](images/sdd2_A1/github_config_rstudio1.png)

- Ensuite, cliquez sur le bouton `Create RSA key...`. La phrase de passe n'est pas nécessaire (il est même préférable de la laisser vide si vous voulez utiliser Github sans rien devoir taper à chaque fois). Cliquez sur le bouton `Create`.

![](images/sdd2_A1/github_config_rstudio2.png)

- Vous obtenez alors une fenêtre similaire à celle ci-dessous (bien sûr avec des données différentes). Ceci confirme que votre clé cryptographique a été créée localement. Fermez cette fenêtre pour revenir à la boite de dialogue de configuration de RStudio Server.

![](images/sdd2_A1/github_config_rstudio3.png)

- Dans la boite de dialogue de configuration de RStudio Server, section `Git/SVN` cliquez sur le lien `View public key` qui apparait une fois la clé créée\ :

![](images/sdd2_A1/github_config_rstudio4.png)

- La clé apparait dans une fenêtre, déjà présélectionnée. Copiez-là dans le presse-papier (`Ctrl-C` ou clic bouton droit et sélection de `Copy` dans le menu contextuel), puis fermez cette fenêtre.

![](images/sdd2_A1/github_config_rstudio5.png)

- Dans votre navigateur web favori, naviguez vers https://github.com, loggez-vous, et accédez aux paramètres de votre compte Github (menu déroulant en haut à droite, entrée `Settings`)\ :

![](images/sdd2_A1/github_config_rstudio6.png)

- Dans les paramètres de votre compte, cliquez sur la rubrique `SSH and GPG keys`, ensuite sur le bouton vert `New SSH key`

![](images/sdd2_A1/github_config_rstudio7.png)

- Collez-y votre clé à partir du presse-papier dans la zone `Key`. Vous pouvez lui donner un nom évocateur dans le champ `Title`. Ensuite, cliquez sur `Add SSH key`.

![](images/sdd2_A1/github_config_rstudio8.png)

- Déloggez, puis reloggez-vous dans RStudio Server pour que les changements soient pris en compte. La prochaine action sur Github depuis RStudio pourrait encore déclencher la demande de votre login et mot de passe, mais ensuite, les opérations devraient se faire directement.

Si vous éprouvez toujours des difficultés à faire collaborer R et RStudio avec Git et Github, voyez [https://happygitwithr.com](Happy Git and Github for the UseR) (en anglais) qui explique les différentes procédures bien plus en détails.
