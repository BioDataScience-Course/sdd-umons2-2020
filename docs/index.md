--- 
title: "A Test Bookdown Book"
author: "Philippe Grosjean & Guyliann Engels"
date: "2019-08-28"
site: bookdown::bookdown_site
output:
  bookdown::gitbook:
    includes:
      after_body: disqus.html
documentclass: book
bibliography: [book.bib, packages.bib]
biblio-style: apalike
link-citations: yes
github-repo: phgrosjean/bookdown-test
url: 'http\://phgrosjean.sciviews.org/bookdown-test/'
description: "Tutoriel interactif pour la science des données avec R & SciViews-R."
cover-image: "images/front-cover.png"
---

# Préambule {-}



Contenu de cet ouvrage en bref...

![](images/front-cover.png)

----

_Le matériel dans cet ouvrage est distribué sous licence [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/deed.fr)._

----

##### System information {-}


```r
sessioninfo::session_info()
```

```
# ─ Session info ──────────────────────────────────────────────────────────
#  setting  value                       
#  version  R version 3.4.4 (2018-03-15)
#  os       macOS  10.14.6              
#  system   x86_64, darwin15.6.0        
#  ui       X11                         
#  language (EN)                        
#  collate  en_US.UTF-8                 
#  ctype    en_US.UTF-8                 
#  tz       Europe/Brussels             
#  date     2019-08-28                  
# 
# ─ Packages ──────────────────────────────────────────────────────────────
#  package     * version date       lib source        
#  assertthat    0.2.1   2019-03-21 [1] CRAN (R 3.4.4)
#  bookdown      0.10    2019-05-10 [1] CRAN (R 3.4.4)
#  cli           1.1.0   2019-03-19 [1] CRAN (R 3.4.4)
#  crayon        1.3.4   2017-09-16 [1] CRAN (R 3.4.1)
#  digest        0.6.20  2019-07-04 [1] CRAN (R 3.4.4)
#  evaluate      0.13    2019-02-12 [1] CRAN (R 3.4.4)
#  htmltools     0.3.6   2017-04-28 [1] CRAN (R 3.4.0)
#  knitr         1.22    2019-03-08 [1] CRAN (R 3.4.4)
#  magrittr      1.5     2014-11-22 [1] CRAN (R 3.4.0)
#  Rcpp          1.0.1   2019-03-17 [1] CRAN (R 3.4.4)
#  rmarkdown     1.12    2019-03-14 [1] CRAN (R 3.4.4)
#  sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 3.4.4)
#  stringi       1.4.3   2019-03-12 [1] CRAN (R 3.4.4)
#  stringr       1.4.0   2019-02-10 [1] CRAN (R 3.4.4)
#  withr         2.1.2   2018-03-15 [1] CRAN (R 3.4.4)
#  xfun          0.7     2019-05-14 [1] CRAN (R 3.4.4)
#  yaml          2.2.0   2018-07-25 [1] CRAN (R 3.4.4)
# 
# [1] /Library/Frameworks/R.framework/Versions/3.4/Resources/library
```
