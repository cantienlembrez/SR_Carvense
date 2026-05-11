# Simulations

`BasicFunctions.py` Contient les fonctions de base : mutations recombinaison et clonalité

`Models.py` contient les différents modèles de reproduction sexuelle et sélection  

`Main.py` contient la boucle principale des simulations. Utilisation :
```
python Main.py MODEL CP N NBSAVE INT NBREPLICATES SEED PARAM_NAME:Value
```
chaque simulations est stocké dans un dossier Simulations/M_NX_Y, X taille de la pop, Y numéro de la répétition

- `MODEL`: modèle de reproduction sexuée. "M1" : Dioïque non biaisé, "M2" : Dioïque biaisé, "M3" : Trioécie (voir le modèle de Nguyen et Pannell 2025).
- `CP` : modèle avec pérennité / clonalité. 0 : générations non chevauchantes, 1 : pérennité / clonalité sans biais,  2 : pérennité / clonalité avec biais
- `N` : taille de la population
- `NBSAVE` : nombre de fois ou l'état de la population est sauvegardé (après la génération max)
- `INT` : intervalle des sauvegardes
- `NBREPLICATES` nombre de répétitions 
- `SEED` seed de l'aléatoire pour l'ensemble des répétitions
- `PARAM_NAME` paramètres :
    - `Gmax` (défaut : 6*N) générations max avant les sauvegarde de l'état de la pop.
    - modèle de reproduction dioïque biaisé
        - `Sm` (défaut : None) survie des mâles dans le 
    - modèle de reproduction trioïque
        - `g` (défaut : None) le nombre de graines produites par une femelle par rapport à la production d'un hermaphrodite
        - `a` (défaut : None) le nombre de grains de pollen produits par un mâle par rapport à la production d'un hermaphrodite
        - `em` (défaut : None) la proportion de grains de pollen perdus du fait de la CMS chez un mâle
        - `s` (défaut : None) probabilité d'autofécondation chez un hermaphrodite
        - `d` (défaut : None) probabilité qu'une graine issue de l'autofécondation soit avortive
    - modèle perrenité clonalité :
        - `p` (défaut : None) probabilité de survie d'une année sur l'autre 
        - `MaxAge` (défaut : None -> pas de limite d'âge) Age à partir duquel un individu ne peut plus survivre l'année suivante
        - `c` (défaut : None) probailité de se cloner après avoir survécu
        - `K` (défaut : None) dans les modèles avec biais facteur diminuant la probabilité de se cloner pour les mâles

# Analyse
`ScriptAnalysis/AnalysisFuncs.py` contient les fonctions de base pour l'analyse (hétéroozygotie, nb génotypes, ...)
`ScriptAnalysis/scriptCompSm.py` analyse modèle dioïque biaisé et comparaison avec modèle null
`ScriptAnalysis/scriptCompTrio.py` analyse modèle trioécie et comparaison avec modèle null
`ScriptAnalysis/scriptCompCp.py` analyse modèle clonal et comparaison avec modèle null
`ScriptAnalysis/scriptCompCpb.py` analyse modèle clonal avec biais et comparaison avec modèle null
`ScriptAnalysis/Compdiff.py` comparaison des différentiations entre males et femlles pour une même taille ($N$) entre tout les modèles (pour voir comment se comportent les indices)

# Autres
`TrioecyFreqs.py` explorations des paramètres sur les équilibres du modèle de Nguyen et Pannell 2025
`Ne_trio.py` application de la formule générale du Ne basé sur la coalescence (Charlesworth et Charlesworth 2010) au cas de la trioécie.

