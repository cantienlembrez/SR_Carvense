# Stage Sex Ratio de la cirse des champs (*Cirsium arvense*)
Stagiaire : Cantien Lembrez

Au sein de l'IDEEV

Encadré par Diala Abu Awad

# Mardi 10/03

## Problème
 - Observations sur les populations de *C. arvense* de l'IDEEV :
 	-  Sex ratio fortement biasé en faveur des femelles
	 - Pas d'observation de fleurs hermaphrodites
	 - Différence de phénologie entre males et femelles : de facon surprenante les femelles fleurissent avant les mâles
  - des hermaphrodites ont été observé dans certaines populations (par exemple Kay 1985 )

## Hypothèses
- (H1) la différence de sex ratio est causé par une mortalité plus importante des mâles 
- (H2) la différence de phénologie suffit à expliquer le biais de sex ratio (comme l'echantillonage est effectué à un instant t)
- (H3) la présence d'un petit nombre d'hérmaphrodites dans la population pourrait introduire un biais de sex ratio m/(m+f) avec un mécanisme tel que décrit par  Nguyen et Pannel 2025 avec un allèle nucléaire qui produit des mâles et une stérilité mâle cytoplasmique (CMS)

**But** : Construire un modèle  avec ces hypothèses (au moins H1 et H3 dans un premier temps) pour voir les conséquences  en terme de diversité notamment pour des loci microsatellites (en cours de séquençage / séquenceur mort) et comparer avec les attendus pour un modèle Wright-Fisher (WF) (H0).

Dans un premier temps je me concentre sur la taille efficace des populations ($N_e$)  (*i.e.* la taille que devrait avoir une pop type WF pour qu'elle ait le même équilibre mutation dérive que la population d'interêt).


**Hypothèses communes des modèles** :
- populations de taille fixe ($N$)
- Génération non chevauchante / reproduction simultanée (pose problème pour H2 on va  pas pouvoir la tester avec cette approche) 

# Mercredi 11/03

Creation du modèle H0 avec SLiM 5 (Haller et al. 2025) avec des loci microsatellites nucléaire et des loci sur le cytotype.
Abandon de SLiM : la façon de le  sexe est implémenté rends difficile la création de modèle plus complexes pour H3. Passage à python 

# Jeudi 12/03

Ecriture du programme python pour H0 + biblio sur les microsatellites :

- estimation de la taille efficace avec microsatellites : 
	- facilement implémentable (Xu et Fu 2004) :
		- méthode  basé sur la variance en nombre de répétition $V_s$, $\hat\theta_{V_s} =2V_s$ avec $\theta = 4N_e\mu$ à l'équilibre mutation dérive, estimateur sans biais mais variance fonction de $\theta^2$
		- méthode basé sur l'homozigotie $\tilde\theta_F =\frac12\left(\frac{1}{F^2} -1\right)$ variance augmente moins vite avec $\theta$ mais estimateur biaisé ($\hat \theta_F = f(\tilde{\theta_F})$ non biaisé construit avec régression) 
- plus difficilement utilisables basé sur maximum de vraissemeblance (par exemple Nikolic et Chevalet 2014) : probablement pas necessaire avec des données simulées
- contexte général microsatellite *à lire* (Alves et al. 2024)

# Vendrdi 13/03

Dans les simulations les estimation de $N_e$ produisent des résultats incohérents systématiquements trops faibles par rapport aux attendus théoriques ; recherche du problème.
 
# Lundi 16/03

Implementation modèle pour H1.
et recherche sur le même problème que jour précédent.

# Mardi 17/03

Resolution du problème des jours précédent.
Début d'implémentation du cas avec CMS.

# Mecredi 18/03

Modèle avec CMS/trioecie implémenté (modèle similaire au modèle sans pollen limitation de Nguyen et Pannel 2025) .
avec de petites simulations tests on obtient les même fréquences à l'équilibre que dans le papier.

**Remarque** : le maintien de la trioecie dans le modèle dépends d'un coefficient de dépression de consanguinité $d$ constant or c'est une propriété dynamique des populations qui dépends du fardeau génétique, qu'on chercherait à modeliser en plus des loci neutres dans un second temps. 
Probleme : peut on maitenir la trioecie en  laissant evoluer dynamiquement la depréssion de consanguinité sans  $d$ fixe explicite.
 
# Jeudi 19/03

Pour l'instant on va s'interesser aux cas avec $d=0$ donc complétement neutre (par rapport au problème du jour précédent).
Por l'instant on a négligé la clonalité dans le cycle de reproduction : c'est complexe à integrer sur une seule génération (necessiterai d'être spatialement explicite et d'avoir des informations comment marche précisément) une solution pour au moins l'integrer d'une génération à l'autre serais d'echantilloner des survivants  d'une génération sur l'autre repartent des racines.

# Vendredi 20/03
Exploration des paramètres du modèle  de Nguyen et Pannel (voir TrioecyFreqs.py).
- On peut maintenir une faible fréquence d'hermaphrodites même avec une depression de consanguinité faible (ou nulle)
- pour avoir un sex ratio biasé vers les femelles il faut avoir le paramètre $em$ grand (proportion de pollen avortifà cause de la CMS chez les mâles plus il est grand moins le restaurateur lié au Y est efficace)

# Lundi 23/03- Mecredi 1/04
## Simulations
Je sais plus précisément quand ont été fait les différentes choses 
paramètres des simulations
- 100 réplicats par jeux de param
- N:50, 200, 1000, 5000 (pour le modèle nul on prends N:Ne des autres modèles) pour avoir pour avoir différentes forces de dérive
- Durée 6N génération pour assurer qu'on ait atteint l'équilibre
- Pour modèle avec mortalité différentielle Sm = 0.5 et Sm = 0.75
- Voir plus tard pour params des autres modèles

Simulations qui ont été lancées au 1er Avril. SEED indique l'argument de np.random.seed indiqué avant chaque ensemble de 100 simulations :

| SEED | Model |   Commit  |   N  |               Params               |
|:----:|:-----:|:---------:|:----:|:----------------------------------:|
|   0  |   M2  |  510d7eb  |  50  |              Sm = 0.5              |
|   1  |   M2  |  510d7eb  |  200 |              Sm = 0.5              |
|   2  |   M2  |  510d7eb  | 1000 |              Sm = 0.5              |
|   3  |   M2  |  510d7eb  | 5000 |              Sm = 0.5              |
|   4  |   M1  |  510d7eb  |  44  |                  /                 |
|   5  |   M1  |  510d7eb  |  178 |                  /                 |
|   6  |   M1  |  510d7eb  |  889 |                  /                 |
|   7  |   M1  |  510d7eb  | 4444 |                  /                 |
|   8  |   M2  |  510d7eb  |  50  |              Sm = 0.75             |
|   9  |   M2  |  510d7eb  |  200 |              Sm = 0.75             |
|  10  |   M2  |  510d7eb  | 1000 |              Sm = 0.75             |
|  11  |   M2  | non lancé | 5000 |              Sm = 0.75             |
|  12  |   M1  | non lancé |  49  |                  /                 |
|  13  |   M1  | non lancé |  196 |                  /                 |
|  14  |   M1  | non lancé |  980 |                  /                 |
|  15  |   M1  | non lancé | 4898 |                  /                 |
|  16  |   M3  |  4f8e646  |  50  |  g=1.3, a=4.3, s=0.5, d=0, em=0.9  |
|  17  |   M3  |  4f8e646  |  200 |  g=1.3, a=4.3, s=0.5, d=0, em=0.9  |
|  18  |   M3  |  4f8e646  | 1000 |  g=1.3, a=4.3, s=0.5, d=0, em=0.9  |
|  19  |   M3  | non lancé | 5000 |  g=1.3, a=4.3, s=0.5, d=0, em=0.9  |
|  20  |   M1  | non lancé |   ?  |                  /                 |
|  21  |   M1  | non lancé |   ?  |                  /                 |
|  22  |   M1  | non lancé |   ?  |                  /                 |
|  23  |   M1  | non lancé |   ?  |                  /                 |
| 24   |   M3  | 4f8e646   |  50  | g=1.1, a=4.8, s=0.8, d=0.3, em=0.9 |
| 25   |   M3  | 4f8e646   |  200 | g=1.1, a=4.8, s=0.8, d=0.3, em=0.9 |
| 26   |   M3  | 4f8e646   | 1000 | g=1.1, a=4.8, s=0.8, d=0.3, em=0.9 |
| 27   |   M3  | 4f8e646   | 5000 | g=1.1, a=4.8, s=0.8, d=0.3, em=0.9 |

## Bibliographie et clarification du contexte
Biais de sex-ratio vers les femelles (dans les pops de l'IDEEV mais aussi dans d'autres populations voir, entre autres, Kay 1985, Lalonde et Roitberg 1994, ou la review de Heimann et Cussans 1996)
On a différentes possibilités pour la cause d'un biais de sexe ratio (au niveau des mechanisme et au niveau des causes ultimes) :
- c'est un effet du aux différences de phénologie males femelles (femelles fleurissent avant les mâles et plus longtemps données des populations de l'IDEEV)
- c'est causé par une différence dans la clonalité entre males et femelles : 
	- On observe un compromis entre allocation repro sexué et repro assexué
	- généralement on suppose que le coût des fonction mâles est plus faible que celui des fonctions femelles donc la clonalité devrait être associé à des sex-ratio biasés en faveur des males. Pour *C. arvense* on observe l'inverse mais on pourrait avoir un cout plus elevées pour les mâles (On peut aussi noter que dans Field et al, Evolution, 2013 ils trouvent qu'être une espèce clonale est plus souvent associé à un sex ratio biasé en faveur des femelles. Même si c'est juste une association statistique que les auteurs interpretent comme une conséquence du fait que la clonalité corrèle avec d'autres traits d'histoire de vie qui pourrait être plus important du point de vue du sexe ratio) 
- Combinaison clonalité et effet de métapopulation (voir Field et al, Annals of Botany, 2013) -> combinaison d'echantillonage des pops lors d'evenement de colonisation et de la clonalité qui ralentit la convergence du sex ratio vers l'equilibre. C'est peu probable étant donné qu'on devrait observeer des populations dans lesquels on a un biais vers les mâles.
- Différence de mortalité (c'est pas encore clair les raisons théoriques liées à cette hypothèse c'est aussi une histoire de différence d'allocation) 
- Trioecie avec une CMS dans (un cas avec peu d'hermaphrodites -> on est proche de la fixation de la CMS, Nguyen et Pannell)
- Autres (différences de fertilisation selon pollen avec X/Y (ou avec Z/W) voir Lloyd 1974 (pas encore lu), *etc*).

parmi ces hypothèses on choisit de s'interesser aux suivantes notamment : 
- Différence de mortalité
- Différence de clonalité
- Trioecie

chaque mécanisme au dessus peut avoir une influence sur la diversité génétique même si ils ne sont pas forcément impliqué dans le sexe ratio. On va donc croiser certains des scénarios. 

## Modélisation de la clonalité
Comme dit précédemment, on a peu d'informations spatiales sur la repro asexué de la cirse.
On choisit de modéliser des populations de champs (comme c'est pour elles qu'on a des données)
entre chque génération :
- une proportion $p$ survit au fauchage modèlise la perénité (nombre d'années max ? on sait pas trop ~2an voir Leathwick et Bourdôt 2012) 
- une proportion $c$ produit un ramet avec ces racines qui se retrouve physiquement séparé du reste du clone modélise la clonalité en tant que telle 
	- voir Solé et Al 2004 et Bodo Slotta et al 2010 pour des estimations des proportions de génotypes uniques ; il semble qu'il y ait une coquille dans l'article de Bodo Slotta proportion >1 dans la table 1 ?

# Point Mardi 7/04

### Clonalité
le modèle finalement retenu est le suivant : étant donné que la clonalité et la pérénité dépendent toutes les deux de l'investissemnt dans les racines après la floraison on procède ainsi:
- une proportion $p$ survit au fauchage modèlise la perénité (nombre d'années max ? on sait pas trop ~2an voir Leathwick et Bourdôt 2012) on sépare males et femelles $p_m = p\times K$ et $p_f = p$ pour faire un modèle ou le biais provient de d'une différence d'investissement dans la pérénité
- une proportion **des individus survivants (i.e. on suppose que si pas assez investi dans racines pas capable de repro vegetative)** $c$ produit un ramet avec ces racines qui se retrouve physiquement séparé du reste du clone modélise la clonalité en tant que telle 

### Simulations & programmes
Pas précisé avant mais on simule 10 loci microsats 

Jusqu'à présent j'utilisais une loi de Poisson $\mathcal{P}(Nlocus\times Ncopies\times \mu)$ (pour approx $\mathcal{B}(Nlocus\times Ncopies, \mu)$) pour tirer le nombre de mutations à la fin de chaque cycle (on peut les tirer tous à la fin étant donné que toutes les mutations sont neutres). Comme on a un taux de mutation élevé pour les loci microsatellites ($\mu = 10^{-3}$), la loi de poisson peut approximer plutôt mal la binomiale.
Maintenant je tire explicitement dans une loi binomiale même si c'est plus couteux (visuellement les resultats tests semblent pas très différents mais idéalement il faudrait refaire les simulations avec le nouveau code).


premiers resultats sans analyses approfondies
- beaucoup d'extinctions de pops de taille 50 (c'est probablement causé par la façon d'initialisé la pop. Avec 4 * 1/4 pour tous les phénotypes sexuels) pas très important cf point 2
- Toutes les populations avec $N=200$ et 20% des pops avec $N=1000$ ont fixé la CMS (c'est logique les paramètres pour avoirs des fréquences des différents phénotypes demandent d'être proches de la fixations pas peu d'hermaphrodites et de males sans CMS)
- Pas de diff visuelles entre modele avec diff de mortalité et modèle nul avec meme taille efficace pour heterozygotie sous Hw, Nb allèles par locus, $N_e$ quelque soit la facon de l'estimer 

# Références
