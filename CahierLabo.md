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

# Point 17/04

## Modele clonalité (definitif ?)

J'avais pas précisé mais on néglige les mutations somatiques (pour les microsats c'est discutable : la majorité des mutations -> slippage de replication (Ellegren 2004) qui pour le coups ce produisent dans lignées somatiques
Mais également une partie est causée par erreurs de recombinaison. Donc on aurait pu faire le choix de prendre un $\mu_{clonal}$ plus faible que $\mu$)

**modification :** le modèle avec biais diffère du modèle sans biais au niveau la probabilité de clonage ($c_f = c$ et $c_m = cK$), au lieu d'une différence de survie (plus raccord avec la litterature)

## Simulations lancées 
j'ai accidentellement kill le job avec la seed 10 (cf tab précédent) relancé sur nouvelle version code 
nouvelles simulatiosn lancées avec code commit a2fcdd3 : 
- Sans seed (à cause d'une fausse manip)
	- Clonalité non biaisé $p=0.1$ $c=1$ AgeMax = 2 N = 1000 (j'ai également pour p=0.1 et c=0.1 par erreur)
		- Modele repro sexué sans biais
		- Modèle cms (g=1.3, a=4.3, s=0.5, d=0, em=0.9)
	- Clonalité biasé : $p=0.1$ $c=1$ AgeMax = 2 K = 0 N = 1000
		
Avec seed cf tableau (Lundi 23/03- Mecredi 1/04) :
	- Seeds : 10 à 15 (et 19 avec commit précédent f16a881)

## Expression Analytique de Ne dans le cas avec trioecie

Idée d'approche pour trouver une expression Analytique de $N_e$ partir de (Charlesworth et Charlesworth 2010) chapitre 5.2 

### Compréhension du Charlesworth et Charlesworth

Ils proposent la formule générale pour $P_c$ la probabilité que deux
allèles pris aléatoirement chez deux individus différents coalescent
d'une génération sur l'autre :
$$P_c = \sum_{rsu}\alpha_r\alpha_s\beta_ {rsu}\Theta_{rsu}\gamma_{rsu}$$
avec :

-   $\alpha_r$ et $\alpha_s$, *the probabilities that they come from
    indivuals of classes(compartements) r and s, respectively*

-   $\beta_ {rsu}$, *the probability that the pair of alleles with
    origins $r$ and $s$ both come from a parent of sex $u$*

-   $\Theta_{rsu}$ *the probability that both alleles come from the same
    individual, given that an $rs$ pair of alleles was derived from a
    parent of sex $u$*

-   $\gamma_{rsu}$ *the probability that the alleles coalesce whthin
    that individual*

**remarque** : les probabilités $\beta_ {rsu}$, $\Theta_{rsu}$ et
$\gamma_{rsu}$ sont des probas conditionnelles qui supposent que tous
les évènements précédents se soient réalisés

Si j'ai bien compris $\alpha_r$ sont \"du point de vue\" de l'allèle
*i.e.* pour un allèle autosomal, $\alpha_f = \alpha_m = \frac12$ sans
dépendre des effectifs. Comprendre $\alpha_r$ : probabilité que l'allèle
choisi soit dans la classe $r$

les $\beta$ pour autosome valent tous $\frac14$ par ce que c'est (je
note $P(r \leftarrow u)$ la proba que un allèle chez un $r$ vienne d'un
parent $u$)
:$$\beta_ {rsu}  =  P(r \leftarrow u) \times P(s \leftarrow u)$$ or
quand chaque sexe a la même contribution à la descendance donc
$P(r \leftarrow f) = P(r \leftarrow m) = \frac12$. (Est ce que c'est
vrai parce que tout les individus ont forcément un père et une mère ? Si
on transpose au modèle Nguyen et Pannell je ne suis pas sur que cette
hypothèse soit respectée : faut il intégrer les coefficients $g$,
$\alpha$ et $s$ de leur modèle à cet endroit ?) . Pour comprendre le
principe, on peut faire pour les allèles liées au X (pour un mâle X
vient forcément de la mère) :

-   $\beta_{mmf}=  P(m \leftarrow f) \times P(m \leftarrow f)  =  1\times 1 = 1$

-   $\beta_{fmf} =  P(f \leftarrow f) \times P(m \leftarrow f)  = \frac 12\times 1 = \frac12$

-   $\beta_{fmm} =  P(f \leftarrow m) \times P(m \leftarrow m)  = \frac 12\times 0=0$

-   $\beta_{mmf} =  P(m \leftarrow f) \times P(m \leftarrow f)  = \frac 12\times\frac 12=\frac 14$

-   etc.

Pour $\gamma_{rsu}$, la proba de coalescence dans un individu c'est
simplement $\frac12 + \frac{F_{is}}{2}$ premier terme identité de
provenance et second : identité de copie.

Pour les coefficients $\Theta_{rsu}$, je vais refaire son développement
de $\Theta_{fff}$ en détaillant plus. On cherche la proba que deux
allèles chez deux filles différentes viennent de la même mère sachant
qu'ils viennent chacun d'une femelle. Notons $d_{fk}$ le nombre de fille
($d$ pour daugther) produit par la k-eme femelle. soit
$\bar d_{f} = \frac{1}{N_f}\sum_{k=0}^{N_f}d_{fk}$ le nombre moyen de
filles par mère. La proba qu'**une** fille soit la fille de la k-eme
mère est de
$$\frac{d_{fk}}{\sum_{k=0}^{N_f}d_{fk}} = \frac{d_{fk}}{N_f\bar d_{f}}$$
donc la proba que **deux filles viennent de la même mère** est de (si la
pop est suffisamment grande pour négliger les -1 au dénominateur) :
$$\frac{d_{fk}(d_{fk} -1)}{(N_f\bar d_{f})^2}$$ Et la proba que **toutes
les deux viennent de la même mère n'importe laquelle** s'obtient en
sommant sur k :
$$\Theta_{fff} = \sum_{k=0}^{N_f}\frac{d_{fk}(d_{fk} -1)}{(N_f\bar d_{f})^2} 
    \label{thetabrut}$$

Le développement suivant permet de simplifier l'expression
([\[thetabrut\]](#thetabrut){reference-type="ref"
reference="thetabrut"}) en utilisant des information sur la distribution
à la place du nombre de descendant de chaque femelle. On veut tomber sur
:
$$\Theta_{fff} = \frac{\Delta V_{ff}/\bar d_f^2 + 1}{N_f} = \frac{1}{(N_f\bar{d}_f)^2}\left(N_f(\Delta V_{ff} + \bar{d}_f^2)\right)
    \label{objectif}$$ Si tous les individus produisent le même nombre
de descendants, alors la distribution du nombre de descendants pour un
individu suit une loi de poisson de paramètre $\bar d_{f}$ (pour le
développement je vais noter autrement que dans le livre, je note $FF$ la
distribution du nb de filles produite par une mère). Si cette hypothèse
n'est respectée on peut calculer la déviation à la loi de poisson :
$$\Delta V_{ff} = Var(FF) - E(FF) = E(FF^2) - E(FF)^2 - E(FF)= \sum \frac {d_{fk}^2}{N_f} - \bar d_f^2 - \bar d_f
    \label{delta}$$

en injectant [\[delta\]](#delta){reference-type="ref" reference="delta"}
dans [\[objectif\]](#objectif){reference-type="ref"
reference="objectif"} : $$\begin{aligned}
        & \frac{1}{(N_f\bar{d}_f)^2}\left[N_f\left( \sum \frac {d_{fk}^2}{N_f} - \bar d_f^2 - \bar d_f + \bar{d}_f^2\right)\right]\\
        & \frac{1}{(N_f\bar{d}_f)^2}\left( \sum  d_{fk}^2- \bar d_fN_f \right) =  \frac{1}{(N_f\bar{d}_f)^2}\left( \sum  d_{fk}^2- N_f\frac{\sum d_{fk}}{N_f} \right) \\
        & \sum_{k=0}^{N_f}\frac{d_{fk}(d_{fk} -1)}{(N_f\bar d_{f})^2}  = \Theta_{fff} 
\end{aligned}$$

ensuite ils proposent une simplification supplémentaire pour n'avoir
plus que $\Delta V_F$ la déviation de la distribution des descendant
produits quelque soit le sexe par une femelle à loi de poisson.

Dans leur tableau des $\Theta$ livre, pour les $\Theta_{rsu}$ avec
$r\neq s$ ils font intervenir $C_{rmf}$ la covariance entre le nombre de
fils et filles produit par le sexe $r$ à la place de $\Delta V_{rf}$ ça
doit venir du même type d'argument. Logiquement on devrait avoir :
$$C_{fmf}= \sum_{k=0}^{N_f}\frac{(\bar d_f - d_{fk})(\bar s_f  - s_{fk})}{N_f}$$
($s_{fk}$ est le nb de fils produit par la k-eme femelle.)

### Application à Nguyen et Pannell

#### Choix approche

Je pense que le meilleur choix est de partir avec 4 classes : mâles
($M$) mâles avec cms ($C$) hermaphrodites ($H$) et femelles ($F$). Dans
le modèle, les effectifs de chaque sexe sont stables, il n'y a pas de
raison pour que ils ait de déviation à la loi de poisson (on a
$\Delta V_{xx} = 0$). Je vois pas non plus de raison pour que les
covariance soient non nulles (sauf peut-être : à cause de
l'autofécondation $C_{HHM}\neq0$ mais normalement non les tirages de
chaque descendant sont indépendants).

On devrait donc avoir $\Theta_{rsu} = \frac{1}{N_u}$ et la complexité
c'est dériver les $\beta_{rsu}$. Il faut trouver les contribution des
parents aux classes des enfants.

Avec 3 classes je pense que les $\beta$ sont pas beaucoup plus simple à
trouver et il faut en plus trouver les déviations aux lois de poisson
pour les mâles avec ou sans cms.

#### Calculs des coefficients $\beta$

revient à calculer 16 $P(r \leftarrow u)$ (rappel la proba que un allèle
chez un $r$ vienne d'un parent $u$ )

##### Cas simple : contribution à M {#cas-simple-contribution-à-m .unnumbered}

Les croisement qui donnent M (XYc) sont XXc x XYc et XXc x XYn. La
contribution génétique des parents (pour allele autosomal) c'est 1/2
pollen et 1/2 ovules.

On trouve facilement que $P(M \leftarrow F) = 0$, et que
$P(M \leftarrow H) = 1\times \frac12$. Pour $M$ et $C$ la contribution
devrait dépendre des fréquences relatives de pollen. Pollen total
produit par $M$ : $\alpha N_M$ (c'est le $\alpha$ du modèle de trioecie
pas celui de l'équation (1)) et pollen total produit par $C$ :
$(1-em)\alpha N_C$.
$P(M \leftarrow M) = \frac{\alpha N_M}{\alpha N_M + (1-em)\alpha N_C }\times \frac12 = \frac{N_M}{2(N_M + (1-em )N_C)}$
et $P(M \leftarrow C) = \frac{(1-em )N_C}{2(N_M + (1-em )N_C)}$.

Les autres $P(r \leftarrow u)$ devraient être trouvables de la même
manière (les hermaphrodites vont complexifier l'affaire pour certaines
classes).

#### Autres coefficients

les coefficient $\alpha_i$ probabilité que l'allèle choisit soit dans la
classe $i$ devraient tous valoir 1/4. Pour $\gamma_{rsu}$, il faudrait
trouver l'expression du $F_{IS}$ à l'équilibre (ça doit se trouver en
considérant le taux d'autofécondation effectif
$\frac {N_H}{N}\times s(1-d)$).


# Point Analyse Simulations

Cette section a été modifié plusiseurs fois pour tracer les changement suivre les commits git.

### Comparaison modèle mortalité et modèle null (scriptCompSm.py)

tableau : en premier la moyenne pour le modèle focal en 2eme celle du modèle null $N_e$ equivalente, puis resultats test de wilcoxon pour la comparaison des deux groupes.

Pour le premier tableau sur les allèles
colone :
- Ne (F) : Ne estimé avec $\hat\theta_F$
- Ne (Vs) : Ne estimé avec $\hat\theta_{V_s}$
- He (iam) : heterozigotie attendue sous panmixie pour loci infinite allele model 
- Ne : estimation Ne à partir colonne precédente 
- Fis global pour les loci microsat $L=10$ (moins les locis non polymorphes) $$\frac1L\times\frac{\sum H_o}{\sum H_e}$$

##### Avec Sm = 0.5
| M surv b | M1 eq    | Ne (F)                                | Ne (Vs)                               | He (iam)                           | Ne (iam)         | Fis                                        | N all                              |
|----------|----------|---------------------------------------|---------------------------------------|------------------------------------|------------------|--------------------------------------------|------------------------------------|
| $N=50$   | $N=44$   | 44.29, 44.93 (w = 2463, p = 0.83)     | 41.66, 40.84 (w = 2458, p = 0.81)     | 0.142, 0.158 (w  = 1730, p = 0.55) | 41.53, 46.08     | -1.5e-3, -3.6e-3 (w = 1895, p = **0.030**) | 1.74, 1.68 (w = 1587, p = 0.064)   |
| $N=200$  | $N=178$  | 174.35, 185.53 (w = 2139, p = 0.18)   | 159.94, 162.25 (w = 2330, p = 0.50)   | 0.402, 0.400 (W = 2403, p = 0.67)  | 167.93, 166.60   | -6.3e-3, -6.0e-3 (w = 2516, p = 0.98)      | 3.14, 3.09 (w = 1530.5, p = 0.19)  |
| $N=1000$ | $N=889$  | 916.18, 891.87 (w = 2161, p = 0.21)   | 895.14, 845.31 (w = 2205, p = 0.27)   | 0.786, 0.789 (w = 2505, p = 0.95)  | 918.35, 934.65   | -1.1e-4, -4.9e-5 (w = 2261, p = 0.66)      | 6.08, 6.05 (w = 2162.5, p = 0.66)  |
| $N=5000$ | $N=4444$ | 4359.67, 4389.08 (w = 2439, p = 0.77) | 4000.02, 4023.73 (w = 2455, p = 0.81) | 0.948, 0.945 (w = 2171, p = 0.22)  | 4586.67, 4329.27 | -1.6e-5, -2.3e-5 (w = 2419, p = 0.71)      | 12.22, 12.22 (w = 2225.5, p =0.71) |




# Références
