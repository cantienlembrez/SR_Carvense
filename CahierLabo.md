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
 

# Références