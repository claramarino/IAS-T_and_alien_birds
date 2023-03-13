# Traits of birds involved in biological invasions
---

Scripts and data for reproducing the results obtained by Marino & Bellard (2023) in the paper "When origin, reproduction ability, and diet define the role of birds in invasions"

## Description of the Data and file structure

### 1. Traits and groups of birds involved in biological invasions

Data file 00_Traits_alien_iast_birds.rds contains the taxonomy, the traits and the impact information of the 850 bird species involved in biological invasions (864 lines x 21 columns).
Lines refer to species (note that 14 species are duplicated because they are in both alien and IAS-threatened groups).

Column description:
- binomial (factor): binomial name of species as in IUCN Red List
- Ecosystem (binary): impact mechanism includes habitat/ecosystem modification (1) or not (0)
- Direct_sp_effect (binary): impact mechanism includes a direct effect on species (1) or not (0)
- Indirect_sp_effect (binary): impact mechanism includes an indirect effect on species (1) or not (0)
- combi (factor): combination of the three impact mechanisms (columns "Ecosystem", "Direct_sp_effect" and "Indirect_sp_effect" pasted together)
- group (factor): species is alien ("EICAT") or IAS-threatened ("IAS-T")
- group2 (factor): species is alien with high impact ("EICAT_imp"), alien with low impact ("EICAT_no_imp"), alien with no information on impact ("EICAT_DD"), IAS-threatened with high impact ("IAS-T_imp") or IAS-threatened with low impact ("IAS-T_no_imp")
- Species1 (factor): binomial name of species as in AVONET
- Hand.Wing.Index (numeric): Hand-wing index (ratio, no unit)
- Tail.Length (numeric): length of the tail (mm)
- hab_sum (numeric): number of different habitats used
- insular_endemic (binary): species is insular endemic (1) or not (0)
- volant (binary): species has the ability to fly (1) or not(0)
- Trophic.Level (factor): trophic level
- Primary.Lifestyle (factor): preferential foraging niche
- ln.Mass (numeric): log-transformed body mass (gram)
- ln.Clutch (numeric): log-transformed clutch size (number of eggs)
- ln.Beak.Depth (numeric): log-transformed beak depth (mm)
- ln.Beak.Length_Nares (numeric): log-transformed beak length (mm)
- bioreg (factor): bioregion of origin 
- insul_level (numeric): species is mainland endemic (1), present on mainland and islands (2), insular endemic (3)


### 2. Traits of worldwide birds

Data file 00_Traits_worldwide_birds.rds contains the taxonomy and the traits of all worldwide birds (10 943 lines x 13 columns).
Lines refer to species.

Column description:
- binomial (character): binomial name of species as in IUCN Red List
- Beak.Depth (numeric): beak depth (mm)
- Hand.Wing.Index (numeric): Hand-wing index (ratio, no unit)
- Tail.Length (numeric): length of the tail (mm)
- Mass (numeric): body mass (gram)
- hab_sum (numeric): number of different habitats used
- Clutch (numeric): clutch size (number of eggs)
- Beak.Length_Nares (numeric): beak length (mm)
- Trophic.Level (character): trophic level
- Primary.Lifestyle (character): preferential foraging niche
- bioreg (factor): bioregion of origin 
- insul_level (numeric): species is mainland endemic (1), present on mainland and islands (2), insular endemic (3)

### How to use the scripts

The first script computes the axes of the multidimensionnal functional space based on the traits of both alien and IAS-threatened birds. The second script performs the statistical analyses based on the functional space and creates the figures 1 and 4. The third script performs the trait-by-trait analysis and creates the figured 2 and 3. The last script describes the steps implemented for evaluating the prediction of impact type based on the position of birds in the functional space.


## Sharing/access Information

Links to other publicly accessible locations of the data:
https://github.com/claramarino/IAS-T_and_alien_birds

Was data derived from another source?
Data came from AVONET (Tobias et al. 2022), IUCN Red List of threatened species (IUCN 2022), Amniote database (Myhrvold et al. 2015), and EICAT database (Evans et al. 2016)

Evans, T. G., Kumschick, S., & Blackburn, T. M. (2016). Application of the Environmental Impact Classification for Alien Taxa (EICAT) to a global assessment of alien bird impacts. Diversity and Distributions, 22(9), 919–931. https://doi.org/10.1111/ddi.12464

IUCN. (2022). The IUCN Red List of Threatened Species. https://www.iucnredlist.org.

Myhrvold, N. P., Baldridge, E., Chan, B., Sivam, D., Freeman, D. L., & Ernest, S. K. M. (2015). An amniote life-history database to perform comparative analyses with birds, mammals, and reptiles. Ecology, 96(11), 3109–000. https://doi.org/10.1890/15-0846r.1

Tobias, J. A., Sheard, C., Pigot, A. L., Devenish, A. J. M., Yang, J., Neate-Clegg, M. H. C., Alioravainen, N., Weeks, T. L., Barber, R. A., Walkden, P. A., MacGregor, H. E. A., Jones, S. E. I., Vincent, C., Phillips, A. G., Marples, N. M., Montaño-Centellas, F., Leandro-Silva, V., Claramunt, S., Darski, B., … Schleunning, M. (2021). AVONET: morphological, ecological and geographical data for all birds. Ecology Letters, April, 1–17. https://doi.org/10.1111/ele.13898




