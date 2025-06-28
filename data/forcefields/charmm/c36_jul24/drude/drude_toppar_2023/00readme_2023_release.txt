Overview of the 2023 polarizable Drude oscillator force field release.
July 2023

!! The Drude LJ terms with zero Rmin previously required to
   simulations with OpenMM are no longer needed.

2023_7:2023 versions of main, carbohydrate, lipid and model compound toppar files created. 

General comments
----------------

2023 version files created motivated by the additional optimizaiton of
the lipid parameters to be compatible with LJ PME. Based on this the
recommnended truncation scheme for the entire Drude FF with periodic
boundardy conditions in CHARMM is as follows. LJPME may also be used
with OpenMM and Gromacs and it is being implemented in NAMD.

Yu, Y., Venable, R.M., Thirman, J., Chatterjee, P., Kumar, A., Pastor,
R.W.,* Roux, R.,* MacKerell, A.D., Jr.,* and Klauda, J.*, “Drude
Polarizable Lipid Force Field with Explicit Treatment of Long-Range
Dispersion: Parametrization and Validation for Saturated and
Mono-unsaturated Zwitterionic Lipids,” Journal of Chemical Theory and
Computation, 19: 2590-2605, 2023, PMC10404126,
10.1021/acs.jctc.3c00203

set 3 12.0  ! cutim
set 4 12.0  ! cutnb
set 5 9.0   ! ctonnb
set 6 9.0   ! ctofnb
set 7 switch
set 8 atom
set 9 vatom

set ord 6      !for both electrostatic and LJPME
set kap 0.40   !electrostatic PME
set ljkap 0.40 !LJ PME
set fftxd 32   !for both electrostatic and LJPME
set fftyd 32   
set fftzd 32

! set nonbond variable
update cdie inbfrq -1 imgfrq -1 ihbfrq 0 -
ewald pmewald kappa @kap fftx @fftxd ffty @fftyd fftz @fftzd order @ord -
ljpme dkappa @ljkap dftx @fftxd dfty @fftyd dftz @fftzd dorder @ord -
@7 @8 @9 vswitch cutimg @3 cutnb @4 ctofnb @6 ctonnb @5 ! bycb !may improve performance on some machines

## main file and proteins (name change from "master" to "main")
toppar_drude_main_protein_2023a.str
see diff_main_2019_2023.txt for changes from 2019h version

!summary
!update of alkene parameters required for lipids associated with 2023 lipid LJ PME update. This impacts a
!significant number of parameters
!atom types added for glycopetide and modified carbohydrates
!removal of special lipid carbonyl atom types and merging of lipid/nucleic acid phosphate atom types and parameters
!sulfate ion added and associated molecular anion parameters
!PRES NNEU corrected and uncommented
!PRES GNNU added ! neutral N-terminus for Glycine; charges from MAM1, parameters added
!PRES HS2 corrected to reproduce RESI HSE
!CMAP terms associated with alkenes deleted
!
!Note: Yiling's ion-biomolecule optimization work will require removal/alteration of 
!associated NBFIX/NBTHOLE parameters
!

## d_amino acids
toppar_drude_d_aminoacids_2023a.str 
!copy cp toppar_drude_main_protein_2023a.str toppar_drude_d_aminoacids_2023a.str followed by
!change residue names (eg ALA to DALA)
!"ATOM CA   CD31C" to "ATOM CA   CD31DC"
!change of IC tables for to D configuration
!NOTE: check the the final CMAPs are based on Jing's 2017 update for L aas.

## carbohydates
toppar_drude_carbohydrate_2023a.str
see diff_carbo_2019_2023.txt for changes from 2019 version
!summary
!addition of N- and O-linked glycosidic linkages to proteins. 2022
!
!corrections to hexapyranose C5, C6 and O5 atom charges
!
!add residues 
!RESI AXYL           0.000  ! alpha-D-xylose (a-glucose w/o exocyclic moiety) og
!RESI BXYL           0.000  ! beta-D-xylose (b-glucose w/o exocyclic moiety) og
!RESI AFUC           0.000  !  alpha-L-Fuctose
!RESI BFUC           0.000  !  alpha-L-Fuctose
!RESI ARHM           0.000  !  alpha-L-rhamnose
!RESI BRHM           0.000  !  beta-L-rhamnose
!RESI AQUI           0.000  ! alpha-D-Quinovose
!RESI BQUI           0.000  ! beta-D-Quinovose
!RESI ARHMNA         0.000  !  alpha-L-rhamnose	with C2	N-acetyl
!RESI BRHMNA         0.000  !  beta-L-rhamnose	with C2	N-acetyl
!RESI ANE5AC        -1.000  !  alpha-Neu5Ac (N-acetyl-alpha-D-neuraminic acid)
!RESI BNE5AC        -1.000  !  beta-Neu5Ac (N-acetyl-alpha-D-neuraminic acid)
!
! patches corrections and clarifications based on the following nomenclature
!Glycosidac patch descriptions
!Initial characters are numbers: glycosidac linkages between pyranoses
!Initial characters are two letters indicating specific classes of monosaccharides being linked
!"p" indicates a pyranose
!"f" indicates a furanose with an exocyclic hydroxyl on C1 (eg. RESI BFRU, beta-Fructofuranose)
!"a" indicates a furanose withOUT an exocyclic hydroxyl on C1 (eg. RESI ADEO, alpha-Deoxy-Ribose)
!
!new patches
!
!PRES 1AOME          0.000 ! adding axial O-methyl to C1 on pyranose
!PRES SGPA          -0.107 ! glycosidic linkage for connecting SER at the alpha position
!PRES SGPB          -0.107 ! glycosidic linkage for connecting SER at the beta position
!PRES TGPA          -0.047 ! glycosidic linkage for connecting THR at the alpha position
!PRES NGLA          -0.075 ! N-glycosidic link at ASN (alpha), previous charge 0.064
!PRES SA23AB        -0.129 ! sialic acid alpha 2->3 equatorial, og
!PRES SA26AT        -0.479 ! (i)2->6(i-1) equatorial at C2 and equatorial at C6
!
!numerous parameters for glycosidic linkages to proteins
!

## model compounds
toppar_drude_model_2023a.str
!see diff_model_2019_2023.txt
!summary
!RESI PRP1 to PRPA
!Alkene model compounds added
!RESI 2BTE          0.000 ! 2-butene electrostatic parameters predicted by DNN
!RESI 2BTD          0.000 ! 2-butene double bond electrostatic parameters based on analogy with 5-decene
!RESI 2PTE          0.000 ! 2-pentene electrostatic parameters predicted by DNN
!RESI 2PTD          0.000 ! 2-pentene double bond electrostatic parameters based on analogy with 5-decene
!RESI 2HEX          0.000 ! 2-hexene electrostatic parameters predicted by DNN
!RESI 2HED          0.000 ! 2-hexene double bond electrostatic parameters based on analogy with 5-decene
!RESI 3HEX          0.000 ! 3-hexene electrostatic parameters predicted by DNN
!RESI 3HED          0.000 ! 3-hexene double bond electrostatic parameters based on analogy with 5-decene
!RESI 5DEC          0.000 ! 5-decene electrostatic parameters predicted by DNN
!RESI 5DED          0.000 ! 5-decene: charges adjused so central -CH2-CH=CH-CH2- is neutral. To neutralize adjacent -CH2- carbon adjusted
!RESI DMPN to DMP including update of atom types and parameters to merge nucleic acid and lipid DMP models
! new anionic model compounds and mulitple updates to all anionic model compounds atom types, electrostatic and bonded parameters.
!RESI SO4          -2.000 ! Sulfate ion
!RESI MSO4         -1.000 ! Methylsulfate
!RESI MSNA         -1.000 ! Methylsulfonate
!RESI NMSM         -1.000 ! N-methyl sulfamate 
!RESI NESM         -1.000 ! N-Ethyl sulfamate
!RESI MEO          -1.000 ! Methoxide
!RESI ETO          -1.000 ! Ethoxide
!misc house cleaning of files

##lipids
toppar_drude_lipid_2023a.str
!update of lipids to account for LJ PME as published in JCTC, Yu et al. 2023, File identical to 2022d version

## nucleic acids
toppar_drude_nucleic_acid_2017d.str
!this version is based on toppar_drude_nucleic_acid_2017c.str with pre2020 atom types
!and parameters added to allow for its use with toppar_drude_main_protein_2023a.str
!due to the removal of 3 atom types associated with the conflict between the lipid
!and nucleic acid ester and phosphate groups completed in 2022

!As of July 2023 additional optimization of the nucleic acids was ongoing.
!for Nucleic acids continue to use 2017c version with master_2023a 





