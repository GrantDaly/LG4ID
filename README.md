# LG4ID
Long-G4 (Guanine Quadruplex) Identifier

Program associated with Williams J.D. et al., 2020.

Jonathan D Williams, Dominika Houserova, Bradley R Johnson, Brad Dyniewski, Alexandra Berroyer, Hannah French, Addison A Barchie, Dakota D Bilbrey, Jeffrey D Demeis, Kanesha R Ghee, Alexandra G Hughes, Naden W Kreitz, Cameron H McInnis, Susanna C Pudner, Monica N Reeves, Ashlyn N Stahly, Ana Turcu, Brianna C Watters, Grant T Daly, Raymond J Langley, Mark N Gillespie, Aishwarya Prakash, Erik D Larson, Mohan V Kasukurthi, Jingshan Huang, Sue Jinks-Robertson, Glen M Borchert, Characterization of long G4-rich enhancer-associated genomic regions engaging in a novel loop:loop ‘G4 Kissing’ interaction, Nucleic Acids Research, Volume 48, Issue 11, 19 June 2020, Pages 5907–5925, [https://doi.org/10.1093/nar/gkaa357](https://doi.org/10.1093/nar/gkaa357)

# LG4ID Summary
>To identify a panel of loci containing extensive G4 sequence motifs, we searched for large genomic stretches significantly enriched for guanine-triplets (G-triplets) instead of focusing on shorter, more rigidly defined G4 motifs. G-triplets were counted because they are the basic sequence necessary for G4 structure formation. Further search parameters (e.g. window size) were based on the immunoglobulin switch region Sμ (Supplementary Information 1), which is a G4-forming recombination site recognized by mismatch repair factors (24,52). The Sμ guanine density of 120 G-triplets/1.5 kb window, which is a much lower density of G-triplets compared to other well-known G4s such as telomeres (65), was used to train our analyses. Modeling our LG4ID search program on these parameters, we identified 301 loci containing a density of at least 80 GGG repeats/kb (Figure 2A). The 301 long G4-capable regions (LG4s) we identified in the human genome ranged from 199 to 4973 bp in length (subset shown in Figure 2A). Although the initial search window was 1.5 kb, several smaller length regions contained a high density of G-triplets surrounding the larger repetitive unit and, therefore, met our minimal G-triplet requirements.

>Python 3-based program for identification of LG4s–LG4ID
In order to identify long G4-capable regions (LG4s) present in the human genome, we wrote a Python 3 program to search a FASTA formatted sequence file for long G-triplet regions likely to form G4. The program identifies LG4 motifs based on the density of G-triplets within 1.5 kb sequence windows sliding one base pair per iteration and does not account for loop length. In order to define the minimal density of G-triplets needed to call a LG4, we modeled our search program after the G4-capable sequence density within the human immunoglobulin mu (Sμ) switch region (11–13,24,36). Sliding windows were applied, with a minimal output threshold of (GGG) X 121 for every 1.5 kb of sequence. To identify G-rich sequences on the + and – stands of the genome, CCC density was also determined as above then combined with the GGG search output data (for more details see Supplementary Information 1). LG4ID was used to identify all LG4s (with parameters stated above) in the human genome (hg38) after which each LG4’s location and identity were confirmed manually.

# Supplemental Methods
>LG4 (inputString, windowSize=1500, mustExceed=120)
LG4ID finds all non-overlapping "GGG" and "CCC" in inputString, and records the index of
each "GGG" and "CCC" into "GGG Boolean List" and "CCC Boolean List", respectively. To
clarify "non-overlapping", "GGGGGG" is two "GGG", not four "GGG". "GGG Boolean List"
and "CCC Boolean List" are both of length(inputString) - 2, and used to reference all "GGG"
and "CCC" in the "inputString". Such Boolean lists are used to obtain linear search-time
complexity.

>Next, LG4ID iterates through "inputString“ (Boolean Lists in this case) in "windowSize"
overlapping bins, advancing by one index at a time. This guarantees that all LG4s will be
found, as opposed to non-overlapping bin approaches. For each bin LG4ID counts “GGG”s
and “CCC”s separately and obtains the maximum value among those, which is called the
"maximum Hits". The "maximum Hits" is calculated by only one strand and is thus a
pessimistic measure. If "maximum Hits" exceeds "mustExceed“ (‘120’) then the bin is a valid
LG4, otherwise it is not. All valid LG4 bins which overlap or neighbor by start and stop
coordinates are joined into contiguous LG4s, which are returned to the user.
