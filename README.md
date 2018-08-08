# Solanum-Plastid-Assembly

Three scripts were used to create a draft assembly and correct any SNPs, short InDels, or Long InDels (soft clips).

Solanum_plastid_assembly_1.sh has a reference plastid and fastq reads as starting inputs. The fasta outputs need to be manually assembled into a draft plastid sequencing using overlapping fragments and orientation determined by NUCmer.

Solanum_plastid_assembly_2.sh has the draft plastid (only 1 IR) and raw reads as the inputs. Also, the length of each region (LSC, IR, SSC) needs to be filled in. The output is a vcf file that contains SNPs and short InDels that are errors in the draft assembly.

Solanum_plastid_assembly_3.sh takes the output vcf file from Solanum_plastid_assembly_2.sh and corrects for the errors. Then it reruns the SNP and short InDel caller, and also checks for softclips. The outputs may be empty, but if there are still mistakes in the draft assembly, they must be corrected for and this script should be re-ran until there are no more mistakes in the assembly. The final product is a complete genome, error free.
