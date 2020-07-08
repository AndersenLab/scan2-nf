#!/usr/bin/env nextflow

// Script to run scan2 and scan2000 permutations for linkagemapping traits

// params.in = riail phenotype for just 1 trait
date = new Date().format( 'yyyyMMdd' )

params.nperm = 1000
params.set = 2
params.cross = "N2xCB4856cross_full"
nperms = params.nperm
params.out = 'Analysis-${date}'

// generate seeds random numbers
seeds = Channel
			.from(1..100000)
			.randomSample(nperms)

// riail phenotype channel
riails = Channel.fromPath("${params.in}")

// do the scan2
process scan2 {
	cpus 4
	memory '20 GB'
	publishDir params.out, mode: 'copy'

	input:
		file("pheno") from riails

	output:
		file("scan2.Rda") into scan2_object
		file("mapcross.Rda") into crossobj
		file("scan2plot.png")


	"""
	#!/usr/bin/env Rscript --vanilla

	library(dplyr)
	library(readr)
	library(qtl)
	library(linkagemapping)

	#insert cross data
	crossobj <- get(linkagemapping::load_cross_obj("${params.cross}"))

	# read RIAIL phenotype data for 1 trait
	pheno <- readr::read_tsv("$pheno")

	#create completed cross object with pheno data set
	mapcross <- linkagemapping::mergepheno(crossobj, pheno, set = ${params.set})

	# save cross object
	save(mapcross, file = "mapcross.Rda")

	# run scan2
	scan2 <- qtl::scantwo(mapcross, pheno.col=3, method="mr")

	# make output into dataframe
	save(scan2, file = "scan2.Rda")

	# plot scan2
	png("scan2plot.png")
	plot(scan2)
	dev.off()

	"""
}

// save mapcross for each seed
crossobj_split = seeds.combine(crossobj)

// do the permutations
process scan2000 {
	cpus 4
	memory '20 GB'
	tag { s }

	input:
		set val("s"), file("mapcross") from crossobj_split

	output:
		file("scan2thousand*.tsv") into scan2000

	"""
	#!/usr/bin/env Rscript --vanilla

	library(dplyr)
	library(readr)
	library(qtl)
	library(linkagemapping)

	# load cross object
	load("$mapcross")

	# run scan2 with 5 permutations
	scan2thousand <- qtl::scantwo(mapcross, n.perm=1, pheno.col=3, method="mr")

	# make output into dataframe
	df <- data.frame(full = scan2thousand\$full[[1]], fv1 = scan2thousand\$fv1[[1]], int = scan2thousand\$int[[1]], 
						add = scan2thousand\$add[[1]], av1 = scan2thousand\$av1[[1]], one = scan2thousand\$one[[1]])

	# save dataframe
	readr::write_tsv(df, paste0("scan2thousand", "_", $s, ".tsv"))

	"""
}

// add all the tsv together 
process concatenate_quantiles {

    publishDir params.out, mode: 'copy'
    
    input:
		val("scans") from scan2000.toSortedList()

    output:
        file("scantwothousand.tsv") into perms

    """
    # use this to only print the header of the first line
	awk 'FNR>1 || NR==1' ${scans.join(" ")} > scantwothousand.tsv
    """

}

// summarize scan2 and add thresholds for significance from permutations
process summarize_scan2 {

	publishDir params.out, mode: "copy"

	input:
	file("scantwothousand") from perms
	file("scan2") from scan2_object

	output:
	file("scan2summary.tsv")


	"""
	#!/usr/bin/env Rscript --vanilla

	library(dplyr)
	library(readr)
	library(linkagemapping)
	library(qtl)

	#insert cross data
	crossobj <- get(linkagemapping::load_cross_obj("${params.cross}"))

	# load scan2
	load("$scan2")

	# load perms
	perms <- readr::read_tsv("$scantwothousand")

	# summarize scantwo and add GWER thresholds defined by permutations
	scan2_summary <- summary(scan2) %>%
	    dplyr::mutate(fv1_thresh = quantile(perms\$fv1, probs = 0.95),
	                  full_thresh = quantile(perms\$full, probs = 0.95),
	                  add_thresh = quantile(perms\$add, probs = 0.95),
	                  av1_thresh = quantile(perms\$av1, probs = 0.95),
	                  int_thresh = quantile(perms\$int, probs = 0.95)) %>%
	    dplyr::select(trait, chr1:int_thresh) %>%
	    dplyr::mutate(pos1f = as.character(pos1f),
	                  pos2f = as.character(pos2f),
	                  pos1a = as.character(pos1a),
	                  pos2a = as.character(pos2a))

	# riail marker conversion
	mappos <- qtl::pull.map(crossobj, as.table = TRUE) %>%
	    dplyr::mutate(marker = rownames(.),
	                  cM = as.character(pos)) %>%
	    dplyr::select(-pos, -chr) %>%
	    dplyr::distinct(cM, .keep_all = T) 

	# convert genetic pos to genomic pos
	scan2_summary <- scan2_summary %>%
	    # pos1f
	    dplyr::left_join(mappos, by = c("pos1f" = "cM")) %>%
	    dplyr::mutate(pos1f = as.numeric(stringr::str_split_fixed(marker, "_", 2)[,2])) %>%
	    dplyr::select(-marker) %>%
	    # pos2f
	    dplyr::left_join(mappos, by = c("pos2f" = "cM")) %>%
	    dplyr::mutate(pos2f = as.numeric(stringr::str_split_fixed(marker, "_", 2)[,2])) %>%
	    dplyr::select(-marker) %>%
	    # pos1a
	    dplyr::left_join(mappos, by = c("pos1a" = "cM")) %>%
	    dplyr::mutate(pos1a = as.numeric(stringr::str_split_fixed(marker, "_", 2)[,2])) %>%
	    dplyr::select(-marker) %>%
	    # pos2a
	    dplyr::left_join(mappos, by = c("pos2a" = "cM")) %>%
	    dplyr::mutate(pos2a = as.numeric(stringr::str_split_fixed(marker, "_", 2)[,2])) %>%
	    dplyr::select(-marker)

	readr::write_tsv(scan2_summary, "scan2summary.tsv")

	"""


}
