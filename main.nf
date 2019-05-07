#!/usr/bin/env nextflow

// Script to run scan2 and scan2000 permutations for linkagemapping traits

// params.in = riail phenotype for just 1 trait
params.nperm = 1000
params.set = 2
params.cross = "N2xCB4856cross_full"
nperms = params.nperm
params.out = '.'

// generate seeds random numbers
seeds = Channel
			.from(1..100000)
			.randomSample(nperms)

// riail phenotype channel
riails = Channel.fromPath("${params.in}")

// do the scan2
process scan2 {
	publishDir params.out, mode: 'copy'

	input:
		file("pheno") from riails

	output:
		file("scan2.Rda")
		file("mapcross.Rda") into crossobj

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

	"""
}

// save mapcross for each seed
crossobj_split = seeds.combine(crossobj)

// do the permutations
process scan2000 {
	cpus 4
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
        file("scantwothousand.tsv")

    """
    # use this to only print the header of the first line
	awk 'FNR>1 || NR==1' ${scans.join(" ")} > scantwothousand.tsv
    """

}
