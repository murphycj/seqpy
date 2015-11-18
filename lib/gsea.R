
parse_gmt <- function(gmt) {
	#
	# parse a MSigDB GMT file
	#

	fin <- file(gmt,open="r")
	pathways = list()
	lines <- readLines(fin)
	for (i in 1:length(lines)) {
		temp = strsplit(lines[i],"\t")
		pathways[[temp[[1]][1]]] = temp[[1]][-c(1,2)]
	}
	close(fin)

	return(pathways)
}
