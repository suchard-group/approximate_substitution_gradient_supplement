formatAlnForXML <- function(char.aln) {
  str <- c(
    "  <!-- The sequence alignment (each sequence refers to a taxon above).         -->\n",
    paste0("  <!-- ntax=",dim(char.aln)[1]," nchar=",dim(char.aln)[2],"                                                    -->\n"),
    paste0('  <alignment id="',"alignment",'" dataType="nucleotide">\n')
  )
  
  for (i in 1:dim(char.aln)[1]) {
    str <- c(str,
      "    <sequence>\n",
      paste0('      <taxon idref="',row.names(char.aln)[i],'"/>\n'),
      paste0("      ",paste0(toupper(char.aln[i,]),collapse=""),"\n"),
      "    </sequence>\n"
    )
  }
  
  str <- c(str,"  </alignment>\n")
  
  return(str)
}
