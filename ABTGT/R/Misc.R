getclass<-function (x, classy) inherits(get(x), classy)

findy<-function (classy){
  return(unique(c(ls("package:ABTGT")[unlist(lapply(ls("package:ABTGT"),getclass, classy = classy))], ls(envir = .GlobalEnv)[unlist(lapply(ls(envir = .GlobalEnv),getclass, classy = classy))])))
}
