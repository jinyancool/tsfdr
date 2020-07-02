installmosek <- function(targetDir = "~/bin"){
  if(!dir.exists(targetDir)){
    dir.create(targetDir)
  }
  if(!"Rmosek" %in% installed.packages()[,"Package"]) {
    install.packages(new.packages)
  }
  switch(Sys.info()[['sysname']],
         Windows= {ostype <- "windows"},
         Linux  = {ostype <- "linux"},
         Darwin = {ostype <- "osx"})
  mosekUrl <- sprintf("https://download.mosek.com/stable/9.2.14/mosektools%s64x86.tar.bz2", ostype)
  download.file(mosekUrl, destfile = sprintf("%s/mosektools%s64x86.tar.bz2", targetDir, ostype))
  cmd <- sprintf("cd %s && tar -xvf mosektools%s64x86.tar.bz2", targetDir, ostype)
  system(cmd)
  Rmosek::mosek_attachbuilder(sprintf("%s/mosek/9.2/tools/platform/%s64x86/bin", targetDir, ostype))
  install.rmosek()
  if(!file.exists("~/mosek/mosek.lic")){
    message("Please visit: https://www.mosek.com/products/academic-licenses/ to get the licsense.
             Check email, put licsense file at: ~/mosek/mosek.lic")
  }
}
