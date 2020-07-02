installmosek <- function(){
  targetDir <- "~/bin"
  if(!dir.exists(targetDir)){
    dir.create(targetDir)
  }
  if(!"Rmosek" %in% installed.packages()[,"Package"]) {
    install.packages(new.packages)
  }
  os <- Sys.info()[['sysname']]
  if (os == "Darwin"){
    mosekUrl <- "https://download.mosek.com/stable/9.2.14/mosektoolsosx64x86.tar.bz2"
    download.file(mosekUrl, destfile = file.path(targetDir,"mosektoolsosx64x86.tar.bz2"))
    cmd <- sprintf("cd ~/bin && tar -xvf mosektoolsosx64x86.tar.bz2")
    system(cmd)
    Rmosek::mosek_attachbuilder("~/bin/mosek/9.2/tools/platform/osx64x86/bin")
  }
  if (os == "Linux"){
    mosekUrl <- "https://download.mosek.com/stable/9.2.14/mosektoolslinux64x86.tar.bz2"
    download.file(mosekUrl, destfile = file.path(targetDir,"mosektoolslinux64x86.tar.bz2"))
    cmd <- sprintf("cd ~/bin && tar -xvf mosektoolslinux64x86.tar.bz2")
    system(cmd)
    Rmosek::mosek_attachbuilder("~/bin/mosek/9.2/tools/platform/linux64x86/bin")
  }
  install.rmosek()
  if(!file.exists("~/mosek/mosek.lic")){
    message("Please visit: https://www.mosek.com/products/academic-licenses/ to get the licsense.
             Check email, put licsense file at: ~/mosek/mosek.lic")
  }
}
