# downloads the 5.7.1 version of ASTRAL to current working directory and
# sets the global 'astralPath' option to point to it
# this avoids the user needing to download it and essentially makes java the
# only external dependency

#' Download Astral
#'
#' @description
#' Retrieves the Astral v5.7.1 .jar file from github and sets a global option
#' called 'astralPath' to a be the location of that file. Note that astralPath
#' can also be passed to the callAstral function manually if you'd prefer to use your
#' own version of Astral.
#'
#' @param install.dir Path to the directory where Astral will be downloaded to. By default, this is the current working directory.
#' @param method Method is an optional argument passed to [utils::download.file()]; see [utils::download.file()] for more info.
#'
#' @returns The function returns nothing (void) and is called entirely to download Astral and set the astralPath global option.
#'
#' @seealso [UCETools::callAstral]
#'
#' @importFrom utils download.file unzip
#'
#' @export

getAstral <- function(install.dir = getwd(), method){
    ## Basic input checks
    if(!dir.exists(install.dir)){
        stop("Install path does not exist or is not a directory.")
    }
    # method is an optional argument for download.file(); see download.file() docs
    if(missing(method)){
        method = "auto"
    }

    ## Get ASTRAL v5.7.1 from github and decompress it
    exitStatus <- download.file(url = "https://github.com/smirarab/ASTRAL/archive/refs/tags/v5.7.1.zip",
                                method = method,
                                destfile = file.path(install.dir,"ASTRAL.zip"))
    # check if download failed
    if(exitStatus != 0){
        stop("Download finished with non-zero exit status. Try again, or download ASTRAL and set astralPath manually.")
    }
    unzip(zipfile = file.path(install.dir,"ASTRAL.zip"),
          overwrite = TRUE,
          exdir = install.dir)
    # If I don't sleep the script for 1 second the file.remove and 2nd unzip don't run #JustRThings
    Sys.sleep(1)
    file.remove(file.path(install.dir,"ASTRAL.zip"))
    # Why is there a 2nd zip file nested in this???
    # Why does it have the actual .jar file that I need???????
    unzip(zipfile = file.path(install.dir,"ASTRAL-5.7.1","Astral.5.7.1.zip"),
          overwrite = TRUE,
          exdir = file.path(install.dir,"ASTRAL-5.7.1"))
    cat(paste("Downloaded ASTRAL v5.7.1 to ", install.dir, "\n\n"))

    ## Set global astralPath var
    # Had to add 'normalizePath()' here because of inconsistency between file.path(), getwd(), and tempdir()
    # This seems to be new since I wrote this two years ago; I think file.path() behaviour changed on Windows? 2024-06-21
    options(astralPath = normalizePath(file.path(install.dir,"ASTRAL-5.7.1","Astral","astral.5.7.1.jar")))
    cat(paste("Global 'astralPath' option set to ", getOption("astralPath"), "\n"))
}
