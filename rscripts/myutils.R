
library(tools)

######################################################################
# getBasename
######################################################################
getBasename <- function( filename ){
    suffix <- file_ext(filename)
    path <- dirname(filename)
    #basename <- basename(gsub( "(.*)\\.(.*)", "\\1", filename))
    basename <- basename( file_path_sans_ext(filename) )
    lst <- c(path, basename, suffix)
    return(lst)
}# END getBasename

######################################################################
# readMultiTable
######################################################################
readMultiTable <- function( filename ){
    x <- readLines( filename )
    x <- x[1:15] # This is a hack! delete if you can!
    start <- grep("^[/]+$", x)
    mark <- vector('integer', length(x))
    mark[start] <- 1
    # determine limits of each table
    mark <- cumsum(mark)
    # split the data for reading
    df <- lapply(split(x, mark), function(.data){
        .input <- read.table(textConnection(.data), skip=2, header=FALSE, sep="\t")
        attr(.input, 'name') <- .data[2]  # save the name
        .input
    })
    # rename the list
    names(df) <- sapply(df, attr, 'name')
    return(df)
} #END readMultiTable

