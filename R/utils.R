#' Calculates the depth of a list
#' @seealso https://stackoverflow.com/a/13433689/1987629 (Function copied from)
#' @param this A list whose depth needs to be calculated
#' @param thisdepth The default value if the list is flat
#' @return An integer specifying the depth of a list
#' @examples
#' depth(list(list(1, 2), list(2, 3))) # 3
#' depth(1:5) # 1
#' depth(list(list(1, 2), list(2, list(3, 4)))) # 4
#' depth(list(1:2, 3:4)) # 2
#' @export
depth <- function(this, thisdepth=1){
  if(!is.list(this)){
    return(thisdepth)
  } else {
    return(max(unlist(
      lapply(this, depth, thisdepth = thisdepth + 1))))
  }
}
