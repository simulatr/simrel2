#' Parse a character string to a list
#' @description These function parsers characters into list and
#'   list into characters. They help to parse the argument and makes the
#'   user input easier and intutive.
#' @param chr A character string that represents a list or a vector
#'   The default is null. See notes for details.
#' @return A list or a vector
#' @note A semicolon in the string separates the elements of a list while
#'   A comma in the string separates the elements of a vector. See example.
#' @examples
#' chr2list("1, 2, 3") # c(1, 2, 3)
#' chr2list("1, 2; 3, 4") # list(c(1, 2), c(3, 4))
#' chr2list("1:2; 3:7") # list(1:2, 3:7)
#' chr2list("1, 2, 3;") # list(c(1, 2, 3))
#' @export
chr2list <- function(chr = NULL) {
  if (is.null(chr)) return(NULL)
  x <- unlist(strsplit(chr, ";"))
  x <- gsub("[[:space:]]", "", x)
  ret <- lapply(x, function(y) {
    if(grepl(":", y)) {
      eval(parse(text = y))
    } else {
      as.numeric(unlist(strsplit(y, ",")))
    }
  })
  if (!grepl(";", chr))
    ret <- ret[[1]]
  return(ret)
}

#' Parse List/Vector to a character string
#' @param list A list or a vector which need to be deparsed as a string
#'   The default is NULL and it will return NULL
#' @return A character string representing the list
#' @note This is the reverse for of @seealso chr2list
#' @examples
#' list2chr(c(1, 2, 3)) # "1; 2; 3"
#' list2chr(list(c(1, 2), c(3, 4))) # "1, 2; 3, 4"
#' list2chr(list(1:2, 3:7)) # "1, 2; 3, 4, 5, 6, 7"
#' @export
list2chr <- function(list = NULL) {
  if (is.null(list)) return(NULL)
  depth <- depth(list)
  if (depth > 2) {
    stop("List should not nest more than 2 level deep.")
  } else {
    inner_chr <- sapply(list, paste0, collapse = ", ")
    ret <- paste(inner_chr, collapse = "; ")
  }
  return(ret)
}

