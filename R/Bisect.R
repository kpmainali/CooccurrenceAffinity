#' Bisections for finding a root of a function
#'
#' Find a root of a function by the method of Bisections
#'
#' @details   This function finds the root of the increasing function ffnc over the scalar interval intrv by the Method of Bisections. The function must be increasing but need not be smooth, and it must have a negative sign (value less than -tol) at the left endpoint of  intrv  and positive sign (value greater than tol) at the right endpoint. The method of Bisection is used in successive iterations to successively halve the width of the interval in which the root lies.
#'
#' @param ffnc an increasing function of a single scalar argument
#' @param intrv an interval over which the root of ffnc is sought
#' @param tol a tolerance determining when the successive bisections of the interval within which the root will lie have become small enough to stop
#'
#' @return This function returns a vector consisting of two numbers. The first named root is an estimate of the root x  solving ffnc(x) = 0, valid within an error of tol.
#' The second output vector element named fval is the value of the function ffnc at root.
#' It should be very close to 0 unless the function happens to jump from a value less than 0 to a value greater than 0 at  root.
#'
#' @author Eric Slud
#'
#' @references to be added
#'
#' @example
#' inst/examples/Bisect_example.R
#'
#' @export

Bisect <-
  function(ffnc, intrv, tol=1e-8) {
    #  root of increasing possibly non-smooth function by bisections
    if(ffnc(intrv[1])>= -tol | ffnc(intrv[2])<= tol)
      return("Same sign at endpts!")
    while(abs(intrv[2]-intrv[1]) > tol) {
      newp = (intrv[2]+intrv[1])/2
      newv = ffnc(newp)
      if(newv<0) intrv[1]=newp else intrv[2]=newp }
    newv=ffnc((intrv[2]+intrv[1])/2)
    c(root=(intrv[1]+intrv[2])/2, fval=newv) }
