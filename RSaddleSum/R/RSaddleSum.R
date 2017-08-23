
saddleSum.init <- function(weights) 
{
	ans <- .Call("saddleSumCreate", as.numeric(weights))
        if (is.null(ans)) 
          stop("Could not create saddleSum object.")
        return(ans)
}

saddleSum.pvalue <- function(saddleSum.data, score, num.hits, cutoff.pval=1e-3, maxiter=50, tol=1e-10)
{
        if (is.null(attr(saddleSum.data,"saddleSum_ptr")))
          stop("Invalid saddleSum object.")
	ans <- .Call("saddleSumPvalue", saddleSum.data, as.numeric(score), as.integer(num.hits), as.numeric(cutoff.pval),
                     as.integer(maxiter), as.numeric(tol))
        if (is.null(ans)) 
          stop("Could not evaluate saddleSum pvalue.")
        return(ans)
}
