#-------------------------------------------------------------------------------
#
#           sim2DSinSurf
#
#-------------------------------------------------------------------------------

#' simulated large spatial data sets from combination of sine wave surfaces with random amplitudes, freqencies, and spatial shifts.
#'
#' simulated large spatial data sets from combination of sine wave surfaces with random amplitudes, freqencies, and spatial shifts.
#'
#' @param x a vector of x-coordinates
#' @param y a vector of y-coordinates
#'
#' @return data.frame of simulated surface (z) column.
#'
#' @author Jay Ver Hoef
#' @export

sim2DSinSurf = function(x,y) 
{
    
	rottheta = runif(1)*pi
	newxy = matrix(c(cos(rottheta),-sin(rottheta), sin(rottheta), 
				cos(rottheta)), nrow = 2) %*% rbind(x,y)
	xr = newxy[1,]
	yr = newxy[2,]
	z = sin(runif(1)*2*pi*(xr + runif(1)*pi)) + 
			sin(runif(1)*2*pi*(yr + runif(1)*pi))
			
	for(i in 1:99) {
			rottheta = runif(1)*pi
			newxy = matrix(c(cos(rottheta),-sin(rottheta), sin(rottheta), 
				cos(rottheta)), nrow = 2) %*% rbind(x,y)
			xr = newxy[1,]
			yr = newxy[2,]
			z = z + runif(1)*(1 - i/100)*(sin(runif(1)*2*i*pi*(xr + runif(1)*pi)) + 
				sin(runif(1)*2*i*pi*(yr + runif(1)*pi)))
	}

	data.frame(x = x,  y = y, z = z)
}

