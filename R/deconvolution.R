# Deconvolution Strategy Pattern
# Ingests a request to run deconvolution, returns a normalized result regardless of strategy
# @author jbard


run_deconvolution <- function(strategy,){
  if(strategy="bayesprism"){
    run_bayesprism()
  } else if (strategy=="music") {
    run_music()
  } else if (strategy=="gedit"){
    run_gedit()
  } else {
    warning("Invalid deconvolution strategy specified")
  }
}

run_music <-function(){}
run_bayesprism <-function(){}
run_gedit <-function(){}



