#' Add migration/gene flow between two populations to a demographic model
#'
#' This function adds the assumption to the model that some individuals
#' 'migrate' from one sub-population to another, i.e. they leave the one
#' and become a member of the other. This is usually used to model ongoing
#' gene flow through hybridisation after the populations separated. 
#'
#' You can enter a time ('time.start') at which the migration is 
#' assumed to start (looking backwards in time). From that time on, a 
#' fixed number of migrants move from population 'pop.from' to
#' population 'pop.to' each generation. This number is given via this 
#' feature's parameter, which equals 4*N0*m,  where m is the 
#' fraction of 'pop.to' that is replaced with migrants each generation. 
#' If 'pop.to' has also size Ne, than this is just the
#' expected number of individuals that migrate each generation.
#' 
#' You can add different mutation rates at different times to your model.
#' Then each rate will be used for the period from its time point to
#' the next. Migration from and to an population always ends with the 
#' speciation event in which the population is created.
#'
#' @param rate  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "M" will use
#'            an parameter with name M that you have previously 
#'            created. You can also use R expression here, so "2*M"
#'            or "5*M+2*tau" (if tau is another parameter) will also
#'            work (also this does not make much sense).
#' @param pop.from The population from which the individuals leave.
#' @param pop.to The population to which the individuals move.
#' @param time.start The time point at which the migration with this rate
#'            starts.
#' @return    The demographic model with migration
#' @export
#'
#' @examples
#' # Asymmetric migration for two populations
#' dm <- dm.createDemographicModel(c(25,25), 100) +
#'   feat_migration(par_const(0.5), 1, 2) +
#'   feat_migration(par_const(0.75), 2, 1)
#'   
#' # Symmetric Migration
#' dm <- dm.createDemographicModel(c(25,25), 100) +
#'   feat_migration(par_range('m', 0.1, 2), symmetric=TRUE)
feat_migration <- function(rate, pop_from, pop_to, 
                           symmetric=FALSE, time_start="0") {
  if (symmetric) {
    return(Feature$new("migration_sym", rate, time_point=time_start))
  } else {
    return(Feature$new("migration", rate, pop_source=pop_from, pop_sink=pop_to, 
                       time_point=time_start))
  }
}