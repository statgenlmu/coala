# Keep a user modifiable list of availible simulation programs in a private
# enviroment
if (!exists("simprograms")) simprograms <- new.env()

create_simprog <- function(name, possible_features, possible_sum_stats,
                             sim_func=NULL, print_cmd_func=NULL,
                             priority=50) {

  # Create the simulation program
  simprograms[[name]] <- list(name=name,
                               possible_features=possible_features,
                               possible_sum_stats=possible_sum_stats,
                               sim_func=sim_func,
                               print_cmd_func=print_cmd_func,
                               priority=priority)
}

get_simprog <- function(name) {
  return(simprograms[[name]])
}
