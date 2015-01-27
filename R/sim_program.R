#---------------------------------------------------------------
# SimProgram.R
#
#
# Author:   Paul R. Staab & Lisha Naduvilezhath
# Email:    staab (at) bio.lmu.de
# Date:     2012-10-05
# Licence:  GPLv3 or later
#----------------------------------------------------------------

# Keep a user modifiable list of availible simulation programs in a private
# enviroment (jaatha's own env is read-only after package load)
if (!exists(".jaatha")) .jaatha <- new.env()
if (!exists("sim_progs", envir=.jaatha)) .jaatha$sim_progs <- list()

createSimProgram <- function(name, possible_features, possible_sum_stats,
                             sim_func=NULL, finalization_func=NULL,
                             print_cmd_func=NULL,
                             priority=50) {

  # Create the simulation program
  .jaatha$sim_progs[[name]] = list(name=name,
                                   possible_features=possible_features,
                                   possible_sum_stats=possible_sum_stats,
                                   sim_func=sim_func,
                                   finalization_func=finalization_func,
                                   print_cmd_func=print_cmd_func,
                                   priority=priority)
}

getSimProgram <- function(name) {
  return(.jaatha$sim_progs[[name]])
}
