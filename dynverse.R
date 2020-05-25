library(tidyverse)
library(dyngen)

params <- simple_params
options(ncores = 1)
model <- invoke(generate_model_from_modulenet, params$model)
plot_net(model)


simulation <- invoke(simulate_multiple, params$simulation, model$system)
plot_simulation_space_time(simulation)

gs <- invoke(extract_goldstandard, params$gs, simulation, model)
experiment <- invoke(run_experiment, params$experiment, simulation, gs)
normalisation <- invoke(dynnormaliser::normalise_filter_counts, 
                        params$normalisation, experiment$counts)
task <- wrap_dyngen_dataset("readme_dataset", params, model, simulation, 
                            gs, experiment, normalisation)


#simulate_table