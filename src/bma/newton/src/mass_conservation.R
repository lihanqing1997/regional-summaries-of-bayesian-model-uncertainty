# The switch permits matched validation without changing the posterior or budget.
cfg$mass_conservation <- TRUE
cfg$mass_conservation_passes <- 4L

ars_conserve_mass <- function(nodes) {
  if(!isTRUE(cfg$mass_conservation))return(nodes)
  started <- proc.time()[["elapsed"]]
  ids <- names(nodes)
  parents <- which(vapply(nodes,function(node)
    length(node$children)>0L && all(node$children %in% ids),logical(1)))
  if(!length(parents))return(nodes)
  children <- lapply(nodes[parents],function(node)match(node$children,ids))
  lower <- vapply(nodes,`[[`,numeric(1),"log_lower")
  upper <- vapply(nodes,`[[`,numeric(1),"log_upper")
  fit <- ars_propagate_log_mass(lower,upper,parents,children,cfg$mass_conservation_passes)
  changed <- which(fit$log_lower>lower | fit$log_upper<upper)
  for(i in changed) {
    nodes[[i]]$log_lower <- fit$log_lower[i]
    nodes[[i]]$log_upper <- fit$log_upper[i]
  }
  if(!is.null(ars_tracker)) {
    previous <- ars_tracker$mass_conservation
    if(is.null(previous))previous <- list(calls=0L,changed_node_updates=0L,seconds=0)
    previous$calls <- previous$calls+1L
    previous$changed_node_updates <- previous$changed_node_updates+length(changed)
    previous$seconds <- previous$seconds+proc.time()[["elapsed"]]-started
    ars_tracker$mass_conservation <- previous
  }
  nodes
}

# Keep the vendor search unchanged. Insert propagation before both state solves
# in a copied function, so the new endpoints remain in its local node history.
ars_add_mass_conservation <- function(engine) {
  replacements <- 0L
  visit <- function(expr) {
    if(!is.call(expr))return(expr)
    if(identical(expr[[1]],as.name("rrs_tree_epsilon_search_state"))) {
      replacements <<- replacements+1L
      return(substitute({nodes <- ars_conserve_mass(nodes); CALL},list(CALL=expr)))
    }
    for(i in seq_along(expr))if(is.call(expr[[i]]))expr[[i]] <- visit(expr[[i]])
    expr
  }
  body(engine) <- visit(body(engine))
  stopifnot(replacements==2L)
  engine
}
