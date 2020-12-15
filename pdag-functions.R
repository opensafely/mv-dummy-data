
library('tidyverse')
library('dagitty')


#"not in" infix function
`%ni%` <- Negate(`%in%`)

#get all functions used in a formula, and exclude variables
all.funs = function(expr){all.names(expr, unique=TRUE)[!(all.names(expr, unique=TRUE) %in% all.vars(expr))]}

node <- function(variable_formula, missing_rate=~0, keep=TRUE, known=FALSE){
  
  stopifnot(rlang::is_formula(variable_formula))
  stopifnot(rlang::is_formula(missing_rate))
  stopifnot(is.logical(keep))
  
  l <- list(
    variable_formula = variable_formula,
    missing_rate = missing_rate,
    keep = keep,
    known = known
  )
  
  class(l) <- append(class(l), "node")
  l 
  
}


pdag_create <- function(list, known_variables=NULL){
  # converts list to data frame which is a bit easier to work with, and embellishes with some useful columns.
  # The function performs a few checks on the the list, for instance to make sure the p-DAG is indeed acyclic
  # and that variables using in the expressions are defined elsewhere.
  
  # The known_variables argument is for passing a character vector of variables names
  # for variables that are already defined externally in a 
  # given dataset, which can be passed to pdag_simulate
  
  # whilst the variable_formula is the variable name itself, this is to help with the pdag_simulate function
  # it doesn't actually lead to self-dependence (eg var depends on var)
  
  
  stopifnot("'list' must be a list where each element is an object of class 'node'" = all(sapply(list, function(x){"node" %in% class(x)})))
  
  
  df <- list %>%
    enframe(name="variable", value="list") %>% unnest_wider("list") %>%
    mutate(
      dependencies = map(variable_formula, ~all.vars(.)),
      missing_formula = map(missing_rate, ~{
        rhs <- deparse1(rlang::f_rhs(.))
        fun <- as.formula(paste0("~rbernoulli(n=1, p=", rhs, ")"))
        fun
      }),
      known = FALSE,
    )
  
  if(!is.null(known_variables)){
    
    df_available <-
      tibble(
        variable = known_variables,
        variable_formula = list(~as.formula(paste0("~",variable))),
        missing_rate = list(~0),
        keep = TRUE,
        dependencies = list(character()),
        missing_formula = list(as.formula("~rbernoulli(n=1, p=0)")),
        known = TRUE,
      )
    
    df <- bind_rows(df_available, df)
  }
  
  df <- df %>% mutate(in_order = row_number())
  
  dagitty <- pdag_to_dagitty(df)
  
  stopifnot("graph is not acyclic" = isAcyclic(dagitty)) 
  stopifnot("not all dependencies are defined" = all(simplify(unique(flatten(df$dependencies))) %in% df$variable))
  stopifnot("variable names are not unique" = length(df$variable) == length(unique(df$variable)))
  
  df
}


pdag_to_dagitty <- function(pdag){
  # convert pdag_df to a dagitty object
  
  pdag1 <- pdag %>% 
    mutate(
      dagitty_str = map2_chr(variable, dependencies, ~paste0(.x, " <- ", "{", paste0(.y, collapse=" ") , "}")),
    )
  dagitty(paste0("dag {", paste0(pdag1$dagitty_str, collapse=" "), "}"))
}


pdag_plot <- function(pdag){
  dagitty <- pdag_to_dagitty(pdag)
  plot(graphLayout(dagitty))
}


pdag_simulate <- function(pdag, data=NULL, pop_size, keep_all=FALSE){
  
  dagitty <- pdag_to_dagitty(pdag)
  
  pdag1 <- pdag %>% mutate(
    pdag_fun = pmap(lst(variable, variable_formula), function(variable, variable_formula){
      
      function(tib){
        row_num <- seq_len(nrow(tib))
        x <- simplify(map(row_num,  ~eval(rlang::f_rhs(variable_formula), tib[.,])))
        tib1 <- tib
        tib1[variable] <- x
        tib1
      }
    }),
    
  )
  
  #reorder based on dependencies so that simulation will create variables in the right order
  pdag_ordered <- dagitty %>%
    topologicalOrdering() %>%
    enframe(name="variable", value='topological_order') %>% 
    unnest(topological_order) %>% 
    left_join(pdag1, ., by='variable') %>% 
    arrange(topological_order)
  
  
  # simulate complete dataset (with a patient ID variable in the initiated dataset)
  if (is.null(data))
    tbl0 <- tibble(ptid = seq_len(pop_size))
  else
    tbl0 <- data
  
  pdag_ordered_unknown <- pdag_ordered %>% filter(!known)
  
  tblsim_complete <- compose(!!!pdag_ordered_unknown$pdag_fun, .dir='forward')(tbl0)
  
  missing_formula <-  setNames(pdag_ordered_unknown$missing_formula, pdag_ordered_unknown$variable)
  
  # make some values missing, according to missing_fun
  tblsim <- pmap_df(
    lst(tblsim_complete = tblsim_complete[pdag_ordered_unknown$variable], missing_formula, simdat=list(tblsim_complete)), 
    function(tblsim_complete, missing_formula, simdat){
      NA_type_ <- NA
      mode(NA_type_) <- typeof(tblsim_complete)
      row_num <- seq_len(nrow(simdat))
      mask <- map_lgl(row_num, ~eval(rlang::f_rhs(missing_formula), simdat[.,]))
      if_else(mask, NA_type_, tblsim_complete)
    }
  )
  
  tblsim <- bind_cols(tbl0, tblsim)
  
  # choose which variables to return
  returnvars <- pdag1 %>% filter(keep | keep_all) %>% pluck("variable")
  
  tblsim %>% select(names(tbl0), all_of(returnvars))
}

