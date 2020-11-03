
library('tidyverse')
library('dagitty')

# # create a new temporary "package" called "project_functions" containing functions and other objects that you don't want to appear in the global environment
# # the parent environment is the "stats" package
# if("project_functions" %in% search()) detach(project_functions)
# 
# with((project_functions <- new.env(parent=as.environment("package:stats"))),
#   {
  #list project-specific functions here:
  
    #"not in" infix function
    `%ni%` <- Negate(`%in%`)
    
    #get all functions used in a formula, and exclude variables
    all.funs = function(expr){all.names(expr, unique=TRUE)[!(all.names(expr, unique=TRUE) %in% all.vars(expr))]}
    
    node <- function(variable_formula, missing_rate=~0, keep=TRUE){
      
      stopifnot(rlang::is_formula(variable_formula))
      stopifnot(rlang::is_formula(missing_rate))
      stopifnot(is.logical(keep))
      
      l <- list(
        variable_formula = variable_formula,
        missing_rate = missing_rate,
        keep = keep
      )
      
      class(l) <- append(class(l), "node")
      l 
      
    }
    
    pdag_create <- function(list){
      # converts list to data frame which is a bit easier to work with, and embellishes with some useful columns.
      # The function performs a few checks on the the list, for instance to make sure the p-DAG is indeed acyclic
      # and that variables using in the expressions are defined elsewhere.
      
      stopifnot("each element in 'list' must be a list of class 'node'" = all(sapply(list, function(x){"node" %in% class(x)})))
      
      df <- list %>%
        enframe(name="variable", value="list") %>% unnest_wider("list") %>%
        mutate(
          in_order = row_number(),
          dependencies = map(variable_formula, ~all.vars(.)),
          missing_formula = map(missing_rate, ~{
            rhs <- deparse1(rlang::f_rhs(.))
            fun <- as.formula(paste0("~rbernoulli(n=1, p=", rhs, ")"))
            fun
          }),
        )
      
      dagitty <- pdag_to_dagitty(df)
      
      stopifnot("graph is not acyclic" = isAcyclic(dagitty)) 
      stopifnot("not all dependencies are defined" = all(simplify(unique(flatten(df$dependencies))) %in% df$variable))
      
      df
    }
    
    
    pdag_to_dagitty <- function(df){
      # convert abn_df to a dagitty object
      
      df1 <- df %>% 
        mutate(
          dagitty_str = map2_chr(variable, dependencies, ~paste0(.x, " <- ", "{", paste0(.y, collapse=" ") , "}")),
        )
      dagitty(paste0("dag {", paste0(df1$dagitty_str, collapse=" "), "}"))
    }
    
    
    pdag_plot <- function(pdag){
      dagitty <- pdag_to_dagitty(pdag)
      plot(graphLayout(dagitty))
    }
    
    
    pdag_simulate <- function(df, pop_size, keep_all=FALSE){
      
      dagitty <- pdag_to_dagitty(df)
      
      df1 <- df %>% mutate(
        df_fun = pmap(lst(variable, variable_formula), function(variable, variable_formula){
          
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
      df_ordered <- dagitty %>%
        topologicalOrdering() %>%
        enframe(name="variable", value='topological_order') %>% 
        unnest(topological_order) %>% 
        left_join(df1, ., by='variable') %>% 
        arrange(topological_order)
      
      
      # simulate complete dataset (with a patient ID variable in the initiated dataset)
      tbl0 <- tibble(ptid = seq_len(pop_size))
      tblsim_complete <- compose(!!!df_ordered$df_fun, .dir='forward')(tbl0)
      
      
      missing_formula <- append(list(ptid =~ FALSE), setNames(df_ordered$missing_formula, df_ordered$variable))
      
      # make some values missing, according to missing_fun
      tblsim <- pmap_df(
        lst(tblsim_complete, missing_formula, simdat=list(tblsim_complete)), 
        function(tblsim_complete, missing_formula, simdat){
          NA_type_ <- NA
          mode(NA_type_) <- typeof(tblsim_complete)
          row_num <- seq_len(nrow(simdat))
          mask <- map_lgl(row_num, ~eval(rlang::f_rhs(missing_formula), simdat[.,]))
          if_else(mask, NA_type_, tblsim_complete)
        }
      )
      
      # choose which variables to return
      returnvars <- df1 %>% filter(keep | keep_all) %>% pluck("variable")
      
      tblsim %>% select(ptid, all_of(returnvars))
    }
    
#   }
# )
# attach(project_functions);rm(project_functions)