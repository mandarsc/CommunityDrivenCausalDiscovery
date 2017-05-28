# Author: Mandar Chaudhary (mschaudh@ncsu.edu)
# PhD student at NC State Univeristy, USA

library(bnlearn)
library(pcalg)
library(igraph)
library(Matrix)
library(corpcor)

# Generate random simulated DAGs and corresponding data sets

num_vars <- c(50, 100, 200, 500)

dag_data <- dags <- amat <- list()
cor_mat <- adj_nodes <- list()
dist_mat <- all_node_dist_mat <- list()


for(i in 1:length(num_vars[i])){
    
    # Select the sparsity level
    sparsity <- c('p', '0.1', '0.15', '0.2', '0.25')
    num_iter <- 20
    pc_cdcd_shd <- pc_shd <- cdcd_shd <- matrix(0, num_iter, 1)
    pc_tests <- matrix(0, num_iter, 20)
    cdcd_tests <- matrix(0, num_iter, 20)

    tpr_results <- matrix(0, num_iter, 2, dimnames=list(c(), c('CCD-TPR', 'PC-TPR')))
    fpr_results <- matrix(0, num_iter, 2, dimnames=list(c(), c('CCD-FPR', 'PC-FPR')))
    tdr_results <- matrix(0, num_iter, 2, dimnames=list(c(), c('CCD-TDR', 'PC-TDR')))

    # Set working directory containing simulated data sets and load the .RData files
    # setwd("/scratch/hdd1/home/mschaudh/Second_Chapter/Simulated_Data/")
    # load(paste('Vars_', num_vars[i], '_Sparsity_', sparsity, '_Data.RData', sep=''))
    
    for(j in 1:num_iter){
        print(cat('n=', num_vars[i], ' j= ', j))  
        test_counter <- 0

        res_obj <- marginalIndepTest(cor(dag_data[[j]]), nrow(dag_data[[j]]), test_counter)
        amat <- res_obj[[1]]
        test_counter <- res_obj[[2]]

        cor_mat <- cor(dag_data[[j]])
        
        p <- ncol(amat)
        n <- nrow(dag_data[[j]])
        vars <- 1:num_vars[i]
        
        pVals <- matrix(-Inf, nrow=p, ncol=p)
        G <- matrix(FALSE, nrow=p, ncol=p)
        idx <- which(amat!=0, arr.ind=TRUE)
        G[idx] <- TRUE

        sepset <- lapply(1:nrow(amat), function(.) vector("list", nrow(amat)))
        idx <- which(amat!=0, arr.ind=TRUE)
        amat[idx] <- cor_mat[idx]
        
        bool_amat <- amat
        bool_amat[bool_amat!=0] <- 1
        row.names(amat) <- colnames(amat) <- row.names(bool_amat) <- colnames(bool_amat) <- 1:num_vars[i]

        num_edge_tests <- array(0, 19)
        pc_num_edge_tests <- array(0, 20)
        edge_list <- list()
        pdsep_list <- vector("list", p)
        cdcd_obj <- ccd(G, pVals, cor(dag_data[[j]]), amat, bool_amat, sepset, edge_list, pdsep_list, 
                    nrow(dag_data[[j]]), 1:num_vars[i], num_edge_tests, test_counter)

        pc_obj <- pc(suffStat=list(C=cor(dag_data[[j]]), n=nrow(dag_data[[j]])), indepTest=gaussCItest, alpha=0.05, p=ncol(dag_data[[j]]), verbose=F, skel.method="stable")


        pc_cdcd_shd[j, ] <- shd(pc_obj, cdcd_obj)
        pc_shd[j, ] <- shd(dags[[j]], pc_obj)
        cdcd_shd[j, ] <- shd(dags[[j]], cdcd_obj)

        pc_num_edge_tests[1:length(pc_obj@n.edgetests)] <- pc_obj@n.edgetests
        pc_tests[j, 1:length(pc_num_edge_tests)] <- pc_num_edge_tests
        cdcd_tests[j, 1:length(cdcd_obj@n.edgetests)] <- cdcd_obj@n.edgetests

        tpr_results[j, 'CCD-TPR'] <- compareGraphs(cdcd_obj@graph, dags[[j]])['tpr']
        fpr_results[j, 'CCD-FPR'] <- compareGraphs(cdcd_obj@graph, dags[[j]])['fpr']
        tdr_results[j, 'CCD-TDR'] <- compareGraphs(cdcd_obj@graph, dags[[j]])['tdr']

        tpr_results[j, 'PC-TPR'] <- compareGraphs(pc_obj@graph, dags[[j]])['tpr']
        fpr_results[j, 'PC-FPR'] <- compareGraphs(pc_obj@graph, dags[[j]])['fpr']
        tdr_results[j, 'PC-TDR'] <- compareGraphs(pc_obj@graph, dags[[j]])['tdr']
    }
    # Set path to store results in a csv file
    # path <- "/scratch/hdd1/home/mschaudh/Second_Chapter/Results_AllComms/"
    # setwd(path)
    dir.create(paste(getwd(), '/Sparsity_', sparsity, sep=''))
    setwd(paste(getwd(), '/Sparsity_', sparsity, sep=''))
    write.csv(cbind(tpr_results, fpr_results, tdr_results), paste(getwd(), '/Metrics_', num_vars[i], '_n_50_21_40', '.csv', sep=''))
    write.csv(cbind(pc_cdcd_shd, pc_shd, cdcd_shd), paste(getwd(), '/SHD_', num_vars[i], '_n_50_21_40','.csv', sep=''))
    write.csv(cbind(pc_tests, cdcd_tests), paste(getwd(), '/CI_Tests_', num_vars[i], '_n_50_21_40','.csv', sep=''))
}

ccd <- function(G, pVals, cor_mat, amat, bool_amat, sepset, edge_list, pdsep_list, n, vars, num_edge_tests, test_counter){
    cl <- match.call()
    # res_obj <- buildSkeleton(G, pVals, cor_mat, amat, bool_amat, sepset, edge_list, pdsep_list, n, vars, num_edge_tests)
    res_obj <- buildSkeletonAllComms(G, pVals, NULL, cor_mat, amat, bool_amat, sepset, edge_list, pdsep_list, n, vars, num_edge_tests)
    G <- res_obj[[1]]
    sepset <- res_obj[[2]]
    pVals <- res_obj[[3]]
    amat <- res_obj[[4]]
    num_edge_tests <- res_obj[[8]]
    non_zero_idx <- which(num_edge_tests!=0)

    Gobject <- if (sum(G) == 0) {
        new("graphNEL", nodes = vars)
    }
    else {
        colnames(G) <- rownames(G) <- vars
        as(G, "graphNEL")
    }
    cdcd_skel <- new("pcAlgo", graph = Gobject, call = cl, n = integer(0), 
        max.ord = as.integer(length(which(num_edge_tests!=0))+1), n.edgetests = c(test_counter, num_edge_tests[non_zero_idx]), 
        sepset = sepset, pMax = pVals, zMin = matrix(NA, 1, 1))
    cdcd_skel@call <- cl
    cdcd_graph <- switch("relaxed", rand=udag2pdag(cdcd_skel), retry=udag2pdagSpecial(cdcd_skel)$pcObj, relaxed = udag2pdagRelaxed(cdcd_skel, verbose = FALSE, 
                solve.confl = FALSE))
    return(cdcd_graph)
}

conditionalIndepTest <- function(G, pVals, sepset, cor_mat, amat_cpy, i, j_var, S_i, edge_list, pdsep_list, N, l, num_edge_tests){
  # Parameters:
  # 1. cor_mat: correlation matrix
  # 2. amat: adjacency matrix with edge weights between two nodes
  # 3. i, j_var: test for conditional independence between i and variables in j_var
  # 4. S_i: conditioning set 
  # 5. N, l: number of samples, size of conditioning set
  sep_set <- c()
  l_sub <- l
  # for(l_sub in 1:l){
    for(j in j_var){
      # print(cat('j: ',j))
      S_sub_i <- S_i[!S_i %in% j]
      # print(j_var)
        if(length(S_sub_i)>=l_sub && G[i, j]){
            # print(l_sub)
            if(l_sub>1){
              # print(S_sub_i)
              power_set <- as.matrix(combn(S_sub_i, l_sub))
              # power_set <- t(apply(power_set, 1, as.numeric))
            }else{
              power_set <- matrix(S_sub_i, nrow=1)
            }
            for(s in 1:ncol(power_set)){
                if(length(edge_list) > 0){
                    # print(edge_list)
                    edge_idx <- Position(function(x) identical(x, c(i, j)), edge_list, nomatch=0)
                    if(edge_idx!=0){
                        edge_pdsep_list <- pdsep_list[[i]][[j]]
                        if(is.null(edge_pdsep_list)){
                            pdsep_list[[i]][[j]] <- list()
                        }
                        else if(length(edge_pdsep_list) < l_sub){
                            edge_pdsep_list[[l_sub]] <- as.matrix(power_set[, s])
                            pdsep_list[[i]][[j]] <- edge_pdsep_list
                        }else{
                            pdsep_mat <- edge_pdsep_list[[l_sub]]
                            pdsep_present <- sum(apply(pdsep_mat, 2, function(x){
                                                length(intersect(x, power_set[, s]))==length(power_set[, s])}))
                            if(pdsep_present==0){
                                edge_pdsep_list[[l_sub]] <- cbind(edge_pdsep_list[[l_sub]], as.matrix(power_set[, s]))
                                pdsep_list[[i]][[j]] <- edge_pdsep_list
                                # print(pdsep_list[[i]][[j]])
                            }else{
                                next
                            }
                        }                        
                    }else{
                        edge_idx <- length(edge_list)
                        edge_list[[edge_idx+1]] <- c(i, j)
                        edge_pdsep_list <- list()
                        edge_pdsep_list[[l_sub]] <- as.matrix(power_set[, s])
                        pdsep_list[[i]][[j]] <- edge_pdsep_list                        
                    }
                }else{
                    edge_idx <- 1
                    edge_list[[edge_idx]] <- c(i, j)
                    edge_pdsep_list <- list()
                    edge_pdsep_list[[l_sub]] <- as.matrix(power_set[, s])
                    pdsep_list[[i]][[j]] <- edge_pdsep_list
                }
                p_val <- gaussCItest(x=i, y=j, S=power_set[, s], suffStat=list(C=cor_mat, n=N))
                # test_counter <- test_counter + 1
                # print(cat('l_sub:', l_sub, ' num_edge_tests[l_sub]:', num_edge_tests[l_sub]))
                num_edge_tests[l_sub] <- num_edge_tests[l_sub] + 1
                # print(cat('i:', i, 'j:', j, 'S:', power_set[, s], 'p-val:', p_val, 'test counter: ', num_edge_tests[l_sub]))
                if(p_val>=0.05){
                    # print(cat('i:', i, 'j:', j, 'S:', power_set[, s], 'p-val:', p_val))
                    G[i, j] <- G[j, i] <- FALSE
                    # print(sum(G))
                    amat_cpy[i, j] <- amat_cpy[j, i] <- 0
                    sepset[[i]][[j]] <- as.vector(power_set[, s])
                    sep_set <- union(sep_set, as.vector(power_set[, s]))
                    break
                }
                if(pVals[i, j] < p_val) pVals[i, j] <- p_val
            }  
        }
    }
  return(list(G, sepset, pVals, amat_cpy, edge_list, pdsep_list, num_edge_tests))
}

findInterCommVars <- function(c1_vars, c2_vars, bool_amat){
    inter_comm_conn <- bool_amat[c1_vars, !c1_vars, drop=FALSE]
    intra_comm_conn <- bool_amat[c1_vars, c1_vars, drop=FALSE]

    n_inter_edges <- apply(inter_comm_conn, 1, sum)
    n_intra_edges <- apply(intra_comm_conn, 1, sum)

    c1_z <- names(which((n_intra_edges>0 & n_inter_edges>0)==TRUE))
    c1_z <- names(which((n_inter_edges>0)==TRUE))

    inter_comm_conn <- bool_amat[c2_vars, c1_vars, drop=FALSE]
    intra_comm_conn <- bool_amat[c2_vars, c2_vars, drop=FALSE]

    n_inter_edges <- apply(inter_comm_conn, 1, sum)
    n_intra_edges <- apply(intra_comm_conn, 1, sum)

    c2_z <- names(which((n_intra_edges>0 & n_inter_edges>0)==TRUE))
    c2_z <- names(which((n_inter_edges>0)==TRUE))
    # print(cat('c1_z:',c1_z))
    # print(cat('c2_z:',c2_z))

    return(list(c1_z, c2_z))
}   

marginalIndepTest <- function(amat, N, test_counter){
  temp_amat <- amat
  vars <- 1:nrow(amat)
  for(i in vars){
    amat[i, i] <- 0
    y <- which(amat[i, ]!=0)
    for(j in y){
      p_val <- gaussCItest(i, j, S=NULL, suffStat=list(C=temp_amat, n=N))
      test_counter <- test_counter + 1
      if(p_val>0.05)
        amat[i, j] <- amat[j, i] <- 0
    }
  }
  return(list(amat, test_counter))  
}

diffCommunity <- function(communities_old, communities_new){
    flag <- TRUE
    num_old_comms <- length(unique(communities_old$membership))
    num_new_comms <- length(unique(communities_new$membership))
    old_comm_vars <- sort(as.numeric(communities_old$names))
    new_comm_vars <- sort(as.numeric(communities_new$names))

    if(num_old_comms==1 & num_new_comms==1 & identical(old_comm_vars, new_comm_vars)){
        print(cat('old_comm_vars: ', old_comm_vars))
        return(!flag)
    }else{
        return(flag)
    }
}

removeInterCommEdges <- function(c_z, c_vars, G, pVals, sepset, vars, bool_amat, cor_mat, amat, 
                        edge_list, pdsep_list, n, max_degree, num_edge_tests){
    # print(c_z)
    # for(l in 1:max_degree){
    l <- max_degree
    bool_amat_cpy <- bool_amat
    S <- list()
    # p <- 1:nrow(bool_amat)
    # p <- p[!p %in% vars]

    for(i in c_z){
        adj_comm <- names(which(bool_amat[i, vars[c_vars]]!=0))
        # adj_non_comm <- names(which(bool_amat[i, ]!=0))
        adj_non_comm <- names(which(bool_amat[i, vars[!c_vars]]!=0))
        S[[i]] <- c(adj_comm, adj_non_comm)
    }

    for(i in c_z){
        # print(cat('c_vars: ', dim(bool_amat), ' c_vars: ', length(c_vars)))
        # c_i_z <- setdiff(c_z, i)
        adj_comm <- names(which(bool_amat[i, vars[c_vars]]!=0))
        adj_non_comm <- names(which(bool_amat[i, vars[!c_vars]]!=0))
        # print(c_i_z)
        # print(vars[vars==c_i_z])
        # S_i <- c(adj_comm, adj_non_comm)
        X_j <- c(adj_non_comm, adj_comm)

        # if(l>2)
        # print('=============================')
        # print(cat('l: ', l, 'S_i: ', S[[i]], 'X_j: ', X_j, ' vars[c_vars]:', vars[c_vars]))
        if(length(X_j)==0 | length(S[[i]])==0 | length(S[[i]]) < l){
            # print('Next')
            next
        }else{
          res_obj <- conditionalIndepTest(G, pVals, sepset, cor_mat, amat, as.numeric(i), as.numeric(unlist(X_j)), as.numeric(S[[i]]), 
                    edge_list, pdsep_list, n, l, num_edge_tests)
          G <- res_obj[[1]]
          sepset <- res_obj[[2]]
          pVals <- res_obj[[3]]
          amat <- res_obj[[4]]
          edge_list <- res_obj[[5]]
          pdsep_list <- res_obj[[6]]
          num_edge_tests <- res_obj[[7]]
        }
        amat_temp <- amat[vars, vars]
        amat_temp[amat_temp!=0] <- 1
        bool_amat[vars, vars] <- amat_temp
    }
    # }
    # print('=====================')
    return(list(G, sepset, pVals, amat, edge_list, pdsep_list, num_edge_tests))
}

buildSkeletonAllComms <- function(G, pVals, communities_old, cor_mat, amat, bool_amat, sepset, edge_list, pdsep_list, n, vars, num_edge_tests){
  # print(cat('vars: ', vars))
    if(length(vars)<=3){
        # print(cat('Returning with vars', vars))
        if(sum(bool_amat[vars, vars])==6){
            res_obj <- removeInterCommEdges(vars, !vector(mode="logical", length(vars)), G, pVals, sepset, vars, bool_amat, cor_mat, amat, 
                        edge_list, pdsep_list, n, max_degree=1, num_edge_tests)
            G <- res_obj[[1]]
            sepset <- res_obj[[2]]
            pVals <- res_obj[[3]]
            amat <- res_obj[[4]]
            edge_list <- res_obj[[5]]
            pdsep_list <- res_obj[[6]]
            num_edge_tests <- res_obj[[7]]
        }
        # print('num_edge_tests: ', num_edge_tests)
        return(list(G, sepset, pVals, amat, bool_amat, edge_list, pdsep_list, num_edge_tests))
    }
  
    cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
    if(is.null(E(cor_graph)$weight))
        return(list(G, sepset, pVals, amat, bool_amat, edge_list, pdsep_list, num_edge_tests))

    # communities <- cluster_fast_greedy(cor_graph, weights=abs(E(cor_graph)$weight))
    # community_membership <- cutat(communities, no=2)

    communities <- cluster_louvain(cor_graph, weights=abs(E(cor_graph)$weight))
    community_membership <- communities$membership
    # plot(cor_graph, vertex.color=community_membership)
    num_edge_tests <- as.integer(num_edge_tests)

    if(!is.null(communities_old)){
        if(!diffCommunity(communities_old, communities)){
          return(list(G, sepset, pVals, amat, bool_amat, edge_list, pdsep_list, num_edge_tests))  
        }
    }
  
  # Perform 0-1 conditional independence tests
    for(l in 1:length(unique(community_membership))){
    # print('0-1 CI tests')
        max_degree <- 1L
        c_l_vars <- community_membership==l
        c_k_vars <- community_membership!=l
        c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
        c_l_z <- c_z[[1]]
        c_k_z <- c_z[[2]]

        if(length(c_l_z)>0 & length(c_k_z)>0){
            res_obj <- removeInterCommEdges(c_l_z, c_l_vars, G, pVals, sepset, vars, bool_amat, cor_mat, amat, 
                                          edge_list, pdsep_list, n, max_degree, num_edge_tests)
            G <- res_obj[[1]]
            sepset <- res_obj[[2]]
            pVals <- res_obj[[3]]
            amat <- res_obj[[4]]
            edge_list <- res_obj[[5]]
            pdsep_list <- res_obj[[6]]
            num_edge_tests <- res_obj[[7]]

            # new_cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
            # plot(new_cor_graph, vertex.color=community_membership)
        }else{
          if(identical(sort(as.numeric(communities$names)), sort(vars))){
            communities <- cluster_fast_greedy(cor_graph, weights=abs(E(cor_graph)$weight))
            community_membership <- cutat(communities, no=2)
            
            c_l_vars <- community_membership==l
            c_k_vars <- community_membership!=l
            c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
            c_l_z <- c_z[[1]]
            c_k_z <- c_z[[2]]
            
            res_obj <- removeInterCommEdges(c_l_z, c_l_vars, G, pVals, sepset, vars, bool_amat, cor_mat, amat, 
                                            edge_list, pdsep_list, n, max_degree, num_edge_tests)
            G <- res_obj[[1]]
            sepset <- res_obj[[2]]
            pVals <- res_obj[[3]]
            amat <- res_obj[[4]]
            edge_list <- res_obj[[5]]
            pdsep_list <- res_obj[[6]]
            num_edge_tests <- res_obj[[7]]
            
            res_obj <- removeInterCommEdges(c_k_z, c_k_vars, G, pVals, sepset, vars, bool_amat, cor_mat, amat, 
                                            edge_list, pdsep_list, n, max_degree, num_edge_tests)
            G <- res_obj[[1]]
            sepset <- res_obj[[2]]
            pVals <- res_obj[[3]]
            amat <- res_obj[[4]]
            edge_list <- res_obj[[5]]
            pdsep_list <- res_obj[[6]]
            num_edge_tests <- res_obj[[7]]
            
            # new_cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
            # plot(new_cor_graph, vertex.color=community_membership)
          }
        }
        amat_temp <- amat[vars, vars]
        amat_temp[amat_temp!=0] <- 1
        bool_amat[vars, vars] <- amat_temp
    }

  
  # Update c_z and then perform build 0-1 CI graph on each community
    for(l in 1:length(unique(community_membership))){
        c_names <- communities$names
        c_l_vars <- community_membership==l
        c_k_vars <- community_membership!=l
        # print(row.names(bool_amat))
        c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
        c_l_z <- c_z[[1]]
        c_k_z <- c_z[[2]]

        # print(cat('c_l_z', c_l_z, ' c_k_z', c_k_z))
        # print('Communities do not have an edge between them')
        res_obj <- buildSkeletonAllComms(G, pVals, communities, cor_mat, amat, bool_amat, sepset, edge_list, pdsep_list, 
                                       n, as.numeric(c_names[c_l_vars]), num_edge_tests)
        G <- res_obj[[1]]
        sepset <- res_obj[[2]]
        pVals <- res_obj[[3]]
        amat <- res_obj[[4]]
        bool_amat <- res_obj[[5]]
        edge_list <-res_obj[[6]]
        pdsep_list <- res_obj[[7]]
        num_edge_tests <- res_obj[[8]]        
    }
  
    # Remove intra-comm & inter-comm edges of higher order i.e., |S|>=2
    max_degree <- max_degree + 1
    max_comm_degree <- 0
    repeat{
        for(l in 1:length(unique(community_membership))){
            c_names <- communities$names
            c_l_vars <- community_membership==l
            c_k_vars <- community_membership!=l
            # print(row.names(bool_amat))
            c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
            c_l_z <- c_z[[1]]
            c_k_z <- c_z[[2]]

            # print(max(apply(bool_amat[c_l_z, , drop=FALSE], 1, function(x) sum(x))))
            if(length(c_l_z)>0 & length(c_k_z)>0){
                max_conn <- max(apply(bool_amat[c_l_z, , drop=FALSE], 1, function(x) sum(x)))
                if(max_degree<=(max_conn-1)){
                    # print('Removing edges with higher-order CI tests')
                    res_obj <- removeInterCommEdges(c_l_z, c_l_vars, G, pVals, sepset, vars, bool_amat, cor_mat, amat, 
                                                  edge_list, pdsep_list, n, max_degree, num_edge_tests)
                    G <- res_obj[[1]]
                    sepset <- res_obj[[2]]
                    pVals <- res_obj[[3]]
                    amat <- res_obj[[4]]
                    edge_list <- res_obj[[5]]
                    pdsep_list <-res_obj[[6]]
                    num_edge_tests <- res_obj[[7]]

                    # new_cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
                    # plot(new_cor_graph, vertex.color=community_membership)            
                    max_comm_degree <- max(max_comm_degree, max_conn)
                }
            }
            amat_temp <- amat[vars, vars]
            amat_temp[amat_temp!=0] <- 1
            bool_amat[vars, vars] <- amat_temp
        }
        max_degree <- max_degree + 1
        if(max_degree > max_comm_degree)   break
    }   
    # print('num_edge_tests: ', num_edge_tests)
    return(list(G, sepset, pVals, amat, bool_amat, edge_list, pdsep_list, num_edge_tests))
}
