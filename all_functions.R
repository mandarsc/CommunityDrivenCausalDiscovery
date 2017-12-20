buildStableSkeleton <- function(G, pVals, communities_old, cor_mat, amat, data, unique_vals, bool_amat, sepset, n, vars, num_edge_tests){
    # Exit split phase 
    if(length(vars)<=3){
        if(sum(bool_amat[vars, vars])>=2){
            # Uncomment lines 8-13 if want to run PC-stable on subgraph
            temp <- matrix(1, nrow=length(vars), ncol=length(vars))
            diag(temp) <- 0
            bool_amat[vars, vars] <- temp
            amat[vars, vars] <- temp
            G[vars, vars] <- TRUE
            diag(G) <- FALSE

            skel_obj <- pcStableSkeleton(bool_amat, amat, G, data, cor_mat, pVals, unique_vals, vars, sepset)
            bool_amat <- skel_obj[[1]]
            amat <- skel_obj[[2]]
            G <- skel_obj[[3]]
            sepset <- skel_obj[[4]]
        }
        return(list(G, sepset, pVals, amat, bool_amat, num_edge_tests))
    }

    # Perform community detection
    cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
    if(is.null(E(cor_graph)$weight)){
        return(list(G, sepset, pVals, amat, bool_amat, num_edge_tests))
    }
    communities <- cluster_louvain(cor_graph, weights=abs(E(cor_graph)$weight))
    community_membership <- communities$membership
    num_comms <- length(unique(community_membership))

    if(num_comms==1){
        communities <- cluster_fast_greedy(cor_graph, weights=abs(E(cor_graph)$weight))
        community_membership <- cutat(communities, no=2)
        communities$membership <- community_membership
        num_comms <- 2
        # plot(cor_graph, vertex.color=community_membership)
    }
    plot(cor_graph, vertex.color=community_membership)        

    # Exit split phase if there is no change in communities
    if(!is.null(communities_old) && !diffCommunity(communities_old$membership, community_membership, communities_old$names, communities$names)){
        print(cat('No diff in comms: ', vars))
        temp <- matrix(1, nrow=length(vars), ncol=length(vars))
        # Uncomment lines 64-68 if want to run PC-stable on subgraph
        diag(temp) <- 0
        bool_amat[vars, vars] <- temp
        amat[vars, vars] <- temp
        G[vars, vars] <- TRUE
        diag(G) <- FALSE

        skel_obj <- pcStableSkeleton(bool_amat, amat, G, data, cor_mat, pVals, unique_vals, vars, sepset)
        bool_amat <- skel_obj[[1]]
        amat <- skel_obj[[2]]
        G <- skel_obj[[3]]
        sepset <- skel_obj[[4]]
        return(list(G, sepset, pVals, amat, bool_amat, num_edge_tests))
    }else{
        # Remove edges between inter- and intra-community variables
        remove_edge_obj <- removeEdges(G, pVals, communities, cor_mat, amat, data, unique_vals, bool_amat, sepset, n, vars, num_edge_tests, TRUE)
        G <- remove_edge_obj[[1]]
        sepset <- remove_edge_obj[[2]]
        pVals <- remove_edge_obj[[3]]
        amat <- remove_edge_obj[[4]]
        bool_amat <- remove_edge_obj[[5]]
    }

    # Find the variable set adjancent to inter-community variables in each community and store them in z_vars 
    z_vars <- var_set <- vector(mode="list", num_comms)
    bool_amat_comm <- bool_amat
    amat_comm <- amat
    G_comm <- G
    for(l in 1:num_comms){
        c_l_vars <- community_membership==l
        c_k_vars <- community_membership!=l
        c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
        z_vars[[l]] <- c_z[[2]]
    }

    # Enter split phase by iterating over each community and splitting with variable set X U Z
    # where X is the set of variables in a community, C_i, and Z is the variable set of inter-community variables
    # not in C_i
    for(l in 1:num_comms){
        c_names <- communities$names
        c_l_vars <- community_membership==l

        new_var_set <- sort(as.numeric(union(c_names[c_l_vars], z_vars[[l]])))
        same_var_set <- 0
        if(l>1){
            same_var_set <- sum(unlist(lapply(var_set, function(x) all(x %in% new_var_set) && all(new_var_set %in% x))))
        }
        if(same_var_set==0){
            var_set[[l]] <- new_var_set
            # Uncomment lines 120-121 if doing overlap, otherwise uncomment lines 122-123
            res_obj <- buildStableSkeleton(G_comm, pVals, communities, cor_mat, amat_comm, data, unique_vals, bool_amat_comm, sepset,  
                                       n, sort(union(as.numeric(c_names[c_l_vars]), z_vars[[l]])), num_edge_tests)
            G <- G & res_obj[[1]]
            sepset <- res_obj[[2]]
            pVals <- res_obj[[3]]
            amat_temp <- res_obj[[4]]
            idx <- which(amat_temp==0, arr.ind=TRUE)
            amat[idx] <- 0
            bool_amat[idx] <- 0
        }
    }
    # Combine all the G_sub into one graph
    return(list(G, sepset, pVals, amat, bool_amat, num_edge_tests))
}

pcStableSkeleton <- function(bool_amat, amat, G, data, cor_mat, pVals, unique_vals, vars, sepset){
    rem_edges <- which(bool_amat[vars, vars]==1, arr.ind=TRUE)
    # print(rem_edges)
    rem_edges[, 1] <-vars[rem_edges[, 1]]
    rem_edges[, 2] <-vars[rem_edges[, 2]]
    print(cat('vars: ', vars))
    # print(rem_edges)
    l <- 0
    S <- vector(mode="list", length=max(vars))
    repeat{
        for(e in 1:nrow(rem_edges)){
            i <- rem_edges[e, 1]
            j <- rem_edges[e, 2]
            # print(cat('i: ', i, ' j: ', j))
            if(G[i, j]){
                if(l==0){
                    ci_obj <- conditionalIndepTest(G, pVals, sepset, cor_mat, amat, data, unique_vals, i, j, NULL, nrow(data), l)
                }else{
                    for(v in vars){
                        S[[v]] <- vars[which(bool_amat[v, vars]==TRUE, arr.ind=TRUE)]
                        # print(cat('S: ', S[[v]], ' l: ', l))
                    }
                    ci_obj <- conditionalIndepTest(G, pVals, sepset, cor_mat, amat, data, unique_vals, i, j, S[[i]], nrow(data), l)
                }
                G <- ci_obj[[1]]
                amat <- ci_obj[[2]]
                pVals <- ci_obj[[3]]
                sepset <- ci_obj[[4]]                   
            }
        }
        idx <- which(G==FALSE, arr.ind=TRUE)
        bool_amat[idx] <- 0
        if(max(apply(bool_amat[vars, vars], 2, sum))-1 < l) break                
        l <- l + 1
    }
    return(list(bool_amat, amat, G, sepset))
}

removeEdges <- function(G, pVals, communities, cor_mat, amat, data, unique_vals, bool_amat, sepset, n, vars, num_edge_tests, diffComm=TRUE){
    community_membership <- communities$membership
    num_comms <- length(unique(community_membership))

    bool_amat_comm <- bool_amat
    c_l_z <- list()    
    # Get a list of all inter-community variables in each community
    for(l in 1:num_comms){
        c_l_vars <- community_membership==l
        c_k_vars <- community_membership!=l
        c_z <- findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
        c_l_z[[l]] <- c_z[[1]]
    }
 
    max_conn <- array(0, num_comms)
    max_comm_degree <- max_degree <- 0

    repeat{
        # Build the conditioning set for all vars
        S <- buildSepset(num_comms, community_membership, bool_amat_comm, vars, diffComm)
        for(l in 1:num_comms){
            c_l_vars <- community_membership==l
            # Update c_l_z[[l]] if a variable is no longer an inter-community variable
            for(c_l in c_l_z[[l]]){
                if(sum(bool_amat_comm[c_l, unlist(c_l_z)])==0)
                    c_l_z[[l]] <- c_l_z[[l]][!c_l_z[[l]] %in% c_l]    
            }
            if(length(c_l_z[[l]])>0){
                max_conn[l] <- getMinMBSize(data, unlist(c_l_z), bool_amat_comm, vars, unique_vals)
                print(cat('c_l_z: ', c_l_z[[l]], ' max_conn: ', max_conn, ' max_degree: ', max_degree))
                if(max_degree <= max_conn[l]){
                    # Test the edges with |S|=max_degree
                    res_obj <- removeInterCommEdges(c_l_z[[l]], c_l_vars, G, pVals, sepset, vars, list(bool_amat, S), cor_mat, amat, data, unique_vals,
                                                  n, max_degree, num_edge_tests, diffComm)
                    G <- res_obj[[1]]
                    sepset <- res_obj[[2]]
                    pVals <- res_obj[[3]]
                    amat_temp <- res_obj[[4]]
                    idx <- which(amat_temp==0, arr.ind=TRUE)
                    amat[idx] <- 0
                    bool_amat[idx] <- 0
                }                
            }
        }
        # print(cat('max_conn: ', max_conn))
        max_comm_degree <- max(max_conn)
        bool_amat_comm[vars, vars] <- bool_amat[vars, vars]
        max_degree <- max_degree + 1
        if(max_degree > max_comm_degree)   break                   
    }
    return(list(G, sepset, pVals, amat, bool_amat))
}

conditionalIndepTest <- function(G, pVals, sepset, cor_mat, amat, data, unique_vals, i, j_var, S_i, N, l){
    # Parameters:
    # 1. cor_mat: correlation matrix
    # 2. amat: adjacency matrix with edge weights between two nodes
    # 3. i, j_var: test for conditional independence between i and variables in j_var
    # 4. S_i: conditioning set 
    # 5. N, l: number of samples, size of conditioning set
    power_set <- matrix(0, nrow=1, ncol=1)
    for(j in j_var){
        # print(cat('j: ',j))
        if(!is.null(S_i)){
            S_sub_i <- S_i[!S_i %in% j]            
            if(length(S_sub_i)>=l){
                # print(l)
                if(l==1){
                    power_set <- matrix(S_sub_i, nrow=1)
                }else{
                    power_set <- as.matrix(combn(S_sub_i, l))
                }
            }else{
                break
            }
        }

        # print(j_var)
        start_time <- proc.time()
        cond_set_size <- ifelse(l==0, 1, ncol(power_set))
        for(s in 1:cond_set_size){
            if(l>0){
                cond_set <- power_set[, s]
            }else{
                cond_set <- NULL
            }
            if(is.null(unique_vals))
                p_val <- gaussCItest(x=i, y=j, S=cond_set, suffStat=list(C=cor_mat, n=N))
            else
                p_val <- disCItest(x=i, y=j, S=cond_set, suffStat=list(dm=data, nlev=unique_vals, adaptDF=FALSE))

            if(p_val>=0.05){
                print(cat('i:', i, 'j:', j, 'S:', cond_set, 'p-val:', p_val))
                # print(l_sub)
                G[i, j] <- G[j, i] <- FALSE
                amat[i, j] <- amat[j, i] <- 0
                # print(sepset)
                print(cat(' length: ', sepset[[i]][[j]]))
                sepset[[i]][[j]] <- as.vector(cond_set)
                
                if(pVals[i, j] < p_val){
                    pVals[i, j] <- p_val
                } 
                break
            }
        }
    }
  return(list(G, amat, pVals, sepset))
}

removeInterCommEdges <- function(c_z, c_vars, G, pVals, sepset, vars, bool_list, cor_mat, amat, data, unique_vals,  
                        n, max_degree, num_edge_tests, diffComm=TRUE){

    l <- max_degree
    bool_amat <- bool_list[[1]]
    S <- bool_list[[2]]

    if(is.null(S)){
        S <- list()
        for(i in c_z){
            adj_comm_idx <- which(bool_amat[i, vars[c_vars]]!=0, arr.ind=TRUE)
            adj_non_comm_idx <- which(bool_amat[i, vars[!c_vars]]!=0, arr.ind=TRUE)
            adj_comm <- colnames(bool_amat)[vars[c_vars]]
            adj_comm <- adj_comm[adj_comm_idx]
            adj_non_comm <- colnames(bool_amat)[vars[!c_vars]]
            adj_non_comm <- adj_non_comm[adj_non_comm_idx]
            S[[i]] <- c(adj_comm, adj_non_comm)
        }        
    }

    # print(cat('c_z: ', c_z))
    for(i in c_z){
        adj_comm_idx <- which(bool_amat[i, vars[c_vars]]!=0)
        adj_non_comm_idx <- which(bool_amat[i, vars[!c_vars]]!=0)
        adj_comm <- colnames(bool_amat)[vars[c_vars]]
        adj_comm <- adj_comm[adj_comm_idx]
        adj_non_comm <- colnames(bool_amat)[vars[!c_vars]]
        adj_non_comm <- adj_non_comm[adj_non_comm_idx]
        # if(!diffComm){
        #     X_j <- vars[which(bool_amat[i, vars]!=0)]
        #     # print(cat('vars:', vars, 'X_j: ', X_j))
        # }else{
            X_j <- adj_non_comm
        # } 

        if(length(X_j)==0 | length(S[[i]])==0 | length(S[[i]]) < l){
            next
        }else{
            start_time <- proc.time()
            # print(cat('i: ', i, ' j: ', X_j, 'S: ', S[[i]]))
            res_obj <- conditionalIndepTest(G, pVals, sepset, cor_mat, amat, data, unique_vals, as.numeric(i), as.numeric(unlist(X_j)), as.numeric(S[[i]]), n, l)
            # print(cat('CI Time: ', proc.time()['elapsed']-start_time['elapsed']))
            G <- res_obj[[1]]
            amat <- res_obj[[2]]
            pVals <- res_obj[[3]]
            sepset <- res_obj[[4]]

            amat_temp <- amat[vars, vars]
            amat_temp[amat_temp!=0] <- 1
            bool_amat[vars, vars] <- amat_temp
        }
    }
    return(list(G, sepset, pVals, amat))
}

getMinMBSize <- function(data, c_z, bool_amat, vars, unique_vals){
    # Iterate over each variable in c_z and get the size of their parent-child set
    # Select the minimum size of the parent-child set
    pc_size <- c()
    data <- data.frame(data)
    # print(cat('c_z:', c_z, ' sum:', sum(c_z)))

    for(v in c_z){
        nbrs <- vars[which(bool_amat[v, vars]!=0)]
        if(length(nbrs)>0){
            # print(c(as.numeric(v), as.numeric(nbrs)))
            data_sub <- data[, c(v, nbrs)]
            if(!is.null(unique_vals)){
                for(v_i in 1:length(c(v, nbrs))){
                    data_sub[, v_i] <- as.factor(data_sub[, v_i])
                }
            }
            mb_obj <- mmhc(data_sub)
            pc_size <- c(pc_size, length(mb_obj$nodes[[1]]$parents)+length(mb_obj$nodes[[1]]$children))
            # pc_size <- c(pc_size, length(mb_obj$nodes[[1]]$mb))
            # print(cat('v:  ', v, ' nbrs: ', nbrs, ' parents: ', length(mb_obj$nodes[[1]]$parents), ' children: ', length(mb_obj$nodes[[1]]$children)))
        }
    }
    # print((min(pc_size)))
    return(min(pc_size))
}

build01Graph <- function(data, disc_flag, alpha, amat=NULL, vars=NULL, edge_list=NULL, pdsep_list=NULL, sepset=NULL){
    if(is.null(vars) & is.null(edge_list) & is.null(pdsep_list)){
        p <- ncol(data)
        seq_p <- seq_len(p)
        vars <- 1:p
        edge_list <- vector("list", 20)
        pdsep_list <- vector("list", p)    
        sepset <- lapply(seq_p, function(.) vector("list", p))        
    }else{
        seq_p <- vars
    }
    sepset_test <- lapply(1:ncol(data), function(.) vector("list", ncol(data)))
    bool_amat <-  matrix(TRUE, nrow=ncol(data), ncol=ncol(data))
    new_amat <- matrix(0, nrow=ncol(data), ncol=ncol(data))

    N <- nrow(data)
    diag(bool_amat) <- FALSE
    amat_lower_tri <- lower.tri(bool_amat[vars, vars])
    num_edge_tests <- array(0, 19)    
    rem_edges <- which(amat_lower_tri==1, arr.ind=T)

    if(!is.null(vars)){
        rem_edges[, 1] <- vars[rem_edges[, 1]]
        rem_edges[, 2] <- vars[rem_edges[, 2]]        
    }

    p_val <- array(0, nrow(rem_edges))

    if(disc_flag){
        cor_amat <- cor(data)
    }else{
        unique_vals <- apply(data, 2, function(x) length(unique(x)))
    }

    # Iterate over each edge X_i -- X_j and perform CI test ind(X_i, X_j, Z) with conditioning set Z
    # where Z= empty_set U X\{X_i, X_j}
    for(e in 1:nrow(rem_edges)){
        p_vals <- c()
        sepset_temp <- c()
        ord <- 0L
        x <- rem_edges[e, 1]
        y <- rem_edges[e, 2]
        nbrs_bool <- bool_amat[vars, x]
        nbrs_bool[which(vars==y)] <- FALSE
        nbrs <- seq_p[nbrs_bool]
        length_nbrs <- length(nbrs)

        # Build 01-graph
        idx <- 1
        edge_pdsep_list <- list()
        while(ord <= 1){
            S <-seq_len(ord)
            repeat{
                num_edge_tests[ord+1] <- num_edge_tests[ord+1] + 1
                if(disc_flag)
                    p_vals <- c(p_vals, gaussCItest(x, y, S=nbrs[S], suffStat=list(C=cor_amat, n=N)))
                else
                    p_vals <- c(p_vals, disCItest(x, y, S=nbrs[S], suffStat=list(dm=data, nlev=unique_vals, adaptDF=FALSE)))

                next_set <- getNextSet(length_nbrs, ord, S)
                # sepset_test[[x]][[y]] <- sepset_test[[y]][[x]] <- nbrs[S]
                if(ord==0){
                    sepset_temp <- c(sepset_temp, 0)
                }else{
                    edge_list[[idx]] <- c(x, y)
                    if(idx > 1)
                        edge_pdsep_list[[ord]] <- cbind(edge_pdsep_list[[ord]], as.matrix(nbrs[S]))
                    else 
                        edge_pdsep_list[[ord]] <- as.matrix(nbrs[S])
                    pdsep_list[[x]][[y]] <- pdsep_list[[y]][[x]] <- edge_pdsep_list
                    idx <- idx + 1
                    sepset_temp <- c(sepset_temp, nbrs[S])
                }
                if(next_set$wasLast)
                    break
                S <- next_set$nextSet
            }
            ord <- ord + 1L
        }
        max_pval_idx <- which.max(p_vals)
        p_val[e] <- p_vals[max_pval_idx]
        sepset_test[[x]][[y]] <- sepset_test[[y]][[x]] <- sepset_temp[max_pval_idx]
    }
    # print(pdsep_list[[8]][[10]])
    p_val_adj <- p.adjust(p_val, method="BH")
    adj_idx <- which(p_val_adj < alpha)
    non_adj_idx <- which(p_val_adj >= alpha)

    # For the edges that remain, build adjacency matrix with variable association as the edge weights
    for(idx in adj_idx){
        x <- rem_edges[idx, 1]
        y <- rem_edges[idx, 2]
        # gk_tau <- GKtau(data[, x], data[, y])
        # new_amat[x, y] <- new_amat[y, x] <- max(gk_tau$tauxy, gk_tau$tauyx)
        if(disc_flag){
            new_amat[x, y] <- new_amat[y, x] <- abs(cor(data[, x], data[, y]))
        }else{
            new_amat[x, y] <- new_amat[y, x] <- cv.test(data[, x], data[, y])                
        }
    }
    for(idx in non_adj_idx){
        x <- rem_edges[idx, 1]
        y <- rem_edges[idx, 2]
        sepset[[x]][[y]] <- sepset[[y]][[x]] <- sepset_test[[x]][[y]]
    }
    return(list(new_amat, num_edge_tests, sepset, edge_list, pdsep_list))
}

buildCorGraph <- function(data, disc_flag, alpha){
    num_vars <- ncol(data)
    bool_amat <- matrix(TRUE, nrow=num_vars, ncol=num_vars)
    amat <- matrix(0, nrow=num_vars, ncol=num_vars)
 
    diag(bool_amat) <- FALSE
    amat_lower_tri <- lower.tri(bool_amat)
    rem_edges <- which(amat_lower_tri==1, arr.ind=T)

    for(e in 1:nrow(rem_edges)){
        x <- rem_edges[e, 1]
        y <- rem_edges[e, 2]
        if(disc_flag){
            cor_obj <- chisq.test(data[, x], data[, y], correct=FALSE)
            cor_val <- sqrt(cor_obj$statistic/(length(data[, x]) * (min(length(unique(data[,x])), length(unique(data[, y])))-1)))            
        }else{
            cor_obj <- cor.test(data[, x], data[, y])
            cor_val <- abs(cor_obj$estimate)
        }    
        if(cor_obj$p.value < alpha){
                amat[x, y] <- amat[y, x] <- cor_val
        }            
    }
    return(amat)
}

cdcd_pc <- function(G, pVals, cor_mat, amat, data, bool_amat, sepset, edge_list, pdsep_list, n, vars, num_edge_tests, test_counter, run_time){
    cl <- match.call()
    if(is.null(cor_mat)){
        unique_vals <- apply(data, 2, function(x) length(unique(x)))
    }

    res_obj <- buildStableSkeleton(G, pVals, NULL, cor_mat, amat, data, unique_vals, bool_amat, sepset, n, vars, num_edge_tests)
    G <- res_obj[[1]]
    sepset <- res_obj[[2]]
    pVals <- res_obj[[3]]
    amat <- res_obj[[4]]
    # num_edge_tests <- res_obj[[8]]
    # run_time <- res_obj[[9]]
    # non_zero_idx <- which(num_edge_tests!=0)

    start_time <- proc.time()
    Gobject <- if (sum(G) == 0) {
        print(vars)
        new("graphNEL", nodes = as.character(vars))
    }
    else {
        colnames(G) <- rownames(G) <- vars
        as(G, "graphNEL")
    }
    cdcd_skel <- new("pcAlgo", graph = Gobject, call = cl, n = integer(0), 
        max.ord = as.integer(4), n.edgetests = c(1,2,3,4), 
        sepset = sepset, pMax = pVals, zMin = matrix(NA, 1, 1))
    cdcd_skel@call <- cl
    cdcd_graph <- switch("relaxed", rand=udag2pdag(cdcd_skel), retry=udag2pdagSpecial(cdcd_skel)$pcObj, relaxed = udag2pdagRelaxed(cdcd_skel, verbose = FALSE, 
                solve.confl = FALSE))
    # run_time$Other <- run_time$Other + proc.time()['elapsed'] - start_time['elapsed']
    return(list(cdcd_graph))
}

findInterCommVars <- function(c1_vars, c2_vars, bool_amat){
    inter_comm_conn <- bool_amat[c1_vars, !c1_vars, drop=FALSE]
    intra_comm_conn <- bool_amat[c1_vars, c1_vars, drop=FALSE]
    n_inter_edges <- apply(inter_comm_conn, 1, sum)
    n_intra_edges <- apply(intra_comm_conn, 1, sum)

    c1_z <- names(which((n_intra_edges>0 & n_inter_edges>0)==TRUE))
    c1_z <- as.numeric(names(which((n_inter_edges>0)==TRUE)))
    # print(cat('c1_z: ', c1_z, ' n_inter_edges: ', n_inter_edges))

    inter_comm_conn <- bool_amat[c2_vars, c1_vars, drop=FALSE]
    intra_comm_conn <- bool_amat[c2_vars, c2_vars, drop=FALSE]

    n_inter_edges <- apply(inter_comm_conn, 1, sum)
    n_intra_edges <- apply(intra_comm_conn, 1, sum)

    c2_z <- names(which((n_intra_edges>0 & n_inter_edges>0)==TRUE))
    c2_z <- as.numeric(names(which((n_inter_edges>0)==TRUE)))

    return(list(c1_z, c2_z))
}   

cv.test <- function(x,y) {
  CV = sqrt(chisq.test(x, y, correct=FALSE)$statistic /
    (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  return(as.numeric(CV))
}

buildSepset <- function(num_comms, community_membership, bool_amat, vars, diffComm){
    S <- vector(mode="list", length=max(vars))
    # S <- list()
    # print('=============')
    for(l in 1:num_comms){
        c_l_vars <- community_membership==l
        c_k_vars <- community_membership!=l
        c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
        c_l_z <- c_z[[1]]
        if(!diffComm){
            c_l_z <- as.numeric(union(c_l_z, vars[c_l_vars]))
            # print(cat('c_l_z: ', c_l_z))
        }
        for(var in c_l_z){
            adj_comm_idx <- which(bool_amat[var, vars[c_l_vars]]!=0, arr.ind=TRUE)
            adj_non_comm_idx <- which(bool_amat[var, vars[!c_l_vars]]!=0, arr.ind=TRUE)
            adj_comm <- colnames(bool_amat)[vars[c_l_vars]]
            adj_comm <- adj_comm[adj_comm_idx]
            adj_non_comm <- colnames(bool_amat)[vars[!c_l_vars]]
            adj_non_comm <- adj_non_comm[adj_non_comm_idx]
            S[[var]] <- c(adj_comm, adj_non_comm)
            # print(cat('var: ', var, 'S[[var]]: ', S[[var]]))
        }        
    }
    return(S)    
}

diffCommunity <- function(communities_old_membership, communities_new_membership, comm_old_names, comm_new_names){
    # Returns TRUE if two communities are different otherwise FALSE
    if(!identical(sort(comm_old_names), sort(comm_new_names))){
        return(TRUE)
    }else{
        num_old_comms <- length(unique(communities_old_membership))
        num_new_comms <- length(unique(communities_new_membership))
        if(num_new_comms==num_old_comms){
            comm_old_list <- lapply(1:num_old_comms, function(x) comm_old_names[communities_old_membership==x])
            comm_new_list <- lapply(1:num_new_comms, function(x) comm_new_names[communities_new_membership==x])
            counter <- 0
            for(c in 1:num_old_comms){
                counter <- counter + sum(unlist(lapply(comm_new_list, function(x) identical(sort(x), sort(comm_old_list[[c]])))))
            }
            if(counter==num_old_comms){
                return(FALSE)
            }else{
                return(TRUE)
            }
        }
        else{
            return(TRUE)
        }
    }
}

# Functions mentioned below are not used
# checkSepset <- function(data, sep_set, i, j_var, amat_cpy, G, sepset){
#     for(j in j_var){
#         if(j %in% unlist(sep_set) & !G[i, j]){
#             print(cat('i: ', i, ' j:', j))
#             G[i, j] <- G[j, i] <- TRUE
#             amat_cpy[i, j] <- amat_cpy[j, i] <- cv.test(data[, i], data[, j])
#             # sepset[[i]][[j]] <- NULL
#         }
#     }
#     return(list(G, amat_cpy, sepset))
# }

# marginalIndepTest <- function(data, N, test_counter, disc_flag, run_time){
#     temp_amat <- amat <-  matrix(1, nrow=ncol(data), ncol=ncol(data))
#     vars <- 1:ncol(data)
#     diag(amat) <- 0
    
#     start_time <- proc.time()
#     if(disc_flag){
#         cor_amat <- cor(data)
#     }else{
#         for(i in vars){
#             y <- which(amat[i, ]!=0)
#             for(j in y){
#                 cv_val <- cv.test(data[, i], data[, j])
#                 amat[i, j] <- amat[j, i] <- cv_val
#             }
#         }
#         # amat <- calcEdgeWgt(data)
#         unique_vals <- apply(data, 2, function(x) length(unique(x)))
#     }
#     run_time$MI <- run_time$MI + proc.time()['elapsed'] - start_time['elapsed']
#     for(i in vars){
#         y <- which(amat[i, ]!=0)
#         for(j in y){
#             if(disc_flag){
#                 p_val <- gaussCItest(i, j, S=NULL, suffStat=list(C=cor_amat, n=N))
#             }
#             else{
#                 p_val <- disCItest(i, j, S=NULL, suffStat=list(dm=data, nlev=unique_vals, adaptDF=FALSE))
#                 # p_val <- disCItest(i, j, S=NULL, suffStat=list(dm=data, nlev=unique_vals, adaptDF=FALSE))                
#             }
#             test_counter[1] <- test_counter[1] + 1
#             if(p_val>0.05)
#                 amat[i, j] <- amat[j, i] <- 0
#         }
#     }
#     run_time$MI <- run_time$MI + proc.time()['elapsed'] - start_time['elapsed']
#     print(run_time$MI)
#     return(list(amat, test_counter, run_time))  
# }

# calcEdgeWgt <- function(data){
#     adj_mat <- matrix(1, nrow=ncol(data), ncol=ncol(data))
#     diag(adj_mat) <- 0
#     adj_list <- lapply(seq_len(ncol(data)), function(x) adj_mat[x, ])
#     start_time <- proc.time()
#     N <- nrow(data)
#     amat_list <- lapply(adj_list, function(x){
#         i <- which(x==0)
#         sapply(1:length(x), function(z) as.numeric(sqrt(chisq.test.mod(table(data[, i], data[, z]), N)/
#             (N*(min(length(unique(data[, i])), length(unique(data[, z])))-1)))))
#     })
#     # print(proc.time()['elapsed']-start.time['elapsed'])
#     amat <- matrix(unlist(amat_list), ncol=ncol(data), byrow=TRUE)
#     return(amat)
# }

# chisq.test.mod <- function (x, n){
#     nr <- as.integer(nrow(x))
#     nc <- as.integer(ncol(x))
#     sr <- rowSums(x)
#     sc <- colSums(x)
#     E <- (sr %o% sc)/n
#     return(sum((abs(x - E) - 0)^2/E))