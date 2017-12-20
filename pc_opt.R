library(bnlearn)
library(pcalg)
library(igraph)
library(Matrix)
library(corpcor)
library(GoodmanKruskal)

# Generate random simulated DAGs and corresponding data sets

num_vars <- c(10, 20, 50, 100, 200, 500, 1000)

data <- dag <- amat <- list()
cor_mat <- adj_nodes <- list()
cdcd_obj <- pc_obj <- pc_stable_obj <- dist_mat <- all_node_dist_mat <- list()
num_edge_tests <- sepset <- edge_list <- pdsep_list <- list()
i <- 6
num_iter <- 50

pc_cdcd_shd <- pc_shd <- pc_stable_shd <- cdcd_shd <- matrix(0, num_iter, 1)
pc_tests <- matrix(0, num_iter, 20)
pc_stable_tests <- matrix(0, num_iter, 20)
cdcd_tests <- matrix(0, num_iter, 20)

tpr <- matrix(0, num_iter, 3, dimnames=list(c(), c('CDCD-TPR', 'PC-TPR', 'PC-STABLE-TPR')))
fpr <- matrix(0, num_iter, 3, dimnames=list(c(), c('CDCD-FPR', 'PC-FPR', 'PC-STABLE-FPR')))
tdr <- matrix(0, num_iter, 3, dimnames=list(c(), c('CDCD-TDR', 'PC-TDR', 'PC-STABLE-TDR')))
perf_metrics <- matrix(0, num_iter, 9, dimnames=list(c(), c('CDCD-TP', 'PC-TP', 'PC-STABLE-TP', 'CDCD-FP', 'PC-FP', 'PC-STABLE-FP', 'CDCD-FN', 'PC-FN', 'PC-STABLE-FN')))
cdcd_rt <- matrix(0, num_iter, 4, dimnames=list(c(), c('MI', 'Order1', 'Order2', 'Other')))
pc_rt <- pc_stable_rt <- array(0, num_iter)

# change the network and sample size for different networks and sample sizes
network <- "alarm"
sample_size <- 1000
r_data_file <- paste(network, '_', sample_size, '.RData', sep='')
load(paste(network, '.rda', sep=''))
load(r_data_file)
# sample_size <- 5000

for(j in 1:10){

    # dag[[j]] <- as.graphNEL(bn)
    print(cat('j= ', j))
    alpha <- 0.05 
    dag[[j]] <- randomDAG(num_vars[i], prob=(2/(num_vars[i]-1)))
    data[[j]] <- rmvDAG(n=sample_size, dag=dag[[j]], errDist="normal")

    test_counter <- 0
    unique_vals <- cor_mat <- NULL

    # Uncomment next two lines if data is continuous
    # real_data <- rbn(bn, n=sample_size)
    # data[[j]] <- real_data

    # Comment line 60 if data is continuous
    # while(TRUE){
    #     real_data <- rbn(bn, n=sample_size)
    #     data[[j]] <- matrix(0, nrow=nrow(real_data), ncol=ncol(real_data))
    #     for(c in 1:ncol(data[[j]])){
    #         data[[j]][, c] <- as.numeric(real_data[, c])
    #         data[[j]][, c] <- data[[j]][, c] - 1
    #     }
        unique_vals <- apply(data[[j]], 2, function(x) length(unique(x)))
    #     if(length(which(unique_vals==1))==0)
    #         break            
    # }

    run_time <- list(MI=0, Order1=0, Order2=0, Other=0)
    
    start_time <- proc.time()
    # Build undirected independence graph
    undir_obj <- undirIndepGraph(data[[j]], is.null(unique_vals), alpha)
        
    amat[[j]] <- buildCorGraph(data[[j]], !is.null(unique_vals), alpha)
    idx <- which(amat[[j]]!=0, arr.ind=TRUE)
    # print(pdsep_list[[j]][[8]][[10]])

    if(is.null(unique_vals)){
        org_amat <- wgtMatrix(dag[[j]], transpose=FALSE)
        org_idx <- which(org_amat!=0, arr.ind=TRUE)
        org_amat[org_idx] <- 1
        cor_mat <- cor(data[[j]])
        amat[[j]][idx] <- cor_mat[idx]
    }else{
        org_amat <- amat(bn)
    }
    
    # Initilialize stuff
    row.names(org_amat) <- colnames(org_amat) <- 1:ncol(data[[j]])
    org_graph <- graph_from_adjacency_matrix(org_amat, mode="directed")
    dag_bn <- as.bn(igraph.to.graphNEL(org_graph))

    p <- ncol(amat[[j]])
    n <- nrow(data[[j]])
    vars <- 1:ncol(data[[j]])
    num_edge_tests <- array(0, 19)
    pVals <- matrix(-Inf, nrow=p, ncol=p)
    G <- matrix(FALSE, nrow=p, ncol=p)
    pc_num_edge_tests <- array(0, 20)
    # edge_list <- vector("list", 20)
    # pdsep_list <- vector("list", p)
    sepset <- lapply(seq_len(p), function(.) vector("list", p))

    G[idx] <- TRUE
    bool_amat_cpy <- bool_amat <- amat_cpy <- amat[[j]]
    bool_amat[bool_amat!=0] <- 1
    row.names(amat[[j]]) <- colnames(amat[[j]]) <- row.names(bool_amat) <- colnames(bool_amat) <- vars

    # Run CDCD
    res_obj <- cdcd_pc(G, pVals, cor_mat, amat[[j]], data[[j]], bool_amat, sepset, edge_list[[j]], pdsep_list[[j]], 
                nrow(data[[j]]), vars, num_edge_tests, test_counter, run_time)
    cdcd_obj[[j]] <- res_obj[[1]]
    cdcd_rt[j,] <- unlist(res_obj[[2]])

    if(is.null(unique_vals)){
        start_time <- proc.time()
        pc_obj[[j]] <- pc(suffStat=list(C=cor(data[[j]]), n=nrow(data[[j]])), indepTest=gaussCItest, alpha=alpha, p=ncol(data[[j]]), verbose=F, skel.method="original")
        pc_rt[j] <- proc.time()['elapsed'] - start_time['elapsed']
        start_time <- proc.time()
        pc_stable_obj[[j]] <- pc(suffStat=list(C=cor(data[[j]]), n=nrow(data[[j]])), indepTest=gaussCItest, alpha=alpha, p=ncol(data[[j]]), verbose=F, skel.method="stable")
        pc_stable_rt <- proc.time()['elapsed'] - start_time['elapsed']
    }else{
        start_time <- proc.time()
        pc_obj[[j]] <- pc(suffStat=list(dm=data[[j]], nlev=unique_vals, adaptDF=FALSE), indepTest=disCItest, alpha=alpha, p=ncol(data[[j]]), skel.method="original", verbose=F)
        pc_rt[j] <- proc.time()['elapsed'] - start_time['elapsed']
        start_time <- proc.time()
        pc_stable_obj[[j]] <- pc(suffStat=list(dm=data[[j]], nlev=unique_vals, adaptDF=FALSE), indepTest=disCItest, alpha=alpha, p=ncol(data[[j]]), skel.method="stable", verbose=F)
        pc_stable_rt[j] <- proc.time()['elapsed'] - start_time['elapsed']
    }            

    print(cat('CDCD run_time: ', cdcd_rt[j, ]))
    print(cat('PC-stable run_time: ', pc_stable_rt[j]))
    print(cat('PC run_time: ', pc_rt[j]))   
     
    pc_cdcd_shd[j, ] <- pcalg::shd(pc_obj[[j]], cdcd_obj[[j]])
    pc_shd[j, ] <- pcalg::shd(dag[[j]], pc_obj[[j]])
    cdcd_shd[j, ] <- pcalg::shd(dag[[j]], cdcd_obj[[j]])
    pc_stable_shd[j, ] <- pcalg::shd(dag[[j]], pc_stable_obj[[j]])

    pc_num_edge_tests[1:length(pc_obj[[j]]@n.edgetests)] <- pc_obj[[j]]@n.edgetests
    pc_tests[j, 1:length(pc_num_edge_tests)] <- pc_num_edge_tests
    cdcd_tests[j, 1:length(cdcd_obj[[j]]@n.edgetests)] <- cdcd_obj[[j]]@n.edgetests
    pc_stable_tests[j, 1:length(pc_stable_obj[[j]]@n.edgetests)] <- pc_stable_obj[[j]]@n.edgetests

    tpr[j, 'CDCD-TPR'] <- compareGraphs(cdcd_obj[[j]]@graph, dag[[j]])['tpr']
    fpr[j, 'CDCD-FPR'] <- compareGraphs(cdcd_obj[[j]]@graph, dag[[j]])['fpr']
    tdr[j, 'CDCD-TDR'] <- compareGraphs(cdcd_obj[[j]]@graph, dag[[j]])['tdr']

    tpr[j, 'PC-TPR'] <- compareGraphs(pc_obj[[j]]@graph, dag[[j]])['tpr']
    fpr[j, 'PC-FPR'] <- compareGraphs(pc_obj[[j]]@graph, dag[[j]])['fpr']
    tdr[j, 'PC-TDR'] <- compareGraphs(pc_obj[[j]]@graph, dag[[j]])['tdr']

    tpr[j, 'PC-STABLE-TPR'] <- compareGraphs(pc_stable_obj[[j]]@graph, dag[[j]])['tpr']
    fpr[j, 'PC-STABLE-FPR'] <- compareGraphs(pc_stable_obj[[j]]@graph, dag[[j]])['fpr']
    tdr[j, 'PC-STABLE-TDR'] <- compareGraphs(pc_stable_obj[[j]]@graph, dag[[j]])['tdr']

    perf_metrics[j, c(1, 4, 7)] <- unlist(bnlearn::compare(dag_bn, as.bn(cdcd_obj[[j]]@graph)))
    perf_metrics[j, c(2, 5, 8)] <- unlist(bnlearn::compare(dag_bn, as.bn(pc_obj[[j]]@graph)))
    perf_metrics[j, c(3, 6, 9)] <- unlist(bnlearn::compare(dag_bn, as.bn(pc_stable_obj[[j]]@graph))) 
}

print('==== CI Tests ====')
print(cat('CDCD: ', mean(apply(cdcd_tests[1:j,],1,sum)), ' PC: ', mean(apply(pc_tests[1:j,],1,sum)), ' PC-Stable: ', mean(apply(pc_stable_tests[1:j,],1,sum))))
print('==== SHD ====')
print(cat('CDCD: ', mean(cdcd_shd[1:j]), ' PC: ', mean(pc_shd[1:j]), ' PC-Stable: ', mean(pc_stable_shd[1:j])))
print('==== Performance Metrics ====')
print(cat('TDR: ', apply(tdr[1:j,],2,mean), ' TPR: ', apply(tpr[1:j,],2,mean), ' FPR: ', apply(fpr[1:j,],2,mean)))
print(apply(perf_metrics[1:j, ], 2, mean))
print(apply(pc_stable_tests[1:j,],2,mean))
print(apply(cdcd_tests[1:j,],2,mean))
# save(dag, amat, num_edge_tests, sepset, edge_list, pdsep_list, cdcd_obj, pc_obj, pc_stable_obj, file="cdcd01_minmb_andes_2000.RData")    
if(!is.null(unique_vals)){
    save(dag, data, amat, num_edge_tests, sepset, edge_list, pdsep_list, cdcd_obj, pc_obj, pc_stable_obj, file=paste("cdcd01_minmb_", network, '_', sample_size, ".RData", sep=''))    
}else{
    save(dag, data, amat, num_edge_tests, sepset, edge_list, pdsep_list, cdcd_obj, pc_obj, pc_stable_obj, file=paste("cdcd01_minmb_", num_vars[i],  '_', sample_size,".RData", sep=''))
}
# path <- "/scratch/hdd1/home/mschaudh/Second_Chapter/Results_AllComms/"
# setwd(path)
# dir.create(paste(getwd(), '/Sparsity_', sparsity, sep=''))
# setwd(paste(getwd(), '/Sparsity_', sparsity, sep=''))
# write.csv(cbind(tpr, fpr, tdr), paste(getwd(), '/Metrics_', num_vars[i], '_n_200_21_40', '.csv', sep=''))
# write.csv(cbind(pc_cdcd_shd, pc_shd, cdcd_shd), paste(getwd(), '/SHD_', num_vars[i], '_n_200_21_40','.csv', sep=''))
# write.csv(cbind(pc_tests, cdcd_tests), paste(getwd(), '/CI_Tests_', num_vars[i], '_n_200_21_40','.csv', sep=''))
# }

conditionalIndepTest <- function(G, pVals, sepset, cor_mat, amat_cpy, data, unique_vals, i, j_var, S_i, edge_list, pdsep_list, N, l, num_edge_tests){
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
                            # print(cat('i: ',i, 'j:', j, 'S:', as.matrix(power_set[, s]), ' pdsep_mat ', pdsep_mat[, 1], ' l_sub:', l_sub, ' edge_pdsep_list:', unlist(edge_pdsep_list)))
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
                if(is.null(unique_vals))
                    p_val <- gaussCItest(x=i, y=j, S=power_set[, s], suffStat=list(C=cor_mat, n=N))
                else
                    p_val <- disCItest(x=i, y=j, S=power_set[, s], suffStat=list(dm=data, nlev=unique_vals, adaptDF=FALSE))
                test_counter <- test_counter + 1
                # print(cat('l_sub:', l_sub, ' num_edge_tests[l_sub]:', num_edge_tests[l_sub]))
                num_edge_tests[l_sub] <- num_edge_tests[l_sub] + 1
                # if(l_sub<=2)
                    # print(cat('i:', i, 'j:', j, 'S:', power_set[, s], 'p-val:', p_val, 'test counter: ', num_edge_tests[l_sub]))
                if(p_val>=0.05){
                    # if(l_sub<=3)
                        # print(cat('i:', i, 'j:', j, 'S:', power_set[, s], 'p-val:', p_val, 'test counter: ', num_edge_tests[l_sub]))
                    G[i, j] <- G[j, i] <- FALSE
                    # print(sum(G))
                    amat_cpy[i, j] <- amat_cpy[j, i] <- 0
                    sep_set <- union(sep_set, as.vector(power_set[, s]))
                    if(pVals[i, j] < p_val){
                        sepset[[i]][[j]] <- as.vector(power_set[, s])
                        pVals[i, j] <- p_val
                    } 
                    # S_i <- S_i[!S_i %in% j]
                    break
                }
            }  
        }
    }
  return(list(G, sepset, pVals, amat_cpy, edge_list, pdsep_list, num_edge_tests))
}

cdcd_pc <- function(G, pVals, cor_mat, amat, data, bool_amat, sepset, edge_list, pdsep_list, n, vars, num_edge_tests, test_counter, run_time){
    cl <- match.call()
    if(is.null(cor_mat)){
        unique_vals <- apply(data, 2, function(x) length(unique(x)))
    }

    # res_obj <- buildSkeleton(G, pVals, cor_mat, amat, bool_amat, sepset, edge_list, pdsep_list, n, vars, num_edge_tests)
    # full_graph <- graph_from_adjacency_matrix(amat, weighted=TRUE, mode="undirected")
    # coords <- layout.davidson.harel(full_graph)
    graph_list <- list()
    graph_list[[1]] <- amat
    # graph_list[[2]] <- coords
    # graph_list[[3]] <- buildSepset(num_comms=1, !vector(mode="logical", length(vars)), bool_amat, vars, FALSE)
    # res_obj <- buildSkeletonAllComms(G, pVals, NULL, cor_mat, graph_list, data, unique_vals, bool_amat, sepset, edge_list, pdsep_list, n, vars, num_edge_tests)
    res_obj <- buildStableSkeleton(G, pVals, NULL, cor_mat, graph_list, data, unique_vals, bool_amat, sepset, edge_list, pdsep_list, n, vars, num_edge_tests, run_time)
    G <- res_obj[[1]]
    sepset <- res_obj[[2]]
    pVals <- res_obj[[3]]
    amat <- res_obj[[4]]
    num_edge_tests <- res_obj[[8]]
    run_time <- res_obj[[9]]
    non_zero_idx <- which(num_edge_tests!=0)

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
        max.ord = as.integer(length(which(num_edge_tests!=0))+1), n.edgetests = c(test_counter, num_edge_tests[non_zero_idx]), 
        sepset = sepset, pMax = pVals, zMin = matrix(NA, 1, 1))
    cdcd_skel@call <- cl
    cdcd_graph <- switch("relaxed", rand=udag2pdag(cdcd_skel), retry=udag2pdagSpecial(cdcd_skel)$pcObj, relaxed = udag2pdagRelaxed(cdcd_skel, verbose = FALSE, 
                solve.confl = FALSE))
    # run_time$Other <- run_time$Other + proc.time()['elapsed'] - start_time['elapsed']
    return(list(cdcd_graph, run_time))
}

mergeSepset <- function(sepset, pVals, sepset_sub, pVal_sub, c_l, c_k_z, data=NULL, unique_vals=NULL){
 
    for(l_z in c_l){
        for(k_z in c_k_z){
            # print(cat('Merge sepset ', ' c_l ', l_z, ' c_k_z ', k_z))
            # print(cat(' sepset_sub ', length(sepset_sub[[l_z]])))
            if(!is.null(sepset_sub[[l_z]][[k_z]])){
                if(pVals[l_z, k_z] < pVal_sub[l_z, k_z]){
                    sepset[[l_z]][[k_z]] <- sepset_sub[[l_z]][[k_z]]   
                    pVals[l_z, k_z] <- pVals[k_z, l_z] <- pVal_sub[l_z, k_z]
                }
            }
            if(!is.null(sepset_sub[[k_z]][[l_z]])){
                # print(cat('sepset ', sepset[[k_z]][[l_z]], ' sepset_sub ', sepset_sub[[k_z]][[l_z]]))
                if(pVals[k_z, l_z] < pVal_sub[k_z, l_z]){
                    sepset[[k_z]][[l_z]] <- sepset_sub[[k_z]][[l_z]]   
                    pVals[k_z, l_z] <- pVals[l_z, k_z] <- pVal_sub[k_z, l_z]
                }
            }
        }
        for(l_l_z in c_l){
            if(l_z == l_l_z)
                next
            # print(cat('l_z:', l_z, ' l_l_z:', l_l_z, ' sepset_sub:', length(sepset_sub[[l_z]])))
            if(!is.null(sepset_sub[[l_z]][[l_l_z]]) && !is.null(sepset_sub[[l_l_z]][[l_z]])){
                if(!(all(sepset_sub[[l_l_z]][[l_z]] %in% sepset_sub[[l_z]][[l_l_z]]) && all(sepset_sub[[l_z]][[l_l_z]] %in% sepset_sub[[l_l_z]][[l_z]]))){
                    if(!is.null(unique_vals)){
                        p_val_1 <- disCItest(x=l_z, y=l_l_z, S=sepset_sub[[l_z]][[l_l_z]], suffStat=list(dm=data, nlev=unique_vals, adaptDF=FALSE))
                        p_val_2 <- disCItest(x=l_z, y=l_l_z, S=sepset_sub[[l_l_z]][[l_z]], suffStat=list(dm=data, nlev=unique_vals, adaptDF=FALSE))
                    }else{
                        p_val_1 <- gaussCItest(x=l_z, y=l_l_z, S=sepset_sub[[l_z]][[l_l_z]], suffStat=list(C=cor(data), n=nrow(data)))
                        p_val_2 <- gaussCItest(x=l_z, y=l_l_z, S=sepset_sub[[l_l_z]][[l_z]], suffStat=list(C=cor(data), n=nrow(data)))                        
                    }

                    if(p_val_1 < p_val_2){
                        sepset[[l_z]][[l_l_z]] <- NULL
                        sepset[[l_l_z]][[l_z]] <- sepset_sub[[l_l_z]][[l_z]]
                        pVals[l_l_z, l_z] <- pVals[l_z, l_l_z] <- p_val_2
                    }else{
                        sepset[[l_z]][[l_l_z]] <- sepset_sub[[l_z]][[l_l_z]]
                        sepset[[l_l_z]][[l_z]] <- NULL
                        pVals[l_l_z, l_z] <- pVals[l_z, l_l_z] <- p_val_1
                    }
                }
            }
            if(!is.null(sepset_sub[[l_z]][[l_l_z]])){
                if(pVals[l_z, l_l_z] <- pVal_sub[l_z, l_l_z]){
                    sepset[[l_z]][[l_l_z]] <- sepset_sub[[l_z]][[l_l_z]] 
                    pVals[l_z, l_l_z] <- pVals[l_l_z, l_z] <- pVal_sub[l_z, l_l_z]   
                }
            }
            if(!is.null(sepset_sub[[l_l_z]][[l_z]])){
                if(pVals[l_l_z, l_z] <- pVal_sub[l_l_z, l_z]){
                    sepset[[l_l_z]][[l_z]] <- sepset_sub[[l_l_z]][[l_z]]
                    pVals[l_l_z, l_z] <- pVals[l_z, l_l_z] <- pVal_sub[l_l_z, l_z]
                }
            }
        }
    }
    # print(cat('sepset ', sepset[[4]][[5]], ' sepset_sub ', sepset_sub[[4]][[5]]))
    # print(cat('sepset ', sepset[[5]][[4]], ' sepset_sub ', sepset_sub[[5]][[4]]))
 
    return(list(sepset, pVals))
}

buildSkeleton <- function(G, pVals, communities_old, cor_mat, graph_list, data, unique_vals, bool_amat, sepset, edge_list, pdsep_list, n, vars, num_edge_tests){
    amat <- graph_list[[1]]
    coords <- graph_list[[2]]

    full_graph <- graph_from_adjacency_matrix(amat, weighted=TRUE, mode="undirected")
    colors <- rep(1, nrow(amat))

    if(length(vars)<=3){
        # print(cat('Returning with vars', vars))
        if(sum(bool_amat[vars, vars])==6){
            res_obj <- removeInterCommEdges(vars, !vector(mode="logical", length(vars)), G, pVals, sepset, vars, list(bool_amat, NULL), cor_mat, amat, data, unique_vals,
                        edge_list, pdsep_list, n, max_degree=1, num_edge_tests)
            G <- res_obj[[1]]
            sepset <- res_obj[[2]]
            pVals <- res_obj[[3]]
            amat <- res_obj[[4]]
            edge_list <- res_obj[[5]]
            pdsep_list <- res_obj[[6]]
            num_edge_tests <- res_obj[[7]]

            amat_temp <- amat[vars, vars]
            amat_temp[amat_temp!=0] <- 1
            bool_amat[vars, vars] <- amat_temp
            sub_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
            colors <- communities_old$membership
            plot(sub_graph, vertex.color=colors, main=paste("Graph with 3 variables", paste(vars, collapse=" ")))
        }
        # print('num_edge_tests: ', num_edge_tests)
        return(list(G, sepset, pVals, amat, bool_amat, edge_list, pdsep_list, num_edge_tests))
    }

    cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
    if(is.null(E(cor_graph)$weight)){
        return(list(G, sepset, pVals, amat, bool_amat, edge_list, pdsep_list, num_edge_tests))
    }

    communities <- cluster_louvain(cor_graph, weights=abs(E(cor_graph)$weight))
    community_membership <- communities$membership
    colors[vars] <- community_membership + 1
    num_comms <- length(unique(community_membership))
    # if(!is.null(communities_old))
    #     print(cat('old: ', communities_old$names))
    # print(cat('new: ', communities$names))

    if(!is.null(communities_old) && !diffCommunity(communities_old, communities)){
        print(cat('Old comm: ', communities_old$membership))
        print(cat('New comm: ', communities$membership))
        S <- buildSepset(num_comms, community_membership, bool_amat, vars)
        if(num_comms==1){
            communities <- cluster_fast_greedy(cor_graph, weights=abs(E(cor_graph)$weight))
            community_membership <- cutat(communities, no=2)            
            num_comms <- 2
        }
        for(l in 1:num_comms){
            c_l_vars <- community_membership==l
            c_k_vars <- community_membership!=l
            c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
            c_l_z <- c_z[[1]]
            c_k_z <- c_z[[2]]
            # print(cat('c_l_z: ', c_l_z))
            # print(cat('c_k_z: ', c_k_z))
            if(length(c_l_z)>0){
                res_obj <- removeInterCommEdges(c_l_z, c_l_vars, G, pVals, sepset, vars, list(bool_amat, S), cor_mat, amat, data, unique_vals,
                                              edge_list, pdsep_list, n, max_degree=1, num_edge_tests)
                G <- res_obj[[1]]
                sepset <- res_obj[[2]]
                pVals <- res_obj[[3]]
                amat <- res_obj[[4]]
                edge_list <- res_obj[[5]]
                pdsep_list <- res_obj[[6]]
                num_edge_tests <- res_obj[[7]]

                # new_cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
                # full_graph <- graph_from_adjacency_matrix(amat, weighted=TRUE, mode="undirected")
                # plot(full_graph, vertex.color=colors, layout=coords, main=paste("Graph with 3 variables", paste(vars, collapse=" ")))
            }
            amat_temp <- amat[vars, vars]
            amat_temp[amat_temp!=0] <- 1
            bool_amat[vars, vars] <- amat_temp
            sub_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
            colors <- communities$membership
            plot(sub_graph, vertex.color=colors, main=paste("Graph with same communities variables", paste(vars, collapse=" ")))            
        }
        return(list(G, sepset, pVals, amat, bool_amat, edge_list, pdsep_list, num_edge_tests))
    }
    # plot(full_graph, vertex.color=colors, layout=coords, main=paste("Full graph"))
    
    num_edge_tests <- as.integer(num_edge_tests)
    # Perform 0-1 conditional independence tests
    S <- list()
    if(num_comms==1){
        max_degree <- 1L    
        communities <- cluster_fast_greedy(cor_graph, weights=abs(E(cor_graph)$weight))
        community_membership <- cutat(communities, no=2)
        comm_id <- unique(community_membership)
        c_l_vars <- community_membership==comm_id[1]
        c_k_vars <- community_membership!=comm_id[1]

        c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
        c_l_z <- c_z[[1]]
        c_k_z <- c_z[[2]]

        S <- buildSepset(num_comms=2, community_membership, bool_amat, vars)
        res_obj <- removeInterCommEdges(c_l_z, c_l_vars, G, pVals, sepset, vars, list(bool_amat, S), cor_mat, amat, data, unique_vals,
                                        edge_list, pdsep_list, n, max_degree, num_edge_tests)
        G <- res_obj[[1]]
        sepset <- res_obj[[2]]
        pVals <- res_obj[[3]]
        amat <- res_obj[[4]]
        edge_list <- res_obj[[5]]
        pdsep_list <- res_obj[[6]]
        num_edge_tests <- res_obj[[7]]

        amat_temp <- amat[vars, vars]
        amat_temp[amat_temp!=0] <- 1
        bool_amat[vars, vars] <- amat_temp
         
        res_obj <- removeInterCommEdges(c_k_z, c_k_vars, G, pVals, sepset, vars, list(bool_amat, S), cor_mat, amat, data, unique_vals,
                                        edge_list, pdsep_list, n, max_degree, num_edge_tests)
        G <- res_obj[[1]]
        sepset <- res_obj[[2]]
        pVals <- res_obj[[3]]
        amat <- res_obj[[4]]
        edge_list <- res_obj[[5]]
        pdsep_list <- res_obj[[6]]
        num_edge_tests <- res_obj[[7]]

        amat_temp <- amat[vars, vars]
        amat_temp[amat_temp!=0] <- 1
        bool_amat[vars, vars] <- amat_temp

    }else{
        S <- buildSepset(num_comms, community_membership, bool_amat, vars)
        for(l in 1:num_comms){
            max_degree <- 1L
            c_l_vars <- community_membership==l
            c_k_vars <- community_membership!=l
            c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
            c_l_z <- c_z[[1]]
            c_k_z <- c_z[[2]]

            if(length(c_l_z)>0){
                res_obj <- removeInterCommEdges(c_l_z, c_l_vars, G, pVals, sepset, vars, list(bool_amat, S), cor_mat, amat, data, unique_vals,
                                              edge_list, pdsep_list, n, max_degree, num_edge_tests)
                G <- res_obj[[1]]
                sepset <- res_obj[[2]]
                pVals <- res_obj[[3]]
                amat <- res_obj[[4]]
                edge_list <- res_obj[[5]]
                pdsep_list <- res_obj[[6]]
                num_edge_tests <- res_obj[[7]]

                # new_cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
                # full_graph <- graph_from_adjacency_matrix(amat, weighted=TRUE, mode="undirected")
                # plot(full_graph, vertex.color=colors, layout=coords, main=paste("Graph with 3 variables", paste(vars, collapse=" ")))
            }
            amat_temp <- amat[vars, vars]
            amat_temp[amat_temp!=0] <- 1
            bool_amat[vars, vars] <- amat_temp
        }        
    }    

    # Build 0-1 CI graph on each community
    for(l in 1:num_comms){
        c_names <- communities$names
        c_l_vars <- community_membership==l
        c_k_vars <- community_membership!=l
        # print(row.names(bool_amat))
        c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
        c_l_z <- c_z[[1]]
        c_k_z <- c_z[[2]]
        # print(row.names(bool_amat))

        # print(cat('c_l_z', c_l_z, ' c_k_z', c_k_z))
        # print('Communities do not have an edge between them')
        res_obj <- buildSkeleton(G, pVals, communities, cor_mat, list(amat, coords), data, unique_vals, bool_amat, sepset, edge_list, pdsep_list, 
                                       n, as.numeric(union(c_names[c_l_vars], c_k_z)), num_edge_tests)
        G <- res_obj[[1]]
        temp_sepset <- res_obj[[2]]
        pVals <- res_obj[[3]]
        amat <- res_obj[[4]]
        temp_bool_amat <- res_obj[[5]]
        edge_list <-res_obj[[6]]
        pdsep_list <- res_obj[[7]]
        num_edge_tests <- res_obj[[8]]

        if(l>1){
            temp <- bool_amat[c(c_l_z, c_k_z), c(c_l_z, c_k_z)] & temp_bool_amat[c(c_l_z, c_k_z), c(c_l_z, c_k_z)]
            bool_amat[c(c_l_z, c_k_z), c(c_l_z, c_k_z)] <- temp
            colnames(bool_amat) <- row.names(bool_amat) <- 1:nrow(bool_amat)
            c_l_z_temp <- as.numeric(c_l_z)
            c_k_z_temp <- as.numeric(c_k_z)
            sepset <- mergeSepset(sepset, temp_sepset, c_l_z_temp, c_k_z_temp)
        }else{
            bool_amat <- temp_bool_amat
        }

        full_graph <- graph_from_adjacency_matrix(amat, weighted=TRUE, mode="undirected")
        # plot(full_graph, vertex.color=colors, layout=coords, main=paste("Returning from recursion", paste(vars, collapse=" ")))      
    }
  
    # Remove intra-comm & inter-comm edges of higher order i.e., |S|>=2
    max_degree <- max_degree + 1
    max_comm_degree <- 0

    repeat{
        S <- buildSepset(num_comms, community_membership, bool_amat, vars)
        for(l in 1:num_comms){
            c_names <- communities$names
            c_l_vars <- community_membership==l
            c_k_vars <- community_membership!=l
            # print(row.names(bool_amat))
            c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
            c_l_z <- c_z[[1]]

            # print(max(apply(bool_amat[c_l_z, , drop=FALSE], 1, function(x) sum(x))))
            if(length(c_l_z)>0){
                max_conn <- max(apply(bool_amat[c_l_z, , drop=FALSE], 1, function(x) sum(x)))
                # max_conn <- max(unlist(lapply(S[c_l_z], function(x) length(x))))
                if(max_degree<=(max_conn-1)){
                    # print('Removing edges with higher-order CI tests')
                    res_obj <- removeInterCommEdges(c_l_z, c_l_vars, G, pVals, sepset, vars, list(bool_amat, NULL), cor_mat, amat, data, unique_vals,
                                                  edge_list, pdsep_list, n, max_degree, num_edge_tests)
                    G <- res_obj[[1]]
                    sepset <- res_obj[[2]]
                    pVals <- res_obj[[3]]
                    amat <- res_obj[[4]]
                    edge_list <- res_obj[[5]]
                    pdsep_list <-res_obj[[6]]
                    num_edge_tests <- res_obj[[7]]

                    new_cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
                    plot(new_cor_graph, vertex.color=community_membership)            
                    max_comm_degree <- max(max_comm_degree, max_conn)
                    
                    # full_graph <- graph_from_adjacency_matrix(amat, weighted=TRUE, mode="undirected")
                    # plot(full_graph, vertex.color=colors, layout=coords, main=paste("Higher Order CI Tests", paste(vars[c_l_vars], collapse=" ")))
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

removeInterCommEdges <- function(c_z, c_vars, G, pVals, sepset, vars, bool_list, cor_mat, amat, data, unique_vals,  
                        edge_list, pdsep_list, n, max_degree, num_edge_tests, diffComm=TRUE){

    l <- max_degree
    bool_amat <- bool_list[[1]]
    S <- bool_list[[2]]
    bool_amat_cpy <- bool_amat

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

    for(i in c_z){
        # print(cat('c_vars: ', dim(bool_amat), ' c_vars: ', length(c_vars)))
        # c_i_z <- setdiff(c_z, i)
        adj_comm_idx <- which(bool_amat[i, vars[c_vars]]!=0, arr.ind=TRUE)
        adj_non_comm_idx <- which(bool_amat[i, vars[!c_vars]]!=0, arr.ind=TRUE)
        adj_comm <- colnames(bool_amat)[vars[c_vars]]
        adj_comm <- adj_comm[adj_comm_idx]
        adj_non_comm <- colnames(bool_amat)[vars[!c_vars]]
        adj_non_comm <- adj_non_comm[adj_non_comm_idx]
        # # print(c_i_z)
        # print(vars[vars==c_i_z])
        # S_i <- c(adj_comm, adj_non_comm)
        X_j <- c(adj_non_comm, adj_comm)
        # if(length(vars)<=3 | !diffComm){
        #     X_j <- adj_comm
        # }else{
        #     X_j <- adj_non_comm
        # }
        # S[[i]] <- c(adj_comm, adj_non_comm)

        # if(l==1){
        #     print('=============================')
            # print(cat('vars: ', vars, 'i: ', i, 'S[[i]]: ', S[[i]], ' X_j: ' , X_j, ' vars[c_vars]:', vars[c_vars]))
        #     # print(bool_amat[i, vars[!c_vars]])
        # }
        if(length(X_j)==0 | length(S[[i]])==0 | length(S[[i]]) < l){
            # print('Next')
            next
        }else{
          res_obj <- conditionalIndepTest(G, pVals, sepset, cor_mat, amat, data, unique_vals, as.numeric(i), as.numeric(unlist(X_j)), as.numeric(S[[i]]), 
                    edge_list, pdsep_list, n, l, num_edge_tests)
          G <- res_obj[[1]]
          sepset <- res_obj[[2]]
          pVals <- res_obj[[3]]
          amat <- res_obj[[4]]
          edge_list <- res_obj[[5]]
          pdsep_list <- res_obj[[6]]
          num_edge_tests <- res_obj[[7]]
          # print(cat('in function: ', pVals[9, 4]))
        }
        amat_temp <- amat[vars, vars]
        amat_temp[amat_temp!=0] <- 1
        bool_amat[vars, vars] <- amat_temp
    }
    # }
    # print('=====================')
    return(list(G, sepset, pVals, amat, edge_list, pdsep_list, num_edge_tests))
}

findInterCommVars <- function(c1_vars, c2_vars, bool_amat){
    inter_comm_conn <- bool_amat[c1_vars, !c1_vars, drop=FALSE]
    intra_comm_conn <- bool_amat[c1_vars, c1_vars, drop=FALSE]

    n_inter_edges <- apply(inter_comm_conn, 1, sum)
    n_intra_edges <- apply(intra_comm_conn, 1, sum)

    c1_z <- names(which((n_intra_edges>0 & n_inter_edges>0)==TRUE))
    c1_z <- as.numeric(names(which((n_inter_edges>0)==TRUE)))

    inter_comm_conn <- bool_amat[c2_vars, c1_vars, drop=FALSE]
    intra_comm_conn <- bool_amat[c2_vars, c2_vars, drop=FALSE]

    n_inter_edges <- apply(inter_comm_conn, 1, sum)
    n_intra_edges <- apply(intra_comm_conn, 1, sum)

    c2_z <- names(which((n_intra_edges>0 & n_inter_edges>0)==TRUE))
    c2_z <- as.numeric(names(which((n_inter_edges>0)==TRUE)))
    # print(cat('c1_z:',c1_z))
    # print(cat('c2_z:',c2_z))

    return(list(c1_z, c2_z))
}   

cv.test <- function(x,y) {
  CV = sqrt(chisq.test(x, y, correct=FALSE)$statistic /
    (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  # print.noquote("Cramér V / Phi:")
  return(as.numeric(CV))
}

marginalIndepTest <- function(data, N, test_counter, disc_flag=FALSE, run_time){
    temp_amat <- amat <- matrix(1, nrow=ncol(data), ncol=ncol(data))
    
    # start_time <- proc.time()
    if(!disc_flag){
        cor_amat <- cor(data)
    }
    else{
        for(i in 1:ncol(data)){
            amat[i, i] <- 0
            y <- which(amat[i, ]!=0)
            for(j in y){
                cv_val <- cv.test(as.numeric(data[, i]), as.numeric(data[, j]))
                if(cv_val==0)
                    amat[i, j] <- amat[j, i] <- 0.001
                else
                    amat[i, j] <- amat[j, i] <- cv_val    
            }
        }
        unique_vals <- apply(data, 2, function(x) length(unique(x)))
    }
    vars <- 1:nrow(amat)
    start_time <- proc.time()
    for(i in vars){
        amat[i, i] <- 0
        y <- which(amat[i, ]!=0)
        for(j in y){
            if(!disc_flag)
                p_val <- gaussCItest(i, j, S=NULL, suffStat=list(C=cor_amat, n=N))
            else
                p_val <- disCItest(i, j, S=NULL, suffStat=list(dm=data, nlev=unique_vals, adaptDF=FALSE))
        test_counter <- test_counter + 1
        if(p_val>0.05)
            amat[i, j] <- amat[j, i] <- 0
        }
    }
    run_time <- proc.time()['elapsed'] - start_time['elapsed']

    return(list(amat, test_counter, run_time))  
}

buildSepset <- function(num_comms, community_membership, bool_amat, vars, diffComm){
    S <- list()
    for(l in 1:num_comms){
        c_l_vars <- community_membership==l
        c_k_vars <- community_membership!=l
        c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
        c_l_z <- c_z[[1]]
        
        if(!diffComm){
            c_l_z <- c(c_l_z, setdiff(vars[c_l_vars], c_l_z))
        }
        for(var in c_l_z){
            adj_comm_idx <- which(bool_amat[var, vars[c_l_vars]]!=0, arr.ind=TRUE)
            adj_non_comm_idx <- which(bool_amat[var, vars[!c_l_vars]]!=0, arr.ind=TRUE)
            adj_comm <- colnames(bool_amat)[vars[c_l_vars]]
            adj_comm <- adj_comm[adj_comm_idx]
            adj_non_comm <- colnames(bool_amat)[vars[!c_l_vars]]
            adj_non_comm <- adj_non_comm[adj_non_comm_idx]
            S[[var]] <- c(adj_comm, adj_non_comm)
            # print(cat('var: ', var, ' S[[var]]', S[[var]]))
        }        
    }
    return(S)    
}

diffCommunity <- function(communities_old_membership, communities_new_membership, comm_old_names, comm_new_names){
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

buildSkeletonAllComms <- function(G, pVals, communities_old, cor_mat, graph_list, data, unique_vals, bool_amat, sepset, edge_list, pdsep_list, n, vars, num_edge_tests){
    # print(cat('vars: ', vars))
    amat <- graph_list[[1]]
    coords <- graph_list[[2]]

    full_graph <- graph_from_adjacency_matrix(amat, weighted=TRUE, mode="undirected")
    colors <- rep(1, nrow(amat))

    if(length(vars)<=3){
        # print(cat('Returning with vars', vars))
        if(sum(bool_amat[vars, vars])>=4){
            res_obj <- removeInterCommEdges(vars, !vector(mode="logical", length(vars)), G, pVals, sepset, vars, list(bool_amat, NULL), cor_mat, amat, data, unique_vals,
                        edge_list, pdsep_list, n, max_degree=1, num_edge_tests)
            G <- res_obj[[1]]
            sepset <- res_obj[[2]]
            pVals <- res_obj[[3]]
            amat <- res_obj[[4]]
            edge_list <- res_obj[[5]]
            pdsep_list <- res_obj[[6]]
            num_edge_tests <- res_obj[[7]]

            amat_temp <- amat[vars, vars]
            amat_temp[amat_temp!=0] <- 1
            bool_amat[vars, vars] <- amat_temp

            plot(full_graph, vertex.color=colors, layout=coords, main=paste("Graph with 3 variables", paste(vars, collapse=" ")))
        }
        # print('num_edge_tests: ', num_edge_tests)
        return(list(G, sepset, pVals, amat, bool_amat, edge_list, pdsep_list, num_edge_tests))
    }
  
    cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
    if(is.null(E(cor_graph)$weight)){
        return(list(G, sepset, pVals, amat, bool_amat, edge_list, pdsep_list, num_edge_tests))
    }

    communities <- cluster_louvain(cor_graph, weights=abs(E(cor_graph)$weight))
    community_membership <- communities$membership
    colors[vars] <- community_membership + 1
    # plot(full_graph, vertex.color=colors, layout=coords, main=paste("Full graph"))
    
    num_edge_tests <- as.integer(num_edge_tests)

    # Perform 0-1 conditional independence tests
    S <- list()
    num_comms <- length(unique(community_membership))
    
    if(num_comms==1){
        max_degree <- 1L    
        communities <- cluster_fast_greedy(cor_graph, weights=abs(E(cor_graph)$weight))
        community_membership <- cutat(communities, no=2)
        comm_id <- unique(community_membership)
        c_l_vars <- community_membership==comm_id[1]
        c_k_vars <- community_membership!=comm_id[1]

        c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
        c_l_z <- c_z[[1]]
        c_k_z <- c_z[[2]]

        # for(var in c_l_z){
        #     S[[var]] <- c(names(which(bool_amat[var, vars[c_l_vars]]!=0)), names(which(bool_amat[var, vars[!c_l_vars]]!=0)))
        # }

        res_obj <- removeInterCommEdges(c_l_z, c_l_vars, G, pVals, sepset, vars, list(bool_amat, NULL), cor_mat, amat, data, unique_vals,
                                        edge_list, pdsep_list, n, max_degree, num_edge_tests)
        G <- res_obj[[1]]
        sepset <- res_obj[[2]]
        pVals <- res_obj[[3]]
        amat <- res_obj[[4]]
        edge_list <- res_obj[[5]]
        pdsep_list <- res_obj[[6]]
        num_edge_tests <- res_obj[[7]]

        amat_temp <- amat[vars, vars]
        amat_temp[amat_temp!=0] <- 1
        bool_amat[vars, vars] <- amat_temp
         
        # for(var in c_k_z){
        #     S[[var]] <- c(names(which(bool_amat[var, vars[c_l_vars]]!=0)), names(which(bool_amat[var, vars[!c_l_vars]]!=0)))
        # }

        res_obj <- removeInterCommEdges(c_k_z, c_k_vars, G, pVals, sepset, vars, list(bool_amat, NULL), cor_mat, amat, data, unique_vals,
                                        edge_list, pdsep_list, n, max_degree, num_edge_tests)
        G <- res_obj[[1]]
        sepset <- res_obj[[2]]
        pVals <- res_obj[[3]]
        amat <- res_obj[[4]]
        edge_list <- res_obj[[5]]
        pdsep_list <- res_obj[[6]]
        num_edge_tests <- res_obj[[7]]

        amat_temp <- amat[vars, vars]
        amat_temp[amat_temp!=0] <- 1
        bool_amat[vars, vars] <- amat_temp

    }else{
        # for(l in 1:num_comms){
        #     c_l_vars <- community_membership==l
        #     c_k_vars <- community_membership!=l
        #     c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
        #     c_l_z <- c_z[[1]]
        #     for(var in c_l_z){
        #         S[[var]] <- c(names(which(bool_amat[var, vars[c_l_vars]]!=0)), names(which(bool_amat[var, vars[!c_l_vars]]!=0)))
        #     }
        # }

        for(l in 1:num_comms){
            max_degree <- 1L
            c_l_vars <- community_membership==l
            c_k_vars <- community_membership!=l
            c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
            c_l_z <- c_z[[1]]
            c_k_z <- c_z[[2]]

            if(length(c_l_z)>0){
                res_obj <- removeInterCommEdges(c_l_z, c_l_vars, G, pVals, sepset, vars, list(bool_amat, NULL), cor_mat, amat, data, unique_vals,
                                              edge_list, pdsep_list, n, max_degree, num_edge_tests)
                G <- res_obj[[1]]
                sepset <- res_obj[[2]]
                pVals <- res_obj[[3]]
                amat <- res_obj[[4]]
                edge_list <- res_obj[[5]]
                pdsep_list <- res_obj[[6]]
                num_edge_tests <- res_obj[[7]]

                # new_cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
                full_graph <- graph_from_adjacency_matrix(amat, weighted=TRUE, mode="undirected")
            }
            amat_temp <- amat[vars, vars]
            amat_temp[amat_temp!=0] <- 1
            bool_amat[vars, vars] <- amat_temp
        }        
    }

    # for(l in 1:length(unique(community_membership))){
    #     # print('0-1 CI tests')

    #     if(length(c_l_z)>0 & length(c_k_z)>0){
    #         for(i in c_l_z){
    #             adj_comm <- names(which(bool_amat[i, vars[c_l_vars]]!=0))
    #             adj_non_comm <- names(which(bool_amat[i, vars[!c_l_vars]]!=0))
        
    #             S[[i]] <- c(adj_comm, adj_non_comm)
    #         }
    #         res_obj <- removeInterCommEdges(list(c_l_z, FALSE), c_l_vars, G, pVals, sepset, vars, list(bool_amat, S), cor_mat, amat, data, unique_vals,
    #                                       edge_list, pdsep_list, n, max_degree, num_edge_tests)
    #         G <- res_obj[[1]]
    #         sepset <- res_obj[[2]]
    #         pVals <- res_obj[[3]]
    #         amat <- res_obj[[4]]
    #         edge_list <- res_obj[[5]]
    #         pdsep_list <- res_obj[[6]]
    #         num_edge_tests <- res_obj[[7]]

    #         # new_cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
    #         full_graph <- graph_from_adjacency_matrix(amat, weighted=TRUE, mode="undirected")
    #         # plot(full_graph, vertex.color=colors, layout=coords, main="Order 1 CI Tests")
    #     }else{
    #         communities <- cluster_fast_greedy(cor_graph, weights=abs(E(cor_graph)$weight))
    #         community_membership <- cutat(communities, no=2)
    #         num_comms <- unique(community_membership)

    #         c_l_vars <- community_membership==num_comms[1]
    #         c_k_vars <- community_membership!=num_comms[1]
    #         c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
    #         c_l_z <- c_z[[1]]
    #         c_k_z <- c_z[[2]]
            
    #         res_obj <- removeInterCommEdges(list(c_l_z, FALSE), c_l_vars, G, pVals, sepset, vars, bool_amat, cor_mat, amat, data, unique_vals,
    #                                         edge_list, pdsep_list, n, max_degree, num_edge_tests)
    #         G <- res_obj[[1]]
    #         sepset <- res_obj[[2]]
    #         pVals <- res_obj[[3]]
    #         amat <- res_obj[[4]]
    #         edge_list <- res_obj[[5]]
    #         pdsep_list <- res_obj[[6]]
    #         num_edge_tests <- res_obj[[7]]
            
    #         res_obj <- removeInterCommEdges(list(c_k_z, FALSE), c_k_vars, G, pVals, sepset, vars, bool_amat, cor_mat, amat, data, unique_vals,
    #                                         edge_list, pdsep_list, n, max_degree, num_edge_tests)
    #         G <- res_obj[[1]]
    #         sepset <- res_obj[[2]]
    #         pVals <- res_obj[[3]]
    #         amat <- res_obj[[4]]
    #         edge_list <- res_obj[[5]]
    #         pdsep_list <- res_obj[[6]]
    #         num_edge_tests <- res_obj[[7]]
            
    #         # new_cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
    #         # plot(new_cor_graph, vertex.color=community_membership)
    #         colors[vars] <- community_membership + 1
    #         full_graph <- graph_from_adjacency_matrix(amat, weighted=TRUE, mode="undirected")
    #         # plot(full_graph, vertex.color=colors, layout=coords, main=paste("Order 1 CI Tests with 2 Communities", paste(vars[c_l_vars], collapse=" ")))
    #     }
    #     amat_temp <- amat[vars, vars]
    #     amat_temp[amat_temp!=0] <- 1
    #     bool_amat[vars, vars] <- amat_temp
    # }

  
    # Build 0-1 CI graph on each community
    S <- list()
    for(l in 1:num_comms){
        c_names <- communities$names
        c_l_vars <- community_membership==l
        # print(row.names(bool_amat))

        # print(cat('c_l_z', c_l_z, ' c_k_z', c_k_z))
        # print('Communities do not have an edge between them')
        res_obj <- buildSkeletonAllComms(G, pVals, communities, cor_mat, list(amat, coords), data, unique_vals, bool_amat, sepset, edge_list, pdsep_list, 
                                       n, as.numeric(c_names[c_l_vars]), num_edge_tests)
        G <- res_obj[[1]]
        sepset <- res_obj[[2]]
        pVals <- res_obj[[3]]
        amat <- res_obj[[4]]
        bool_amat <- res_obj[[5]]
        edge_list <-res_obj[[6]]
        pdsep_list <- res_obj[[7]]
        num_edge_tests <- res_obj[[8]]

        full_graph <- graph_from_adjacency_matrix(amat, weighted=TRUE, mode="undirected")
        # plot(full_graph, vertex.color=colors, layout=coords, main=paste("Returning from recursion", paste(vars, collapse=" ")))      
    }
  
    # Remove intra-comm & inter-comm edges of higher order i.e., |S|>=2
    max_degree <- max_degree + 1
    max_comm_degree <- 0

    repeat{
        # for(l in 1:num_comms){
        #     c_l_vars <- community_membership==l
        #     c_k_vars <- community_membership!=l
        #     c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
        #     c_l_z <- c_z[[1]]
        #     for(var in c_l_z){
        #             S[[var]] <- c(names(which(bool_amat[var, vars[c_l_vars]]!=0)), names(which(bool_amat[var, vars[!c_l_vars]]!=0)))
        #     }        
        # }

        for(l in 1:num_comms){
            c_names <- communities$names
            c_l_vars <- community_membership==l
            c_k_vars <- community_membership!=l
            # print(row.names(bool_amat))
            c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
            c_l_z <- c_z[[1]]

            # print(max(apply(bool_amat[c_l_z, , drop=FALSE], 1, function(x) sum(x))))
            if(length(c_l_z)>0){
                max_conn <- max(apply(bool_amat[c_l_z, , drop=FALSE], 1, function(x) sum(x)))
                # max_conn <- max(unlist(lapply(S[c_l_z], function(x) length(x))))
                if(max_degree<=(max_conn-1)){
                    # print('Removing edges with higher-order CI tests')
                    res_obj <- removeInterCommEdges(c_l_z, c_l_vars, G, pVals, sepset, vars, list(bool_amat, NULL), cor_mat, amat, data, unique_vals,
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
                    
                    full_graph <- graph_from_adjacency_matrix(amat, weighted=TRUE, mode="undirected")
                    plot(full_graph, vertex.color=colors, layout=coords, main=paste("Higher Order CI Tests", paste(vars[c_l_vars], collapse=" ")))
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

buildSkeletonOld <- function(G, pVals, cor_mat, amat, bool_amat, sepset, edge_list, pdsep_list, n, vars, num_edge_tests){
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
    
    communities <- cluster_fast_greedy(cor_graph, weights=abs(E(cor_graph)$weight))
    community_membership <- cutat(communities, no=2)

    # communities <- cluster_louvain(cor_graph, weights=abs(E(cor_graph)$weight))
    # community_membership <- communities$membership
    # plot(cor_graph, vertex.color=community_membership)
    num_edge_tests <- as.integer(num_edge_tests)
    
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

            amat_temp <- amat[vars, vars]
            amat_temp[amat_temp!=0] <- 1
            bool_amat[vars, vars] <- amat_temp        
        }
    }
    # bool_amat <- amat
    # bool_amat[bool_amat!=0] <- 1

    # Update c_z and then build 0-1 CI graph on each community
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
        res_obj <- buildSkeleton(G, pVals, cor_mat, amat, bool_amat, sepset, edge_list, pdsep_list, 
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
                max_conn <- max(apply(bool_amat[c_l_z, vars, drop=FALSE], 1, function(x) sum(x)))
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
        # bool_amat <- amat
        # bool_amat[bool_amat!=0] <- 1 
        max_degree <- max_degree + 1
        if(max_degree > max_comm_degree)   break
    }   
    # print('num_edge_tests: ', num_edge_tests)
    return(list(G, sepset, pVals, amat, bool_amat, edge_list, pdsep_list, num_edge_tests))
}

removeInterCommEdgesV1 <- function(c_z, c_vars, G, pVals, sepset, vars, bool_amat, cor_mat, amat, 
                        edge_list, pdsep_list, n, max_degree, num_edge_tests){
    # print(c_z)
    # for(l in 1:max_degree){
    l <- max_degree
    bool_amat_cpy <- bool_amat
    S <- list()

    # Splitting graph by removing edges with |S|=1
    if(l==1){
        for(i in c_z){
            adj_comm <- names(which(bool_amat[i, vars[c_vars]]!=0))
            adj_non_comm <- names(which(bool_amat[i, vars[!c_vars]]!=0))
            S[[i]] <- c(adj_comm, adj_non_comm)
        }
        for(i in c_z){
            print('=============================')
            adj_comm <- names(which(bool_amat[i, vars[c_vars]]!=0))
            adj_non_comm <- names(which(bool_amat[i, vars[!c_vars]]!=0))
            # print(c_i_z)
            # print(vars[vars==c_i_z])
            # S_i <- c(adj_comm, adj_non_comm)
            X_j <- c(adj_non_comm, adj_comm)
            print(cat('l: ', l, 'S_i: ', S[[i]], 'X_j: ', X_j, ' vars[c_vars]:', vars[c_vars]))
            if(length(X_j)==0 | length(S[[i]]) < l){
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
    }else{
        # Going bottom-up with higher-order CI tests
        s_adj_comm <- s_adj_non_comm <- s_adj_p <- list()
        p <- 1:nrow(bool_amat)
        p <- p[!p %in% vars]
        
        for(i in c_z){
            s_adj_comm[[i]] <- names(which(bool_amat[i, vars[c_vars]]!=0))
            s_adj_non_comm[[i]] <- names(which(bool_amat[i, vars[!c_vars]]!=0))
            s_adj_p[[i]] <- names(which(bool_amat[i, p]!=0))
        }
        for(i in c_z){
            adj_comm <- names(which(bool_amat[i, vars[c_vars]]!=0))
            adj_non_comm <- names(which(bool_amat[i, vars[!c_vars]]!=0))
            adj_p <- names(which(bool_amat[i, p]!=0))

            print('=============================')

            for(j in adj_non_comm){
                S[[i]] <- c(s_adj_comm[[i]], s_adj_non_comm[[i]][!s_adj_non_comm[[i]] %in% j])
                if(length(S[[i]]) < l){
                    # print('Next')
                    next
                }else{
                    print(cat('l: ', l, 'S_i: ', S[[i]], 'X_j: ', j, ' vars[c_vars]:', vars[c_vars]))
                    res_obj <- conditionalIndepTest(G, pVals, sepset, cor_mat, amat, as.numeric(i), as.numeric(j), as.numeric(S[[i]]), 
                        edge_list, pdsep_list, n, l, num_edge_tests)
                    G <- res_obj[[1]]
                    sepset <- res_obj[[2]]
                    pVals <- res_obj[[3]]
                    amat <- res_obj[[4]]
                    edge_list <- res_obj[[5]]
                    pdsep_list <- res_obj[[6]]
                    num_edge_tests <- res_obj[[7]]
                }
            }
            for(j in adj_p){
                S[[i]] <- c(s_adj_comm[[i]], s_adj_p[[i]][!s_adj_p[[i]] %in% j])
                if(length(S[[i]]) < l){
                    # print('Next')
                    next
                }else{
                    print(cat('l: ', l, 'S_i: ', S[[i]], 'X_j: ', j, ' vars[c_vars]:', vars[c_vars]))
                    res_obj <- conditionalIndepTest(G, pVals, sepset, cor_mat, amat, as.numeric(i), as.numeric(j), as.numeric(S[[i]]), 
                        edge_list, pdsep_list, n, l, num_edge_tests)
                    G <- res_obj[[1]]
                    sepset <- res_obj[[2]]
                    pVals <- res_obj[[3]]
                    amat <- res_obj[[4]]
                    edge_list <- res_obj[[5]]
                    pdsep_list <- res_obj[[6]]
                    num_edge_tests <- res_obj[[7]]
                }
            }

            for(j in adj_comm){
                S[[i]] <- c(s_adj_comm[[i]][!s_adj_comm[[i]] %in% j], s_adj_non_comm[[i]], s_adj_p[[i]])
                if(length(S[[i]]) < l){
                    # print('Next')
                    next
                }else{
                    print(cat('l: ', l, 'S_i: ', S[[i]], 'X_j: ', j, ' vars[c_vars]:', vars[c_vars]))
                    res_obj <- conditionalIndepTest(G, pVals, sepset, cor_mat, amat, as.numeric(i), as.numeric(j), as.numeric(S[[i]]), 
                        edge_list, pdsep_list, n, l, num_edge_tests)
                    G <- res_obj[[1]]
                    sepset <- res_obj[[2]]
                    pVals <- res_obj[[3]]
                    amat <- res_obj[[4]]
                    edge_list <- res_obj[[5]]
                    pdsep_list <- res_obj[[6]]
                    num_edge_tests <- res_obj[[7]]
                }
            }
            amat_temp <- amat[c(vars, adj_p), c(vars, adj_p)]
            amat_temp[amat_temp!=0] <- 1
            bool_amat[c(vars, adj_p), c(vars, adj_p)] <- amat_temp
        }
    }
    # print('=====================')
    return(list(G, sepset, pVals, amat, edge_list, pdsep_list, num_edge_tests))
}

buildSkeletonBkp <- function(G, pVals, cor_mat, amat, bool_amat, sepset, edge_list, pdsep_list, n, vars, num_edge_tests){
    # print(cat('vars: ', vars))
    if(length(vars)<=3){
        # print(cat('Returning with vars', vars))
        # print('num_edge_tests: ', num_edge_tests)
        return(list(G, sepset, pVals, amat, bool_amat, edge_list, pdsep_list, num_edge_tests))
    }

    cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
    if(is.null(E(cor_graph)$weight))
        return(list(G, sepset, pVals, amat, bool_amat, edge_list, pdsep_list, num_edge_tests))
    
    communities <- cluster_fast_greedy(cor_graph, weights=abs(E(cor_graph)$weight))
    community_membership <- cutat(communities, no=2)
    # plot(cor_graph, vertex.color=community_membership)
    # communities <- cluster_louvain(cor_graph, weights=abs(E(cor_graph)$weight))
    # community_membership <- communities$membership
    num_edge_tests <- as.integer(num_edge_tests)
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

            new_cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
            # plot(new_cor_graph, vertex.color=community_membership)
        }
    }
    bool_amat <- amat
    bool_amat[bool_amat!=0] <- 1

    # Update c_z and then build 0-1 CI graph on each community
    for(l in 1:length(unique(community_membership))){
        c_names <- communities$names
        c_l_vars <- community_membership==l
        c_k_vars <- community_membership!=l
        # print(row.names(bool_amat))
        c_z <-findInterCommVars(c_l_vars, c_k_vars, bool_amat[vars, vars])
        c_l_z <- c_z[[1]]
        c_k_z <- c_z[[2]]
        
        # If communities do not share any edges 
        if(length(c_l_z)==0 | length(c_k_z)==0){
            # print(cat('c_l_z', c_l_z, ' c_k_z', c_k_z))
            # print('Communities do not have an edge between them')
            res_obj <- buildSkeleton(G, pVals, cor_mat, amat, bool_amat, sepset, edge_list, pdsep_list, 
                        n, as.numeric(c_names[c_l_vars]), num_edge_tests)
            G <- res_obj[[1]]
            sepset <- res_obj[[2]]
            pVals <- res_obj[[3]]
            amat <- res_obj[[4]]
            bool_amat <- res_obj[[5]]
            edge_list <-res_obj[[6]]
            pdsep_list <- res_obj[[7]]
            num_edge_tests <- res_obj[[8]]        

            # res_obj <- buildSkeleton(G, pVals, cor_mat, amat, bool_amat, sepset, edge_list, pdsep_list, 
            #             n, as.numeric(c_names[c_k_vars]), num_edge_tests)
            # G <- res_obj[[1]]
            # sepset <- res_obj[[2]]
            # pVals <- res_obj[[3]]
            # amat <- res_obj[[4]]
            # bool_amat <- res_obj[[5]]
            # edge_list <-res_obj[[6]]
            # pdsep_list <- res_obj[[7]]
            # num_edge_tests <- res_obj[[8]]        
        }else{
            # print('Communities have edge between them')
            # z_flag <- ifelse(length(c_l_z) < length(c_k_z), 1, 0)
            # if(z_flag){
            #     new_c_vars <- union(as.numeric(c_names[c_l_vars]), c_l_z)
            # }else{
            #     new_c_vars <- union(as.numeric(c_names[c_l_vars]), c_k_z)
            # }
            new_c_vars <- as.numeric(c_names[c_l_vars])
            res_obj <- buildSkeleton(G, pVals, cor_mat, amat, bool_amat, sepset, edge_list, pdsep_list, 
                        n, new_c_vars, num_edge_tests)
            G <- res_obj[[1]]
            sepset <- res_obj[[2]]
            pVals <- res_obj[[3]]
            amat <- res_obj[[4]]
            bool_amat <- res_obj[[5]]
            edge_list <-res_obj[[6]]
            pdsep_list <- res_obj[[7]]
            num_edge_tests <- res_obj[[8]]        
        }
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

                    new_cor_graph <- graph_from_adjacency_matrix(amat[vars, vars], weighted=TRUE, mode="undirected")
                    # plot(new_cor_graph, vertex.color=community_membership)            
                    max_comm_degree <- max(max_comm_degree, max_conn)
                }
            }
        }
        bool_amat <- amat
        bool_amat[bool_amat!=0] <- 1 
        max_degree <- max_degree + 1
        if(max_degree > max_comm_degree)   break
    }   
    # print('num_edge_tests: ', num_edge_tests)
    return(list(G, sepset, pVals, amat, bool_amat, edge_list, pdsep_list, num_edge_tests))
}