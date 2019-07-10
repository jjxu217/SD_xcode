/* The subroutine adds the newly formed cut to the cutsType structure. For the optimality cuts, if there is no room in the cutsType structure
 * then the reduceCuts() subroutine is invoked to remove the 'loose' and 'old' cuts. For the feasibility cuts, we check if the new cut is a
 * duplicate of existing cut before it is added to the pool. */

if ( dropCut(cell, cell->iCutIdx) ){
    errMsg("algorithm", "reduceCuts", "failed to drop a cut", 0);
    return -1;
}
cell->cuts->vals[cell->iCutIdx] = cut;

return cell->iCutIdx;





int addCut2Master(oneProblem *master, oneCut *cut, vector vectX, int lenX) {
    intvec     indices;
    int     cnt;
    static int cummCutNum = 0;
    double new_beta[2 * lenX + 1];
    
    /* Cut right-hand side */
    if ( config.MASTER_TYPE == PROB_QP || config.MASTER_TYPE == PROB_MILP)
        cut->alphaIncumb = cut->alpha - vXv(cut->beta, vectX, NULL, lenX);
    
    /* Set up the cut name */
    sprintf(cut->name, "cut_%04d", cummCutNum++);
    
    /* add the cut to the cell cuts structure as well as on the solver */
    if (config.MASTER_TYPE == PROB_MILP){
        /* Set up indices */
        if (!(indices = (intvec) arr_alloc(2 * lenX + 1, int)))
            errMsg("Allocation", "addcut2Master", "fail to allocate memory to coefficients of beta",0);
        for (cnt = 1; cnt <= lenX; cnt++){
            indices[cnt] = cnt - 1;
            indices[lenX + cnt] = lenX + cnt;
            new_beta[cnt] = cut->beta[cnt];
            new_beta[lenX + cnt] = -cut->beta[cnt];
        }
        indices[0] = lenX;
        new_beta[0] = cut->beta[0];
        
        if ( addRow(master->lp, 2 * lenX + 1, cut->alphaIncumb, GE, 0, indices, new_beta, cut->name) ) {
            errMsg("solver", "addcut2Master", "failed to add new row to problem in solver", 0);
            return 1;
        }
    }
    else{
        /* Set up indices */
        if (!(indices = (intvec) arr_alloc(lenX + 1, int)))
            errMsg("Allocation", "addcut2Master", "fail to allocate memory to coefficients of beta",0);
        for (cnt = 1; cnt <= lenX; cnt++)
            indices[cnt] = cnt - 1;
        indices[0] = lenX;
        
        if ( addRow(master->lp, lenX + 1, cut->alphaIncumb, GE, 0, indices, cut->beta, cut->name) ) {
            errMsg("solver", "addcut2Master", "failed to add new row to problem in solver", 0);
            return 1;
        }
    }
    
    cut->rowNum = master->mar++;
    
#ifdef CUT_CHECK
    writeProblem(master->lp,"master_wNewCut.lp");
#endif
    
    mem_free(indices);
    return 0;
}//END addCuts2Master()




int addL1Norm (probType *prob, batchSummary *batch) {
    int i, j, cnt;
    int status;
    int NumABS;
    intvec lb, ub;
    vector obj, rhs, rmatbeg, rmatind, rmatval;
    string sense;
    
    NumABS = prob->sp->macsz * batch->cnt;
    
    obj = (vector) arr_alloc(NumABS, double);
    lb = (intvec) arr_alloc(NumABS, int);
    ub = (intvec) arr_alloc(NumABS, int);
    
    for (i = 0; i < NumABS; i++){
        obj[i] = batch->quadScalar;
        lb[i] = 0;
        ub[i] = 1;
    }
    /*add new variables: d_i = abs(X_i - incumbX_i)*/
    status = CPXnewcols (env, batch->sp->lp, NumABS, obj, lb, ub, NULL, NULL);
    if(status)
        errMsg("addL1Norm", "CPXnewcols", "failed to add new cols", 0);
    
    /*add the constraints: d_i >= X_i - incumbX_i, d_i >= -X_i +incumbX_i, which are
     -x_i +d_i >= -incumbX_i,
     x_i +d_i >= incumbX_i,*/
    rhs = (vector) arr_alloc(2*NumABS, double);
    sense = (string) arr_alloc(2*NumABS, char);
    
    rmatbeg = (vector) arr_alloc(2*NumABS, double);
    rmatind = (vector) arr_alloc(4*NumABS, double);
    rmatval = (vector) arr_alloc(4*NumABS, double);
    
    for (i = 0; i < batch->cnt; i++){
        for (j = 0; j < prob->sp->macsz; j++){
            cnt = i * batch->cnt + j;
            rhs[2 * cnt] = -batch->incumbX[j+1], rhs[2 * cnt+ 1] = batch->incumbX[j+1];
            sense[2 * cnt] = 'G', sense[2 * cnt + 1] = 'G';
            rmatbeg[2 * cnt] = 4 * cnt, rmatbeg[2 * cnt + 1] = 4 * cnt + 2;
            
            rmatind[4 * cnt] = rmatind[4 * cnt + 2] = i * (prob->sp->macsz + 1) + j;
            rmatind[4 * cnt + 1] = rmatind[4 * cnt + 3] = batch->cnt * (prob->sp->macsz + 1) + cnt
            
            rmatval[4 * cnt] = -1;
            rmatval[4 * cnt + 1] = rmatval[4 * cnt + 2] = rmatval[4 * cnt + 3] = 1
            
        }
    }
    
    status = CPXaddrows (env, batch->sp->lp, 0, 2*NumABS, 6*NumABS, rhs,
                         sense, rmatbeg, rmatind, rmatval,
                         NULL, NULL);
}



char solname[BLOCKSIZE];
sprintf(solname, "%s%s/%s.lp", inputDir, probName, probName);
