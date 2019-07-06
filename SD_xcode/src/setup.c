/*
 * setup.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 
 Modified by Jiajun on Feb, 27,2019
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "twoSD.h"

extern configType config;

int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell, batchSummary **batch, vector *meanSol) {
	vector	lb = NULL;
	int 	t;

	/* setup mean value problem which will act as reference for all future computations */
	(*meanSol) = meanProblem(orig, stoc);
	if ( meanSol == NULL ) {
		errMsg("setup", "setupAlgo", "failed to setup and solve mean value problem", 0);
		return 1;
	}

	/* calculate lower bounds for each stage */
	lb = calcLowerBound(orig, tim, stoc);  //Jiajun Check, new LB in IP?
	if ( lb == NULL )  {
		errMsg("setup", "setupAlgo", "failed to compute lower bounds on stage problem", 0);
		mem_free(lb); return 1;
	}
	/* decompose the problem into master and subproblem */
	(*prob) = newProb(orig, stoc, tim, lb, config.TOLERANCE);
	if ( (*prob) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to update probType with elements specific to algorithm", 0);
		mem_free(lb); return 1;
	}

#ifdef DECOMPOSE_CHECK
	printDecomposeSummary(stdout, orig->name, tim, (*prob));
#endif

	/* ensure that we have a linear programs at all stages.*/
	t = 0;
	while ( t < tim->numStages ) {
		if ( (*prob)[t]->sp->type  != PROB_LP )
			printf("Note :: Stage-%d problem is a mixed-integer program.\n", t);
        t++;
	}
    //Jiajun TODO: Will solve first stage as a MIP problem.check newCell
	/* create the cells which will be used in the algorithms */
	(*cell) = newCell(stoc, (*prob), (*meanSol));
	if ( (*cell) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to create the necessary cell structure", 0);
		return 1;
	}

	if ( config.NUM_REPS > 1 )
		(*batch)  = newBatchSummary((*prob)[0], config.NUM_REPS);

	mem_free(lb);
	return 0;
}//END setupAlgo()

/* This function is used to create cells used in the algorithm */
cellType *newCell(stocType *stoc, probType **prob, vector xk) {
	cellType    *cell = NULL;
	int			length;

	/* allocate memory to all cells used in the algorithm. The first cell belongs to the master problem, while the rest correspond to each of the
	 * sub-agents in the problem.  */
	if (!(cell = (cellType *) mem_malloc(sizeof(cellType))) )
		errMsg("Memory allocation", "newCell", "failed to allocate memory to cell",0);
	cell->master = cell->subprob = NULL;
	cell->candidX = cell->incumbX = NULL;
	cell->piM = cell->djM = NULL;
	cell->cuts = cell->fcuts = NULL;
	cell->lambda = NULL; cell->sigma = NULL; cell->delta = NULL; cell->omega = NULL;
	cell->pi_ratio = NULL;

	/* setup the master problem */
	cell->master = newMaster(prob[0]->sp, prob[0]->lb);
	if ( cell->master == NULL ) {
		errMsg("setup", "newCell", "failed to setup the master problem", 0);
		return NULL;
	}
	/* setup the subproblem */
	cell->subprob = newSubprob(prob[1]->sp);

	/* -+-+-+-+-+-+-+-+-+-+-+ Allocating memory to other variables that belongs to master mcell +-+-+-+-+-+-+-+-+-+- */
	cell->k 	= 0;
	cell->LPcnt = 0;
	if (prob[0]->lb == 0)
		cell->lbType = TRIVIAL;
	else
		cell->lbType = NONTRIVIAL;
	cell->lb = prob[0]->lb;

	/* candidate solution and estimates */
	cell->candidX 			= duplicVector(xk, prob[0]->num->cols);
	cell->candidEst 		= prob[0]->lb + vXvSparse(cell->candidX, prob[0]->dBar);

	/* incumbent solution and estimates */
    //Jiajun: add master type MIP
	if (config.MASTER_TYPE == PROB_QP || config.MASTER_TYPE == PROB_MIQP || config.MASTER_TYPE == PROB_MILP) {
		cell->incumbX   = duplicVector(xk, prob[0]->num->cols);
		cell->incumbEst = cell->candidEst;
		cell->quadScalar= config.MIN_QUAD_SCALAR;     						/* The quadratic scalar, 'sigma'*/
		cell->iCutIdx   = 0;
		cell->iCutUpdt  = 0;
		cell->incumbChg = TRUE;
	}
	else {
		cell->incumbX   = NULL;
		cell->incumbEst = 0.0;
		cell->quadScalar= 0.0;
		cell->iCutIdx   = -1;
		cell->iCutUpdt  = -1;
		cell->incumbChg = FALSE;
	}
	cell->gamma 			= 0.0;
	cell->normDk_1 			= 0.0;
	cell->normDk 			= 0.0;
    
    cell->RepeatedTime = 0;
    if ( !(cell->argmax_best_candid = (vector) arr_alloc(config.MAX_ITER, double)) )
        errMsg("allocation", "newCell", "cell->argmax_best_candid", 0);
    if ( !(cell->pi_best_candid  = (intvec) arr_alloc(config.MAX_ITER, int)) )
        errMsg("allocation", "newCell", "cell->pi_best_candid", 0);
    
    if ( !(cell->argmax_best_incumb = (vector) arr_alloc(config.MAX_ITER, double)) )
        errMsg("allocation", "newCell", "cell->argmax_best_incumb", 0);
    if ( !(cell->pi_best_incumb  = (intvec) arr_alloc(config.MAX_ITER, int)) )
        errMsg("allocation", "newCell", "cell->pi_best_incumb", 0);
    
    
	/* lower bounding approximations held in cuts structure */
	cell->maxCuts = config.CUT_MULT * prob[0]->num->cols + 3;
	cell->cuts 	  = newCuts(cell->maxCuts);

	/* solution parts of the cell */
	if ( !(cell->djM = (vector) arr_alloc(prob[0]->num->cols + 2, double)) )
		errMsg("allocation", "newCell", "cell->di", 0);
	if ( !(cell->piM = (vector) arr_alloc(prob[0]->num->rows + cell->maxCuts + 1, double)) )
		errMsg("allocation", "newCell", "cell->piM", 0);

	/* stochastic elements: we need more room to store basis information when the cost coefficients are random. */
	if ( prob[1]->num->rvdOmCnt > 0 )
		length = prob[1]->num->rvdOmCnt*config.MAX_ITER + config.MAX_ITER / config.TAU + 1;
	else
		length = config.MAX_ITER + config.MAX_ITER / config.TAU + 1;
    //Jiajun BUG:
    cell->basis  = newBasisType(2*config.MAX_ITER, prob[1]->num->cols, prob[1]->num->rows, WORDLENGTH);
    //Jiajun: if the RHS doesn't contain RV, we don't need lambda/delta 
    if (prob[1]->num->rvRowCnt > 0){
        cell->lambda = newLambda(length, 0, prob[1]->num->rvRowCnt);
        cell->delta  = newDelta(length);
    }
	cell->sigma  = newSigma(length, prob[1]->num->cntCcols, 0);
	cell->omega  = newOmega(prob[1]->num->numRV, config.MAX_ITER);

	cell->optFlag 			= FALSE;

	/* Dual stability test is disabled DUAL_STABILITY is false. */
	if ( !config.DUAL_STABILITY ) {
		cell->dualStableFlag = TRUE;
		cell->pi_ratio = NULL;
	}
	else {
		cell->dualStableFlag 	= FALSE;
		if ( !(cell->pi_ratio = (vector) arr_alloc(config.SCAN_LEN, double)) )
			errMsg("allocation", "newCell", "cell->pi_ratio", 0);
	}

	cell->spFeasFlag = TRUE;
	cell->fcuts		= newCuts(cell->maxCuts);
	cell->fcutsPool = newCuts(cell->maxCuts);
	cell->feasCnt 		= 0;
	cell->infeasIncumb 	= FALSE;
	cell->fUpdt[0] = cell->fUpdt[1] = 0;

	cell->time.repTime = cell->time.iterTime = cell->time.masterIter = cell->time.subprobIter = cell->time.optTestIter = cell->time.argmaxIter = 0.0;
	cell->time.iterAccumTime = cell->time.masterAccumTime = cell->time.subprobAccumTime = cell->time.optTestAccumTime = cell->time.argmaxAccumTime = 0.0;

    
	/* construct the QP using the current incumbent */
	if ( config.MASTER_TYPE == PROB_QP ) {
		if ( constructQP(prob[0], cell, cell->incumbX, cell->quadScalar) ) {
			errMsg("setup", "newCell", "failed to change the right-hand side after incumbent change", 0);
			return NULL;
		}
		cell->incumbChg = FALSE;
#if defined(SETUP_CHECK)
		if ( writeProblem(cell->master->lp, "newQPMaster.lp") ) {
			errMsg("write problem", "new_master", "failed to write master problem to file",0);
			return NULL;
		}
#endif
	}
    else if ( config.MASTER_TYPE == PROB_MILP ) {
        if ( changeMILPwithL1(cell->master->lp, prob[0]->sp, prob[0]->num->cols) ) {
            errMsg("algorithm", "algoIntSD", "failed to change the proximal term", 0);
            return NULL;
        }
        
        if ( constructMILP(prob[0], cell, cell->incumbX, cell->quadScalar) ) {
            errMsg("setup", "newCell", "failed to change the right-hand side after incumbent change", 0);
            return NULL;
        }
        cell->incumbChg = FALSE;
#if defined(SETUP_CHECK)
        if ( writeProblem(cell->master->lp, "newMILPMaster.lp") ) {
            errMsg("write problem", "new_master", "failed to write master problem to file",0);
            return NULL;
        }
#endif
    }

	return cell;
}//END newCell()

void freeConfig() {

	if (config.RUN_SEED) mem_free(config.RUN_SEED);
	if (config.EVAL_SEED) mem_free(config.EVAL_SEED);

}//END freeConfig()

int cleanCellType(cellType *cell, probType *prob, vector xk) {
	int cnt;

	/* constants and arrays */
	cell->k = 0;
	cell->LPcnt = 0;
	cell->optFlag 		 = FALSE;
	cell->spFeasFlag 	 = TRUE;
	if ( config.DUAL_STABILITY )
		cell->dualStableFlag 	= FALSE;

	copyVector(xk, cell->candidX, prob->num->cols, TRUE);
	cell->candidEst	= prob->lb + vXvSparse(cell->candidX, prob->dBar);

	if (config.MASTER_TYPE == PROB_QP) {
		copyVector(xk, cell->incumbX, prob->num->cols, TRUE);
		cell->incumbEst = cell->candidEst;
		cell->quadScalar= config.MIN_QUAD_SCALAR;
		cell->iCutIdx   = 0;
		cell->iCutUpdt  = 0;
		cell->incumbChg = TRUE;
	}
    cell->RepeatedTime = 0;
	cell->gamma 	= 0.0;
	cell->normDk_1 	= 0.0;
	cell->normDk 	= 0.0;

	/* oneProblem structures and solver elements */
	for ( cnt = prob->num->rows+cell->cuts->cnt+cell->fcuts->cnt-1; cnt >= prob->num->rows; cnt-- )
		if (  removeRow(cell->master->lp, cnt, cnt) ) {
			errMsg("solver", "cleanCellType", "failed to remove a row from master problem", 0);
			return 1;
		}
	cell->master->mar = prob->num->rows;
	
    if ( config.MASTER_TYPE == PROB_QP ) {
        if( changeQPproximal(cell->master->lp, prob->num->cols, cell->quadScalar)) {
            errMsg("algorithm", "cleanCellType", "failed to change the proximal term", 0);
            return 1;
        }
    }
    else if (config.MASTER_TYPE == PROB_MILP){
        if( changeMILPproximal(cell->master->lp, prob->sp->objx, prob->num->cols, cell->quadScalar)) {
            errMsg("algorithm", "cleanCellType", "failed to change the proximal term", 0);
            return 1;
        }
    }

	/* cuts */
	if (cell->cuts) freeCutsType(cell->cuts, TRUE);
	if (cell->fcuts) freeCutsType(cell->fcuts, TRUE);
	if (cell->fcutsPool) freeCutsType(cell->fcutsPool, TRUE);
	cell->feasCnt 		= 0;
	cell->infeasIncumb 	= FALSE;
	cell->fUpdt[0] = cell->fUpdt[1] = 0;
    
    /*Jiajun, added on July 5th, 2019, cuts extra*/
    if(cell->argmax_best_candid)mem_free(cell->argmax_best_candid);
    if(cell->argmax_best_incumb)mem_free(cell->argmax_best_incumb);
    if(cell->pi_best_candid)mem_free(cell->pi_best_candid);
    if(cell->pi_best_incumb)mem_free(cell->pi_best_incumb);

	/* stochastic components */
	if (cell->basis) freeBasisType(cell->basis, TRUE);
    //Jiajun Check
    if(prob->num->rvRowCnt > 0){
        if (cell->delta) freeDeltaType(cell->delta, cell->lambda->cnt, cell->omega->cnt, TRUE);
        if (cell->lambda) freeLambdaType(cell->lambda, TRUE);
    }
	if (cell->sigma) freeSigmaType(cell->sigma, TRUE);
	if (cell->omega) freeOmegaType(cell->omega, TRUE);

	/* reset all the clocks */
	cell->time.repTime = cell->time.iterTime = cell->time.masterIter = cell->time.subprobIter = cell->time.optTestIter = cell->time.argmaxIter = 0.0;
	cell->time.iterAccumTime = cell->time.masterAccumTime = cell->time.subprobAccumTime = cell->time.optTestAccumTime = cell->time.argmaxAccumTime = 0.0;

	if ( config.MASTER_TYPE == PROB_QP ) {
		if ( constructQP(prob, cell, cell->incumbX, cell->quadScalar) ) {
			errMsg("setup", "newCell", "failed to change the right-hand side after incumbent change", 0);
			return 1;
		}

		cell->incumbChg = FALSE;
#if defined(SETUP_CHECK)
		if ( writeProblem(cell->master->lp, "cleanedQPMaster.lp") ) {
			errMsg("write problem", "new_master", "failed to write master problem to file",0);
			return 1;
		}
#endif
	}
    else if ( config.MASTER_TYPE == PROB_MILP ) {
        if ( constructMILP(prob, cell, cell->incumbX, cell->quadScalar) ) {
            errMsg("setup", "newCell", "failed to change the right-hand side after incumbent change", 0);
            return 1;
        }
        
        cell->incumbChg = FALSE;
#if defined(SETUP_CHECK)
        if ( writeProblem(cell->master->lp, "cleanedMILPMaster.lp") ) {
            errMsg("write problem", "new_master", "failed to write master problem to file",0);
            return 1;
        }
#endif
    }

	return 0;
}//END cleanCellType()


void freeCellType(cellType *cell) {

	if ( cell ) {
		if (cell->master) freeOneProblem(cell->master);
		if (cell->candidX) mem_free(cell->candidX);
		if (cell->incumbX) mem_free(cell->incumbX);
		if (cell->piM) mem_free(cell->piM);
		if (cell->djM) mem_free(cell->djM);
		if (cell->cuts) freeCutsType(cell->cuts, FALSE);
		if (cell->fcuts) freeCutsType(cell->fcuts, FALSE);
		if (cell->fcutsPool) freeCutsType(cell->fcutsPool, FALSE);
		if (cell->omega) freeOmegaType(cell->omega, FALSE);
        //if (cell->delta) freeDeltaType(cell->delta, cell->lambda->cnt, cell->omega->cnt, FALSE);
        //if (cell->lambda) freeLambdaType(cell->lambda, FALSE);
//#ifdef DEBUG
//        printf("sigma val in freeCellType():");
//        printVector(cell->sigma->vals[0].piC, 3, NULL);
//#endif
        if (cell->sigma) freeSigmaType(cell->sigma, FALSE);
		if (cell->basis) freeBasisType(cell->basis, FALSE);
		if (cell->pi_ratio) mem_free(cell->pi_ratio);
		mem_free(cell);
	}

}//END freeCellType()
