/*--------------------------------------------------------------------------
 * Ox DySco  Package
 * DySco-Pack pack is fully operational but the author does not provide any liability,
 * including liability for financial losses.
 * The package may be used freely for non-commercial use.
 *
 * The handbook is available at  web addrressXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 *
 * version: 		2.1 (July 2013)
 * copright: 		Philipp Andres, University of Cambrigde, k.philipp.andres@gmail.com
 * 
 *
 *--------------------------------------------------------------------------*/

#include <oxstd.h>
#include <arma.h>
#import <modelbase>
#include <oxprob.h>
#include <oxdraw.h>
#include <oxfloat.h>
#import <maximize>
			
// vModelstruct = [ARLags, MALags, Leverage, AR2comp,MA2comp]
// m_vStruct = [ARLags, MALags, vShapeParas, Levergae, AR2comp, MA2comp]
// vP = [const, AR, MA, Shape, Lev, ARSr, MASr]

// Lists
enum{Y_VAR,X_VAR, ID_VAR};
enum{EX, GA, WBL, LNORM, GGA, LLOG, BURR, F, GB2,DAG};
enum{BFGS, NR, SCORING, BHHH}
enum{SCORE, MEM}

 // m_vStruct = [ARLags[0], MALags[1], vShapeParas[2], Levergae[3], AR2comp[4], MA2comp[5]]
class DySco : Modelbase
{
// member globals
decl m_cDist; decl m_vStruct; decl m_iLev;
decl m_iModelClass;	decl m_iScore; decl m_iRoutine;
decl m_vParStart; decl m_cElapsedTime;
decl m_iResult; decl m_dLogLik;	decl m_asParNames;
decl m_iTwoComp; decl m_mLev; decl m_vP;  //decl m_cT;
decl m_cP; decl m_cMaxIter;	decl m_iAutoStartValues; decl m_iPrint;
decl m_cBins; decl m_vEps; decl m_cTAll;

decl m_vTruePar; decl m_vY;
decl m_cRep; decl m_cBurnIn;
decl m_vLambda; decl m_vU;decl m_vScore; decl m_vMu;
// member functions


AutoStartValues();
APGT(const cd, const ng, const np);

BurrLambdaFilter(const vLambda, const vY, const vP, const vStruct, const iLev);


Citation();
CheckPara(const vPar);
DAGLambdaFilter(vLambda,vU, const vY, const vP, const vStruct);
DAGLik(const vP, const adFunc, const avScore, const amHessian);
DAGAnaCovMatrix(mAnaHess, const vO, const cT, const cMaxpq);
DAGMEMLik(const vP, const adFunc, const avScore, const amHessian);
DAGLikVal(const vP, const cStart);
DySco();
DCSOutputPar();
DCSOutputParLaTex();
DoEstimation(vPar);
DoStaticEstimation();
DoSimDlg();
DoSettingsDlg();
DoEstimateDlg();
DoTestGraphicsDlg();
DoTestDlg();
DrawEps(const iPos);
DrawPitEpsIS(const iPos);
DrawEpsOoS(const iPos);
DrawPredict(const iPos);
DrawPredictOoS(const iPos);
DrawLannePit(const iPos);
DrawPitEps(const iPos, const vEps);
DrawPitEpsOoS(const iPos, const vEps);
DrawPitEpsOoS2(const iPos);
DrawPitScore(const iPos, const vScore);
DrawPitScoreOoS(const iPos, const vScore);
DrawEpsACF(const iPos);
DrawEpsACFOoS(const iPos);
DrawScoreACF(const iPos);
DrawScoreACFOoS(const iPos);
DrawEpsSqACF(const iPos);
DrawEpsSqACFOoS(const iPos);
DrawScoreSqACFOoS(const iPos);
DrawScoreSqACF(const iPos);
DrawQQPlotPits2(const iPos);
Draw(const iPos,const sFileName, const mCoefVal);



Beta(const x, const y);
BURRAnaCovMatrix(mAnaHess, const vO, const cT, const cMaxpq);
BURRLambdaFilter(vLambda,vU, const vY, const vP, const vStruct);
BURRLik(const vP, const adFunc, const avScore, const amHessian);
BURRLikVal(const vP, const cStart);
BURRScore(vScore, const vEps,const vLambda, const vP, const cMaxpq);
BURRMEMScores(vScore,const vY,const vLambda, const vP, const cMaxpq);
BURRMEMLik(const vP, const adFunc, const avScore, const amHessian);
BURRSimulate(const mCoefVal, const vElapsedTime, const cT, const cRep, const vPtrue, const sFilename, const sFileNameData, const iPos);



EXLambdaFilter(vLambda,vU, const vY, const vP, const vStruct);
EXAnaCovMatrix(const mAnaHess, const vO, const cT, const cMaxpq);
EXLik(const vP, const adFunc, const avScore, const amHessian) ;
EXLikVal(LikVal, const vP, const cStart) ;
EXMEMLik(const vP, const adFunc, const avScore, const amHessian) ;
EXScore(vScore, const vEps,const vLambda, const vP,const cMaxpq);
EXMEMScore(vScore, const vEps,const vLambda,const vY, const vP, const cMaxpq);


Estimate();
FEM(const forc, const obs, const bOoS);
FEM2();
FLambdaFilter(vLambda,vU, const vY, const vP, const vStruct);
FLik(const vP, const adFunc, const avScore, const amHessian);
FLikVal(const vP, const cStart);
FScore(vScore, const eps,const vLambda, const vP, const cMaxpq);
FMEMLik(const vP, const adFunc, const avScore, const amHessian);
FMEMScores(vScore,const vY,const vLambda, const vP, const cMaxpq);
FAnaCovMatrix(const mAnaHess, const vO, const cT, const cMaxpq);
FSimulate(const mCoefVal, const vElapsedTime, const cT, const cRep, const vPtrue, const sFilename, const sFileNameData, const iPos);

GALambdaFilter(vLambda,vU, const vY, const vP, const vStruct);
GAMEMScore(vScore, const vEps,const vLambda,const vY,const vP,const cMaxpq, const cT);
GAMEMLik(const vP, const adFunc, const avScore, const amHessian);
GAAnaCovMatrix(mAnaHess,const vO, const cT, const cMaxpq);
GAScore(vScore, const vEps,const vLambda,const vP,const cMaxpq, const cT);
GALik(const vP, const adFunc, const avScore, const amHessian);
GALikVal(const vP, const cStart);
GGALambdaFilter(vLambda,vU, const vY, const vP, const vStruct, const cT, const maxpq);
GGALik(const vP, const adFunc, const avScore, const amHessian);
GGALikVal(const vP, const cStart);
GGAAnaCovMatrix(mAnaHess, const vO, const cT, const cMaxpq);
GGAMEMLik(const vP, const adFunc, const avScore, const amHessian) ;
GASimulate(mCoefVal, const vElapsedTime, const cT, const cRep, const vPtrue, const sFilename, const sFileNameData, const iPos);
GB2Lik(const vP, const adFunc, const avScore, const amHessian);
GB2LambdaFilter(const vLambda, const vY, const vP, const cT);
GB2AnaCovMatrix(const mAnaHess, const vO);
GB2Scores(const vScore, const eps,const vLambda, const vP, const cT, const cMaxpq);
GB2LikVal(const vP, const cStart);
GB2MEMLik(const vP, const adFunc, const avScore, const amHessian) ;

GetPackageName();
GetPackageVersion();
GetModelClassString();
GetDistriString();
GetNumAlgoString();
GetStartParams();
GetLevAux(vLevaux, const vLev);
GetNumStdErr(vNumStdErr,const vP, const iDist);
GetMEMAnaStdErr(vAnaStdErr, const vP, const iDist, const cT, const cMaxpq);
GetStaticNumStdErr(vNumStdErr,const vP, const iDist);
GetAnaStdErr(vAnaStdErr, const vP, const iDist, const cT, const cMaxpq);
GetStaticLL(const vPara, const adFunc, const avScore, const amHessian);
GetNumberParams(const vStruct);
 GetNumStdErrMEM(vNumStdErr, const vP, const iDist);
GetShapePara(const vP, const vStruct, const cDist);
GetEps();
GetPredict();
GetPit(const vEps, const vShapePara);
GetScores();
GetPitScore(const vScore, const vShapePara);
GetPitConfUpper(const avF, const cUpper);
GetPitConfLower(const avF, const cUpper);
 
InitData();
InitPar();
InitGlobals();

LannePit(const vY, const mu1, const mu2, const gamma1, const gamma2, const ProbPi);
LambdaFilter(vLambda, vU, const vY, const vP, const vStruct, const cDist);
LambdaFilter2comp(vLambda, vU, const vY, const vP, const vStruct, const cDist);
LBQTest(const ma, const vLag);
LBQTest2(const vLag);
LBQTestScores(const vLag);
LLOGLik(const vP, const adFunc, const avScore, const amHessian);
LLOGLikVal(const vP, const cStart);
LNORMLikVal(const vP, const cStart);
LLOGLambdaFilter(vLambda,vU, const vY, const vP, const vStruct);
LLOGScore(vScore, const vEps,const vLambda, const vP, const cMaxpq);
LLOGAnaCovMatrix(const mAnaHess, const vO, const cT, const cMaxpq);
LLOGMEMLik(const vP, const adFunc, const avScore, const amHessian);
LLOGMEMScore(const vScore, const vY,const vLambda, const vP, const cMaxpq);
LLOGSimulate(const mCoefVal, const vElapsedTime, const cT, const cRep, const vPtrue, const sFilename, const sFileNameData, const iPos);


LNORMLik(const vP, const adFunc, const avScore, const amHessian);
LogNormalLikStatic(const vY,  const vPara, const cMaxpq, const cT);
LNORMLambdaFilter(vLambda,vU, const vY, const vP, const vStruct, const cT, const maxpq);
LNORMScore(vScore, const vEps,const vLambda, const vP, const cMaxpq, const cT);
LNORMAnaCovMatrix(mAnaHess,const vO, const cT,const cMaxpq, const vEps);
LNORMMEMLik(const vP, const adFunc, const avScore, const amHessian);
LNORMMEMScore(vScore, const vEps,const vLambda,const vY , const vP,const cMaxpq, const cT);
MEMLambdaFilter(const vLambda, const vY, const vP, const vStruct);
MEMLNORMAnaCovMatrix(mAnaHess,const vO, const cT,const cMaxpq, const vEps);
MEMLambdaFilter2comp(vLambda, const vY, const vP, const vStruct);
MEMWBLAnaCovMatrix(const mAnaHess,const vO);
MEMGGAnaCovMatrix(const mAnaHess, const vO, const cT, const cMaxpq);
MEMGGScore(const vScore, const eps,const vLambda, const vP, const cT, const cMaxpq);
MEMBurrAnaCovMatrix(const mAnaHess, const vO, const cT, const cMaxpq);
MEMLLOGAnaCovMatrix(const mAnaHess, const vO, const cT, const cMaxpq);
MEMFAnaCovMatrix(const mAnaHess, const vO, const cMaxpq);


MEMOutputPar();

LBQStatPitEps(const vEps, const vLag);

NormTest1(const vX) ;


OutputHeader(const sTitle);
Output();
ParNames();
PrintLikStatsLatex();
probgengamma(const x, const nu , const beta, const alpha);
PrintLikStats();
probloglogistic(const x, const scale, const beta);
probburr(const x, const scale, const beta, const alpha);
probdag(const x, const scale, const nu, const xi);
probgb2(const x,const scale, const nu, const xi,const varsig);
PrintSampleStats();

QQPlot(const vX, const cDist, const vShapePara, const iPos);
QQPlotDrawEps(const iPos, const vEps);
QQPlotDrawEpsOoS(const iPos, const vEps);
QQPlotDrawScore(const iPos, const vU);
QQPlotDrawScoreOoS(const iPos, const vU);
QQPlotScores(const iPos, const cDist, const vShapePara, const vScores);
QQPlotPits2(const iPos);
QQPlotPits(const iPos,const vScore);


ReceiveData();														 
ReceiveModel();
ReceiveMenuChoice(const sMenu);
Report(const sMethod, const cElapsedTime, const mCoefVal, const vP, const cRep, const cT,const cRepFailed);


SetStruct(vStruct, const vModelStruct, const cDist);
SplitPara(amPar, const vP, const vStruct);
SetModel(const vModelStruct);
SetDistribution(const cDist);
SetModelClass(const iModelClass);
SetScore(const iScore);
SetRoutine(const iRoutine);
SetStartParams(const vP);
SetTwoComp();
SetcT();
SetMaxIter(const cMaxIter);
SendSpecials();
SendVarStatus();
SendMenu(const sMenu);
SetAutoStartVals(const iAutoStartVals);
SetNumberBins(const cBins);
SetTruePar(const vBeta);
SetInSamplePeriod(const cStart, const cEnd);

UserStartValues(const vInitPar);


WBLMEMLik(const vP, const adFunc, const avScore, const amHessian);
WBLAnaCovMatrix(mAnaHess,const vO, const cT,const cMaxpq);
WBLLambdaFilter(vLambda,vU, const vY, const vP, const vStruct, const cT, const maxpq);
WBLMEMScore(vScore, const vEps,const vLambda,const vY , const vP,const cMaxpq, const cT);
WBLLik(const vP, const adFunc, const avScore, const amHessian);
WBLLikVal(const vP, const cStart);
WBLScore(vScore, const vEps,const vLambda, const vP,const cMaxpq, const cT);
WBLSimulate(const mCoefVal, const vElapsedTime, const cT, const cRep, const vPtrue, const sFilename, const sFileNameData, const iPos);

}

DySco::DySco()
{
	Modelbase();
	println("---- ", GetPackageName(), " ", GetPackageVersion(),
	" session started at ", time(), " on ", date(), " ----");
	InitGlobals();	
}

DySco::InitGlobals()
{
 //m_iPrint = 1;
 m_cBins = 100;
 m_cMaxIter = 1000;
 m_iScore = 1;
}


DySco::GetPackageName()
{
 return " DySco ";
}
DySco::GetPackageVersion()
{
 return "2.0";
}

DySco::Citation()
{
println("");
println(" -----  Please cite DySco-Pack whenever used ------- ");
println("");
	println("@TECHREPORT{Dynamic Score Models for Positive Variables (DySco-Pos), OX Pack Manual v2.0,\n"
            " AUTHOR = {Philipp K. Andres},\n",
            " TITLE = {Dynamic Score Models for Positive Variables (DySco-Pos): An Object-Oriented Package for Ox},\n",
            " INSTITUTION = {Faculty of Economics, Cambridge University},\n",
            " YEAR = {2012},\n",
            " type = {Manual},\n",
            " number = {v210},\n",
            " address = {Cambridge, UK},\n",
            " month = {September},\n",
            "}");
println("---------------------------------------------------  ");
}

// Initalize Data
DySco::InitData() 
{
//	Database::Info();
//	m_cTAll  = sizer(m_mData);

	m_vY = GetGroup(Y_VAR);
	m_mLev = GetGroup(X_VAR);

//	//m_mLev = m_mData[][1];//GetGroup(X_VAR);
//	if (rows(m_mData)!=rows(m_mLev)) {oxwarning("Leverage Series not of same length as data series. Please correct."); exit(0);}


//	println(m_vY);
	//println(m_vY);
	m_cTAll  = sizer(m_mData);
//		println(m_mData[0:10][]);
		
	//println("m_cTAll",m_cTAll);

	if (m_mData == <>) {oxwarning("Please select the Y variable"); exit(0);}

	decl cY = columns(m_vY);
	if (cY != 1){ oxwarning("Only univariate models allowed"); exit(0); }

	decl cL = columns(m_mLev);
//	println(m_mLev[1:10][]);
//	println("cL :",cL);
	if (cL > 1){ oxwarning("Only one leverage series allowed"); exit(0); }

	
	if (ismissing(log(m_vY))){oxwarning("Dataset contains zero or negative observations. Please remove them");exit(0);}

    return TRUE;
}

DySco::SetInSamplePeriod(const cStart, const cEnd)
{  SetSelSampleByIndex(cStart, cEnd);
 GetSelSample();
 m_cT = m_iT2sel - m_iT1sel + 1;
}

//Initialize Parameters
DySco::InitPar()
{  Modelbase::InitPar();
   SetParCount(m_cPar);
   
   if (m_vStruct[3] == 1 && m_mLev == <>){
		oxwarning("Leverage selected but no return series choosen"); exit(0);}
	//println("aux", m_vStruct[0]~m_vStruct[1]~m_iScore);

//	if(m_iScore==0);// || (m_vStruct[0]!=1 || m_vStruct[1]!=1 |\ m_iRoutine==SCORING) || (m_vStruct[0]!=1 || m_vStruct[1]!=1 || m_iRoutine==BHHH) )
//	{oxwarning("Method of Scoring and analytical derivatives feasible for DySco(1,1) and log-ACD(1,1) model only"); exit(0);}

   return TRUE;
}

////////////////////////////////////////////////////////////////////////////////
////////////////// Helper Functions ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
DySco::SetAutoStartVals(const iAutoStartVals)
{ m_iAutoStartValues = iAutoStartVals;}

DySco::SetTwoComp() 
{if (m_vStruct[4]!=0 && m_vStruct[5]!=0) m_iTwoComp = 1;
 else if(m_vStruct[4]==0 && m_vStruct[5]==0) m_iTwoComp = 0;
 else if (m_vStruct[4]!=0 && m_vStruct[5]==0) {oxwarning("in 2 comp. moodel, cant set only short run AR elements. Short run MA terms required"); exit(0);}
}
DySco::GetLevAux(vLevaux, const vLev)
{// Purpose: from vectro vLev, generate series of 0 and 1's depending on size of vX
  vLevaux =  vLev .< 0 .? 1 .: 0;
  return vLevaux;
}

DySco::CheckPara(const vPar)
{decl Y;
	Y = (sumc(vPar[1:m_vStruct[0]]).> 1);
	if (Y > 0){
		println("Sum of choosen AR parameters exceeds one. You have entered ", vPar[1:m_vStruct[0]], "\n",
		"Optimisation terminated prematurely. Try using different starting values");
		exit(0);
		}	
	if (m_iTwoComp == 2){
	decl Y2 = (sumc(vPar[m_vStruct[0]+m_vStruct[1]+m_vStruct[2]+m_vStruct[3]:m_vStruct[0]+m_vStruct[1]+m_vStruct[2]+m_vStruct[3]+m_vStruct[4]]).>1);
		if (Y2 > 0){
			println("Sum of choosen short run AR parameters exceeds one. You have entered ", vPar[m_vStruct[0]+m_vStruct[1]+m_vStruct[2]+m_vStruct[3]:m_vStruct[0]+m_vStruct[1]+m_vStruct[2]+m_vStruct[3]+m_vStruct[4]], "\n",
			"Optimisation terminated prematurely. Try using different starting values");
			exit(0);
		}
	}
// long run AR parameters exceed short run AR parameters
	decl Z = 0;
	if (m_iTwoComp == 1){
	decl cAR = vPar[1:m_vStruct[0]];
	decl cAR2comp = vPar[1+m_vStruct[0]+m_vStruct[1]+m_vStruct[2]+m_vStruct[3]:m_vStruct[0]+m_vStruct[1]+m_vStruct[2]+m_vStruct[3]+m_vStruct[4]];
//	println(m_vStruct);
//	println("ARs", cAR);
//	println(cAR2comp);
	decl minp = min(m_vStruct[0], m_vStruct[4]); // smallest lag genth in long and short run MA parameters			
		for (decl i=0; i<minp;i++){
		 decl 	aux = (cAR2comp[i]>cAR[i]); // if returns 0
		 		Z += aux;
		}
	if (Z > 0){
		oxwarning(println("Choosen long run AR parameters dont all exceed short run AR parameters \n",
		"Optimisation terminated prematurely. Found parameter vector is:, ", vPar',
		"\n Try using different starting values"));
	//	println("Parameter Vecor: ", m_vPar');
	//	exit(0);
		}
	}
// long run MA parameters smaller than short run MA parameters
	decl ZZ = 0;
	if (m_iTwoComp == 1){
	decl cMA = vPar[1+m_vStruct[0]:m_vStruct[0]+m_vStruct[1]];
	decl cMA2comp = vPar[1+m_vStruct[0]+m_vStruct[1]+m_vStruct[2]+m_vStruct[3]+m_vStruct[4]:m_vStruct[0]+m_vStruct[1]+m_vStruct[2]+m_vStruct[3]+m_vStruct[4]+m_vStruct[5]];

	decl minq = min(m_vStruct[1], m_vStruct[5]); // smallest lag genth in long and short run MA parameters			
		for (decl i=0; i<minq;i++){
		 decl 	aux = (cMA2comp[i]<cMA[i]); // if returns 0
				ZZ += aux;
		}
	if (ZZ > 0){
	 	oxwarning(println("Choosen long run MA parameters ", m_vP', " not all smaller than short run MA parameters \n",
		"Optimisation terminated prematurely. Found parameter vector is:, ", m_vP',
		"\n Try using different starting values"));
		exit(0);
		}
	}

// Static Parametes exceed 0
	decl YY=0;
	for (decl j=1; j<=m_vStruct[2]; j++){
	 	if (vPar[m_vStruct[0]+m_vStruct[1]+j]<=0) YY+=1;
		}
	if(YY>0){ oxwarning(println("Found parameter estimates for the shape parameters are negative. Please correct by e.g. trying different starting values"));}
	return 1; 
}


DySco::SetDistribution(const cDist)
{if (cDist != EX && cDist != GA &&  cDist != WBL && cDist != LNORM && cDist != LLOG && cDist != BURR
	&& cDist != F && cDist !=GGA && cDist != GB2 && cDist != DAG){println("Inserted Distribution not valid, choose amogst EX, GA, IGA, WBL, LNORM, LLOG, BURR, F, GGA, GB2");
	exit(0);}
m_cDist = cDist;}

DySco::GetDistriString()
{if (m_cDist==EX) return "Exponential";
 else if (m_cDist==GA) return "Gamma";
 else if (m_cDist==WBL) return "Weibull";
 else if (m_cDist==LNORM) return "Log-Normal";
 else if (m_cDist==GGA) return "Gen Gamma";
 else if (m_cDist==LLOG) return "Log-Logistic";
 else if (m_cDist==BURR) return "Burr";
 else if (m_cDist==F) return "F";
 else if (m_cDist==DAG) return "Dagum";
 else if (m_cDist==GB2) return "Generalized Beta Two";}

DySco::SetModelClass(const iModelClass)
{ m_iModelClass = iModelClass;}

DySco::GetModelClassString()
{ if (m_iModelClass==SCORE) return "DySco";
  else if (m_iModelClass==MEM) return "log-ACD";}

DySco::SetScore(const iScore)
{if (iScore!=0 && iScore!=1) {oxwarning(println("SetScore(const iScore) requires 0 or 1 input")); exit(0);}
m_iScore = iScore;}
DySco::SetRoutine(const iRoutine)
{ if (iRoutine!=NR && iRoutine!=BFGS && iRoutine!=BHHH  && iRoutine!=SCORING) {oxwarning("SetRoutine(const iRoutine) requires either NR, BGFS or BHHH as input"); exit(0);}
m_iRoutine = iRoutine;}
  
DySco::GetNumAlgoString()
{ if (m_iScore==0 && m_iRoutine==BFGS) return "BFGS, analytical derivatives";
  else if (m_iScore==1 && m_iRoutine==BFGS) return "BFGS, numerical derivatives";
  else if (m_iRoutine==NR) return "Newton Raphson, approx derivatives & Hessian";
  else if (m_iRoutine==SCORING) return "Method of Scoring";
  else if (m_iRoutine==BHHH) return "BHHH";
  else {oxwarning("An invalid combination of m_iScore and m_iRoutine was choosen"); exit(0);}
}

//DySco::SetcT()
//{ m_cT = m_iT2sel - m_iT1sel +1;
//}
//
DySco::SetStruct(vStruct, const vModelStruct, const cDist)
{// vModelstruct = [ARLags, MALags, Leverage, AR2comp,MA2comp]
// m_vStruct = [ARLags, MALags, vShapeParas, Levergae, AR2comp, MA2comp]
 // rule out negative AR/MA lags
 if (minc(vModelStruct).<0){ oxwarning("Specified AR or MA lags cannot be negative. Leverage must take binary values.  Please correct");
 						exit(0);}
 vStruct[0:1] = vModelStruct[0:1];						
 if (cDist == EX || cDist == LNORM) vStruct[2][0] = 0;
 if (cDist == GA || cDist == WBL || cDist == LLOG) vStruct[2][0] = 1;
 if (cDist == GGA || cDist == BURR|| cDist==F|| cDist==DAG) vStruct[2][0] = 2;
 if (cDist == GB2) vStruct[2][0] = 3;
 vStruct[3:] = vModelStruct[2:];
 return vStruct;
}

DySco::SetMaxIter(const cMaxIter)
{m_cMaxIter = cMaxIter;}

DySco::GetNumberParams(const vStruct)
{m_cP = sumc(vStruct)+1;
//println("m_cP",m_cP); println(vStruct); println(sumc(vStruct));
}

DySco::SetModel(const vModelStruct)
{ decl vStruct = new matrix[6][1];
 m_vStruct = SetStruct(vStruct, vModelStruct, m_cDist);
 GetNumberParams(m_vStruct);
 SetTwoComp();
}

DySco::SetStartParams(const vP)
{ m_vParStart = vP;}

DySco::GetStartParams()
{return m_vParStart;}

DySco::GetShapePara(const vP, const vStruct, const cDist)
{ decl aux = sumc(vStruct[0:1]);
decl vShapePara;
//if (cDist==EX) oxwarning("Exponential Distribution: No shape parameter");
if (cDist==GA || cDist==WBL || cDist==LLOG || cDist==LNORM) vShapePara = vP[aux+1];
if (cDist==BURR || cDist==F || cDist==GGA|| cDist==DAG) vShapePara = vP[aux+1:aux+2];
if (cDist==GB2) vShapePara = vP[aux+1:aux+3];
return vShapePara;
}

DySco::SplitPara(amPar, const vP, const vStruct) 
{// purpose: take vP and model structure to array of parameters
//vStruct = [#AR, #MA,#ShapeParas,#Leverage#,#AR2comp, #MA2comp]
// ampar=[ConstPar, ARPar, MAPar, ShapePara, LeveragePar, Ar2Par, Ma2Par]
//if (rows(vP)!=m_cP) {
//	oxwarning("Number of entered parameters incorrect. Please correct"); exit(0);
//	}
decl nxt = 0;
amPar[0][0] = {vP[nxt]}; // constant
nxt += 1;
amPar[0][1]= {vP[nxt:nxt+vStruct[0]-1]};
nxt += vStruct[0];
amPar[0][2]= {vP[nxt:nxt+vStruct[1]-1]};
nxt += vStruct[1];	
if (vStruct[2]==0) { amPar[0][3]={}; }
else {amPar[0][3]={vP[nxt:nxt+vStruct[2]-1]};}
nxt += vStruct[2];
if (vStruct[3]==1) {amPar[0][4]={vP[nxt]}; nxt+=+1;}
else amPar[0][4]={};
if(vStruct[4]!=0){ amPar[0][5]= {vP[nxt: nxt+vStruct[4] -1]};  nxt += vStruct[4];}
else if (vStruct[4]==0) {amPar[0][5]={};}
if(vStruct[5]!=0){ amPar[0][6]= {vP[nxt: nxt+vStruct[5] -1]};  nxt +=  vStruct[5];}
else if (vStruct[5]==0) {amPar[0][6]={};}
return 1; 
}

DySco::ParNames()
{decl aNames = {}; aNames = {"Constant"};decl nxt = 1;
  if (m_vStruct[0]!=0){
 	 for (decl i=1; i<=m_vStruct[0]; i++){
		aNames ~= {sprint("AR(",i,")")};
		nxt+=1;} }
  if (m_vStruct[1]!=0){
   for (decl i=1; i<=m_vStruct[1]; i++){
		aNames ~= {sprint("MA(",i,")")};
		nxt +=1;}}
 if (m_vStruct[2]!=0){
	for (decl i = 1; i<=m_vStruct[2]; i++){
		aNames ~= {sprint("Shape Para ",i)};
		nxt +=nxt+1;
	}
  }
  if (m_cDist == LNORM)  aNames ~=  {"sigma^(2)"};
  if (m_vStruct[3]==1){
	aNames ~= {sprint("Leverage")};
  }

  if (m_vStruct[4]!=0){
	for (decl i=1; i<=m_vStruct[4]; i++){
		aNames ~= {sprint("AR short run(",i,")")};
		nxt +=1;}}
	if (m_vStruct[5]!=0){
	for (decl i=1; i<=m_vStruct[5]; i++){
		aNames ~= {sprint("MA short run(",i,")")};
		nxt +=1;}}

return m_asParNames = aNames;
}

 
///////////////////////////////////////////////////////////////////////////////////
////////////////////// General Dynamics //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
DySco::MEMLambdaFilter(const vLambda, const vY, const vP, const vStruct)
{ decl cT=rows(vY);		decl p = vStruct[0]; decl q = vStruct[1];
decl maxpq = max(p, q);
decl vLevAux=new matrix[cT][1];
if (vStruct[3]==1) {vLevAux=GetLevAux(&vLevAux, m_mLev);}
// starting val for first round is log(meanc(vY))
vLambda[0][0:maxpq-1] = log(meanr(vY));
decl amPar =  new array[7];
SplitPara(&amPar,vP,vStruct);
for(decl i=maxpq; i<cT;i++){
	vLambda[0][i] = amPar[0][0]*(1-sumc(amPar[1][0]))+amPar[1][0]'*reversec(vLambda[0][i-p:i-1])+amPar[2][0]'*reversec(log(vY[i-q:i-1]./exp(vLambda[0][i-q:i-1])));
	if (vStruct[3]==1){vLambda[0][i] += amPar[4][0]*vLevAux[i-1][0];}
	}
decl cStart = meanc(vLambda[0]);
vLambda[0][0:maxpq-1] = cStart;
// redo- lambda recursion but with starting value set to average of first lambda recursion
for(decl i=maxpq; i<cT;i++){
	vLambda[0][i] = amPar[0][0]*(1-sumc(amPar[1][0]))+amPar[1][0]'*reversec(vLambda[0][i-p:i-1])+amPar[2][0]'*reversec(log(vY[i-q:i-1]./exp(vLambda[0][i-q:i-1])));
	if (vStruct[3]==1){vLambda[0][i] += amPar[4][0]*vLevAux[i-1][0];}
	}
}

DySco::MEMLambdaFilter2comp(vLambda, const vY, const vP, const vStruct)
{decl cT = rows(vY); decl p = vStruct[0]; decl q = vStruct[1]; decl maxpq = max(p, q);
decl pS = vStruct[4]; decl qS = vStruct[5]; decl maxpqSR = max(pS, qS);
decl vLevAux=new matrix[cT][1];
if (vStruct[3]==1) {vLevAux=GetLevAux(&vLevAux, m_mLev);}
vLambda[0][0:maxpq-1] = log(meanr(vY));
decl vLambdaSR= new matrix[1][cT];
decl vLambdaLR= new matrix[1][cT];
decl amPar =  new array[7];
SplitPara(&amPar,vP,vStruct);
for(decl i=max(maxpq, maxpqSR); i<cT;i++){
	vLambdaSR[0][i] = amPar[5][0]'*reversec(vLambdaSR[0][i-pS:i-1])' + amPar[6][0]'*reversec(log(vY[i-qS:i-1]./exp(vLambdaSR[0][i-qS:i-1])'));
	if (vStruct[3]==1){vLambdaSR[0][i] += amPar[4][0]*vLevAux[i-1];}
	vLambdaLR[0][i] = amPar[1][0]'*reversec(vLambdaLR[0][i-p:i-1])' +amPar[2][0]'*reversec(log(vY[i-q:i-1]./exp(vLambdaLR[0][i-q:i-1])'));
	vLambda[0][i] = amPar[0][0]*(1-sumc(amPar[1][0]))+ vLambdaLR[0][i] + vLambdaSR[0][i];
	}
}
DySco::LambdaFilter(vLambda, vU, const vY, const vP, const vStruct, const cDist)
{decl cT = rows(vY); decl p = vStruct[0]; decl q = vStruct[1]; decl maxpq = max(p, q);
decl vLevAux=new matrix[cT][1];
if (vStruct[3]==1) {vLevAux=GetLevAux(&vLevAux, m_mData[][1]);}
if (m_iModelClass == SCORE) vLambda[0][0:maxpq-1] = vP[0];
decl amPar =  new array[7];
SplitPara(&amPar,vP,vStruct);
for(decl i=maxpq; i<cT;i++){
	vLambda[0][i] = amPar[0][0]*(1-sumc(amPar[1][0]))+amPar[1][0]'*reversec(vLambda[0][i-p:i-1])+amPar[2][0]'*reversec(vU[0][i-q:i-1]);
	if (vStruct[3]==1)  {vLambda[0][i] += amPar[4][0]*vLevAux[i-1][0];}

	if (cDist==EX) 				vU[0][i] = -1+vY[i]*exp(-vLambda[0][i]);
	else if (cDist==GA)			vU[0][i] = -amPar[3][0]+vY[i]*exp(-vLambda[0][i]);
	else if (cDist==WBL)		vU[0][i] = -1+(vY[i]*exp(-vLambda[0][i]))^(amPar[3][0]);
	else if (cDist==GGA)		vU[0][i] = (vY[i]*exp(-vLambda[0][i])).^(amPar[3][0][0]) - amPar[3][0][1];
	else if (cDist==LNORM)		vU[0][i] = log(vY[i])-vLambda[0][i];
	else if (cDist==LLOG)		vU[0][i] = -1+2*vY[i]^(amPar[3][0])*exp(-amPar[3][0]*vLambda[0][i])/(1+vY[i]^(amPar[3][0])*exp(-amPar[3][0]*vLambda[0][i]));
	else if (cDist==BURR)		vU[0][i] = -1+(amPar[3][0][1]+1)*vY[i]^(amPar[3][0][0])*exp(-amPar[3][0][0]*vLambda[0][i])/(1+vY[i]^(amPar[3][0][0])*exp(-amPar[3][0][0]*vLambda[0][i]));
	else if (cDist==F)			vU[0][i] = -0.5*amPar[3][0][0]+ 0.5*(amPar[3][0][0]+amPar[3][0][1])*amPar[3][0][0]*vY[i-1]*exp(-vLambda[0][i-1])/(amPar[3][0][1]+amPar[3][0][0]*vY[i-1]*exp(-vLambda[0][i-1]));
	else if (cDist==DAG)		vU[0][i] = -amPar[3][0][1]+ (amPar[3][0][1]+1)*vY[i]^(amPar[3][0][0])*exp(-amPar[3][0][0]*vLambda[0][i])/(1+vY[i]^(amPar[3][0][0])*exp(-amPar[3][0][0]*vLambda[0][i]));
	else if (cDist==GB2)			{ decl eps = vY[i]/exp(vLambda[0][i]);
								vU[0][i] = -amPar[3][0][1]+(amPar[3][0][1]+amPar[3][0][2])*eps^(amPar[3][0][0])/(1+eps^(amPar[3][0][0]));}
	}

}

DySco::LambdaFilter2comp(vLambda, vU, const vY, const vP, const vStruct, const cDist)
{decl cT = rows(vY); decl p = vStruct[0]; decl q = vStruct[1]; decl maxpq = max(p, q);
decl pS = vStruct[4]; decl qS = vStruct[5]; decl maxpqSR = max(pS, qS);
decl vLevAux=new matrix[cT][1];
if (vStruct[3]==1) {vLevAux=GetLevAux(&vLevAux, m_mData[][1]);}
vLambda[0][0:maxpq-1] = vP[0];
decl amPar =  new array[7];
SplitPara(&amPar,vP,vStruct);
decl vLambdaSR = new matrix[1][cT];
decl vLambdaLR = new matrix[1][cT];
for(decl i=max(maxpq, maxpqSR); i<cT;i++){
	vLambdaSR[0][i] = amPar[5][0]'*reversec(vLambdaSR[0][i-pS:i-1]) +amPar[6][0]'*reversec(vU[0][i-qS:i-1]);
	if (vStruct[3]==1){vLambdaSR[0][i] += amPar[4][0]*vLevAux[i-1][0];}
	vLambdaLR[0][i] = amPar[1][0]'reversec(vLambdaLR[0][i-p:i-1])+amPar[2][0]'*reversec(vU[0][i-q:i-1]);
	vLambda[0][i] = amPar[0][0] + vLambdaSR[0][i] + vLambdaLR[0][i];
	if (cDist==EX) 				vU[0][i] = -1+vY[i]*exp(-vLambda[0][i]);
	else if (cDist==GA)			vU[0][i] = -amPar[3][0]+vY[i]*exp(-vLambda[0][i]);
	else if (cDist==WBL)		vU[0][i] = -1+(vY[i]*exp(-vLambda[0][i]))^(amPar[3][0]);
	else if (cDist==GGA)		vU[0][i] = (vY[i]*exp(-vLambda[0][i])).^(amPar[3][0][0]) - amPar[3][0][1];
	else if (cDist==LNORM)		vU[0][i] = log(vY[i])-vLambda[0][i];
	else if (cDist==LLOG)		vU[0][i] = -1+2*vY[i]^(amPar[3][0])*exp(-amPar[3][0]*vLambda[0][i])/(1+vY[i]^(amPar[3][0])*exp(-amPar[3][0]*vLambda[0][i]));
	else if (cDist==BURR)		vU[0][i] = -1+(amPar[3][0][1]+1)*vY[i]^(amPar[3][0][0])*exp(-amPar[3][0][0]*vLambda[0][i])/(1+vY[i]^(amPar[3][0][0])*exp(-amPar[3][0][0]*vLambda[0][i]));
	else if (cDist==F)			vU[0][i] = -0.5*amPar[3][0][0]+ 0.5*(amPar[3][0][0]+amPar[3][0][1])*amPar[3][0][0]*vY[i-1]*exp(-vLambda[0][i-1])/(amPar[3][0][1]+amPar[3][0][0]*vY[i-1]*exp(-vLambda[0][i-1]));
	else if (cDist==DAG)		vU[0][i] = -amPar[3][0][1]+ (amPar[3][0][1]+1)*vY[i]^(amPar[3][0][0])*exp(-amPar[3][0][0]*vLambda[0][i])/(1+vY[i]^(amPar[3][0][0])*exp(-amPar[3][0][0]*vLambda[0][i]));
	else if (cDist==GB2)			{ decl eps = vY[i]/exp(vLambda[0][i]);
								vU[0][i] = -amPar[3][0][1]+(amPar[3][0][1]+amPar[3][0][2])*eps^(amPar[3][0][0])/(1+eps^(amPar[3][0][0]));}

	}

}


////////////////////////////////////////////////////////////////////////////////
///////////////////  Ex Distribution ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
DySco::EXLambdaFilter(vLambda,vU, const vY, const vP, const vStruct)
{decl cT = rows(vY); decl p = vStruct[0]; decl q = vStruct[1]; decl maxpq = max(p, q);
decl vLevAux=new matrix[cT][1];
if (vStruct[3]==1) {vLevAux=GetLevAux(&vLevAux, m_mLev);}
decl u = 0;
vLambda[0][0:maxpq-1] = vP[0];//log(vY[0]);
decl amPar =  new array[7];
SplitPara(&amPar,vP,vStruct);
for(decl i=maxpq; i<cT;i++){
	vLambda[0][i] = amPar[0][0]*(1-sumc(amPar[1][0]))+amPar[1][0]'*reversec(vLambda[0][i-p:i-1])+amPar[2][0]'*reversec(vU[0][i-q:i-1]);
	if (vStruct[3]==1){vLambda[0][i] += amPar[4][0]*vLevAux[i-1][0];}
	vU[0][i] = -1+vY[i]*exp(-vLambda[0][i]);
	}
}

DySco::EXLik(const vP, const adFunc, const avScore, const amHessian) 
{ decl cT = m_iT2sel-m_iT1sel+1;
decl vLambda = 	new matrix[m_iT2sel-m_iT1sel+1][1];
decl vU = 		new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) 	EXLambdaFilter(&vLambda,&vU, m_vY, vP, m_vStruct);
else  				LambdaFilter2comp(&vLambda, &vU, m_vY, vP, m_vStruct, EX);
adFunc[0] = double((-sumc(vLambda[1:])- sumc(m_vY[1:]./exp(vLambda[1:])))/cT);
if (isarray(avScore)){  // if analytical score should be used in maximisation
	decl vScore = new matrix[3][1];
	decl eps =m_vY./exp(vLambda);
	EXScore(&vScore, eps,vLambda, vP,max(m_vStruct[1],m_vStruct[2]));
	//println("vScore", vScore');
	(avScore[0])[0] = vScore[0];
	(avScore[0])[1] = vScore[1];
	(avScore[0])[2] = vScore[2];
}
if (amHessian){
	decl mAnaHess = new matrix[3][3];
	EXAnaCovMatrix(&mAnaHess,vP, (m_iT1sel-m_iT2sel+1), 1);
    (amHessian[0])= mAnaHess;
}
return 1;
}

DySco::EXLikVal(LikVal, const vP, const cStart) 
{ decl cT =rows(m_mData[][0]);
decl vLambda = 	new matrix[cT][1];
decl vU = 		new matrix[cT][1];
if (m_iTwoComp==0) 	EXLambdaFilter(&vLambda,&vU, m_mData[][0], vP, m_vStruct);
else  				LambdaFilter2comp(&vLambda, &vU, m_mData[][0], vP, m_vStruct, EX);
decl vLambda2  = vLambda[cStart:];
LikVal = double(-sumc(vLambda2)- sumc(m_mData[cStart:][0]./exp(vLambda2)));
return 1;
}


DySco::EXScore(vScore, const vEps,const vLambda, const vP,const cMaxpq)
{ decl cT = rows(vEps);
if (cT!=rows(vLambda)) oxwarning(println("lenght vEps != vLambda"));
decl mDeriv = new matrix[cT][3];
mDeriv[0][0] = (1-vP[1]); ///(vP[1]-vP[2]*vP[3]);
mDeriv[0][1] = mDeriv[0][2] = 0;
mDeriv[1:][0] = constant(1-vP[1],cT-1,1);
mDeriv[1:][1] = -ones(cT-1,1)*vP[0]+vLambda[:cT-2];
mDeriv[1:][2] = vEps[:cT-2]-ones(cT-1,1);
decl x;
for (decl i=1; i<cT; i++){
x = vP[1]-vP[2]*vEps[i-1];
mDeriv[i][] += x*mDeriv[i-1][];}
mDeriv = (vEps[1:]-ones(cT-1,1)).*mDeriv[1:][];
vScore[0][0:2] = sumc(mDeriv)/(cT-cMaxpq);
delete mDeriv;
}
DySco::EXAnaCovMatrix(const mAnaHess, const vO, const cT, const cMaxpq)
{decl a,b,c;
	a = vO[1]-vO[2];
	b = vO[1]^(2)-2*vO[1]*vO[2]+vO[2]^(2)*(1+1);
	c = -vO[2];	
	decl fact = 1/(1-b); // factor in front of dynamic elements of hessian
	// 1) dynamic parameters
	mAnaHess[0][0][0] =  fact*(1-vO[1])^(2)*(1+a)/(1-a) ;
	mAnaHess[0][1][1] =  fact*vO[2]^(2)*(1+a*vO[1])/((1-vO[1]^(2))*(1-a*vO[1]));
	mAnaHess[0][2][2] =  fact;
	mAnaHess[0][0][1] =  fact*a*c*vO[2]*(1-vO[1])/((1-a)*(1-a*vO[1]));
	mAnaHess[0][1][0] =  mAnaHess[0][0][1];
	mAnaHess[0][0][2] =  fact*c*(1-vO[1])/(1-a);
	mAnaHess[0][2][0] =  mAnaHess[0][0][2];
	mAnaHess[0][1][2] =  fact*a*vO[2]/(1-a*vO[1]);
	mAnaHess[0][2][1] =  mAnaHess[0][1][2];	
	mAnaHess[0] = -mAnaHess[0];
	return mAnaHess;
}
////////////////////// MEM stuff //////////////////////////////////////////////
DySco::EXMEMLik(const vP, const adFunc, const avScore, const amHessian)
{ decl cT = m_iT2sel-m_iT1sel+1;
decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) MEMLambdaFilter(&vLambda, m_vY, vP, m_vStruct);
else 	MEMLambdaFilter2comp(&vLambda, m_vY, vP, m_vStruct);
adFunc[0] = double((-sumc(vLambda[1:])- sumc(m_vY[1:]./exp(vLambda[1:])))/cT);
if (isarray(avScore)){  // if analytical score should be used in maximisation
	decl vScore = new matrix[3][1];
	decl eps = m_vY./exp(vLambda);
	EXMEMScore(&vScore, eps,vLambda,m_vY, vP,max(m_vStruct[1],m_vStruct[2]));
	//println("vScore", vScore');
	(avScore[0])[0] = vScore[0];
	(avScore[0])[1] = vScore[1];
	(avScore[0])[2] = vScore[2];
}return 1;
}

DySco::EXMEMScore(vScore, const vEps,const vLambda,const vY, const vP, const cMaxpq)
{decl cT = rows(vEps);
decl mDeriv = new matrix[cT][3];
mDeriv[0][0] = 1;
mDeriv[0][1] = meanc(-vP[0]+vLambda)/(1-vP[1]);
mDeriv[0][2] = meanc(vY)/(1-vP[1]);
mDeriv[1:][0] = constant(1-vP[1],cT,1);
mDeriv[1:][1] = -vP[0] + vLambda;
mDeriv[1:][2] = vY[0:cT-2];
for (decl i=cMaxpq; i<cT; i++){
mDeriv[i][1:] += vP[1]*mDeriv[i-1][1:];
}
mDeriv = dropr(mDeriv,0);
decl vU = -1 + vY[1:].*exp(-vLambda[1:]);
mDeriv = vU.*mDeriv;
vScore[0][0:2] = sumc(mDeriv)/(cT-cMaxpq);
delete mDeriv; return vScore ;
}

////////////////////////////////////////////////////////////////////////////////
/////////////////// GA Distribution ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
DySco::GALambdaFilter(vLambda,vU, const vY, const vP, const vStruct)
{decl cT = rows(vY); decl p = vStruct[0]; decl q = vStruct[1]; decl maxpq = max(p, q);
decl vLevAux=new matrix[cT][1];
if (vStruct[3]==1) { vLevAux=GetLevAux(&vLevAux, m_mLev);}
decl u = 0;
vLambda[0][0:maxpq-1] = vP[0];//log(vY[0]);
decl amPar =  new array[7];
SplitPara(&amPar,vP,vStruct);
for(decl i=maxpq; i<cT;i++){
	vLambda[0][i] = amPar[0][0]*(1-sumc(amPar[1][0]))+amPar[1][0]'*reversec(vLambda[0][i-p:i-1])+amPar[2][0]'*reversec(vU[0][i-q:i-1]);
	if (vStruct[3]==1){vLambda[0][i] += amPar[4][0]*vLevAux[i-1][0];}
	vU[0][i] = -amPar[3][0]+vY[i]*exp(-vLambda[0][i]);
	}
}

DySco::GALik(const vP, const adFunc, const avScore, const amHessian) 
{ decl cT = m_iT2sel-m_iT1sel+1;  decl p = m_vStruct[0]; decl q = m_vStruct[1]; decl maxpq = max(p, q);
decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
decl vU = new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) GALambdaFilter(&vLambda,&vU, m_vY, vP, m_vStruct);
else  LambdaFilter2comp(&vLambda, &vU, m_vY, vP, m_vStruct, GA);
				 
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
adFunc[0] = double( -loggamma(amPar[3][0]) +(-amPar[3][0]*sumc(vLambda[maxpq:])+(amPar[3][0]-1)*sumc(log(m_vY[maxpq:]))-sumc(m_vY[maxpq:].*exp(-vLambda[maxpq:])) )/(cT) );

if (isarray(avScore)){  // if analytical score should be used in maximisation
	decl vScore = new matrix[4][1];
	decl eps = m_vY./exp(vLambda);
	GAScore(&vScore, eps,vLambda, vP,max(m_vStruct[1],m_vStruct[2]), cT);
//	println("vScore", vScore');
	(avScore[0])[0] = vScore[0];
	(avScore[0])[1] = vScore[1];
	(avScore[0])[2] = vScore[2];
	(avScore[0])[3] = vScore[3];
}
if (amHessian){
	decl mAnaHess = new matrix[4][4];
	GAAnaCovMatrix(&mAnaHess,vP, (m_iT1sel-m_iT2sel+1), maxpq);
    (amHessian[0])= mAnaHess;
}
return 1;
}

DySco::GALikVal(const vP, const cStart)
{ decl cT =rows(m_mData[][0]);	decl cTIn = cT-cStart;
decl vLambda = 		new matrix[cT][1];
decl vU = 			new matrix[cT][1];
if (m_iTwoComp==0) 	GALambdaFilter(&vLambda,&vU, m_mData[][0], vP, m_vStruct);
else  				LambdaFilter2comp(&vLambda, &vU, m_mData[][0], vP, m_vStruct, GA);
//DrawTMatrix(0,vLambda',{"Lambda"});
//ShowDrawWindow();
decl vLambda2  = vLambda[cStart:];
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
decl vLikContr = new matrix[cTIn][1];
for(decl i=1;i<cTIn;i++){
vLikContr[i] = -loggamma(amPar[3][0])-amPar[3][0]*vLambda[cStart+i] + (amPar[3][0]-1)*log(m_mData[cStart+i][0]) - m_mData[cStart+i][0]*exp(-vLambda[cStart+i]);
}
decl LikVal = sumc(vLikContr);
return LikVal;
}

DySco::GAScore(vScore, const vEps,const vLambda,const vP,const cMaxpq, const cT)
{decl mDeriv = new matrix[cT][3];
mDeriv[0][0] = (1-vP[1]); ///(vP[1]-vP[2]*vP[3]);
mDeriv[0][1] = mDeriv[0][2] = 0;
mDeriv[1:][0] = constant(1-vP[1],cT-1,1);
mDeriv[1:][1] = -ones(cT-1,1)*vP[0]+vLambda[:cT-2];
mDeriv[1:][2] = vEps[:cT-2]-vP[3]*ones(cT-1,1);
decl x;
for (decl i=1; i<cT; i++){
x = vP[1]-vP[2]*vEps[i-1];
mDeriv[i][] += x*mDeriv[i-1][];
}
mDeriv = (vEps[1:]-vP[3]*ones(cT-1,1)).*mDeriv[1:][];
vScore[0][0:2] = sumc(mDeriv)/(cT-cMaxpq);
delete mDeriv;
//BFGS//
//vScore[0][3] = (-(cT-cMaxpq)*polygamma(vP[3],0) + sumc(log(vEps[1:])))/(100*(cT-cMaxpq));//sumc(log(m_mData[1:]))-sumc(vLambda[1:]))/(cT-cMaxpq);
vScore[0][3] = (-(cT-cMaxpq)*polygamma(vP[3],0) + sumc(log(vEps[cMaxpq:])))/((cT-cMaxpq));//sumc(log(m_mData[1:]))-sumc(vLambda[1:]))/(cT-cMaxpq);
//println(vScore[0]');
}

DySco::GAAnaCovMatrix(mAnaHess,const vO, const cT, const cMaxpq)
{decl a,b,c;
a = vO[1]-vO[2]*vO[3];
b = vO[1]^(2)-2*vO[1]*vO[2]*vO[3]+vO[2]^(2)*(1+vO[3])*vO[3];
c = -vO[2]*vO[3];	
decl fact = vO[3]/(1-b); // factor in front of dynamic elements of hessian
	// 1) dynamic parameters
mAnaHess[0][0][0] =  fact*(1-vO[1])^(2)*(1+a)/(1-a) ;
mAnaHess[0][1][1] =  fact*vO[2]^(2)*vO[3]*(1+a*vO[1])/((1-vO[1]^(2))*(1-a*vO[1]));
mAnaHess[0][2][2] =  fact*vO[3];
mAnaHess[0][0][1] =  fact*a*c*vO[2]*(1-vO[1])/((1-a)*(1-a*vO[1]));
mAnaHess[0][1][0] =  mAnaHess[0][0][1];
mAnaHess[0][0][2] =  fact*c*(1-vO[1])/(1-a);
mAnaHess[0][2][0] =  mAnaHess[0][0][2];
mAnaHess[0][1][2] =  fact*a*vO[2]*vO[3]/(1-a*vO[1]);
mAnaHess[0][2][1] =  mAnaHess[0][1][2];
mAnaHess[0][3][3] = polygamma(vO[3], 1)/(cT-cMaxpq); // old version: did not divide by #observations, this version is correct but gives much smaller std deviation than numerical std errors
mAnaHess[0][3][0] =(1-vO[1])/(1-a);
mAnaHess[0][0][3] = mAnaHess[0][3][0];
mAnaHess[0][3][1] = mAnaHess[0][3][2]= mAnaHess[0][1][3]=mAnaHess[0][3][2]=0;
mAnaHess[0] = -mAnaHess[0];
}

///////// MEM STUFF ////////////////////////////////////////////////////////////
DySco::GAMEMLik(const vP, const adFunc, const avScore, const amHessian)
{ decl cT = m_iT2sel-m_iT1sel+1;  decl p = m_vStruct[0]; decl q = m_vStruct[1]; decl maxpq = max(p, q);
decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) MEMLambdaFilter(&vLambda, m_vY, vP, m_vStruct);
else 	MEMLambdaFilter2comp(&vLambda, m_vY, vP, m_vStruct);
//println("meam lambda", meanc(vLambda));
//println("vP", vP);
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
adFunc[0] = double( -loggamma(amPar[3][0])+(-amPar[3][0]*sumc(vLambda[maxpq:])+(amPar[3][0]-1)*sumc(log(m_vY[maxpq:]))-sumc(m_vY[maxpq:].*exp(-vLambda[maxpq:])) )/(cT) );
if (isarray(avScore)){  // if analytical score should be used in maximisation
	decl vScore = new matrix[4][1];
	decl eps = m_vY./exp(vLambda);
	GAMEMScore(&vScore, eps,vLambda,m_vY, vP,max(m_vStruct[1],m_vStruct[2]),cT);
	//println("vScore", vScore');
	(avScore[0])[0] = vScore[0];
	(avScore[0])[1] = vScore[1];
	(avScore[0])[2] = vScore[2];
	(avScore[0])[3] = vScore[3];
}return 1;
}

DySco::GAMEMScore(vScore, const vEps,const vLambda,const vY,const vP,const cMaxpq, const cT)
{decl mDeriv = new matrix[cT][3];
mDeriv[0][0] = 1;
mDeriv[0][1] = meanc(-vP[0]+vLambda)/(1-vP[1]);
mDeriv[0][2] = meanc(vY)/(1-vP[1]);
mDeriv[1:][0] = constant(1-vP[1],cT-1,1);
mDeriv[1:][1] = -vP[0] + vLambda;
mDeriv[1:][2] = vY[0:cT-2];
for (decl i=cMaxpq; i<cT; i++){
mDeriv[i][1:] += vP[1]*mDeriv[i-1][1:];
}
mDeriv = dropr(mDeriv,0);
decl vU = -vP[3] + vY[1:].*exp(-vLambda[1:]);
mDeriv = vU.*mDeriv;
vScore[0][0:2] = sumc(mDeriv)/(cT-cMaxpq);
delete mDeriv;
vScore[0][3] = (-(cT-cMaxpq)*polygamma(vP[3],0) +sumc(log(vY[cMaxpq:cT-1]))-sumc(vLambda[cMaxpq:cT-1]))/(cT-cMaxpq);
}
////////////////////////////////////////////////////////////////////////////////
///////////////////  WBL Distribution ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
DySco::WBLLambdaFilter(vLambda,vU, const vY, const vP, const vStruct, const cT, const maxpq)
{decl vLev= new matrix[cT][1]; decl p = m_vStruct[0]; decl q = m_vStruct[1];
decl vLevAux=new matrix[cT][1];
if (vStruct[3]==1) {vLevAux=GetLevAux(&vLevAux, m_mLev);}
decl u = 0;
vLambda[0][0:maxpq-1] = vP[0];//log(vY[0]);
decl amPar =  new array[7];
SplitPara(&amPar,vP,vStruct);
for(decl i=maxpq; i<cT;i++){
	vLambda[0][i] = amPar[0][0]*(1-sumc(amPar[1][0]))+amPar[1][0]'*reversec(vLambda[0][i-p:i-1])+amPar[2][0]'*reversec(vU[0][i-q:i-1]);
	if (vStruct[3]==1){vLambda[0][i] += amPar[4][0]*vLevAux[i-1][0];}
	vU[0][i] = -1+(vY[i]*exp(-vLambda[0][i]))^(amPar[3][0]);
	}
}

DySco::WBLLik(const vP, const adFunc, const avScore, const amHessian) 
{ decl cT = m_iT2sel-m_iT1sel+1;  decl p = m_vStruct[0]; decl q = m_vStruct[1]; decl maxpq = max(p, q);
decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
decl vU = new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) WBLLambdaFilter(&vLambda,&vU, m_vY, vP, m_vStruct, cT, maxpq);
else  LambdaFilter2comp(&vLambda, &vU, m_vY, vP, m_vStruct, WBL);
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
adFunc[0] = double(log(amPar[3][0])+(-amPar[3][0]*sumc(vLambda[1:]) + (amPar[3][0]-1)*sumc(log(m_vY[1:]))-sumc((m_vY[1:]./exp(vLambda[1:])).^(amPar[3][0])))/(cT-1) );

if (isarray(avScore)){  // if analytical score should be used in maximisation
	decl vScore = new matrix[4][1];
	decl eps = m_vY./exp(vLambda);
	WBLScore(&vScore, eps,vLambda, vP,max(m_vStruct[1],m_vStruct[2]), cT);
	//println("vScore", vScore');
	(avScore[0])[0] = vScore[0];
	(avScore[0])[1] = vScore[1];
	(avScore[0])[2] = vScore[2];
	(avScore[0])[3] = vScore[3];
}
if (amHessian){
	decl mAnaHess = new matrix[4][4];
	WBLAnaCovMatrix(&mAnaHess,vP, (m_iT1sel-m_iT2sel+1), maxpq);
    (amHessian[0])= mAnaHess;
}
return 1;
}

DySco::WBLLikVal(const vP, const cStart)
{ decl cT =rows(m_mData[][0]);	decl cTIn = cT-cStart;
decl vLambda = 		new matrix[cT][1];
decl vU = 			new matrix[cT][1];
decl p = m_vStruct[0]; decl q = m_vStruct[1]; decl maxpq = max(p, q);
if (m_iTwoComp==0) 	WBLLambdaFilter(&vLambda,&vU, m_mData[][0], vP, m_vStruct, cT, maxpq);
else  				LambdaFilter2comp(&vLambda, &vU, m_mData[][0], vP, m_vStruct, WBL);
//DrawTMatrix(0,vLambda',{"Lambda"});
//ShowDrawWindow();
decl vLambda2  = vLambda[cStart:];
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
decl vLikContr = new matrix[cTIn][1];
for(decl i=1;i<cTIn;i++){
vLikContr[i] = log(amPar[3][0]) - amPar[3][0]*sumc(vLambda[cStart+i]) + (amPar[3][0]-1)*log(m_mData[cStart+1][0])-(m_mData[cStart+i][0]/exp(vLambda[cStart+i])).^(amPar[3][0]) ;
}

//decl LikTest= double(log(amPar[3][0])+(-amPar[3][0]*sumc(vLambda[1:]) + (amPar[3][0]-1)*sumc(log(m_vY[1:]))-sumc((m_vY[1:]./exp(vLambda[1:])).^(amPar[3][0])))/(cT-1) );

decl LikVal = sumc(vLikContr);
return LikVal;
}

DySco::WBLScore(vScore, const vEps,const vLambda, const vP,const cMaxpq, const cT)
{ decl mDeriv = new matrix[cT][3];
mDeriv[0][0] = (1-vP[1]); ///(vP[1]-vP[2]*vP[3]);
mDeriv[0][1] = mDeriv[0][2] = 0;
mDeriv[1:][0] = constant(1-vP[1],cT-1,1);
mDeriv[1:][1] = -ones(cT-1,1)*vP[0]+vLambda[:cT-2];
mDeriv[1:][2] = vEps[:cT-2].^(vP[3])-ones(cT-1,1);
decl x;
for (decl i=1; i<cT; i++){
x = vP[1]-vP[2]*vP[3]*vEps[i-1]^(vP[3]);
mDeriv[i][] += x*mDeriv[i-1][];
}
mDeriv = vP[3]*(vEps[1:].^(vP[3])-ones(cT-1,1)).*mDeriv[1:][];
vScore[0][0:2] = sumc(mDeriv)/(cT-cMaxpq);
delete mDeriv;
vScore[0][3] = (1/vP[3] + cT^(-1)*sumc(log(vEps[1:])) - cT^(-1)*sumc(vEps[1:].^(vP[3]).*log(vEps[1:])))/100;
}

DySco::WBLAnaCovMatrix(mAnaHess,const vO, const cT,const cMaxpq)
{decl a,b,c;
a = vO[1]-vO[2]*vO[3];
b = vO[1]^(2)-2*vO[1]*vO[2]*vO[3]+2*vO[2]^(2)*vO[3]^(2);
c = -vO[2]*vO[3];	
if (b>=1){oxwarning("b>1"); println("b is:" ,b ~vO');continue;}
decl fact = vO[3]^(2)/(1-b); // factor in front of dynamic elements of hessian
// 1) dynamic parameters
mAnaHess[0][0][0] =  fact*(1-vO[1])^(2)*(1+a)/(1-a) ;
mAnaHess[0][1][1] =	 fact*vO[2]^(2)*1*(1+a*vO[1])/((1-vO[1]^(2))*(1-a*vO[1]));
mAnaHess[0][2][2] =  fact*1;
mAnaHess[0][0][1] =  fact*a*c*vO[2]*(1-vO[1])/((1-a)*(1-a*vO[1]));
mAnaHess[0][1][0] =  mAnaHess[0][0][1];
mAnaHess[0][0][2] =  fact*c*(1-vO[1])/(1-a);
mAnaHess[0][2][0] =  mAnaHess[0][0][2];
mAnaHess[0][1][2] =  fact*a*vO[2]*1/(1-a*vO[1]);
mAnaHess[0][2][1] =  mAnaHess[0][1][2];												

mAnaHess[0][3][3] = ((polygamma(1,0)+1)^(2)+M_PI^(2)/6)*vO[3]^(-2);///(cT-cMaxpq);
mAnaHess[0][3][0] =	 (-1+M_EULER)*(1-vO[1])/(1-a);  //-2*(1-vO[1])/(1-a);
mAnaHess[0][0][3] = mAnaHess[0][3][0];
mAnaHess[0][3][1] = mAnaHess[0][3][2]= mAnaHess[0][1][3]=mAnaHess[0][3][2]=0;

mAnaHess[0] = -mAnaHess[0];
}

/////// WBL MEM STUFF \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
DySco::WBLMEMLik(const vP, const adFunc, const avScore, const amHessian)
{ decl cT = m_iT2sel-m_iT1sel+1;  decl p = m_vStruct[0]; decl q = m_vStruct[1]; decl maxpq = max(p, q);
decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) MEMLambdaFilter(&vLambda, m_vY, vP, m_vStruct);
else 	MEMLambdaFilter2comp(&vLambda, m_vY, vP, m_vStruct);
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
adFunc[0] = double(log(amPar[3][0])+(-amPar[3][0]*sumc(vLambda[1:]) + (amPar[3][0]-1)*sumc(log(m_vY[1:]))-sumc((m_vY[1:]./exp(vLambda[1:])).^(amPar[3][0])))/(cT-1) );

if (isarray(avScore)){  // if analytical score should be used in maximisation
	decl vScore = new matrix[4][1];
	decl eps = m_vY./exp(vLambda);
	WBLMEMScore(&vScore, eps,vLambda,m_vY, vP,max(m_vStruct[1],m_vStruct[2]),cT);
	//println("vScore", vScore');
	(avScore[0])[0] = vScore[0];
	(avScore[0])[1] = vScore[1];
	(avScore[0])[2] = vScore[2];
	(avScore[0])[3] = vScore[3];
}return 1;
}

DySco::WBLMEMScore(vScore, const vEps,const vLambda,const vY , const vP,const cMaxpq, const cT)
{decl mDeriv = new matrix[cT][3];
mDeriv[0][0] = 1;
mDeriv[0][1] = meanc(-vP[0]+vLambda)/(1-vP[1]);
mDeriv[0][2] = meanc(vY)/(1-vP[1]);
mDeriv[1:][0] = constant(1-vP[1],cT-1,1);
mDeriv[1:][1] = -vP[0] + vLambda[:cT-2];
mDeriv[1:][2] = vY[:cT-2];
for (decl i=cMaxpq; i<cT; i++){
mDeriv[i][1:] += vP[1]*mDeriv[i-1][1:];
}
mDeriv = dropr(mDeriv,0);
decl aux = (vY./exp(vLambda)).^(vP[3]);
decl vU = -vP[3] + vP[3]*aux[1:].^(vP[3]);
mDeriv = vU.*mDeriv;
vScore[0][0:2] = sumc(mDeriv)/(cT-cMaxpq);
delete mDeriv;
vScore[0][3] = 1/vP[3] - cT^(-1)*sumc(vLambda[cMaxpq:]) + cT^(-1)*sumc(log(vY[cMaxpq:])) - cT^(-1)*sumc(aux[cMaxpq:].^(vP[3]).*log(aux[cMaxpq:]));
}

DySco::MEMWBLAnaCovMatrix(const mAnaHess,const vO)
{ decl a,b,c,d,e,f,g,sigma_e2,fact;
	f = gammafact(1+1/vO[3]);
	g = gammafact(1+2/vO[3]);
	sigma_e2 = g-f^(2);
	
	a = vO[1] - vO[2]*f;
	b = vO[1]^(2)-2*vO[1]*vO[2]*f + vO[2]^(2)*g;
	c = vO[1]*f -vO[2]*g;
	d = vO[2]*f/(1-vO[1]);
	e = vO[2]*sigma_e2/(1-vO[1]^(2))+d^(2);
	fact = vO[3]^(2)/(1-b);
//	println("a ",a,"b ",b,"c ",c,"d ",d,"e ",e,"f ",f);
//	println("sigma ",sigma_e2);
		if (b>=1){oxwarning("b>1"); println("b is:" ,b ~vO');continue;}
	//1) Dynamic parameters
	mAnaHess[0][0][0] = fact*(1-vO[1])^(2)*(1+a)/((1-a));
	mAnaHess[0][1][1] =	fact*(e+ 2*a*(1-vO[1]*a)^(-1)*(vO[2]^(2)*c*f*((1-vO[1])*(1-a))^(-1)+vO[1]*e+vO[2]*d*f));
	mAnaHess[0][2][2] =	fact*(g*(1-a)+2*c*f  )/(1-a);
	mAnaHess[0][0][1] =	fact*(1-vO[1])*(d+ vO[2]*a*f/((1-vO[1])*(1-a)) +a*(vO[1]*d+vO[2]*f+ vO[2]*c/(1-a))/(1-vO[1]*a)   );
	mAnaHess[0][1][0] =	mAnaHess[0][0][1]; 
	mAnaHess[0][0][2] =	fact*(1-vO[1])*(f+c)/(1-a);
	mAnaHess[0][2][0] =	mAnaHess[0][0][2];
	mAnaHess[0][1][2] =	fact*(d*f+ a*(vO[1]*d*f+vO[2]*g+vO[2]*c*f*(1-a)^(-1))*(1-vO[1]*a)^(-1) +vO[2]*c*f/((1-vO[1])*(1-a))  );
	mAnaHess[0][2][1] =	 mAnaHess[0][1][2];
	
	// 2) cross derivatives
    mAnaHess[0][0][3] = (-1+M_EULER)*(1-vO[1])/(1-a);
	mAnaHess[0][1][3] = (-1+M_EULER)*vO[2]*f/((1-a)*(1-vO[1])) ;
	mAnaHess[0][2][3] = (-1+M_EULER)*f/(1-a);
    mAnaHess[0][3][0] =  mAnaHess[0][0][3]; mAnaHess[0][3][1] = mAnaHess[0][1][3]; 	mAnaHess[0][3][2] = mAnaHess[0][2][3];
	// 3) static parameters	
	mAnaHess[0][3][3] = ((polygamma(1,0)+1)^(2)+M_PI^(2)/6)*vO[3]^(-2);///(cT-cMaxpq);



	mAnaHess[0] = -mAnaHess[0];
return mAnaHess[0];
}

////////////////////////////////////////////////////////////////////////////////
///////////////////  LNORM Distribution ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
DySco::LNORMLambdaFilter(vLambda,vU, const vY, const vP, const vStruct, const cT, const maxpq)
{decl vLev= new matrix[cT][1]; decl p = m_vStruct[0]; decl q = m_vStruct[1];
decl vLevAux=new matrix[cT][1];
if (vStruct[3]==1) {vLevAux=GetLevAux(&vLevAux, m_mLev);}
decl u = 0;
vLambda[0][0:maxpq-1] = vP[0];//log(vY[0]);
decl amPar =  new array[7];
SplitPara(&amPar,vP,vStruct);
for(decl i=maxpq; i<cT;i++){
	vLambda[0][i] = amPar[0][0]*(1-sumc(amPar[1][0]))+amPar[1][0]'*reversec(vLambda[0][i-p:i-1])+amPar[2][0]'*reversec(vU[0][i-q:i-1]);
	if (vStruct[3]==1){vLambda[0][i] += amPar[4][0]*vLevAux[i-1][0];}
	vU[0][i] = log(vY[i])-vLambda[0][i];
	}
}

DySco::LNORMLik(const vP, const adFunc, const avScore, const amHessian) 
{ decl cT = rows(m_vY);  decl p = m_vStruct[0]; decl q = m_vStruct[1]; decl maxpq = max(p, q);
decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
decl vU = new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) 		LNORMLambdaFilter(&vLambda,&vU, m_vY, vP, m_vStruct, cT, maxpq);
else  					LambdaFilter2comp(&vLambda, &vU, m_vY, vP, m_vStruct, LNORM);
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
decl sigma2_parameter = sumc((log(m_vY[1:])-vLambda[maxpq:]).^(2))/(cT-maxpq) ;
adFunc[0] = double( -0.5*log(M_2PI) - 0.5*log(sigma2_parameter)- 0.5*(cT-maxpq)^(-1)*sumc((log(m_vY[1:])-vLambda[maxpq:]).^(2))/(sigma2_parameter) - (cT-maxpq)^(-1)*sumc(log(m_vY[1:]))   );
decl eps = m_vY./exp(vLambda);

if (isarray(avScore)){  // if analytical score should be used in maximisation
	decl vScore = new matrix[3][1];
	LNORMScore(&vScore, eps,vLambda, vP,max(m_vStruct[1],m_vStruct[2]), cT);
	//println("vScore", vScore');
	(avScore[0])[0] = vScore[0];
	(avScore[0])[1] = vScore[1];
	(avScore[0])[2] = vScore[2];

}
if (amHessian){
	decl mAnaHess = new matrix[4][4];
	LNORMAnaCovMatrix(&mAnaHess,vP, (m_iT1sel-m_iT2sel+1), maxpq,eps );
    (amHessian[0])= mAnaHess;
}
return 1;
}

DySco::LNORMLikVal(const vP, const cStart)
{ decl cT =rows(m_mData[][0]);	decl cTIn = cT-cStart;
decl vLambda = 		new matrix[cT][1];
decl vU = 			new matrix[cT][1];
decl p = m_vStruct[0]; decl q = m_vStruct[1]; decl maxpq = max(p, q);
if (m_iTwoComp==0) 	LNORMLambdaFilter(&vLambda,&vU, m_mData[][0], vP, m_vStruct, cT, maxpq);
else  				LambdaFilter2comp(&vLambda, &vU, m_mData[][0], vP, m_vStruct, LNORM);
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
decl vLikContr = new matrix[cTIn][1];
decl sigma2_parameter = sumc((log(m_mData[:cStart][0])-vLambda[:cStart]).^(2))/(cStart) ;
 for(decl i=1;i<cTIn;i++){
vLikContr[i] =  -0.5*log(M_2PI) - 0.5*log(sigma2_parameter)- 0.5*((log(m_mData[cStart+i][0])-vLambda[cStart+i])^(2))/(sigma2_parameter) - log(m_mData[cStart+i][0]);
}
decl LikVal = sumc(vLikContr);
return LikVal;
}



DySco::LNORMScore(vScore, const vEps,const vLambda, const vP,const cMaxpq, const cT)
{ decl mDeriv = new matrix[cT][3];
mDeriv[0][0] = (1-vP[1]); ///(vP[1]-vP[2]*vP[3]);
mDeriv[0][1] = mDeriv[0][2] = 0;
mDeriv[1:][0] = constant(1-vP[1],cT-1,1);
mDeriv[1:][1] = -ones(cT-1,1)*vP[0]+vLambda[:cT-2];
mDeriv[1:][2] = log(vEps[:cT-2]) ;
decl x;
for (decl i=1; i<cT; i++){
x = vP[1]-vP[2];
mDeriv[i][] += x*mDeriv[i-1][];
}
decl sigma2_parameter = sumc(log(vEps[cMaxpq:]).^(2))/(cT-cMaxpq) ;
mDeriv = -sigma2_parameter^(-1)*log(vEps[1:]).*mDeriv[1:][];
vScore[0][0:2] = sumc(mDeriv)/(cT-cMaxpq);
}

DySco::LNORMAnaCovMatrix(mAnaHess,const vO, const cT,const cMaxpq, const vEps)
{decl a,b,c;
a = vO[1] -vO[2];
b = vO[1]^(2)- 2*vO[1]*vO[2]+vO[2]^(2);
c = 0;
if (b>=1){oxwarning("b>1"); println("b is:" ,b ~vO');continue;}
decl sigma_u_2 = sumc(log(vEps[cMaxpq:]).^(2))/(cT-cMaxpq);
decl fact = 1/((1-b)*sigma_u_2); // factor in front of dynamic elements of hessian
// 1) dynamic parameters
mAnaHess[0][0][0] =  fact*(1-vO[1])^(2)*(1+a)/(1-a) ;
mAnaHess[0][1][1] =	 fact*vO[2]^(2)*1*(1+a*vO[1])/((1-vO[1]^(2))*(1-a*vO[1]));
mAnaHess[0][2][2] =  fact*1;
mAnaHess[0][0][1] =  fact*a*c*vO[2]*(1-vO[1])/((1-a)*(1-a*vO[1]));
mAnaHess[0][1][0] =  mAnaHess[0][0][1];
mAnaHess[0][0][2] =  fact*c*(1-vO[1])/(1-a);
mAnaHess[0][2][0] =  mAnaHess[0][0][2];
mAnaHess[0][1][2] =  fact*a*vO[2]*1/(1-a*vO[1]);
mAnaHess[0][2][1] =  mAnaHess[0][1][2];												

//mAnaHess[0][3][3] = 2*vO[3]^(2);
//mAnaHess[0][3][0] =	 0; 
//mAnaHess[0][0][3] = mAnaHess[0][3][0];
//mAnaHess[0][3][1] = mAnaHess[0][3][2]= mAnaHess[0][1][3]=mAnaHess[0][3][2]=0;

mAnaHess[0] = -mAnaHess[0];
}

/////// LNORM MEM STUFF \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
DySco::LNORMMEMLik(const vP, const adFunc, const avScore, const amHessian)
{ decl  cT = rows(m_vY);  decl p = m_vStruct[0]; decl q = m_vStruct[1]; decl maxpq = max(p, q);
decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) MEMLambdaFilter(&vLambda, m_vY, vP, m_vStruct);
else 	MEMLambdaFilter2comp(&vLambda, m_vY, vP, m_vStruct);
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
decl sigma2_parameter = sumc((log(m_vY[maxpq:])-vLambda[maxpq:]).^(2))/(cT-maxpq) ;
adFunc[0] = double( -0.5*log(M_2PI) - 0.5*log(sigma2_parameter)- 0.5*(cT-maxpq)^(-1)*sumc((log(m_vY[maxpq:])-vLambda[maxpq:]).^(2))/(sigma2_parameter) - (cT-maxpq)^(-1)*sumc(log(m_vY[maxpq:]))   );

if (isarray(avScore)){  // if analytical score should be used in maximisation
	decl vScore = new matrix[3][1];
	decl eps = m_vY./exp(vLambda);
	LNORMMEMScore(&vScore, eps,vLambda,m_vY, vP,max(m_vStruct[1],m_vStruct[2]),cT);
	//println("vScore", vScore');
	(avScore[0])[0] = vScore[0];
	(avScore[0])[1] = vScore[1];
	(avScore[0])[2] = vScore[2];
}return 1;
}

DySco::LNORMMEMScore(vScore, const vEps,const vLambda,const vY , const vP,const cMaxpq, const cT)
{decl mDeriv = new matrix[cT][3];
mDeriv[0][0] = 1;
mDeriv[0][1] = meanc(-vP[0]+vLambda)/(1-vP[1]);
mDeriv[0][2] = meanc(vY)/(1-vP[1]);
mDeriv[1:][0] = constant(1-vP[1],cT-1,1);
mDeriv[1:][1] = -vP[0] + vLambda[:cT-2];
mDeriv[1:][2] = vY[:cT-2];
for (decl i=cMaxpq; i<cT; i++){
mDeriv[i][1:] += vP[1]*mDeriv[i-1][1:];
}
mDeriv = dropr(mDeriv,0);
decl sigma2_parameter = sumc((log(vY[cMaxpq:])-vLambda[cMaxpq:]).^(2))/(cT-cMaxpq) ;
decl vU = log(vEps[1:])/sigma2_parameter;
mDeriv = vU.*mDeriv;
vScore[0][0:2] = sumc(mDeriv)/(cT-cMaxpq);
}

DySco::MEMLNORMAnaCovMatrix(mAnaHess,const vO, const cT,const cMaxpq, const vEps)
{  decl a,b,c,d,e,f,g,fact,sigma_e2;
	f = exp(0.5);
	g = exp(2);
	a =vO[1]-vO[2]*f;
	b = vO[1]^(2)-2*vO[1]*vO[2]*f+vO[2]*g;
	c = vO[1]*f-vO[2]*g;
	d = vO[2]*f/(1-vO[1]);
	e = vO[2]^(2)*g/(1-vO[1]) + d^(2);
	sigma_e2 = (exp(1)-1)*exp(1);
	decl sigma_u_2 = sumc(log(vEps[cMaxpq:]).^(2))/(cT-cMaxpq);
	fact = 1/((1-b)*sigma_u_2);
 
		if (b>=1){oxwarning("b>1"); println("b is:" ,b ~vO');continue;}
	//1) Dynamic parameters
	mAnaHess[0][0][0] = fact*(1-vO[1])^(2)*(1+a)/((1-a));
	mAnaHess[0][1][1] =	fact*(e+ 2*a*(1-vO[1]*a)^(-1)*(vO[2]^(2)*c*f*((1-vO[1])*(1-a))^(-1)+vO[1]*e+vO[2]*d*f));
	mAnaHess[0][2][2] =	fact*(g*(1-a)+2*c*f  )/(1-a);
	mAnaHess[0][0][1] =	fact*(1-vO[1])*(d+ vO[2]*a*f/((1-vO[1])*(1-a)) +a*(vO[1]*d+vO[2]*f+ vO[2]*c/(1-a))/(1-vO[1]*a)   );
	mAnaHess[0][1][0] =	mAnaHess[0][0][1]; 
	mAnaHess[0][0][2] =	fact*(1-vO[1])*(f+c)/(1-a);
	mAnaHess[0][2][0] =	mAnaHess[0][0][2];
	mAnaHess[0][1][2] =	fact*(d*f+ a*(vO[1]*d*f+vO[2]*g+vO[2]*c*f*(1-a)^(-1))*(1-vO[1]*a)^(-1) +vO[2]*c*f/((1-vO[1])*(1-a))  );
	mAnaHess[0][2][1] =	 mAnaHess[0][1][2];

mAnaHess[0] = -mAnaHess[0];
return mAnaHess[0];
}



////////////////////////////////////////////////////////////////////////////////
/////////////////// GGA  Distribution ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
DySco::GGALambdaFilter(vLambda,vU, const vY, const vP, const vStruct, const cT, const maxpq)
{decl p = vStruct[0]; decl q = vStruct[1];
decl vLevAux=new matrix[cT][1];
if (vStruct[3]==1) {vLevAux=GetLevAux(&vLevAux, m_mLev);}
decl u = 0;
vLambda[0][0:maxpq-1] = vP[0];//log(vY[0]);
decl amPar =  new array[7];
SplitPara(&amPar,vP,vStruct);
for(decl i=maxpq; i<cT;i++){
	vLambda[0][i] = amPar[0][0]*(1-sumc(amPar[1][0]))+amPar[1][0]'*reversec(vLambda[0][i-p:i-1])+amPar[2][0]'*reversec(vU[0][i-q:i-1]);
	if (vStruct[3]==1){vLambda[0][i] += amPar[4][0]*vLevAux[i-1][0];}
	vU[0][i] = (vY[i]*exp(-vLambda[0][i])).^(amPar[3][0][0]) - amPar[3][0][1];
	}
}
DySco::GGALik(const vP, const adFunc, const avScore, const amHessian) 
{ decl cT = m_iT2sel-m_iT1sel+1;  decl p = m_vStruct[0]; decl q = m_vStruct[1]; decl cMaxpq = max(p, q);
decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
decl vU = new matrix[m_iT2sel-m_iT1sel+1][1];
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
if (m_iTwoComp==0) GGALambdaFilter(&vLambda,&vU, m_vY, vP, m_vStruct, cT, cMaxpq);
else  LambdaFilter2comp(&vLambda, &vU, m_vY, vP, m_vStruct, GGA);

adFunc[0] = double(( (cT-cMaxpq)*(log(amPar[3][0][0])-loggamma(amPar[3][0][1])) -amPar[3][0][0]*amPar[3][0][1]*sumc(vLambda[cMaxpq:]) + (amPar[3][0][0]*amPar[3][0][1]-1)*sumc(log(m_vY[cMaxpq:]))- sumc( (m_vY[cMaxpq:]./exp(vLambda[cMaxpq:])).^(amPar[3][0][0]) ))/cT );

//if (isarray(avScore)){  // if analytical score should be used in maximisation
//	decl vScore = new matrix[5][1];
//	decl eps = m_mData[m_iT1sel:m_iT2sel]./exp(vLambda);
//	decl cT = rows(m_mData[m_iT1sel:m_iT2sel]);
//	GGAScore(&vScore, eps,vLambda, vP,max(m_vStruct[1],m_vStruct[2]), cT);
//	println("vScore", vScore');
//	(avScore[0])[0] = vScore[0];
//	(avScore[0])[1] = vScore[1];
//	(avScore[0])[2] = vScore[2];
//	(avScore[0])[3] = vScore[3];
//	(avScore[0])[4] = vScore[4];
//}
if (amHessian){
	decl mAnaHess = new matrix[5][5];
	GGAAnaCovMatrix(&mAnaHess,vP, (m_iT1sel-m_iT2sel+1), 1);
    (amHessian[0])= mAnaHess;
}
return 1;
}


DySco::GGALikVal(const vP, const cStart)
{ decl cT =rows(m_mData[][0]);	decl cTIn = cT-cStart;
decl p = m_vStruct[0]; decl q = m_vStruct[1]; decl maxpq = max(p, q);
decl vLambda = 		new matrix[cT][1];
if(m_iModelClass == SCORE){ 
decl vU = 			new matrix[cT][1];
if (m_iTwoComp==0) 	GGALambdaFilter(&vLambda,&vU, m_mData[][0], vP, m_vStruct, cT, maxpq);
else  				LambdaFilter2comp(&vLambda, &vU, m_mData[][0], vP, m_vStruct, GGA);
}
else if (m_iModelClass == MEM){
	MEMLambdaFilter(&vLambda, m_mData[][0], vP, m_vStruct);
}
//DrawTMatrix(0,vLambda', {"Lambda"});
//ShowDrawWindow();
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
decl vLikContr = new matrix[cTIn][1];
 for(decl i=1;i<cTIn;i++){
vLikContr[i] = log(amPar[3][0][0])-loggamma(amPar[3][0][1]) -amPar[3][0][0]*amPar[3][0][1]*vLambda[cStart+i] + (amPar[3][0][0]*amPar[3][0][1]-1)*log(m_mData[cStart+i][0]) -  (m_mData[cStart+i][0]/exp(vLambda[cStart+i])).^(amPar[3][0][0]) ;
}
decl LikVal = sumc(vLikContr);
return LikVal;
}

DySco::GGAAnaCovMatrix(mAnaHess, const vO, const cT, const cMaxpq)
// Fix for burr was to change Ana Cov Matrix!
{decl a,b,c;
	a = vO[1]-vO[2]*vO[3]*vO[4];
	b = vO[1]^(2)-2*vO[1]*vO[2]*vO[3]*vO[4]+vO[2]^(2)*(1+vO[3])*vO[3]*vO[4]^(2);
	c = -vO[2]*vO[3]*vO[4];	
	decl fact = vO[3]*vO[4]^(2)/(1-b); //vO[3]/(1-b);//vO[3]*vO[4]^(2)/(1-b); // factor in front of dynamic elements of hessian
	decl sigma_u2 = vO[3];
	// 1) dynamic parameters
	mAnaHess[0][0][0] =  fact*(1-vO[1])^(2)*(1+a)/(1-a) ;
	mAnaHess[0][1][1] =  fact*vO[2]^(2)*sigma_u2*(1+a*vO[1])/((1-vO[1]^(2))*(1-a*vO[1]));
	mAnaHess[0][2][2] =  fact*sigma_u2;
	mAnaHess[0][0][1] =  fact*a*c*vO[2]*(1-vO[1])/((1-a)*(1-a*vO[1]));
	mAnaHess[0][1][0] =  mAnaHess[0][0][1];
	mAnaHess[0][0][2] =  fact*c*(1-vO[1])/(1-a);
	mAnaHess[0][2][0] =  mAnaHess[0][0][2];
	mAnaHess[0][1][2] =  fact*a*vO[2]*sigma_u2/(1-a*vO[1]);
	mAnaHess[0][2][1] =  mAnaHess[0][1][2];


	mAnaHess[0][3][3] = polygamma(vO[3], 1);
	//mAnaHess[0][4][4] =  vO[4]^(-2)*(1- polygamma(vO[3],0)*(2+polygamma(vO[3],0))+vO[3]*polygamma(vO[3],1) );
	mAnaHess[0][4][4] = vO[4]^(-2)*(1+vO[3]*(polygamma(vO[3]+1,0)^(2)+polygamma(vO[3]+1,1) ) );
//	mAnaHess[0][4][4] = -vO[4]^(-2)*(1-(polygamma(vO[3]+1,0)^(2)+polygamma(vO[3]+1,1) ) );
	mAnaHess[0][3][4] = polygamma(vO[3],0)/vO[4];
	mAnaHess[0][4][3] = mAnaHess[0][3][4];
	 // 3) correlation dynamic/scale parameter
	mAnaHess[0][0][4] =	-vO[3]*polygamma(vO[3]+1,0)*(1-vO[1])/(vO[4]*(1-a));
//	mAnaHess[0][0][4] = -1*(vO[3]*polygamma(vO[3],0)+1)*(1-vO[1])/(1-a);
//	mAnaHess[0][0][4] = -vO[4]*gammafact(vO[3])*polygamma(1+vO[3],0)*(1-vO[1])/(1-a);
	mAnaHess[0][4][0] = mAnaHess[0][0][4];
	mAnaHess[0][0][3] = vO[4]*(1-vO[1])/(1-a);
	mAnaHess[0][3][0] = mAnaHess[0][0][3];
	
	mAnaHess[0] = -mAnaHess[0];

//	println(mAnaHess[0]);
}


////////////////////////// MEM  stuff ///////////////////////////////
DySco::GGAMEMLik(const vP, const adFunc, const avScore, const amHessian) 
{decl cT = rows(m_vY); decl p = m_vStruct[0]; decl q = m_vStruct[1]; decl cMaxpq = max(p, q);
decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) MEMLambdaFilter(&vLambda, m_vY, vP, m_vStruct);
else 	MEMLambdaFilter2comp(&vLambda, m_vY, vP, m_vStruct);
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
adFunc[0] = double(( (cT-cMaxpq)*(log(amPar[3][0][0])-loggamma(amPar[3][0][1])) -amPar[3][0][0]*amPar[3][0][1]*sumc(vLambda[cMaxpq:]) + (amPar[3][0][0]*amPar[3][0][1]-1)*sumc(log(m_vY[cMaxpq:]))- sumc( (m_vY[cMaxpq:]./exp(vLambda[cMaxpq:])).^(amPar[3][0][0]) ))/cT );

//if (isarray(avScore)){  // if analytical score should be used in maximisation
//	decl vScore = new matrix[5][1];
//	GGAMEMScore(&vScore,m_vY,vLambda, vP, max(m_vStruct[1], m_vStruct[2]));
////	println(vScore');
//	(avScore[0])[0] = vScore[0];
//	(avScore[0])[1] = vScore[1];
//	(avScore[0])[2] = vScore[2];
//	(avScore[0])[3] = vScore[3];
//	(avScore[0])[3] = vScore[4];
//}
return 1;
}

DySco::MEMGGAnaCovMatrix(const mAnaHess, const vO, const cT, const cMaxpq){
// Fix for burr was to change Ana Cov Matrix!
decl a,b,c,d,e,f,g,sigma_e2,fact;
	f = gammafact(vO[3]+1/vO[4])/gammafact(vO[3]);
	g = gammafact(vO[3]+2/vO[4])/gammafact(vO[3]);
	sigma_e2 = g-f^(2);
		
	a = vO[1] - vO[2]*f;
	b = vO[1]^(2)-2*vO[1]*vO[2]*f + vO[2]^(2)*g;
	c = vO[1]*f -vO[2]*g;
	d = vO[2]*f/(1-vO[1]);
	e = vO[2]*sigma_e2/(1-vO[1]^(2))+d^(2);

	fact = vO[3]*vO[4]^(2)/(1-b);
	
//	println("a ",a,"b ",b,"c ",c,"d ",d,"e ",e,"f ",f);
//	println("sigma ",sigma_e2);
		if (b>=1){oxwarning("b>1"); println("b is:" ,b ~vO');continue;}
	//1) Dynamic parameters
	mAnaHess[0][0][0] = fact*(1-vO[1])^(2)*(1+a)/((1-a));
	mAnaHess[0][1][1] =	fact*(e+ 2*a*(1-vO[1]*a)^(-1)*(vO[2]^(2)*c*f*((1-vO[1])*(1-a))^(-1)+vO[1]*e+vO[2]*d*f));
	mAnaHess[0][2][2] =	fact*(g*(1-a)+2*c*f  )/(1-a);
	mAnaHess[0][0][1] =	fact*(1-vO[1])*(d+ vO[2]*a*f/((1-vO[1])*(1-a)) +a*(vO[1]*d+vO[2]*f+ vO[2]*c/(1-a))/(1-vO[1]*a)   );
	mAnaHess[0][1][0] =	mAnaHess[0][0][1]; 
	mAnaHess[0][0][2] =	fact*(1-vO[1])*(f+c)/(1-a);
	mAnaHess[0][2][0] =	mAnaHess[0][0][2];
	mAnaHess[0][1][2] =	fact*(d*f+ a*(vO[1]*d*f+vO[2]*g+vO[2]*c*f*(1-a)^(-1))*(1-vO[1]*a)^(-1) +vO[2]*c*f/((1-vO[1])*(1-a))  );
	mAnaHess[0][2][1] =	 mAnaHess[0][1][2];

	//2) cross derivatives
	mAnaHess[0][3][3] = polygamma(vO[3], 1);
	mAnaHess[0][4][4] = vO[4]^(-2)*(1-(polygamma(vO[3]+1,0)^(2)-polygamma(vO[3]+1,1) ) );

	mAnaHess[0][3][4] = -polygamma(vO[3],0)/vO[4];
	mAnaHess[0][4][3] = mAnaHess[0][3][4];

	 // 3) correlation dynamic/scale parameter
	mAnaHess[0][0][4] = vO[4]*gammafact(vO[3])*polygamma(1+vO[3],0)*(1-vO[1])/(1-a);
	mAnaHess[0][4][0] = mAnaHess[0][0][4];
	mAnaHess[0][0][3] = vO[4]*(1-vO[1])/(1-a);
	mAnaHess[0][3][0] = mAnaHess[0][0][3];


	mAnaHess[0][0][4] = -(vO[3]*polygamma(vO[3],0))*(1-vO[1])/(1-a);
	mAnaHess[0][0][3] = vO[4]*(1-vO[1])/(1-a);
	mAnaHess[0][4][0] = mAnaHess[0][0][4];
	mAnaHess[0][3][0] = mAnaHess[0][0][3];
	
	mAnaHess[0][1][4] = -(vO[3]*polygamma(vO[3],0))*vO[2]*f/((1-vO[1])*(1-a));
	mAnaHess[0][1][3] = vO[4]*vO[2]*f/((1-vO[1])*(1-a));
	mAnaHess[0][4][1] = mAnaHess[0][1][4];
	mAnaHess[0][3][1] = mAnaHess[0][1][3];

	
	mAnaHess[0][2][4] = -(vO[3]*polygamma(vO[3],0))*f/(1-a);
	mAnaHess[0][2][3] = vO[4]*f/(1-a);
	mAnaHess[0][4][2] = mAnaHess[0][2][4];
	mAnaHess[0][3][2] = mAnaHess[0][2][3];
	
	mAnaHess[0] = -mAnaHess[0];

}

DySco::MEMGGScore(const vScore, const eps,const vLambda, const vP, const cT, const cMaxpq)
{
decl mDeriv = new matrix[cT][3];
decl f = gammafact(vP[3]+1/vP[4])/gammafact(vP[3]);
decl a = vP[1]-vP[2]*f;
mDeriv[0][0] = (1-vP[1])/(1-a);///(vP[1]-vP[2]*vP[3]);
mDeriv[0][1] = vP[2]*f/((1-vP[1])*(1-a));
mDeriv[0][2] = f/(1-a);
mDeriv[1:][0] = constant(1-vP[1],cT-1,1);
mDeriv[1:][1] = -ones(cT-1,1)*vP[0]+vLambda[:cT-2];
mDeriv[1:][2] = eps[:cT-2];
decl x;
for (decl i=1;i<cT;i++){
 x = vP[1]-vP[2]*eps[i-1];
 mDeriv[i][0] +=x*mDeriv[i-1][0];
}
mDeriv = (eps[1:].^(vP[4])-vP[3]*ones(cT-1,1) ).*mDeriv[1:][];
vScore[0][0:2] = meanc(mDeriv);
delete mDeriv;
vScore[0][3] = -polygamma(vP[3],0) + vP[4]*meanc(log(eps[1:]));
vScore[0][4] = vP[4]^(-1) + meanc(log(eps[1:]).*(vP[3]-eps[1:].^(vP[4])) );
//println(vScore);
return vScore;

}



////////////////////////////////////////////////////////////////////////////////
///////////////////  LLOG Distribution ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
DySco::LLOGLambdaFilter(vLambda,vU, const vY, const vP, const vStruct)
{decl cT = rows(vY); decl p = vStruct[0]; decl q = vStruct[1]; decl maxpq = max(p, q);
decl vLevAux=new matrix[cT][1];
if (vStruct[3]==1) {vLevAux=GetLevAux(&vLevAux, m_mLev);}
decl u = 0;
vLambda[0][0:maxpq-1] = vP[0];//log(vY[0]);
decl amPar =  new array[7];
SplitPara(&amPar,vP,vStruct);
for(decl i=maxpq; i<cT;i++){
	vLambda[0][i] = amPar[0][0]*(1-sumc(amPar[1][0]))+amPar[1][0]'*reversec(vLambda[0][i-p:i-1])+amPar[2][0]'*reversec(vU[0][i-q:i-1]);
	if (vStruct[3]==1){vLambda[0][i] += amPar[4][0]*vLevAux[i-1][0];}
	vU[0][i] = -1+2*vY[i]^(amPar[3][0])*exp(-amPar[3][0]*vLambda[0][i])/(1+vY[i]^(amPar[3][0])*exp(-amPar[3][0]*vLambda[0][i]));
	}
}

DySco::LLOGLik(const vP, const adFunc, const avScore, const amHessian) 
{
decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
decl vU = new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) LLOGLambdaFilter(&vLambda,&vU, m_vY, vP, m_vStruct);
else  LambdaFilter2comp(&vLambda, &vU, m_vY, vP, m_vStruct, LLOG);

decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
adFunc[0] = double(log(amPar[3][0]) - amPar[3][0]*(m_iT2sel-m_iT1sel)^(-1)*sumc(vLambda[1:]) + (amPar[3][0]-1)*(m_iT2sel-m_iT1sel)^(-1)*sumc(log(m_vY[1:])) -2*(m_iT2sel-m_iT1sel)^(-1)*sumc(log(1+(m_vY[1:]./exp(vLambda[1:])).^(amPar[3][0]))) );

if (isarray(avScore)){  // if analytical score should be used in maximisation
	decl vScore = new matrix[4][1];
	decl eps = m_vY./exp(vLambda);
	LLOGScore(&vScore, eps,vLambda, vP,max(m_vStruct[1],m_vStruct[2]));
	//println("vScore", vScore');
	(avScore[0])[0] = vScore[0];
	(avScore[0])[1] = vScore[1];
	(avScore[0])[2] = vScore[2];
	(avScore[0])[3] = vScore[3];
}
if (amHessian){
	decl mAnaHess = new matrix[4][4];
	LLOGAnaCovMatrix(&mAnaHess,vP, (m_iT1sel-m_iT2sel+1), 1);
    (amHessian[0])= mAnaHess;
}
return 1;
}

DySco::LLOGLikVal(const vP, const cStart)
{ decl cT =rows(m_mData[][0]);	decl cTIn = cT-cStart;
decl vLambda = 		new matrix[cT][1];
decl vU = 			new matrix[cT][1];
decl p = m_vStruct[0]; decl q = m_vStruct[1]; decl maxpq = max(p, q);
if (m_iTwoComp==0) 	LLOGLambdaFilter(&vLambda,&vU, m_mData[][0], vP, m_vStruct);
else  				LambdaFilter2comp(&vLambda, &vU, m_mData[][0], vP, m_vStruct, LLOG);
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
decl vLikContr = new matrix[cTIn][1];
 for(decl i=1;i<cTIn;i++){
vLikContr[i] = log(amPar[3][0]) - amPar[3][0]*vLambda[cStart+i] + (amPar[3][0]-1)*log(m_mData[cStart+i][0]) -2*log(1+(m_mData[cStart+i][0]/exp(vLambda[cStart+i])).^(amPar[3][0]));

}
decl LikVal = sumc(vLikContr);
return LikVal;
}


DySco::LLOGScore(vScore, const vEps,const vLambda, const vP, const cMaxpq)
{decl cT = rows(vEps);
if (rows(vLambda)!=cT || rows(vEps)!=cT) {oxwarning("length of eps and lambda different. Pls correct"); exit(0);}
decl mDeriv = new matrix[cT][3];
mDeriv[0][0] = (1-vP[1]); ///(vP[1]-vP[2]*vP[3]);
mDeriv[0][1] = mDeriv[0][2] = 0;
mDeriv[1:][0] = constant(1-vP[1],cT-1,1);
mDeriv[1:][1] = -ones(cT-1,1)*vP[0]+vLambda[:cT-2];
mDeriv[1:][2] = -ones(cT-1,1)+2*vEps[:cT-2].^(vP[3])./(1+vEps[:cT-2].^(vP[3]));
decl x;
for (decl i=1; i<cT; i++){
x = vP[1]-2*vP[2]*vP[3]*vEps[i-1]^(vP[3])*(1+vEps[i-1]^(vP[3]))^(-2);
mDeriv[i][] += x*mDeriv[i-1][];
}
mDeriv = vP[3]*(-ones(cT-1,1)+2*vEps[1:].^(vP[3])./(1+vEps[1:].^(vP[3]))).*mDeriv[1:][];
vScore[0][0:2] = sumc(mDeriv)/(cT-1);
delete mDeriv;
vScore[0][3] = 1/vP[3] + (sumc(log(vEps[1:]))-2*sumc(vEps[1:].^(vP[3]).*log(vEps[1:])./(1+vEps[1:].^(vP[3]))))/(cT-1);
}

DySco::LLOGAnaCovMatrix(const mAnaHess, const vO, const cT, const cMaxpq)
{
decl a,b,c;
	a = vO[1]-vO[2]*vO[3]/3;
	b = vO[1]^(2)-2*vO[1]*vO[2]*vO[3]/3 + 2*vO[2]^(2)*vO[3]^(2)/15;
	c = 0;	
	decl fact = vO[3]^(2)/(3*(1-b)); // factor in front of dynamic elements of hessian	
//	decl aux = new matrix[4][4]; // omega notation, Andrews version with nu in score (taken from Andres&Harvey*(2012))
	// 1) dynamic parameters
	mAnaHess[0][0][0] =  fact*(1-vO[1])^(2)*(1+a)/(1-a) ;
	mAnaHess[0][1][1] =	 fact*vO[2]^(2)*3^(-1)*(1+a*vO[1])/((1-vO[1]^(2))*(1-a*vO[1]));
	mAnaHess[0][2][2] =  fact*3^(-1);
	mAnaHess[0][0][1] =  fact*a*c*vO[2]*(1-vO[1])/((1-a)*(1-a*vO[1]));
	mAnaHess[0][1][0] =  mAnaHess[0][0][1];
	mAnaHess[0][0][2] =  fact*c*(1-vO[1])/(1-a);
	mAnaHess[0][2][0] =  mAnaHess[0][0][2];
	mAnaHess[0][1][2] =  fact*a*vO[2]*3^(-1)/(1-a*vO[1]);
	mAnaHess[0][2][1] =  mAnaHess[0][1][2];
												
	mAnaHess[0][3][3] = vO[3]^(-2)*(1+M_PI^(2)/3)/3;
	mAnaHess[0][3][0] =	mAnaHess[0][0][3] = mAnaHess[0][3][1] = mAnaHess[0][3][2]= mAnaHess[0][1][3]=mAnaHess[0][3][2]=0;

	mAnaHess[0] = -mAnaHess[0];
return mAnaHess[0];
}

////////////////// MEM Stuff //////////////////////////////////////////////
DySco::LLOGMEMLik(const vP, const adFunc, const avScore, const amHessian) 
{decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) MEMLambdaFilter(&vLambda, m_vY, vP, m_vStruct);
else 	MEMLambdaFilter2comp(&vLambda, m_vY, vP, m_vStruct);
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
adFunc[0] = double(log(amPar[3][0]) - amPar[3][0]*(m_iT2sel-m_iT1sel)^(-1)*sumc(vLambda[1:]) + (amPar[3][0]-1)*(m_iT2sel-m_iT1sel)^(-1)*sumc(log(m_vY[1:])) -2*(m_iT2sel-m_iT1sel)^(-1)*sumc(log(1+(m_vY[1:]./exp(vLambda[1:])).^(amPar[3][0]))) );
if (isarray(avScore)){  // if analytical score should be used in maximisation
	decl vScore = new matrix[4][1];
	LLOGMEMScore(&vScore, m_vY,vLambda, vP, max(m_vStruct[1], m_vStruct[2]));
//	println(vScore');
	(avScore[0])[0] = vScore[0];
	(avScore[0])[1] = vScore[1];
	(avScore[0])[2] = vScore[2];
	(avScore[0])[3] = vScore[3];
}
return 1;
}

DySco::LLOGMEMScore(const vScore, const vY,const vLambda, const vP, const cMaxpq)
{  decl cT = rows(vY);
decl mDeriv = new matrix[cT][3];
mDeriv[0][0] = 1;
mDeriv[0][1] = meanc(-vP[0]+vLambda)/(1-vP[1]);
mDeriv[0][2] = meanc(vY)/(1-vP[1]);
mDeriv[1:][0] = constant(1-vP[1],cT,1);
mDeriv[1:][1] = -vP[0] + vLambda[:cT-2];
mDeriv[1:][2] = vY[:cT-2]./exp(vLambda[:cT-2]);
for (decl i=cMaxpq; i<cT; i++){
mDeriv[i][1:] += vP[1]*mDeriv[i-1][1:];
}
//DrawTMatrix(0,mDeriv[1:][1]',{"M1"});
//ShowDrawWindow();
decl eps = vY./exp(vLambda);
mDeriv = dropr(mDeriv,0);
decl vU = vP[3]*(-1+ 2*eps.^(vP[3])./(1+eps.^(vP[3])));
//println("mean u", meanc(vU));
mDeriv = vU[1:].*mDeriv;
vScore[0][0:2] = sumc(mDeriv)/(cT-cMaxpq);
delete mDeriv;
vScore[0][3] = 1/vP[3] + (-sumc(vLambda[1:])+ sumc(log(vY[1:])) -2*sumc(eps[1:].^(vP[3]).*log(eps[1:])./(1+eps[1:].^(vP[3]))) )/cT;
return vScore ;
}

DySco::MEMLLOGAnaCovMatrix(const mAnaHess, const vO, const cT, const cMaxpq)
{decl a,b,c,d,e,f,g,sigma_e2,fact;
	f = gammafact(1+1/vO[3])*gammafact(1-1/vO[3])/gammafact(1);
	g = gammafact(1+2/vO[3])*gammafact(1-2/vO[3])/gammafact(1);
	sigma_e2 = g-f^(2);
		
	a = vO[1] - vO[2]*f;
	b = vO[1]^(2)-2*vO[1]*vO[2]*f + vO[2]^(2)*g;
	c = vO[1]*f -vO[2]*g;
	d = vO[2]*f/(1-vO[1]);
	e = vO[2]*sigma_e2/(1-vO[1]^(2))+d^(2);

	fact = vO[3]^(2)*1/((1+2)*(1-b));
	
//	println("a ",a,"b ",b,"c ",c,"d ",d,"e ",e,"f ",f);
//	println("sigma ",sigma_e2);
		if (b>=1){oxwarning("b>1"); println("b is:" ,b ~vO');continue;}
	//1) Dynamic parameters
	mAnaHess[0][0][0] = fact*(1-vO[1])^(2)*(1+a)/((1-a));
	mAnaHess[0][1][1] =	fact*(e+ 2*a*(1-vO[1]*a)^(-1)*(vO[2]^(2)*c*f*((1-vO[1])*(1-a))^(-1)+vO[1]*e+vO[2]*d*f));
	mAnaHess[0][2][2] =	fact*(g*(1-a)+2*c*f  )/(1-a);
	mAnaHess[0][0][1] =	fact*(1-vO[1])*(d+ vO[2]*a*f/((1-vO[1])*(1-a)) +a*(vO[1]*d+vO[2]*f+ vO[2]*c/(1-a))/(1-vO[1]*a)   );
	mAnaHess[0][1][0] =	mAnaHess[0][0][1]; 
	mAnaHess[0][0][2] =	fact*(1-vO[1])*(f+c)/(1-a);
	mAnaHess[0][2][0] =	mAnaHess[0][0][2];
	mAnaHess[0][1][2] =	fact*(d*f+ a*(vO[1]*d*f+vO[2]*g+vO[2]*c*f*(1-a)^(-1))*(1-vO[1]*a)^(-1) +vO[2]*c*f/((1-vO[1])*(1-a))  );
	mAnaHess[0][2][1] =	 mAnaHess[0][1][2];

	mAnaHess[0][3][3] =	(vO[3]^(2)*(2+1))^(-1)*(2+1+1*(polygamma(1,1)+polygamma(1,1) + (polygamma(1,0)-polygamma(1,0) + (1-1)*1^(-1))^(2) - (1+1^(2))*1^(-2) ) );

	mAnaHess[0][0][3] = (1-1+1*(polygamma(1,0)-polygamma(1,0)))*(1-vO[1])/((2+1)*(1-a));
	mAnaHess[0][3][0] = mAnaHess[0][0][3];
	mAnaHess[0][0][4] = -vO[3]*(1-vO[1])/( (1+1)*(1-a) );
	mAnaHess[0][4][0] = mAnaHess[0][0][4];

	mAnaHess[0][1][3] = (1-1+1*(polygamma(1,0)-polygamma(1,0)))/((2+1))   *(vO[2]*f/((1-vO[1])*(1-a)));
	mAnaHess[0][3][1] = mAnaHess[0][1][3];

	
	mAnaHess[0][2][3] = (1-1+1*(polygamma(1,0)-polygamma(1,0)))/((2+1))*(f/(1-a));
	mAnaHess[0][3][2] = mAnaHess[0][2][3];
	
	mAnaHess[0] = -mAnaHess[0];
return mAnaHess[0];
}


////////////////////////////////////////////////////////////////////////////////
/////////////////// Burr Distribution ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
DySco::BURRLambdaFilter(vLambda,vU, const vY, const vP, const vStruct)
{ decl cT=rows(vY);	 decl p = vStruct[0]; decl q = vStruct[1];decl maxpq = max(p, q);
decl vLevAux=new matrix[cT][1];
if (vStruct[3]==1) {vLevAux=GetLevAux(&vLevAux, m_mLev);}
decl u = 0; 
vLambda[0][0:maxpq-1] = vP[0];
decl amPar =  new array[7];
SplitPara(&amPar,vP,vStruct);
//decl vU = new matrix[cT][1];
//println(amPar);
//println(amPar[3][0][0]~	amPar[3][0][1]);
//println( amPar[0][0]	 ~amPar[1][0] ~		 amPar[2][0]);
for(decl i=maxpq; i<cT;i++){
	vLambda[0][i] = amPar[0][0]*(1-sumc(amPar[1][0]))+amPar[1][0]'*reversec(vLambda[0][i-p:i-1])+amPar[2][0]'*reversec(vU[0][i-q:i-1]);
	if (vStruct[3]==1){vLambda[0][i] += amPar[4][0]*vLevAux[i-1][0];}
	vU[0][i] = -1+(amPar[3][0][1]+1)*vY[i]^(amPar[3][0][0])*exp(-amPar[3][0][0]*vLambda[0][i])/(1+vY[i]^(amPar[3][0][0])*exp(-amPar[3][0][0]*vLambda[0][i]));
	}
//return vLambda[0]~vU;
}
 
DySco::BURRLik(const vP, const adFunc, const avScore, const amHessian)
{decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
decl vU = new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) BURRLambdaFilter(&vLambda,&vU,m_vY, vP, m_vStruct);
else  LambdaFilter2comp(&vLambda, &vU, m_vY, vP, m_vStruct, BURR);

decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
//println(amPar[3][0][0]~amPar[3][0][1]);
//println(m_iT2sel-m_iT1sel);
adFunc[0] = double(log(amPar[3][0][1]) + log(amPar[3][0][0]) + (m_iT2sel-m_iT1sel)^(-1)*(amPar[3][0][0]-1)*sumc(log(m_vY[1:]))- (m_iT2sel-m_iT1sel)^(-1)*amPar[3][0][0]*sumc(vLambda[1:]) - (m_iT2sel-m_iT1sel)^(-1)*(amPar[3][0][1]+1)*sumc(log((m_vY[1:]./exp(vLambda[1:])).^(amPar[3][0][0])+1)));

//adFunc[0] =	double( log(amPar[3][0][0])+log(amPar[3][0][1]) + (m_iT2sel-m_iT1sel)^(-1)*(amPar[3][0][0]-1)*sumc(log(m_mData[m_iT1sel+1:m_iT2sel])) -(m_iT2sel-m_iT1sel)^(-1)*amPar[3][0][0]*sumc(vLambda[1:])
//					-(m_iT2sel-m_iT1sel)^(-1)*(amPar[3][0][1]+1)*sumc(log(m_mData[m_iT1sel+1:m_iT2sel]./exp(vLambda[1:]).^(amPar[3][0][0])+1)) );

//	return double( ((cT-cMaxpq)*(log(vShapePara[0])+log(vShapePara[1])) + (vShapePara[0]-1)*sumc(log(vY_))
//					- vShapePara[0]*sumc(vLambda_)-(vShapePara[1]+1)*sumc(log( (vY_./exp(vLambda_)).^(vShapePara[0])+1)))/m_cT)
//

if (isarray(avScore)){  // if analytical score should be used in maximisation
	decl vScore = new matrix[5][1];
	decl eps = m_vY./exp(vLambda);
	BURRScore(&vScore, eps,vLambda, vP,max(m_vStruct[0], m_vStruct[1]));
	//println("vScore ", vScore');
//    decl vScore = Burr2(vP, m_cT, 1);
	(avScore[0])[0] = vScore[0];
	(avScore[0])[1] = vScore[1];
	(avScore[0])[2] = vScore[2];
	(avScore[0])[3] = vScore[3];
	(avScore[0])[4] = vScore[4];
}
if (amHessian){
	decl mAnaHess = new matrix[5][5];
	BURRAnaCovMatrix(&mAnaHess,vP, (m_iT2sel-m_iT1sel+1), max(m_vStruct[0], m_vStruct[1]));
    (amHessian[0])= mAnaHess;
}
return 1;
}


DySco::BURRLikVal(const vP, const cStart)
{ decl cT =rows(m_mData[][0]);	decl cTIn = cT-cStart;
decl vLambda = 		new matrix[cT][1];
if(m_iModelClass == SCORE){ 
decl vU = 			new matrix[cT][1];
if (m_iTwoComp==0) 	BURRLambdaFilter(&vLambda,&vU, m_mData[][0], vP, m_vStruct);
else  				LambdaFilter2comp(&vLambda, &vU, m_mData[][0], vP, m_vStruct, BURR);
}
else if (m_iModelClass == MEM){
	MEMLambdaFilter(&vLambda, m_mData[][0], vP, m_vStruct);
}
//DrawTMatrix(0,vLambda', {"Lambda"});
//ShowDrawWindow();
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
decl vLikContr = new matrix[cTIn][1];
 for(decl i=1;i<cTIn;i++){
vLikContr[i] = log(amPar[3][0][1]) + log(amPar[3][0][0]) + (amPar[3][0][0]-1)*log(m_mData[cStart+i][0]) - amPar[3][0][0]*vLambda[cStart+i] - (amPar[3][0][1]+1)*(log((m_mData[cStart+i][0]./exp(vLambda[cStart+i])).^(amPar[3][0][0])+1));

}
decl LikVal = sumc(vLikContr);
return LikVal;
}

DySco::BURRScore(vScore, const vEps,const vLambda, const vP, const cMaxpq)
{decl cT = rows(vEps);
if (rows(vLambda)!=cT) {oxwarning("length of eps and lambda different. Pls correct"); exit(0);}
decl mDeriv = new matrix[cT][3];
mDeriv[0][0] = (1-vP[1]); ///(vP[1]-vP[2]*vP[3]);
mDeriv[0][1] = mDeriv[0][2] = 0;
mDeriv[1:][0] = constant(1-vP[1],cT-1,1);
mDeriv[1:][1] = -ones(cT-1,1)*vP[0]+vLambda[:cT-2];
mDeriv[1:][2] =	(1+vP[4])*(vEps[:cT-2].^(vP[3])./(1+vEps[:cT-2].^(vP[3])))-1;
decl x;
for (decl i=1; i<cT; i++){
x = vP[1]-vP[2]*vP[3]*(1+vP[4])*vEps[i-1]^(vP[3])*(1+vEps[i-1]^(vP[3]))^(-2);
mDeriv[i][] += x*mDeriv[i-1][];
}
mDeriv =vP[3]*((1+vP[4])*vEps[1:].^(vP[3])./(1+vEps[1:].^(vP[3]))-1).*mDeriv[1:][];
vScore[0][0:2] = sumc(mDeriv)/(cT-cMaxpq);
delete mDeriv;
vScore[0][3] = 1/vP[3] + (sumc(log(vEps[1:])))/(cT-cMaxpq) - (vP[4]+1)*(sumc(log(vEps[1:]).*vEps[1:].^(vP[3])./(1+vEps[1:].^(vP[3]))) )/(cT-cMaxpq);
vScore[0][4] = 1/vP[4] - sumc(log(vEps[1:].^(vP[3])+1))/(cT-cMaxpq);
}


DySco::Beta(const x, const y)
{// complete beta function
 decl beta= new matrix[1][1]; beta = gammafact(x)*gammafact(y)/gammafact(x+y);
 return beta;
}

DySco::BURRAnaCovMatrix(mAnaHess, const vO, const cT, const cMaxpq)
{decl a,b,c,sigma_u_2;
	a = vO[1]-vO[2]*vO[3]*vO[4]/(2+vO[4]);
	b = vO[1]^(2)-2*vO[1]*vO[2]*vO[3]*vO[4]/(2+vO[4]) + vO[2]^(2)*vO[3]^(2)*(vO[4]+1)^(2)*Beta(3,vO[4]+2)/(Beta(1,vO[4]));
	c = vO[2]*vO[3]*(vO[4]+1)*(Beta(1,vO[4]))^(-1)*(Beta(2,vO[4]+1) -(vO[4]+1)*Beta(3,vO[4]+1));
	sigma_u_2 = vO[4]/(2+vO[4]);
	decl fact = vO[3]^(2)*vO[4]/((vO[4]+2)*(1-b)); // factor in front of dynamic elements of hessian
	
	decl aux = new matrix[5][5]; // omega notation, Andrews version with nu in score (taken from Andres&Harvey*(2012))
	// 1) dynamic parameters
	mAnaHess[0][0][0] =  fact*(1-vO[1])^(2)*(1+a)/(1-a) ;
	mAnaHess[0][1][1] =	 fact*vO[2]^(2)*sigma_u_2*(1+a*vO[1])/((1-vO[1]^(2))*(1-a*vO[1]));
	mAnaHess[0][2][2] =  fact*sigma_u_2;
	mAnaHess[0][0][1] =  fact*a*c*vO[2]*(1-vO[1])/((1-a)*(1-a*vO[1]));
	mAnaHess[0][1][0] =  mAnaHess[0][0][1];
	mAnaHess[0][0][2] =  fact*c*(1-vO[1])/(1-a);
	mAnaHess[0][2][0] =  mAnaHess[0][0][2];
	mAnaHess[0][1][2] =  fact*a*vO[2]*sigma_u_2/(1-a*vO[1]);
	mAnaHess[0][2][1] =  mAnaHess[0][1][2];
												
	mAnaHess[0][3][3] = vO[3]^(-2) + vO[4]*vO[3]^(-2)*(2+vO[4])^(-1)*(M_PI^(2)/6 + polygamma(vO[4],1)+ (polygamma(vO[4],0)+M_EULER+ (1-vO[4])*vO[4]^(-1))^(2)-(1+vO[4]^(2))*vO[4]^(-2));
	//-vO[3]^(-2) - (vO[4]+1)*vO[3]^(-2)*Beta(2,vO[4]+1)*((polygamma(2,0)-polygamma(vO[4]+1,0))^(2)+polygamma(2,1)-polygamma(vO[4]+1,1));
	//mAnaHess[0][4][4] =	polygamma(vO[4],1)-polygamma(1+vO[4],1);//vO[4]^(-2);//vO[4]^(-2); //polygamma(vO[4],1)-polygamma(1+vO[4],1);

	//mAnaHess[0][4][4] = -vO[4]^(-2) + 4*vO[2]*vO[3]/((vO[4]+2)*(1+vO[4])) -vO[2]^(2)*vO[3]^(2)*(vO[4]+1)*Beta(3,1+vO[4])/Beta(1,vO[4]);
	mAnaHess[0][4][4] = -(-vO[4]^(-2)+2*vO[2]*vO[3]*(1+vO[4])^(-2)-2*vO[2]^(2)*vO[3]^(2)*(2+vO[4])^(-2)*Beta(2,1+vO[4])/Beta(1,vO[4])-2*vO[2]^(2)*vO[3]^(3)*(1+vO[4])^(-1)*Beta(2,1+vO[4])/Beta(1,vO[4]) );
	//mAnaHess[0][4][4] = -vO[4]^(-2);//-2*vO[2]^(2)*vO[3]^(2)*(2+vO[4])^(-2)*Beta(2,1+vO[4])/Beta(1,vO[4])+2*vO[2]^(2)*vO[3]^(3)*(1+vO[4])^(-1)*Beta(2,1+vO[4])/Beta(1,vO[4]);

	mAnaHess[0][3][4] =mAnaHess[0][4][3] =(-M_EULER-polygamma(vO[4],0)+1)*(vO[3]*(1+vO[4]))^(-1); //- (vO[4]*(-M_EULER-polygamma(vO[4],0))-1)*(vO[3]*(1+vO[4]))^(-1);
	//(1+vO[4]*(1-vO[4])+vO[4]*(1+vO[4])*(M_EULER+polygamma(vO[4]+2,0)))*(vO[3]*vO[4]*(1+vO[4])^(2))^(-1) ;

	mAnaHess[0][0][3] = mAnaHess[0][3][0] = -(vO[4]*(M_EULER-1+polygamma(vO[4],0))+1)*(2+vO[4])^(-1)*(1-vO[1])*(1-a)^(-1) ;//-(vO[4]*(polygamma(1,0)-polygamma(vO[4],0)-1)+1)*(1-vO[1])*(2+vO[4])^(-1)*(1-a)^(-1);
	//-(vO[4]*(M_EULER-1+polygamma(vO[4],0))+1)*(2+vO[4])^(-1)*(1-vO[1])*(1-a)^(-1) ;
	mAnaHess[0][1][3] = mAnaHess[0][3][1] = mAnaHess[0][2][3] = mAnaHess[0][3][2] = 0  ;		   


	//mAnaHess[0][0][4] = mAnaHess[0][4][0] = -vO[3]*(vO[4]+1)^(-1)*(1-vO[1])*(1-a)^(-1);
	mAnaHess[0][0][4] = mAnaHess[0][4][0] = -(vO[3]*(1+vO[4])^(-1)-vO[2]*vO[3]^(2)*Beta(2,1+vO[4])/Beta(1,vO[4]))*(1-vO[1])*(1-a)^(-1);
	//println((vO[2]*vO[3]^(2)*Beta(2,1+vO[4])/Beta(1,vO[4]))*(1-vO[1])*(1-a)^(-1));
	//-vO[3]*(vO[4]+1)^(-1)*(1-vO[1])*(1-a)^(-1); 
	mAnaHess[0][1][4] = mAnaHess[0][4][1] =	mAnaHess[0][2][4] = mAnaHess[0][4][2] =	0;

	mAnaHess[0] = -mAnaHess[0];
	//if(choleski(mAnaHess[0])==0) {println("Choleski failed"); continue;}
	//println("determinant", determinant(mAnaHess[0])); //println("inverse", invertgen(mAnaHess[0],30)); 
//return mAnaHess[0];
}
////////////////////////// BURR MEM  STUFF ////////////////////////////////////////////////////////////
DySco::BURRMEMLik(const vP, const adFunc, const avScore, const amHessian)
{ decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) MEMLambdaFilter(&vLambda, m_vY, vP, m_vStruct);
else 	MEMLambdaFilter2comp(&vLambda, m_vY, vP, m_vStruct);
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
adFunc[0] = double(log(amPar[3][0][1]) + log(amPar[3][0][0]) + (m_iT2sel-m_iT1sel)^(-1)*(amPar[3][0][0]-1)*sumc(log(m_vY[1:]))- (m_iT2sel-m_iT1sel)^(-1)*amPar[3][0][0]*sumc(vLambda[1:]) - (m_iT2sel-m_iT1sel)^(-1)*(amPar[3][0][1]+1)*sumc(log((m_vY[1:]./exp(vLambda[1:])).^(amPar[3][0][0])+1)));


if (isarray(avScore)){  // if analytical score should be used in maximisation
	decl vScore = new matrix[5][1];
	decl aux = m_vY./exp(vLambda);
	BURRMEMScores(&vScore,m_vY,vLambda, vP, max(m_vStruct[1], m_vStruct[2]));
//	println(vScore');
	(avScore[0])[0] = vScore[0];
	(avScore[0])[1] = vScore[1];
	(avScore[0])[2] = vScore[2];
	(avScore[0])[3] = vScore[3];
	(avScore[0])[4] = vScore[4];
}
return 1;
}

DySco::BURRMEMScores(vScore,const vY,const vLambda, const vP, const cMaxpq)
{decl cT = rows(vY); decl mDeriv = new matrix[cT][3];
mDeriv[0][0] = 1;
mDeriv[0][1] = meanc(-vP[0]+vLambda)/(1-vP[1]);
mDeriv[0][2] = meanc(vY)/(1-vP[1]);
mDeriv[1:][0] = constant(1-vP[1],cT-1,1);
mDeriv[1:][1] = -vP[0] + vLambda[:cT-2];
mDeriv[1:][2] = vY[:cT-2];
for (decl i=cMaxpq; i<cT; i++){
mDeriv[i][1:] += vP[1]*mDeriv[i-1][1:];
}
mDeriv = dropr(mDeriv,0);
decl vEps = vY./exp(vLambda);
decl vU = vP[3].*(-1+ (vP[4]+1)*vEps[1:].^(vP[3])./(1+vEps[1:].^(vP[3])));
//vU = vP[3].*(-1+2*aux.^(vP[3])./(1+aux.^(vP[3])));
mDeriv = vU.*mDeriv;
vScore[0][0:2] = sumc(mDeriv)/(cT-cMaxpq);
delete mDeriv;
vScore[0][3] = 1/vP[3] + (sumc(log(vEps[1:])))/(cT-cMaxpq) - (vP[4]+1)*(sumc(log(vEps[1:]).*vEps[1:].^(vP[3])./(1+vEps[1:].^(vP[3]))) )/(cT-cMaxpq);
vScore[0][4] = 1/vP[4] - sumc(log(vEps[1:].^(vP[3])+1))/(cT-cMaxpq);
return vScore ;

}

DySco::MEMBurrAnaCovMatrix(const mAnaHess, const vO, const cT, const cMaxpq)
{
decl a,b,c,d,e,f,g,sigma_e2,fact;
	f = gammafact(1+1/vO[3])*gammafact(vO[4]-1/vO[3])/gammafact(vO[4]);
	g = gammafact(1+2/vO[3])*gammafact(vO[4]-2/vO[3])/gammafact(vO[4]);
	sigma_e2 = g-f^(2);
		
	a = vO[1] - vO[2]*f;
	b = vO[1]^(2)-2*vO[1]*vO[2]*f + vO[2]^(2)*g;
	c = vO[1]*f -vO[2]*g;
	d = vO[2]*f/(1-vO[1]);
	e = vO[2]*sigma_e2/(1-vO[1]^(2))+d^(2);

	fact = vO[3]^(2)*vO[4]/((vO[4]+2)*(1-b));

	
//	println("a ",a,"b ",b,"c ",c,"d ",d,"e ",e,"f ",f);
//	println("sigma ",sigma_e2);
		if (b>=1){oxwarning("b>1"); println("b is:" ,b ~vO');continue;}
	//1) Dynamic parameters
	mAnaHess[0][0][0] = fact*(1-vO[1])^(2)*(1+a)/((1-a));
	mAnaHess[0][1][1] =	fact*(e+ 2*a*(1-vO[1]*a)^(-1)*(vO[2]^(2)*c*f*((1-vO[1])*(1-a))^(-1)+vO[1]*e+vO[2]*d*f));
	mAnaHess[0][2][2] =	fact*(g*(1-a)+2*c*f  )/(1-a);
	mAnaHess[0][0][1] =	fact*(1-vO[1])*(d+ vO[2]*a*f/((1-vO[1])*(1-a)) +a*(vO[1]*d+vO[2]*f+ vO[2]*c/(1-a))/(1-vO[1]*a)   );
	mAnaHess[0][1][0] =	mAnaHess[0][0][1]; 
	mAnaHess[0][0][2] =	fact*(1-vO[1])*(f+c)/(1-a);
	mAnaHess[0][2][0] =	mAnaHess[0][0][2];
	mAnaHess[0][1][2] =	fact*(d*f+ a*(vO[1]*d*f+vO[2]*g+vO[2]*c*f*(1-a)^(-1))*(1-vO[1]*a)^(-1) +vO[2]*c*f/((1-vO[1])*(1-a))  );
	mAnaHess[0][2][1] =	 mAnaHess[0][1][2];

	mAnaHess[0][3][3] =	(vO[3]^(2)*(2+vO[4]))^(-1)*(2+vO[4]+vO[4]*(polygamma(vO[4],1)+polygamma(1,1) + (polygamma(vO[4],0)-polygamma(1,0) + (1-vO[4])*vO[4]^(-1))^(2) - (1+vO[4]^(2))*vO[4]^(-2) ) );
 	mAnaHess[0][4][4] = vO[4]^(-2);
	mAnaHess[0][3][4] = -(vO[3]*(1+vO[4]))^(-1)*(polygamma(vO[4],0)-polygamma(1,0)-1 );
	mAnaHess[0][4][3] = mAnaHess[0][3][4];

	mAnaHess[0][0][3] = (1-vO[4]+vO[4]*(polygamma(vO[4],0)-polygamma(1,0)))*(1-vO[1])/((2+vO[4])*(1-a));
	mAnaHess[0][3][0] = mAnaHess[0][0][3];
	mAnaHess[0][0][4] = -vO[3]*(1-vO[1])/( (1+vO[4])*(1-a) );
	mAnaHess[0][4][0] = mAnaHess[0][0][4];

	mAnaHess[0][1][3] = (1-vO[4]+vO[4]*(polygamma(vO[4],0)-polygamma(1,0)))/((2+vO[4]))   *(vO[2]*f/((1-vO[1])*(1-a)));
	mAnaHess[0][1][4] =  -vO[3]*vO[2]*f/( (1+vO[4])*(1-vO[1])*(1-a));
	mAnaHess[0][4][1] = mAnaHess[0][1][4];
	mAnaHess[0][3][1] = mAnaHess[0][1][3];

	
	mAnaHess[0][2][3] = (1-vO[4]+vO[4]*(polygamma(vO[4],0)-polygamma(1,0)))/((2+vO[4]))*(f/(1-a));
	mAnaHess[0][2][4] = -vO[3]*f/((1+vO[4])*(1-a));
	mAnaHess[0][4][2] = mAnaHess[0][2][4];
	mAnaHess[0][3][2] = mAnaHess[0][2][3];

	
	mAnaHess[0] = -mAnaHess[0];
return mAnaHess[0];
}




////////////////////////////////////////////////////////////////////////////////
/////////////////// Dagum Distribution ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
DySco::DAGLambdaFilter(vLambda,vU, const vY, const vP, const vStruct)
{ decl cT=rows(vY);	 decl p = vStruct[0]; decl q = vStruct[1];decl maxpq = max(p, q);
decl vLevAux=new matrix[cT][1];
if (vStruct[3]==1) {vLevAux=GetLevAux(&vLevAux, m_mLev);}
decl u = 0; 
vLambda[0][0:maxpq-1] = vP[0];
decl amPar =  new array[7];
SplitPara(&amPar,vP,vStruct);

for(decl i=maxpq; i<cT;i++){
	vLambda[0][i] = amPar[0][0]*(1-sumc(amPar[1][0]))+amPar[1][0]'*reversec(vLambda[0][i-p:i-1])+amPar[2][0]'*reversec(vU[0][i-q:i-1]);
	if (vStruct[3]==1){vLambda[0][i] += amPar[4][0]*vLevAux[i-1][0];}
	vU[0][i] = -amPar[3][0][1]+ (amPar[3][0][1]+1)*vY[i]^(amPar[3][0][0])*exp(-amPar[3][0][0]*vLambda[0][i])/(1+vY[i]^(amPar[3][0][0])*exp(-amPar[3][0][0]*vLambda[0][i]));
	// Dagum type II
//		vU[0][i] = 1 - (amPar[3][0][1]+1)*vY[i]^(-amPar[3][0][0])*exp(amPar[3][0][0]*vLambda[0][i])/(1+vY[i]^(-amPar[3][0][0])*exp(amPar[3][0][0]*vLambda[0][i]));	
	}
//return vLambda[0]~vU;
}
 
DySco::DAGLik(const vP, const adFunc, const avScore, const amHessian)
{decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
decl vU = new matrix[m_iT2sel-m_iT1sel+1][1];

if (m_iTwoComp==0) DAGLambdaFilter(&vLambda,&vU,m_vY, vP, m_vStruct);
else  LambdaFilter2comp(&vLambda, &vU, m_vY, vP, m_vStruct, DAG);
//	   println("meanc(vU)",meanc(vU));
//	   println("vP : ", vP');  
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
adFunc[0] = double(log(amPar[3][0][0])-log(Beta(amPar[3][0][1],1))+ (m_iT2sel-m_iT1sel)^(-1)*(amPar[3][0][0]*amPar[3][0][1]-1)*sumc(log(m_vY[1:]))-(m_iT2sel-m_iT1sel)^(-1)*amPar[3][0][0]*amPar[3][0][1]*sumc(vLambda[1:]) - (m_iT2sel-m_iT1sel)^(-1)*(amPar[3][0][1]+1)*sumc(log((m_vY[1:]./exp(vLambda[1:])).^(amPar[3][0][0])+1))  );
// Dagum type II
//adFunc[0] = double(log(amPar[3][0][0])+(log(1-amPar[3][0][0]))+ log(amPar[3][0][1]) - (m_iT2sel-m_iT1sel)^(-1)*(amPar[3][0][0]+1)*sumc(log(m_vY[1:]))+(m_iT2sel-m_iT1sel)^(-1)*amPar[3][0][0]*sumc(vLambda[1:]) + (m_iT2sel-m_iT1sel)^(-1)*(-amPar[3][0][1]-1)*sumc(log((m_vY[1:]./exp(vLambda[1:])).^(-amPar[3][0][0])+1))  );

if (isarray(avScore)){  // if analytical score should be used in maximisation
//	(avScore[0])[0] = vScore[0];
//	(avScore[0])[1] = vScore[1];
//	(avScore[0])[2] = vScore[2];
//	(avScore[0])[3] = vScore[3];
//	(avScore[0])[4] = vScore[4];
}
if (amHessian){
//    (amHessian[0])= mAnaHess;
}
return 1;
}


DySco::DAGLikVal(const vP, const cStart)
{ decl cT =rows(m_mData[][0]);	decl cTIn = cT-cStart;
decl vLambda = 		new matrix[cT][1];
if(m_iModelClass == SCORE){ 
decl vU = 			new matrix[cT][1];
if (m_iTwoComp==0) 	DAGLambdaFilter(&vLambda,&vU, m_mData[][0], vP, m_vStruct);
else  				LambdaFilter2comp(&vLambda, &vU, m_mData[][0], vP, m_vStruct, DAG);
}
else if (m_iModelClass == MEM){
	MEMLambdaFilter(&vLambda, m_mData[][0], vP, m_vStruct);
}
//DrawTMatrix(0,vLambda', {"Lambda"});
//ShowDrawWindow();
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
decl vLikContr = new matrix[cTIn][1];
 for(decl i=1;i<cTIn;i++){
vLikContr[i] = log(amPar[3][0][0]) + log(amPar[3][0][1]) + (amPar[3][0][0]*amPar[3][0][1]-1)*log(m_mData[cStart+i][0]) - amPar[3][0][1]*amPar[3][0][0]*vLambda[cStart+i] - (amPar[3][0][1]+1)*(log((m_mData[cStart+i][0]./exp(vLambda[cStart+i])).^(amPar[3][0][0])+1));

}
decl LikVal = sumc(vLikContr);
return LikVal;
}

//DySco::DAGScore(vScore, const vEps,const vLambda, const vP, const cMaxpq)
//{decl cT = rows(vEps);
//if (rows(vLambda)!=cT) {oxwarning("length of eps and lambda different. Pls correct"); exit(0);}
//decl mDeriv = new matrix[cT][3];
//mDeriv[0][0] = (1-vP[1]); ///(vP[1]-vP[2]*vP[3]);
//mDeriv[0][1] = mDeriv[0][2] = 0;
//mDeriv[1:][0] = constant(1-vP[1],cT-1,1);
//mDeriv[1:][1] = -ones(cT-1,1)*vP[0]+vLambda[:cT-2];
//mDeriv[1:][2] =	(1+vP[4])*(vEps[:cT-2].^(vP[3])./(1+vEps[:cT-2].^(vP[3])))-1;
//decl x;
//for (decl i=1; i<cT; i++){
//x = vP[1]-vP[2]*vP[3]*(1+vP[4])*vEps[i-1]^(vP[3])*(1+vEps[i-1]^(vP[3]))^(-2);
//mDeriv[i][] += x*mDeriv[i-1][];
//}
//mDeriv =vP[3]*((1+vP[4])*vEps[1:].^(vP[3])./(1+vEps[1:].^(vP[3]))-1).*mDeriv[1:][];
//vScore[0][0:2] = sumc(mDeriv)/(cT-cMaxpq);
//delete mDeriv;
//vScore[0][3] = 1/vP[3] + (sumc(log(vEps[1:])))/(cT-cMaxpq) - (vP[4]+1)*(sumc(log(vEps[1:]).*vEps[1:].^(vP[3])./(1+vEps[1:].^(vP[3]))) )/(cT-cMaxpq);
//vScore[0][4] = 1/vP[4] - sumc(log(vEps[1:].^(vP[3])+1))/(cT-cMaxpq);
//}
//



DySco::DAGAnaCovMatrix(mAnaHess, const vO, const cT, const cMaxpq)
{decl a,b,c,sigma_u_2;
	a = vO[1]-vO[2]*vO[3]*vO[4]*(vO[4]+1)*gammafact(vO[4]+1)/gammafact(vO[4]+3);
	b = vO[1]^(2)-2*vO[1]*vO[2]*vO[3]*vO[4]*(vO[4]+1)*gammafact(vO[4]+1)/gammafact(vO[4]+3)+2*vO[2]^(2)*vO[3]^(2)*(vO[4]+1)^(3)*vO[4]*gammafact(vO[1]+1)/gammafact(vO[4]+5);
	c = -vO[2]*vO[3]*vO[4]*(vO[4]+1)*gammafact(vO[4]+1)/gammafact(vO[4]+3);

	sigma_u_2 = vO[4]/(2+vO[4]);
	decl fact = vO[3]^(2)*vO[4]/((vO[4]+2)*(1-b)); // factor in front of dynamic elements of hessian
	
	decl aux = new matrix[5][5]; // omega notation, Andrews version with nu in score (taken from Andres&Harvey*(2012))
	// 1) dynamic parameters
	mAnaHess[0][0][0] =  fact*(1-vO[1])^(2)*(1+a)/(1-a) ;
	mAnaHess[0][1][1] =	 fact*vO[2]^(2)*sigma_u_2*(1+a*vO[1])/((1-vO[1]^(2))*(1-a*vO[1]));
	mAnaHess[0][2][2] =  fact*sigma_u_2;
	mAnaHess[0][0][1] =  fact*a*c*vO[2]*(1-vO[1])/((1-a)*(1-a*vO[1]));
	mAnaHess[0][1][0] =  mAnaHess[0][0][1];
	mAnaHess[0][0][2] =  fact*c*(1-vO[1])/(1-a);
	mAnaHess[0][2][0] =  mAnaHess[0][0][2];
	mAnaHess[0][1][2] =  fact*a*vO[2]*sigma_u_2/(1-a*vO[1]);
	mAnaHess[0][2][1] =  mAnaHess[0][1][2];

	mAnaHess[0][3][3] =  vO[3]^(-2)*(2+vO[4])^(-1)*(2+vO[4]+vO[4]*(polygamma(vO[4],0)+polygamma(1,0))+(gammafact(1)-gammafact(vO[4])+(vO[4]+1)*vO[4]^(-1))^(2)-(vO[4]^(2)+1)*vO[4]^(-2));
	mAnaHess[0][4][4] = polygamma(vO[4],0)-polygamma(vO[4]+1,0);

	mAnaHess[0][3][4] = mAnaHess[0][4][3] = - (vO[3]*(vO[4]+1))^(-1)*(gammafact(vO[4])-gammafact(1)-1);
	
	mAnaHess[0][0][3] = mAnaHess[0][3][0] = (2+vO[4])^(-1)*(vO[4]-1 -vO[4]*(gammafact(vO[4])-gammafact(1)) );
	mAnaHess[0][1][3] = mAnaHess[0][3][1] = mAnaHess[0][2][3] = mAnaHess[0][3][2] = 0  ;

	mAnaHess[0][0][4] = mAnaHess[0][4][0] = vO[3]/(1+vO[4]);
	mAnaHess[0][1][4] = mAnaHess[0][4][1] =	mAnaHess[0][2][4] = mAnaHess[0][4][2] =	0;

mAnaHess[0] = -mAnaHess[0];
return mAnaHess[0];
}
////////////////////////// DAG MEM  STUFF ////////////////////////////////////////////////////////////
DySco::DAGMEMLik(const vP, const adFunc, const avScore, const amHessian)
{ decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) MEMLambdaFilter(&vLambda, m_vY, vP, m_vStruct);
else 	MEMLambdaFilter2comp(&vLambda, m_vY, vP, m_vStruct);
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
adFunc[0] = double(log(amPar[3][0][1])-log(Beta(amPar[3][0][1],1))+ (m_iT2sel-m_iT1sel)^(-1)*(amPar[3][0][0]*amPar[3][0][1]-1)*sumc(log(m_vY[1:]))-(m_iT2sel-m_iT1sel)^(-1)*amPar[3][0][0]*amPar[3][0][1]*sumc(vLambda[1:]) - (m_iT2sel-m_iT1sel)^(-1)*(amPar[3][0][1]+1)*sumc(log((m_vY[1:]./exp(vLambda[1:])).^(amPar[3][0][0])+1))  );


if (isarray(avScore)){  // if analytical score should be used in maximisation
//	decl vScore = new matrix[5][1];
//	decl aux = m_vY./exp(vLambda);
//	DAGMEMScores(&vScore,m_vY,vLambda, vP, max(m_vStruct[1], m_vStruct[2]));
////	println(vScore');
//	(avScore[0])[0] = vScore[0];
//	(avScore[0])[1] = vScore[1];
//	(avScore[0])[2] = vScore[2];
//	(avScore[0])[3] = vScore[3];
//	(avScore[0])[4] = vScore[4];
}
return 1;
}

//DySco::DAGMEMScores(vScore,const vY,const vLambda, const vP, const cMaxpq)
//{decl cT = rows(vY); decl mDeriv = new matrix[cT][3];
//
//}



/////////////////////////////////////////////////////////////////
/////////// F Distribution //////////////////////////////////////
/////////////////////////////////////////////////////////////////
DySco::FLambdaFilter(vLambda,vU, const vY, const vP, const vStruct)
{ decl cT=rows(vY);	 decl p = vStruct[0]; decl q = vStruct[1];decl maxpq = max(p, q);
decl vLevAux=new matrix[cT][1];
if (vStruct[3]==1) {vLevAux=GetLevAux(&vLevAux, m_mLev);}
decl u = 0,b,eps; 
vLambda[0][0:maxpq-1] = vP[0];
decl amPar =  new array[7];
SplitPara(&amPar,vP,vStruct);
for(decl i=maxpq; i<cT;i++){
	vLambda[0][i] = amPar[0][0]*(1-sumc(amPar[1][0]))+amPar[1][0]'*reversec(vLambda[0][i-p:i-1])+amPar[2][0]'*reversec(vU[0][i-q:i-1]);
	if (vStruct[3]==1){vLambda[0][i] += amPar[4][0]*vLevAux[i-1][0];}
	eps = vY[i-1]*exp(-vLambda[0][i-1]);
	b = amPar[3][0][0]*eps/(amPar[3][0][0]*eps + amPar[3][0][1]);
	vU[0][i] = -0.5*amPar[3][0][0]+ 0.5*(amPar[3][0][0]+amPar[3][0][1])*b;
	}
//	println( amPar[3][0][0]~amPar[3][0][1]);
//DrawTMatrix(0,vU[0]',{"score"});
//DrawTMatrix(1,vLambda[0]',{"lambda"});
//ShowDrawWindow();
//return vLambda[0]~vU;
}


DySco::FLik(const vP, const adFunc, const avScore, const amHessian)
{decl cT = m_iT2sel-m_iT1sel+1;
decl vLambda = new matrix[cT][1];
decl vU = new matrix[cT][1];
if (m_iTwoComp==0) FLambdaFilter(&vLambda,&vU,m_vY, vP, m_vStruct);
else  LambdaFilter2comp(&vLambda, &vU, m_vY, vP, m_vStruct, F);
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
adFunc[0] = 0.5*amPar[3][0][0]*log(amPar[3][0][0]) +0.5*amPar[3][0][1]*log(amPar[3][0][1])-	 0.5*amPar[3][0][0]*meanc(vLambda[1:]) + (0.5*amPar[3][0][0]-1)*meanc(log(m_vY[1:]))-0.5*(amPar[3][0][0]+amPar[3][0][1])*meanc(log(amPar[3][0][0]*m_vY[1:]./exp(vLambda[1:]) + amPar[3][0][1]) )-log(Beta(0.5*amPar[3][0][0],0.5*amPar[3][0][1]));

if (isarray(avScore)){  // if analytical score should be used in maximisation
	decl vScore = new matrix[5][1];
	decl eps = m_vY./exp(vLambda);
	FScore(&vScore, eps,vLambda, vP,max(m_vStruct[0], m_vStruct[1]));
	(avScore[0])[0] = vScore[0];
	(avScore[0])[1] = vScore[1];
	(avScore[0])[2] = vScore[2];
	(avScore[0])[3] = vScore[3];
	(avScore[0])[4] = vScore[4];
}
if (amHessian){
	decl mAnaHess = new matrix[5][5];
	FAnaCovMatrix(&mAnaHess,vP, (m_iT2sel-m_iT1sel+1), max(m_vStruct[0], m_vStruct[1]));
    (amHessian[0])= mAnaHess;
}
return 1;
}

DySco::FLikVal(const vP, const cStart)
{ decl cT =rows(m_mData[][0]);	decl cTIn = cT-cStart;
decl vLambda = 		new matrix[cT][1];
decl vU = 			new matrix[cT][1];
decl p = m_vStruct[0]; decl q = m_vStruct[1]; decl maxpq = max(p, q);
if (m_iTwoComp==0) 	FLambdaFilter(&vLambda,&vU, m_mData[][0], vP, m_vStruct);
else  				LambdaFilter2comp(&vLambda, &vU, m_mData[][0], vP, m_vStruct, F);
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
decl vLikContr = new matrix[cTIn][1];
 for(decl i=1;i<cTIn;i++){
vLikContr[i] = 0.5*amPar[3][0][0]*log(amPar[3][0][0]) +0.5*amPar[3][0][1]*log(amPar[3][0][1])-	 0.5*amPar[3][0][0]*vLambda[cStart+i] + (0.5*amPar[3][0][0]-1)*log(m_mData[cStart+i][0])-0.5*(amPar[3][0][0]+amPar[3][0][1])*log(amPar[3][0][0]*m_mData[cStart+i][0]./exp(vLambda[cStart+i]) + amPar[3][0][1]) -log(Beta(0.5*amPar[3][0][0],0.5*amPar[3][0][1]));
}
decl LikVal = sumc(vLikContr);
return LikVal;
}


DySco::FScore(vScore, const eps, const vLambda, const vP, const cMaxpq)
{ decl cT=rows(eps);  if (rows(vLambda)!=cT) {oxwarning("length of eps and lambda different. Pls correct"); exit(0);}
decl mDeriv = new matrix[cT][3];
mDeriv[0][0] = (1-vP[1]); ///(vP[1]-vP[2]*vP[3]);
mDeriv[0][1] = mDeriv[0][2] = 0;
mDeriv[1:][0] = constant(1-vP[1],cT-1,1);
mDeriv[1:][1] = -ones(cT-1,1)*vP[0]+vLambda[:cT-2];
mDeriv[1:][2] = -0.5*vP[3]*ones(cT-1,1)+ 0.5*(vP[3]+vP[4])*vP[3]*eps[:cT-2]./(vP[4]+vP[3]*eps[:cT-2]);
decl x;
for (decl i=1; i<cT; i++){
x = vP[1]-0.5*vP[2]*(vP[3]+vP[4])*vP[3]*vP[4]*eps[i-1]*(vP[4]+vP[3]*eps[i-1])^(-2);
mDeriv[i][] += x*mDeriv[i-1][];
}
mDeriv = (-0.5*vP[3]*ones(cT-1,1)+ 0.5*(vP[3]+vP[4])*vP[3]*eps[1:]./(vP[4]+vP[3]*eps[1:])).*mDeriv[1:][];
vScore[0][0:2] = sumc(mDeriv)/(cT-cMaxpq);
delete mDeriv;
//decl aux2 = eps[1:].^(vP[3]).*log(eps[1:])./(1+eps[1:].^(vP[3]));
vScore[0][3] =0.5+0.5*log(vP[3])+0.5*cT^(-1)*sumc(log(eps[1:]))-0.5*cT^(-1)*sumc(log(vP[3].*eps[1:] +vP[4] ))-0.5*cT^(-1)*(vP[3]+vP[4])*sumc(eps[1:]./(vP[3].*eps[1:]+vP[4]))-0.5*(polygamma(0.5*vP[3],0)-polygamma(0.5*vP[3]+0.5*vP[4],0) ) ;
vScore[0][4] =0.5+0.5*log(vP[4])-0.5*cT^(-1)*sumc(log(vP[3]*eps[1:]+vP[4])) - 0.5*cT^(-1)*(vP[3]+vP[4])*sumc((vP[3].*eps[1:]+vP[4]).^(-1))	-0.5*(polygamma(0.5*vP[4],0)-polygamma(0.5*vP[3]+0.5*vP[4],0));
return vScore ;
}


DySco::FAnaCovMatrix(const mAnaHess, const vO, const cT, const cMaxpq)
{
decl a,b,c, sigma_u_2;
	a = vO[1]-0.25*vO[2]*vO[3]*vO[4]/(1+0.5*vO[3]+0.5*vO[4]);
	//b = vO[1]^(2)-2*vO[1]*vO[2]*vO[3]*vO[4]/(2+vO[4]) + 2*vO[2]^(2)*vO[3]^(2)*vO[4]*(vO[4]+1)^(2)/((vO[4]+2)*(vO[4]+3)*(vO[4]+4));
	b = vO[1]^(2)-0.5*vO[1]*vO[2]*vO[3]*vO[4]/(1+0.5*vO[3]+0.5*vO[4]) +	0.25*vO[2]^(2)*(vO[3]+vO[4])^(2)*Beta(0.5*vO[3]+2,0.5*vO[4]+2)/Beta(0.5*vO[3],0.5*vO[4]);

//	c = vO[2]*vO[3]*vO[4]*
	c =  0.125*vO[2]*vO[3]*vO[4]*(vO[3]-(vO[3]+vO[4])*(0.5*vO[3]+1)/(2+0.5*vO[3]+0.5*vO[4]))/(1+0.5*vO[3]+0.5*vO[4]);
   
	sigma_u_2 = vO[3]*vO[4]/(4*(1+0.5*vO[3]+0.5*vO[4]));
	decl fact = sigma_u_2/(1-b); // factor in front of dynamic elements of hessian
	
	decl aux = new matrix[5][5]; // omega notation, Andrews version with nu in score (taken from Andres&Harvey*(2012))
	// 1) dynamic parameters
	mAnaHess[0][0][0] =  fact*(1-vO[1])^(2)*(1+a)/(1-a) ;
	mAnaHess[0][1][1] =	 fact*vO[2]^(2)*sigma_u_2*(1+a*vO[1])/((1-vO[1]^(2))*(1-a*vO[1]));
	mAnaHess[0][2][2] =  fact*sigma_u_2;
	mAnaHess[0][0][1] =  fact*a*c*vO[2]*(1-vO[1])/((1-a)*(1-a*vO[1]));
	mAnaHess[0][1][0] =  mAnaHess[0][0][1];
	mAnaHess[0][0][2] =  fact*c*(1-vO[1])/(1-a);
	mAnaHess[0][2][0] =  mAnaHess[0][0][2];
	mAnaHess[0][1][2] =  fact*a*vO[2]*sigma_u_2/(1-a*vO[1]);
	mAnaHess[0][2][1] =  mAnaHess[0][1][2];
												
	mAnaHess[0][3][3] = -(0.5/vO[3]-1/(vO[3]+vO[4])+0.5*(0.5*vO[3]+1)/(vO[3]*(1+0.5*vO[3]+0.5*vO[4]))-0.25*(polygamma(0.5*vO[3],1) -polygamma(0.5*vO[3]+0.5*vO[4],1))  );
	mAnaHess[0][4][4] = -(-0.5/vO[4]+vO[3]/(vO[4]*(vO[3]+vO[4]))+0.5*(0.5*vO[4]+1)/(vO[4]*(1+0.5*vO[3]+0.5*vO[4]))-0.25*(polygamma(0.5*vO[4],1) -polygamma(0.5*vO[3]+0.5*vO[4],1))  );
//	mAnaHess[0][3][4] = mAnaHess[0][4][3] = -(-0.5 +1/(1+0.5*vO[3]+0.5*vO[4])+0.25*polygamma(0.5*vO[3]+0.5*vO[4],1)   );
	//mAnaHess[0][3][3] = -(0.5/vO[3]-1/(vO[3]+vO[4])+0.25/(1+0.5*vO[3]+0.5*vO[4])-0.25*(polygamma(0.5*vO[3],1)-0.25*polygamma(0.5*vO[3]+0.5*vO[4],1)));
	//mAnaHess[0][4][4] =	-(0.5/vO[4]-1/(vO[3]+vO[4])-0.5*(0.5*vO[4]+1)/(vO[4]*(1+0.5*vO[3]+0.5*vO[4]))-0.25*(polygamma(0.5*vO[4],1)-polygamma(0.5*vO[3]+0.5*vO[4],1)));

	
	mAnaHess[0][3][4] =mAnaHess[0][4][3] = -(-1/(vO[3]+vO[4])+0.25/(1+0.5*vO[3]+0.5*vO[4]) +0.25*polygamma(0.5*vO[3]+0.5*vO[4],1));
	//(1+vO[4]*(1-vO[4])+vO[4]*(1+vO[4])*(M_EULER+polygamma(vO[4]+2,0)))*(vO[3]*vO[4]*(1+vO[4])^(2))^(-1) ;

	mAnaHess[0][0][3] = mAnaHess[0][3][0] = -(0.5*vO[3]/(vO[3]+vO[4])+0.25*vO[4]/(1+0.5*vO[3]+0.5*vO[4])-0.5)*(1-vO[1])*(1-a)^(-1) ;
	mAnaHess[0][1][3] = mAnaHess[0][3][1] = mAnaHess[0][2][3] = mAnaHess[0][3][2] = 0  ;		   


	mAnaHess[0][0][4] = mAnaHess[0][4][0] =-(0.5*vO[3]/(vO[3]+vO[4])-0.25*vO[3]/(1+0.5*vO[3]+0.5*vO[4]))*(1-vO[1])*(1-a)^(-1); 
	mAnaHess[0][1][4] = mAnaHess[0][4][1] =	mAnaHess[0][2][4] = mAnaHess[0][4][2] =	0;

	mAnaHess[0] = -mAnaHess[0];
return mAnaHess;
}

///////////////// F MEM Stuff /////////////////////
DySco::FMEMLik(const vP, const adFunc, const avScore, const amHessian)
{decl cT = m_iT2sel-m_iT1sel+1;
decl vLambda = new matrix[cT][1];
if (m_iTwoComp==0) MEMLambdaFilter(&vLambda, m_vY, vP, m_vStruct);
else 	MEMLambdaFilter2comp(&vLambda, m_vY, vP, m_vStruct);
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
adFunc[0] = 0.5*amPar[3][0][0]*log(amPar[3][0][1]) +0.5*amPar[3][0][1]*log(amPar[3][0][1])-	 0.5*amPar[3][0][0]*cT^(-1)*sumc(vLambda[1:]) + (0.5*amPar[3][0][0]-1)*cT^(-1)*sumc(log(m_vY[1:]))-0.5*(amPar[3][0][0]+amPar[3][0][1])*cT^(-1)*sumc(log(amPar[3][0][1]*m_vY[1:]./exp(vLambda[1:]) + amPar[3][0][1]) )-log(Beta(0.5*amPar[3][0][0],0.5*amPar[3][0][1]));
if (isarray(avScore)){  // if analytical score should be used in maximisation
	decl vScore = new matrix[5][1];
	decl aux =  m_vY./exp(vLambda);
	FMEMScores(&vScore,m_vY, vLambda, vP, max(m_vStruct[0], m_vStruct[1]));
//	println(vScore');
	(avScore[0])[0] = vScore[0];
	(avScore[0])[1] = vScore[1];
	(avScore[0])[2] = vScore[2];
	(avScore[0])[3] = vScore[3];
	(avScore[0])[4] = vScore[4];
}
return 1;
}

DySco::FMEMScores(vScore,const vY,const vLambda, const vP, const cMaxpq)
{decl cT = m_iT2sel-m_iT1sel+1;	
decl mDeriv = new matrix[cT][3];
mDeriv[0][0] = 1;
mDeriv[0][1] = meanc(-vP[0]+vLambda)/(1-vP[1]);
mDeriv[0][2] = meanc(vY)/(1-vP[1]);
mDeriv[1:][0] = constant(1-vP[1],cT-1,1);
mDeriv[1:][1] = -vP[0] + vLambda[:cT-2];
mDeriv[1:][2] = vY[:cT-2];
for (decl i=cMaxpq; i<cT; i++){
mDeriv[i][1:] += vP[1]*mDeriv[i-1][1:];
}
mDeriv = dropr(mDeriv,0);
decl aux = vY./exp(vLambda);
decl vU = -0.5*vP[3] + 0.5*(vP[3]+vP[4])*vP[3].*aux[1:]./(vP[4]+vP[3].*aux[1:]);
//vU = vP[3].*(-1+2*aux.^(vP[3])./(1+aux.^(vP[3])));
mDeriv = vU.*mDeriv;
vScore[0][0:2] = sumc(mDeriv)/(cT-cMaxpq);
delete mDeriv;
vScore[0][3] =0.5+0.5*log(vP[3])+0.5*m_cT^(-1)*sumc(log(aux[1:]))-0.5*m_cT^(-1)*sumc(log(vP[3].*aux[1:] +vP[4] ))-0.5*m_cT^(-1)*(vP[3]+vP[4])*sumc(aux[1:]./(vP[3].*aux[1:]+vP[4]))-0.5*(polygamma(0.5*vP[3],0)-polygamma(0.5*vP[3]+0.5*vP[4],0) ) ;
vScore[0][4] =0.5+0.5*log(vP[4])-0.5*m_cT^(-1)*sumc(log(vP[3]*aux[1:]+vP[4])) - 0.5*m_cT^(-1)*(vP[3]+vP[4])*sumc((vP[3].*aux[1:]+vP[4]).^(-1))	-0.5*(polygamma(0.5*vP[4],0)-polygamma(0.5*vP[3]+0.5*vP[4],0));
return vScore ;

}

DySco::MEMFAnaCovMatrix(const mAnaHess, const vO, const cMaxpq)
{decl a,b,c,d,e,f,g,sigma_e2,fact;
	f = vO[4]/(vO[4]-2);
	g = gammafact(0.5*vO[3]+2)*gammafact(0.5*vO[4]-2)*(gammafact(0.5*vO[3])*gammafact(0.5*vO[4]))^(-1);
	sigma_e2 = g-f^(2);
		
	a = vO[1] - vO[2]*f;
	b = vO[1]^(2)-2*vO[1]*vO[2]*f + vO[2]^(2)*g;
	c = vO[1]*f -vO[2]*g;
	d = vO[2]*f/(1-vO[1]);
	e = vO[2]*sigma_e2/(1-vO[1]^(2))+d^(2);

	decl sigma_u_2 = vO[3]*vO[4]/(4*(1+0.5*vO[3]+0.5*vO[4]));
	fact = sigma_u_2/(1-b);

	
	println("a ",a,"b ",b,"c ",c,"d ",d,"e ",e,"f ",f);
	println("sigma ",sigma_e2);
		if (b>=1){oxwarning("b>1"); println("b is:" ,b ~vO');continue;}
	//1) Dynamic parameters
	mAnaHess[0][0][0] = fact*(1-vO[1])^(2)*(1+a)/((1-a));
	mAnaHess[0][1][1] =	fact*(e+ 2*a*(1-vO[1]*a)^(-1)*(vO[2]^(2)*c*f*((1-vO[1])*(1-a))^(-1)+vO[1]*e+vO[2]*d*f));
	mAnaHess[0][2][2] =	fact*(g*(1-a)+2*c*f  )/(1-a);
	mAnaHess[0][0][1] =	fact*(1-vO[1])*(d+ vO[2]*a*f/((1-vO[1])*(1-a)) +a*(vO[1]*d+vO[2]*f+ vO[2]*c/(1-a))/(1-vO[1]*a)   );
	mAnaHess[0][1][0] =	mAnaHess[0][0][1]; 
	mAnaHess[0][0][2] =	fact*(1-vO[1])*(f+c)/(1-a);
	mAnaHess[0][2][0] =	mAnaHess[0][0][2];
	mAnaHess[0][1][2] =	fact*(d*f+ a*(vO[1]*d*f+vO[2]*g+vO[2]*c*f*(1-a)^(-1))*(1-vO[1]*a)^(-1) +vO[2]*c*f/((1-vO[1])*(1-a))  );
	mAnaHess[0][2][1] =	 mAnaHess[0][1][2];

												
	mAnaHess[0][3][3] = -(0.5/vO[3]-1/(vO[3]+vO[4])+0.5*(0.5*vO[3]+1)/(vO[3]*(1+0.5*vO[3]+0.5*vO[4]))-0.25*(polygamma(0.5*vO[3],1) -polygamma(0.5*vO[3]+0.5*vO[4],1))  );
	mAnaHess[0][4][4] = -(-0.5/vO[4]+vO[3]/(vO[4]*(vO[3]+vO[4]))+0.5*(0.5*vO[4]+1)/(vO[4]*(1+0.5*vO[3]+0.5*vO[4]))-0.25*(polygamma(0.5*vO[4],1) -polygamma(0.5*vO[3]+0.5*vO[4],1))  );
	mAnaHess[0][3][4] =mAnaHess[0][4][3] = -(-1/(vO[3]+vO[4])+0.25/(1+0.5*vO[3]+0.5*vO[4]) +0.25*polygamma(0.5*vO[3]+0.5*vO[4],1));
//	//(1+vO[4]*(1-vO[4])+vO[4]*(1+vO[4])*(M_EULER+polygamma(vO[4]+2,0)))*(vO[3]*vO[4]*(1+vO[4])^(2))^(-1) ;
//
	mAnaHess[0][0][3] = mAnaHess[0][3][0] = -(0.5*vO[3]/(vO[3]+vO[4])+0.25*vO[4]/(1+0.5*vO[3]+0.5*vO[4])-0.5)*(1-vO[1])*(1-a)^(-1) ;
	mAnaHess[0][0][4] = mAnaHess[0][4][0] =-(0.5*vO[3]/(vO[3]+vO[4])-0.25*vO[3]/(1+0.5*vO[3]+0.5*vO[4]))*(1-vO[1])*(1-a)^(-1); 

	mAnaHess[0][1][3] = mAnaHess[0][3][1] = -(0.5*vO[3]/(vO[3]+vO[4])+0.25*vO[4]/(1+0.5*vO[3]+0.5*vO[4])-0.5)*vO[2]*f*((1-vO[1])*(1-a))^(-1);
	mAnaHess[0][1][4] = mAnaHess[0][4][1] =-(0.5*vO[3]/(vO[3]+vO[4])-0.25*vO[3]/(1+0.5*vO[3]+0.5*vO[4]))*vO[2]*f*((1-vO[1])*(1-a))^(-1);

	mAnaHess[0][2][3] = mAnaHess[0][3][2] = -(0.5*vO[3]/(vO[3]+vO[4])+0.25*vO[4]/(1+0.5*vO[3]+0.5*vO[4])-0.5)*f*(1-a)^(-1);
	mAnaHess[0][2][4] = mAnaHess[0][4][2] =-(0.5*vO[3]/(vO[3]+vO[4])-0.25*vO[3]/(1+0.5*vO[3]+0.5*vO[4]))*f*(1-a)^(-1); 


	mAnaHess[0] = -mAnaHess[0];
return mAnaHess;
}


/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// GB2 Distribution ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
DySco::GB2LambdaFilter(const vLambda, const vY, const vP, const cT)
{decl u = 0; decl eps;
	vLambda[0][0] = vP[0];				 
	for(decl i=1; i<cT;i++){
		eps = vY[i-1][0]/exp(vLambda[0][i-1]);
		u = -vP[4]+(vP[4]+vP[5])*eps^(vP[3])/(1+eps^(vP[3]));
		vLambda[0][i] = vP[0]*(1-vP[1])+vP[1]*vLambda[0][i-1]+vP[2]*u;
	 }
}

DySco::GB2Lik(const vP, const adFunc, const avScore, const amHessian) 
{
decl vLambda = new matrix[m_cT][1];
GB2LambdaFilter(&vLambda, m_vY, vP, m_cT) ;
adFunc[0] = double(log(vP[3])-log(Beta(vP[4],vP[5])) + (m_cT-1)^(-1)*(vP[3]*vP[4]-1)*sumc(log(m_vY[1:]))-(m_cT-1)^(-1)*vP[3]*vP[4]*sumc(vLambda[1:]) - (m_cT-1)^(-1)*(vP[4]+vP[5])*sumc(log(1+(m_vY[1:]./exp(vLambda[1:])).^(vP[3])  ) )     );
if (isarray(avScore)){  // if analytical score should be used in maximisation
	decl vScore = new matrix[6][1];
	decl eps = m_vY./exp(vLambda);
	GB2Scores(&vScore, eps,vLambda, vP, m_cT,1);
//    decl vScore = Burr2(vP, m_cT, 1);
//println(vScore');
	vScore = vScore/10;
	(avScore[0])[0] = vScore[0];
	(avScore[0])[1] = vScore[1];
	(avScore[0])[2] = vScore[2];
	(avScore[0])[3] = vScore[3];
	(avScore[0])[4] = vScore[4];
	(avScore[0])[5] = vScore[5];	
}
if (amHessian){
	decl mAnaHess = new matrix[6][6];
	GB2AnaCovMatrix(&mAnaHess,vP);
    (amHessian[0])= mAnaHess;
}

return 1;
}

DySco::GB2LikVal(const vP, const cStart)
{ decl cT =rows(m_mData[][0]);	decl cTIn = cT-cStart;
decl vLambda = 		new matrix[cT][1];
decl vU = 			new matrix[cT][1];
decl p = m_vStruct[0]; decl q = m_vStruct[1]; decl maxpq = max(p, q);
if (m_iTwoComp==0)	LambdaFilter(&vLambda, &vU,m_mData[][0],vP,m_vStruct,GB2);
else  				LambdaFilter2comp(&vLambda, &vU, m_mData[][0], vP, m_vStruct, GB2);
decl amPar =  new array[7]; SplitPara(&amPar,vP,m_vStruct);
decl vLikContr = new matrix[cTIn][1];
 for(decl i=1;i<cTIn;i++){
vLikContr[i] = log(amPar[3][0][0])-log(Beta(amPar[3][0][1],amPar[3][0][2])) + (amPar[3][0][0]*amPar[3][0][1]-1)*log(m_mData[cStart+i][0])-amPar[3][0][0]*amPar[3][0][1]*vLambda[cStart+i] - (amPar[3][0][1]+amPar[3][0][2])*log(1+(m_mData[cStart+i][0]./exp(vLambda[cStart+i])).^(amPar[3][0][0]) );     
}
decl LikVal = sumc(vLikContr);
return LikVal;
}



DySco::GB2Scores(const vScore, const eps,const vLambda, const vP, const cT, const cMaxpq)
{ 
decl mDeriv = new matrix[cT][3];
mDeriv[0][0] = (1-vP[1]);///(1-(vP[1]-vP[2]*vP[3]*vP[4]/(2+vP[4]))); ///(vP[1]-vP[2]*vP[3]);
mDeriv[0][1] = mDeriv[0][2] = 0;
mDeriv[1:][0] = constant(1-vP[1],cT-1,1);
mDeriv[1:][1] = -ones(cT-1,1)*vP[0]+vLambda[:cT-2];
mDeriv[1:][2] =	(vP[4]+vP[5])*eps[:cT-2].^(vP[3])./(1+eps[:cT-2].^(vP[3]))-vP[4];
decl x;
for (decl i=1; i<cT; i++){
x = vP[1]-vP[2]*vP[3]*(vP[4]+vP[5])*eps[i-1]^(vP[3])*(1+eps[i-1]^(vP[3]))^(-2);
mDeriv[i][] += x*mDeriv[i-1][];
}
mDeriv =vP[3]*((vP[4]+vP[5])*eps[1:].^(vP[3])./(1+eps[1:].^(vP[3]))-vP[4]).*mDeriv[1:][];
vScore[0][0:2] = sumc(mDeriv)/(cT-cMaxpq);
delete mDeriv;
decl vB = eps.^(vP[3])./(1+eps.^(vP[3]));
vScore[0][3] = 1/vP[3] + vP[4]*meanc(log(eps[1:])) - (vP[4]+vP[5])*meanc(log(eps[1:]).*vB[1:]); 
vScore[0][4] = polygamma(vP[4]+vP[5],0)-polygamma(vP[4],0)+ vP[3]*meanc(log(eps[1:]))- meanc(log(1+eps[1:].^(vP[3]) ) );
vScore[0][5] = polygamma(vP[4]+vP[5],0)-polygamma(vP[5],0) - meanc(log(eps[1:].^(vP[3])+1));

return vScore ;
}


DySco::GB2AnaCovMatrix(const mAnaHess, const vO)
{
decl a,b,c,sigma_u_2;
	a = vO[1]-vO[2]*vO[3]*(vO[4]+vO[5])*Beta(1+vO[4],1+vO[5])/Beta(vO[4],vO[5]);
	b = vO[1]^(2) -2*vO[1]*vO[2]*vO[3]*(vO[4]+vO[5])*Beta(1+vO[4],1+vO[5])/Beta(vO[4],vO[5]) + vO[2]^(2)*vO[3]^(2)*(vO[4]+vO[5])^(2)*Beta(vO[4]+2,vO[5]+2)/Beta(vO[4],vO[5]);
	c =	vO[2]*vO[3]*vO[4]*vO[5]*(vO[4]-vO[5])/((1+vO[4]+vO[5])*(2+vO[4]+vO[5]));//vO[2]*vO[3]*(vO[4]+vO[5])*(vO[4]*Beta(vO[4]+1,vO[5])-(vO[4]+vO[5])*Beta(vO[4]+2,vO[5]) )/Beta(vO[4],vO[5]);
	//vO[2]*vO[3]*vO[4]*vO[5]*(vO[4]-vO[5])/((1+vO[4]+vO[5])*(2+vO[4]+vO[5]));
	decl fact = vO[3]^(2)*vO[4]*vO[5]/((1+vO[4]+vO[5])*(1-b));
	decl sigma_u2 = vO[4]*vO[5]/(1+vO[4]+vO[5]);

	decl aux = new matrix[6][6]; // omega notation, Andrews version with nu in score (taken from Andres&Harvey*(2012))
	// 1) dynamic parameters
	mAnaHess[0][0][0] =  fact*(1-vO[1])^(2)*(1+a)/(1-a) ;
	mAnaHess[0][1][1] =  fact*vO[2]^(2)*sigma_u2*(1+a*vO[1])/((1-vO[1]^(2))*(1-a*vO[1]));
	mAnaHess[0][2][2] =  fact*sigma_u2;
	mAnaHess[0][0][1] =  fact*a*c*vO[2]*(1-vO[1])/((1-a)*(1-a*vO[1]));
	mAnaHess[0][1][0] =  mAnaHess[0][0][1];
	mAnaHess[0][0][2] =  fact*c*(1-vO[1])/(1-a);
	mAnaHess[0][2][0] =  mAnaHess[0][0][2];
	mAnaHess[0][1][2] =  fact*a*vO[2]*sigma_u2/(1-a*vO[1]);
	mAnaHess[0][2][1] =  mAnaHess[0][1][2];
	// 2) shape parameters
	mAnaHess[0][3][3] = (vO[3]^(2)*(1+vO[4]+vO[5]))^(-1)*(1+vO[4]+vO[5]+vO[4]*vO[5]*(polygamma(vO[4],1)+polygamma(vO[5],1)+ (polygamma(vO[5],0)-polygamma(vO[4],0)+(vO[4]-vO[5])*(vO[4]*vO[5])^(-1))^(2) - (vO[4]^(2)+vO[5]^(2))*vO[4]^(-2)*vO[5]^(-2) ) );
	// wrong but Inf Matrix becomes positive s.d. 
	//mAnaHess[0][3][3] = (vO[3]^(2)*(1+vO[4]+vO[5]))^(-1)*(1+vO[4]+vO[5]+vO[4]*vO[5]*(polygamma(vO[4],1)+polygamma(vO[5],1)+ (polygamma(vO[5],0)-polygamma(vO[4],0)+(vO[4]+vO[5])*(vO[4]*vO[5])^(-1))^(2) - (vO[4]^(2)+vO[5]^(2))*vO[4]^(-2)*vO[5]^(-2) ) );
	mAnaHess[0][4][4] = polygamma(vO[4],1)-polygamma(vO[4]+vO[5],1);
	mAnaHess[0][5][5] = polygamma(vO[5],1)-polygamma(vO[4]+vO[5],1);
	mAnaHess[0][3][4] = -(vO[5]*(polygamma(vO[4],0)-polygamma(vO[5],1))-1)/(vO[3]*(vO[4]+vO[5]));
	mAnaHess[0][4][3] = mAnaHess[0][3][4];
	mAnaHess[0][3][5] = -(vO[4]*(polygamma(vO[5],0)-polygamma(vO[4],0))-1)/(vO[3]*(vO[4]+vO[5]));
	mAnaHess[0][5][3] = mAnaHess[0][3][5];
	mAnaHess[0][4][5] = -polygamma(vO[4]+vO[5],1);
	mAnaHess[0][5][4] = mAnaHess[0][4][5];
	//3) cross derivatives
	mAnaHess[0][0][3] =	 -vO[3]^(-1)*vO[4]*vO[5]*(1+vO[4]+vO[5])^(-1)*(polygamma(vO[4]+1,0)-polygamma(vO[5]+1,0))*(1-vO[1])*(1-a)^(-1);
	mAnaHess[0][3][0] = mAnaHess[0][0][3];
	mAnaHess[0][0][4] =	-(-1+vO[4]/(vO[4]+vO[5]))*(1-vO[1])/(1-a);
	mAnaHess[0][4][0] = mAnaHess[0][0][4];
	mAnaHess[0][0][5] = -vO[4]*(1-vO[1])/((vO[4]+vO[5])*(1-a));
	mAnaHess[0][5][0] = mAnaHess[0][0][5];
																  

	mAnaHess[0] = -mAnaHess[0];
	//if(choleski(mAnaHess[0])==0) {println("Choleski failed"); continue;}
	//println("determinant", determinant(mAnaHess[0])); //println("inverse", invertgen(mAnaHess[0],30)); 
return mAnaHess[0];
}
/////////////////////////////////// MEM  Stuff //////////////////////////////////////////////////////
DySco::GB2MEMLik(const vP, const adFunc, const avScore, const amHessian) 
{decl cT = m_iT2sel-m_iT1sel+1;
decl vLambda = new matrix[cT][1];
if (m_iTwoComp==0) MEMLambdaFilter(&vLambda, m_vY, vP, m_vStruct);
else 	MEMLambdaFilter2comp(&vLambda, m_vY, vP, m_vStruct);
decl amPar =  new array[7];
SplitPara(&amPar,vP,m_vStruct);
adFunc[0] = double(log(vP[3])-log(Beta(vP[4],vP[5])) + (m_cT-1)^(-1)*(vP[3]*vP[4]-1)*sumc(log(m_vY[1:]))-(m_cT-1)^(-1)*vP[3]*vP[4]*sumc(vLambda[1:]) - (m_cT-1)^(-1)*(vP[4]+vP[5])*sumc(log(1+(m_vY[1:]./exp(vLambda[1:])).^(vP[3])  ) )     );
//if (isarray(avScore)){  // if analytical score should be used in maximisation
//	decl vScore = new matrix[5][1];
//	decl aux =  m_vY./exp(vLambda);
//	FMEMScores(&vScore,m_vY, vLambda, vP, max(m_vStruct[0], m_vStruct[1]));
////	println(vScore');
//	(avScore[0])[0] = vScore[0];
//	(avScore[0])[1] = vScore[1];
//	(avScore[0])[2] = vScore[2];
//	(avScore[0])[3] = vScore[3];
//	(avScore[0])[4] = vScore[4];
//}
return 1;
}

DySco::Estimate()
{
//if (m_vParStart==<> || m_vParStart==0){oxwarning("Startingvalues not defined"); exit(0);}
//println("m_vParStart", m_vParStart);

decl vpstart = GetFreePar();  // map pars to estimation format

DoEstimation(1); // do the estimation

decl vpfree = isarray(m_vP) ? m_vP[0] : m_vP;
SetFreePar(vpfree);// map estimated pars to normal format
SetResult(m_iResult);

//println(m_vPar);
//println(MaxConvergenceMsg(m_iResult));

if (m_iResult == 0 || m_iResult == 1)
        m_iModelStatus = MS_ESTIMATED;
    else {
        m_iModelStatus = MS_EST_FAILED;
		oxwarning("Model did not converge. Try e.g. different starting values");println(MaxConvergenceMsg(m_iResult));
		exit(0);
		}

Output();

return m_iModelStatus == MS_ESTIMATED;
}
	 

DySco::DoEstimation(vPa)
{	 ClearModel();
decl ir; decl dFunc; decl time = 0; decl aux;

if (!InitData()) { print("Failed to load data\n"); 
      return FALSE; 
	  }  

if (!InitPar()) { print("Failed to load parameters \n"); 
      return FALSE; 
	   }  
	
MaxControl(m_cMaxIter, 0);

decl vPar = new matrix[m_cP][1];
vPar = GetFreePar();
FreePar(-1);
SetPar(vPar);		   
SetParCount(m_cP); 
SetFreePar(vPar);	 

InitData();
InitPar(); 
	
if(m_vStruct[0]==0 && m_vStruct[1]==0){ // static model estimated
	decl aux2  = DoStaticEstimation();
	vPar = aux2[0]; dFunc = aux2[1]; ir =aux2[2]; time = aux2[3]; aux = aux2[4];
}
else { // non-static model
	if (m_iAutoStartValues==1){
		AutoStartValues(); vPar = m_vParStart;
		}
	else if (m_iAutoStartValues==0)	vPar = m_vParStart;
	else{
		oxwarning(println("starting values not defined \n"));exit(0);
		}
if (m_vStruct[3]==1 && m_iScore==0) {oxwarning("Exact derivatives in numerical optimisation only for models wothout leverage parameter allowed");exit(0);}
if (m_vStruct[3]==1 && m_iRoutine==SCORING) {oxwarning("Exact derivatives in numerical optimisation only for models wothout leverage parameter allowed");exit(0);}


if (m_iModelClass==SCORE) {
    if (m_iRoutine==BFGS){
	time = timer();
	if (m_cDist == EX) 				ir = MaxBFGS(EXLik, &vPar, &dFunc, 0, m_iScore);
	else if (m_cDist==GA)			ir = MaxBFGS(GALik, &vPar, &dFunc, 0, m_iScore);
	else if (m_cDist==WBL)			ir = MaxBFGS(WBLLik, &vPar, &dFunc, 0, m_iScore);
	else if (m_cDist==GGA)			ir = MaxBFGS(GGALik, &vPar, &dFunc, 0, m_iScore);
	else if(m_cDist==LLOG)		    ir = MaxBFGS(LLOGLik, &vPar, &dFunc, 0, m_iScore);
	else if(m_cDist==LNORM)		    ir = MaxBFGS(LNORMLik, &vPar, &dFunc, 0, m_iScore);
	else if (m_cDist==BURR) 		ir = MaxBFGS(BURRLik, &vPar, &dFunc, 0, m_iScore);
	else if (m_cDist==DAG) 			ir = MaxBFGS(DAGLik, &vPar, &dFunc, 0, m_iScore);	
	else if (m_cDist==F)			ir = MaxBFGS(FLik, &vPar, &dFunc, 0, m_iScore);
	else if (m_cDist==GB2)			ir = MaxBFGS(GB2Lik, &vPar, &dFunc, 0, m_iScore);
	time = timer()-time;
	}
  if (m_iRoutine==NR){decl mHess;
	time = timer();
	if (m_cDist == EX) 				ir = MaxNewton(EXLik, &vPar, &dFunc,&mHess,1);
	else if (m_cDist==GA)			ir = MaxNewton(GALik, &vPar, &dFunc,&mHess,1);
	else if (m_cDist==WBL)			ir = MaxNewton(WBLLik, &vPar, &dFunc,&mHess,1);
	else if (m_cDist==GGA)			ir = MaxNewton(GGALik, &vPar, &dFunc,&mHess,1);
	else if(m_cDist==LNORM)			ir = MaxNewton(LNORMLik, &vPar, &dFunc,&mHess,1);
	else if (m_cDist==BURR) 		ir = MaxNewton(BURRLik, &vPar, &dFunc,&mHess,1);
	else if (m_cDist==DAG) 			ir = MaxNewton(DAGLik, &vPar, &dFunc, &mHess, 1);
	else if (m_cDist==F)			ir = MaxNewton(FLik, &vPar, &dFunc,&mHess,1);
	else if (m_cDist==GB2)			ir = MaxNewton(GB2Lik, &vPar, &dFunc,&mHess,1);
	time = timer()-time;
	}
  if (m_iRoutine==SCORING){decl mHess;
	time = timer();
	if (m_cDist == EX)				ir = MaxNewton(EXLik, &vPar, &dFunc,&mHess,0);
	else if (m_cDist==GA)			ir = MaxNewton(GALik, &vPar, &dFunc,&mHess,0);
	else if (m_cDist==WBL)			ir = MaxNewton(WBLLik, &vPar, &dFunc,&mHess,0);
	else if (m_cDist==GGA)			ir = MaxNewton(GGALik, &vPar, &dFunc,&mHess,0);
	else if (m_cDist==LNORM)		ir = MaxNewton(LNORMLik, &vPar, &dFunc,&mHess,0);
  	else if(m_cDist==LLOG) 			ir = MaxNewton(LLOGLik, &vPar, &dFunc,&mHess,0);
	else if (m_cDist==BURR)			ir = MaxNewton(BURRLik, &vPar, &dFunc,&mHess,0);
	else if (m_cDist==DAG) 			ir = MaxNewton(DAGLik, &vPar, &dFunc, &mHess,0);
	else if (m_cDist==F)			ir = MaxNewton(FLik, &vPar, &dFunc,&mHess,0);
	else if (m_cDist==GB2)			ir = MaxNewton(GB2Lik, &vPar, &dFunc,&mHess,0);
	time = timer()-time;
	}
  }
  else if (m_iModelClass==MEM){
	time = timer();
	if (m_cDist ==EX)				ir = MaxBFGS(EXMEMLik, &vPar, &dFunc,0,m_iScore);
	else if (m_cDist==GA)			ir = MaxBFGS(GAMEMLik, &vPar, &dFunc,0,m_iScore);
	else if (m_cDist==WBL)			ir = MaxBFGS(WBLMEMLik, &vPar, &dFunc,0,m_iScore);
	else if (m_cDist==GGA)			ir = MaxBFGS(GGAMEMLik, &vPar, &dFunc,0,m_iScore);
	else if (m_cDist==LNORM)		ir = MaxBFGS(LNORMMEMLik, &vPar, &dFunc,0,m_iScore);
  	else if(m_cDist==LLOG)			ir = MaxBFGS(LLOGMEMLik, &vPar, &dFunc,0,m_iScore);
	else if (m_cDist ==BURR)		ir = MaxBFGS(BURRMEMLik, &vPar, &dFunc,0,m_iScore);
	else if (m_cDist ==DAG)			ir = MaxBFGS(DAGMEMLik, &vPar, &dFunc,0,m_iScore);
	else if (m_cDist == F)			ir = MaxBFGS(FMEMLik, &vPar, &dFunc,0,m_iScore);
	else if (m_cDist == GB2)		ir = MaxBFGS(GB2MEMLik, &vPar, &dFunc,0,m_iScore);
	time = timer()-time;
 }

}

// LogNormal variance parameter estimates seperately
if (m_cDist==LNORM && m_iModelClass==SCORE) { decl maxpq = max(m_vStruct[0], m_vStruct[1]);
decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
decl vU = new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) 	LNORMLambdaFilter(&vLambda,&vU, m_mData[m_iT1sel:m_iT2sel][0], vPar, m_vStruct, m_cT, maxpq);
else  				LambdaFilter2comp(&vLambda, &vU, m_mData[m_iT1sel:m_iT2sel][0], vPar, m_vStruct, LNORM);
decl vParaux;
if (m_iTwoComp!=1 && m_vStruct[3]!=1){
	vParaux = vPar[:m_vStruct[0]+m_vStruct[1]];
	vParaux|= sumc((log(m_mData[m_iT1sel+maxpq:m_iT2sel][0])-vLambda[maxpq:]).^(2))/(m_cT-maxpq);
	}
else {
	vParaux = vPar[:m_vStruct[0]+m_vStruct[1]];
	vParaux|= sumc((log(m_mData[m_iT1sel+maxpq:m_iT2sel][0])-vLambda[maxpq:]).^(2))/(m_cT-maxpq);
	vParaux|= vPar[m_vStruct[0]+m_vStruct[1]+1:];
	}
delete vPar;
vPar = vParaux;
}

if (m_cDist==LNORM && m_iModelClass==MEM) {decl maxpq = max(m_vStruct[0], m_vStruct[1]);
decl vLambda = new matrix[m_iT2sel-m_iT1sel+1][1];
decl vU = new matrix[m_iT2sel-m_iT1sel+1][1];
if (m_iTwoComp==0) 	MEMLambdaFilter(&vLambda, m_mData[m_iT1sel:m_iT2sel][0], vPar, m_vStruct);
else  				MEMLambdaFilter2comp(&vLambda,m_mData[m_iT1sel:m_iT2sel][0], vPar, m_vStruct);
decl vParaux;// =vPar;
if (m_iTwoComp!=1 && m_vStruct[3]!=1){
	vParaux = vPar[:m_vStruct[0]+m_vStruct[1]];
	vParaux|= sumc((log(m_mData[m_iT1sel+maxpq:m_iT2sel][0])-vLambda[maxpq:]).^(2))/(m_cT-maxpq);
	}
else {
	vParaux = vPar[:m_vStruct[0]+m_vStruct[1]];
	vParaux|= sumc((log(m_mData[m_iT1sel+maxpq:m_iT2sel][0])-vLambda[maxpq:]).^(2))/(m_cT-maxpq);
	vParaux|= vPar[m_vStruct[0]+m_vStruct[1]+1:];
	}
delete vPar;
vPar = vParaux;
//check parameter restrictions hit?
CheckPara(vPar);
}
// define globals
m_vP = vPar;
m_iResult = ir;
m_cElapsedTime = time/100;
m_dLogLik = m_cT*dFunc;

//decl aux;
if (m_iRoutine==SCORING || m_iScore==0) aux = TRUE;
else aux= FALSE;
return {vPar, dFunc, ir, time, aux};
}

DySco::GetNumStdErr(vNumStdErr, const vP, const iDist)
{ //decl cPara = rows(vP);
decl mNumHess; //= new matrix[cPara][cPara];
if (iDist==EX) 			Num2Derivative(EXLik, vP, &mNumHess);
else if (iDist==GA) 	Num2Derivative(GALik, vP, &mNumHess);
else if (iDist==WBL)	Num2Derivative(WBLLik, vP, &mNumHess);
else if (iDist==GGA)	Num2Derivative(GGALik, vP, &mNumHess);
else if (iDist==LLOG)	Num2Derivative(LLOGLik, vP, &mNumHess);
else if (iDist==BURR)	Num2Derivative(BURRLik, vP, &mNumHess);
else if (iDist==DAG)	Num2Derivative(DAGLik, vP, &mNumHess);
else if (iDist==F) 		Num2Derivative(FLik, vP, &mNumHess);
else if (iDist==GB2) 	Num2Derivative(GB2Lik, vP, &mNumHess);

else if (iDist==LNORM){
	Num2Derivative(LNORMLik, vP, &mNumHess);
	//println("mNumHess", mNumHess);
	if(m_iTwoComp == 0 && m_vStruct[3]==0){
//	println("sum", m_vStruct[0]+m_vStruct[1]);
//	println(vP);
//	println(vP[1+m_vStruct[0]+m_vStruct[1]]);
	mNumHess = mNumHess[0:m_vStruct[0]+m_vStruct[1]][0:m_vStruct[0]+m_vStruct[1]]~zeros(1+m_vStruct[0]+m_vStruct[1],1)|zeros(1,1+m_vStruct[0]+m_vStruct[1])~(m_iT2sel-m_iT1sel)*0.5*vP[1+m_vStruct[0]+m_vStruct[1]]^(-2);
	}
	else {decl nxt = m_vStruct[0]+m_vStruct[1],	aux;
	aux = mNumHess[0:nxt][0:nxt]~zeros(nxt+1,1)~mNumHess[0:nxt][nxt+1:nxt+m_vStruct[3]+m_vStruct[4]+m_vStruct[5]];
	aux |=	zeros(1,nxt+1)~(m_iT2sel-m_iT1sel)*0.5*vP[1+m_vStruct[0]+m_vStruct[1]]^(-2)~zeros(1,m_vStruct[3]+m_vStruct[4]+m_vStruct[5]);
	aux |= mNumHess[nxt+1:nxt+m_vStruct[3]+m_vStruct[4]+m_vStruct[5]][0:nxt] ~zeros(m_vStruct[3]+m_vStruct[4]+m_vStruct[5],1)~mNumHess[nxt+1:nxt+m_vStruct[3]+m_vStruct[4]+m_vStruct[5]][nxt+1:nxt+m_vStruct[3]+m_vStruct[4]+m_vStruct[5]];
	//mNumHess[nxt+1:m_vLev+m_cAR2comp+m_cMA2comp-1][0:nxt] ~zeros(m_vLev+m_cAR2comp+m_cMA2comp,1)~ mNumHess[nxt+1:m_vLev+m_cAR2comp+m_cMA2comp-1][nxt+1:m_vLev+m_cAR2comp+m_cMA2comp];
	//		println(aux);
	mNumHess =aux;
	}
}
//println("mNumHess", mNumHess);
decl cMaxpq = max(m_vStruct[0], m_vStruct[1]);
//mNumHess = sqrt(m_cT-cMaxpq)*mNumHess;
decl mInvNumHessian =  invertgen(-mNumHess,3)/(m_cT-cMaxpq);
vNumStdErr  = sqrt(diagonal(mInvNumHessian));
return vNumStdErr;
}


DySco::GetStaticNumStdErr(vNumStdErr, const vP, const iDist)
{ //decl cPara = rows(vP);
decl mNumHess; //= new matrix[cPara][cPara];
Num2Derivative(GetStaticLL, vP, &mNumHess);
decl mInvNumHessian =  invertgen(-mNumHess,3)/(m_cT);
vNumStdErr  = sqrt(diagonal(mInvNumHessian));
return vNumStdErr;
}




DySco::GetNumStdErrMEM(vNumStdErr, const vP, const iDist)
{decl mNumHess; //= new matrix[cPara][cPara];
if (iDist==EX) 			Num2Derivative(EXMEMLik, vP, &mNumHess);
else if (iDist==GA) 	Num2Derivative(GAMEMLik, vP, &mNumHess);
else if (iDist==WBL)	Num2Derivative(WBLMEMLik, vP, &mNumHess);
else if (iDist==GGA)	Num2Derivative(GGAMEMLik, vP, &mNumHess);
else if (iDist==LLOG)	Num2Derivative(LLOGMEMLik, vP, &mNumHess);
else if (iDist==BURR)	Num2Derivative(BURRMEMLik, vP, &mNumHess);
else if (iDist==DAG)	Num2Derivative(DAGMEMLik, vP, &mNumHess);
else if (iDist==F) 		Num2Derivative(FMEMLik, vP, &mNumHess);
else if (iDist==GB2)	Num2Derivative(GB2MEMLik, vP, &mNumHess);
else if (iDist==LNORM){
	Num2Derivative(LNORMMEMLik, vP, &mNumHess);
	//println("mNumHess", mNumHess);
	if(m_iTwoComp == 0 && m_vStruct[3]==0){
//	println("in hessian");
//	println("sum", m_vStruct[0]+m_vStruct[1]);
//	println(vP);
//	println(vP[1+m_vStruct[0]+m_vStruct[1]]);
	mNumHess = mNumHess[0:m_vStruct[0]+m_vStruct[1]][0:m_vStruct[0]+m_vStruct[1]]~zeros(1+m_vStruct[0]+m_vStruct[1],1)|zeros(1,1+m_vStruct[0]+m_vStruct[1])~(m_iT2sel-m_iT1sel)*0.5*vP[1+m_vStruct[0]+m_vStruct[1]]^(-2);
	}
	else {decl nxt = m_vStruct[0]+m_vStruct[1],	aux;
	aux = mNumHess[0:nxt][0:nxt]~zeros(nxt+1,1)~mNumHess[0:nxt][nxt+1:nxt+m_vStruct[3]+m_vStruct[4]+m_vStruct[5]];
	aux |=	zeros(1,nxt+1)~(m_iT2sel-m_iT1sel)*0.5*vP[1+m_vStruct[0]+m_vStruct[1]]^(-2)~zeros(1,m_vStruct[3]+m_vStruct[4]+m_vStruct[5]);
	aux |= mNumHess[nxt+1:nxt+m_vStruct[3]+m_vStruct[4]+m_vStruct[5]][0:nxt] ~zeros(m_vStruct[3]+m_vStruct[4]+m_vStruct[5],1)~mNumHess[nxt+1:nxt+m_vStruct[3]+m_vStruct[4]+m_vStruct[5]][nxt+1:nxt+m_vStruct[3]+m_vStruct[4]+m_vStruct[5]];
	//mNumHess[nxt+1:m_vLev+m_cAR2comp+m_cMA2comp-1][0:nxt] ~zeros(m_vLev+m_cAR2comp+m_cMA2comp,1)~ mNumHess[nxt+1:m_vLev+m_cAR2comp+m_cMA2comp-1][nxt+1:m_vLev+m_cAR2comp+m_cMA2comp];
	//		println(aux);
	mNumHess =aux;
	}
}
//println("mNumHess", mNumHess);
decl cMaxpq = max(m_vStruct[0], m_vStruct[1]);
//mNumHess = sqrt(m_cT-cMaxpq)*mNumHess;

decl mInvNumHessian =  invertgen(-mNumHess,3)/(m_cT-cMaxpq);
vNumStdErr  = sqrt(diagonal(mInvNumHessian));
return vNumStdErr;

}


DySco::GetAnaStdErr(vAnaStdErr, const vP, const iDist, const cT, const cMaxpq)
{ decl mAnaHess = new matrix[m_cP][m_cP];
if (iDist==EX) 			EXAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq);
else if (iDist==GA) 	GAAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq);
else if (iDist==WBL) 	WBLAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq);
else if (iDist==GGA) 	GGAAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq);
else if (iDist==LLOG) 	LLOGAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq);
else if (iDist==BURR) 	BURRAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq);
else if (iDist==DAG) 	DAGAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq);
else if (iDist==F) 		FAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq);
else if (iDist==GB2)	GB2AnaCovMatrix(&mAnaHess, vP);
else if (iDist==LNORM) 	{
	decl vLambda = new matrix[cT][1];
	decl vU = new matrix[cT][1];
	LNORMLambdaFilter(&vLambda,&vU, m_vY, vP, m_vStruct,  cT, cMaxpq);
	decl vEps = new matrix[cT][1];
	vEps = m_vY./exp(vLambda);
	LNORMAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq, vEps);
	mAnaHess = mAnaHess[0:m_vStruct[0]+m_vStruct[1]][0:m_vStruct[0]+m_vStruct[1]]~zeros(1+m_vStruct[0]+m_vStruct[1],1)|zeros(1,1+m_vStruct[0]+m_vStruct[1])~(m_iT2sel-m_iT1sel)*0.5*vP[1+m_vStruct[0]+m_vStruct[1]]^(-2);
	}
//mAnaHess = (cT-cMaxpq)*mAnaHess;
decl mInvAnaHessian =  invertgen(-mAnaHess,3)/(m_cT-cMaxpq);
vAnaStdErr  = sqrt(diagonal(mInvAnaHessian));
return vAnaStdErr;
}


DySco::GetMEMAnaStdErr(vAnaStdErr, const vP, const iDist, const cT, const cMaxpq)
{ decl mAnaHess = new matrix[m_cP][m_cP];
if (iDist==EX) 			EXAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq);
//else if (iDist==GA) 	GAAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq);
//else if (iDist==WBL) 	MEMWBLAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq);
//else if (iDist==GGA) 	MEMGGAAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq);
//else if (iDist==LLOG) 	MEMLLOGAnaCovMatr1258ix(&mAnaHess, vP, cT, cMaxpq);
//else if (iDist==BURR) 	MEMBURRAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq);
//else if (iDist==DAG) 	MEMDAGAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq);
//else if (iDist==F) 		MEMFAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq);
else if (iDist==LNORM) 	{
	decl vLambda = new matrix[cT][1];
	decl vU = new matrix[cT][1];
	MEMLambdaFilter(&vLambda, m_vY, vP,m_vStruct);
	decl vEps = new matrix[cT][1];
	vEps = m_vY./exp(vLambda);
	MEMLNORMAnaCovMatrix(&mAnaHess, vP, cT, cMaxpq, vEps);
	mAnaHess = mAnaHess[0:m_vStruct[0]+m_vStruct[1]][0:m_vStruct[0]+m_vStruct[1]]~zeros(1+m_vStruct[0]+m_vStruct[1],1)|zeros(1,1+m_vStruct[0]+m_vStruct[1])~(m_iT2sel-m_iT1sel)*0.5*vP[1+m_vStruct[0]+m_vStruct[1]]^(-2);
	}
//mAnaHess = (cT-cMaxpq)*mAnaHess;
decl mInvAnaHessian =  invertgen(-mAnaHess,3)/(m_cT-cMaxpq);
vAnaStdErr  = sqrt(diagonal(mInvAnaHessian));
return vAnaStdErr;
}

 
////////////////// Static Model Estimation /////////////////////////////
DySco::LogNormalLikStatic(const vY,  const vPara, const cMaxpq, const cT)
{// To estimate Static LogNormal parameters // vPara= <location para; sigma^(2)>
	return double (-0.5*(cT-cMaxpq)*(log(M_2PI)+log(vPara[1])) - sumc(log(vY)) -0.5*vPara[1]^(-1)*sumc((log(vY)-vPara[0]).^(2)) ); }

DySco::GetStaticLL(const vPara, const adFunc, const avScore, const amHessian) 
{// Purpose: obtain static estimate of model // e.g. for automatic startingvals and DCS(0,0) model
// vPara= [scale para; shape para(s)]
decl maxpq = 1; // no degrees of freedeom reduction necessary
decl cT = m_cT;
decl vLambda = new matrix[cT][1];
vLambda = vPara[0]*ones(cT,1); // set m_lambdaStat to a constant	
// distribution specific stuff
if (m_cDist == EX){
	adFunc[0] = double((-sumc(vLambda[1:])- sumc(m_mData[m_iT1sel+1:m_iT2sel][0]./exp(vLambda[1:])))/cT);
	}	
if (m_cDist == GA){
	adFunc[0] = double( -loggamma(vPara[1])+(-vPara[1]*sumc(vLambda[maxpq:])+(vPara[1]-1)*sumc(log(m_mData[m_iT1sel+maxpq:m_iT2sel][0]))-sumc(m_mData[m_iT1sel+maxpq:m_iT2sel][0].*exp(-vLambda[maxpq:])) )/(cT) );
	}
if (m_cDist == WBL){
	adFunc[0] = double(log(vPara[1])+(-vPara[1]*sumc(vLambda[1:]) + (vPara[1]-1)*sumc(log(m_mData[m_iT1sel+1:m_iT2sel][0]))-sumc((m_mData[m_iT1sel+1:m_iT2sel][0]./exp(vLambda[1:])).^(vPara[1])))/(cT-1) );
}
if (m_cDist == LNORM){
	adFunc[0] = LogNormalLikStatic(m_mData[m_iT1sel+1:m_iT2sel][0], vPara, maxpq, cT);
}
if (m_cDist == LLOG){
	adFunc[0] = double(log(vPara[1]) - vPara[1]*(m_iT2sel-m_iT1sel)^(-1)*sumc(vLambda[1:]) + (vPara[1]-1)*(m_iT2sel-m_iT1sel)^(-1)*sumc(log(m_mData[m_iT1sel+1:m_iT2sel][0])) -2*(m_iT2sel-m_iT1sel)^(-1)*sumc(log(1+(m_mData[m_iT1sel+1:m_iT2sel][0]./exp(vLambda[1:])).^(vPara[1]))) );
}
if (m_cDist == BURR){
	adFunc[0] = double(log(vPara[1]) + log(vPara[2]) + (m_iT2sel-m_iT1sel)^(-1)*(vPara[1]-1)*sumc(log(m_mData[m_iT1sel+1:m_iT2sel][0]))- (m_iT2sel-m_iT1sel)^(-1)*vPara[1]*sumc(vLambda[1:]) - (m_iT2sel-m_iT1sel)^(-1)*(vPara[2]+1)*sumc(log((m_mData[m_iT1sel+1:m_iT2sel][0]./exp(vLambda[1:])).^(vPara[1])+1)));
}
if (m_cDist == DAG){
	adFunc[0] = double(log(vPara[1])+  log(vPara[2]) + (m_iT2sel-m_iT1sel)^(-1)*(vPara[1]*vPara[2]-1)*sumc(log(m_mData[m_iT1sel+1:m_iT2sel][0]))- (m_iT2sel-m_iT1sel)^(-1)*vPara[1]*vPara[2]*sumc(vLambda[1:]) - (m_iT2sel-m_iT1sel)^(-1)*(vPara[2]+1)*sumc(log((m_mData[m_iT1sel+1:m_iT2sel][0]./exp(vLambda[1:])).^(vPara[1])+1)));
	// Dagum type II
//	adFunc[0] = double(log(vPara[1]) + log(-vPara[1]+1) + log(vPara[2]) - (m_iT2sel-m_iT1sel)^(-1)*(vPara[1]+1)*sumc(log(m_mData[m_iT1sel+1:m_iT2sel][0]))+ (m_iT2sel-m_iT1sel)^(-1)*vPara[1]*sumc(vLambda[1:]) - (m_iT2sel-m_iT1sel)^(-1)*(vPara[2]+1)*sumc(log((m_mData[m_iT1sel+1:m_iT2sel][0]./exp(vLambda[1:])).^(-vPara[1])+1)));
}  
if (m_cDist == F){
//	adFunc[0] = double(0.5*vPara[1]*log(vPara[1]) +0.5*vPara[2]*log(vPara[2])-0.5*vPara[1]*(cT)^(-1)*sumc(vLambda[1:]) + (0.5*vPara[1]-1)*(cT)^(-1)*sumc(log(m_mData[m_iT1sel+1:m_iT2sel][0]))-0.5*(vPara[1]+vPara[2])*(cT)^(-1)*sumc(log(vPara[1])*m_mData[m_iT1sel+1:m_iT2sel]./exp(vLambda[1:]) + vPara[2]) )-log(Beta(0.5*vPara[1],0.5*vPara[2])));
	adFunc[0] = double(0.5*vPara[1]*log(vPara[1]) +0.5*vPara[2]*log(vPara[2])-	 0.5*vPara[1]*(cT)^(-1)*sumc(vLambda[1:]) + (0.5*vPara[1]-1)*(cT)^(-1)*sumc(log(m_mData[m_iT1sel+1:m_iT2sel][0]))-0.5*(vPara[1]+vPara[2])*(cT)^(-1)*sumc(log(vPara[1]*m_mData[m_iT1sel+1:m_iT2sel][0]./exp(vLambda[1:]) + vPara[2]) )-log(Beta(0.5*vPara[1],0.5*vPara[2])));
}
if (m_cDist == GGA){
	adFunc[0] = double(( (cT-maxpq)*(log(vPara[1])-loggamma(vPara[2])) -vPara[1]*vPara[2]*sumc(vLambda[maxpq:]) + (vPara[1]*vPara[2]-1)*sumc(log(m_mData[m_iT1sel+maxpq:m_iT2sel][0]))- sumc( (m_mData[m_iT1sel+maxpq:m_iT2sel][0]./exp(vLambda[maxpq:])).^(vPara[1]) ))/cT );

}
if (m_cDist == GB2){
	adFunc[0] = double( (cT-maxpq)*(log(vPara[1])-log(Beta(vPara[2],vPara[3])) ) - (vPara[1]*vPara[2])*sumc(vLambda[maxpq:]) + (vPara[1]*vPara[2]-1)*sumc(log(m_mData[m_iT1sel+maxpq:m_iT2sel][0]))-	(vPara[2]+vPara[3])*sumc(log(1+(m_mData[m_iT1sel+maxpq:m_iT2sel][0]./exp(vLambda[maxpq:])).^(vPara[1]) ) )/cT);

}
		
return 1;                           
}





DySco::DoStaticEstimation()
{// Purpose: Estimate static Model
	ClearModel();
	if (oxversion() < 320)
	{println("\n");println("Requires Ox 3.2 or newer!");
     println("Please update your version of Ox ", oxversion()/100,"\n");
     exit(0);}
// Set Starting Values
	decl cMaxLL, vP;
	if (m_cDist==GA || m_cDist == WBL || m_cDist==LLOG) vP= <1;2>;
	else if (m_cDist==LNORM) vP=<1;0.05>;
	else if (m_cDist==BURR|| m_cDist==GGA)	 vP=<1;2;2>;
	else if (m_cDist==DAG)	   vP=<1;0.8;1.8>;
	else if (m_cDist==F) vP=<1;22;22>;
	else if (m_cDist==GB2) vP=<1;3.4;0.98;1>;
		
// 	if (!InitData()){ 
//	   	print("Failed to load data\n"); 
//       	return; 
//	}  
	
	MaxControl(m_cMaxIter, 0);
	decl ir, time;
	time = timer();
	ir= MaxBFGS(GetStaticLL, &vP, &cMaxLL, 0, TRUE);	
	time = (timer()-time)/100;

	decl auxtest =dropr(vP,0);
    decl X = sumc(auxtest .< 0);
	if (X > 0){
		println("None of the shape parameters can be negative. You have choosen ", vP);exit(0);}
	if (isnan(vP)==1) {oxwarning("DoStaticEstimation(): obtained parameter values cant be NaN \n Try using manual starting values."); exit(0);}
	//m_vStaticPar = vP;			
//	println("vP", vP);
decl aux;
if (m_iRoutine==SCORING || m_iScore==0) aux = TRUE;
else aux= FALSE;

return {vP, cMaxLL, ir, time, aux};
 
}



DySco::AutoStartValues()
{// Purpose: gets automatic starting vals
// Requires: 	m_cP
//				m_mData
//				m_iT1sel m_iT2sel
print("Determining starting values ...\t");	

decl vStartVals = new matrix[m_cP][1];
// 1) Shape Parameters
//println("vStartVals", vStartVals);
decl vStaticPar;
if (m_cDist != EX){
decl aux = DoStaticEstimation();
vStaticPar = aux[0];}
//Fprintln("vStaticPar", vStaticPar);
// ==> m_vStaticPar
// 2) Intercept = log(mean(data))-log(mean(varepsilon))
if (m_cDist==EX || m_cDist==LNORM) vStartVals[0] = log(meanc(m_mData[m_iT1sel:m_iT2sel][0]));
else if (m_cDist==GA) 		vStartVals[0] = log(meanc(m_mData[m_iT1sel:m_iT2sel][0]));//-log(vStaticPar[1]);
else if (m_cDist==WBL) 		vStartVals[0] =	log(meanc(m_mData[m_iT1sel:m_iT2sel][0]));//-log(gammafact(1+1/vStaticPar[1]));
//needs to be finished for GB2 and GG distributon
else 		vStartVals[0] = log(meanc(m_mData[m_iT1sel:m_iT2sel][0])) ;

 
// For DCS(0,0) model: only need fixed starting val for shape parameter
if (m_vStruct[0]>0){
// 3) AR values: use PACF values as starting vals
// 3.1) AR(1) one comp model  ==> set AR(1) fixed to 0.98
if (m_vStruct[0]==1)  vStartVals[1] = 0.98;
// 3.2) AR(p) cone comp model ==> use PACF's: normalize so that they sum to 0.98 and scale all PACF up to lag length m_vStruct[1]
else if (m_vStruct[0]>1){
	 decl vAcf = acf(log(m_mData[m_iT1sel:m_iT2sel][0]),m_vStruct[0]);
	 decl vPacf = pacf(vAcf);
	 delete vAcf;
	 decl cScale =0.98/sumc(vPacf[1:m_vStruct[0]]);
	 vPacf = cScale.*vPacf;
	 vStartVals[1:m_vStruct[0]] = vPacf[1:m_vStruct[0]];
	 }
}

// 4) MA values:
decl nxt = m_vStruct[0]+1;
if (m_vStruct[1]>0){  // Exclude DCS(0,0) model
// 4.1) 1 Comp Model: set (arbitrarily) to 0.1*0.7^(MA lag length-1)
if (m_iTwoComp == 0) vStartVals[m_vStruct[0]+1] = 0.1;
// 4.2) 2 Comp Model: set (arbitrarily) to 0.05*0.7^(MA lag length-1)
else  vStartVals[m_vStruct[0]+1] = 0.05;
for (decl i=1;i<m_vStruct[1];i++){
	nxt += 1;
	vStartVals[nxt]  = vStartVals[nxt-1]*0.7;
	}
nxt+=1;
}

if(m_vStruct[0]!=0 && m_vStruct[1]!=0){
//println("in gererer\ae");
// 5) Shape Parameters
//if (m_cDist != EX)// nxt += 1;
//else{
//{DoStaticEstimation();
if (m_cDist == GA || m_cDist==WBL || m_cDist==LLOG ){//}//|| m_cDist==LNORM){
   vStartVals[nxt] = vStaticPar[1];
   nxt += 1;
   }
else if (m_cDist==BURR || m_cDist==F || m_cDist==GGA|| m_cDist==DAG){
//	println(m_vStaticPar);
	vStartVals[nxt:nxt+1] = vStaticPar[1:2];
	nxt += 2;
	}
else if (m_cDist==GB2){
	vStartVals[nxt:nxt+2] = vStaticPar[1:3];
	nxt += 3;
	}
}
//}
else{
// Use Starting vals from DoStaticEstimation
	if (m_cDist==GA || m_cDist == WBL || m_cDist==LLOG|| m_cDist==LNORM) vStartVals[1]=2;
	else if (m_cDist==BURR|| m_cDist==F || m_cDist==GGA) vStartVals[1:2]=<2;2>;
	else if (m_cDist==GB2) vStartVals = <2;2;2>;
}

//if (m_cDist==GA){vStartVals[0] = log(meanr(m_mData[m_iT1sel:m_iT2sel][0]))-log(vStartVals[3]);}

// 4) Leverage Parameter: set to 0.05 per default
if (m_vStruct[3]==1){
	vStartVals[nxt] = 0.05;
	nxt += 1;
	}

// 5) 2 Comp Model AR Parameters:
// // See above: but aggreate the PACF to 0.8 only
if (m_vStruct[4]==1 && m_iTwoComp==1){
	vStartVals[nxt] = 0.95;
	nxt +=1;
	}
else if (m_vStruct[4]!=1 && m_iTwoComp==1){
	 decl vAcf2 = acf(log(m_mData[m_iT1sel:m_iT2sel][0]),m_vStruct[4]);
	 decl vPacf2 = pacf(vAcf2);
	 decl cScale2 =0.8/sumc(vPacf2[1:m_vStruct[4]]);
	 vPacf2 = cScale2.*vPacf2;
	 vStartVals[nxt:nxt+m_vStruct[4]-1] = vPacf2[1:m_vStruct[4]];
	 nxt += m_vStruct[4];
	 }

// 6) 2 Comp Model MA Parameters:
if (m_iTwoComp == 1){
vStartVals[nxt] = 0.1;
// 3.2) 2 Comp Model: set (arbitrarily) to 0.05*0.7^(MA lag length-1)
for (decl i=1;i<m_vStruct[5];i++){
	nxt += 1;
	vStartVals[nxt]  = vStartVals[nxt-1]*0.7;
	}
}

		
//println("starting vals: ", vStartVals);
print("Done! \t Starting Estimation ... \n");
m_vParStart = vStartVals;
m_iAutoStartValues = 1;
}

//}

DySco::UserStartValues(const vInitPar)
{m_iAutoStartValues=0;
return m_vParStart = vInitPar;
println("Starting Estimation ...");
} 



///////////////////////////////////////////////////////////////////////////
////////////////////////// Output tables //////////////////////////////////
///////////////////////////////////////////////////////////////////////////
DySco::OutputHeader(const sTitle)
{	println("\n---------     ",sTitle , "    ---------");
	println("            -- Philipp Andres (2013) -- ");
    println("%-35s","Estimated model:",GetModelClassString());
	println("%-35s", "Distribution:",GetDistriString());
	println("%-35s", "Maximisation Routine:", GetNumAlgoString() );
	println("%-35s", "Used Dataset:",GetDbName());
	println("%-35s", "# In-Sample Observations:", GetSelSample() );
	println("%-35s", "# Out-of Sample Observations:",m_cTAll-m_iT2sel -1);
	if (m_vStruct[0]==0 && m_vStruct[1]==0) // i.e. static model
	{println("");}
	else  println("%-35s", "Starting Values:", "%-6.2g", GetStartParams()');
	println("%-35s", "Time to convergence: ", m_cElapsedTime," seconds");
	println("%-35s", "Convergence:", MaxConvergenceMsg(m_iResult));
	println("\n ------------------------------------------");
}

DySco::PrintLikStats()
{
  println("----------------------------------");
 println("log-likelihood: ", m_dLogLik,"\n AIC: ", -2*m_dLogLik+2*m_cP, "\n BIC: ", -2*m_dLogLik+m_cP*log(m_iT2sel) );
 println("----------------------------------");
}

DySco::PrintLikStatsLatex()
{
 println(" \\hline");
 println("log-likelihood: &", m_dLogLik,"& AIC:& ", -2*m_dLogLik+2*m_cP, "&BIC:& ", -2*m_dLogLik+2*log(m_cP), "\\\\" );
}

DySco::Output()					
{
	//SetModelStatus(MS_ESTIMATED);
	if (m_iModelStatus!= MS_ESTIMATED) {oxwarning("Need to estimate model first");exit(0);}
	OutputHeader("Package: " ~ GetPackageName()~" , Version: "~GetPackageVersion());
	if (m_iModelClass == SCORE)   DCSOutputPar();
	else if (m_iModelClass == MEM) MEMOutputPar();
	else { exit(0);}
}

DySco::DCSOutputPar()
{
//	println("m_cP", m_cP);
//	println("m_vP", m_vP');
	
	decl aux;
	if (m_cDist!=LNORM) aux=m_cP;
	else aux = m_cP+1;
	decl vNumStdErr = new matrix[1][aux];
    if (m_vStruct[0]==0 && m_vStruct[1]==0) vNumStdErr = GetStaticNumStdErr(&vNumStdErr, m_vP, m_cDist);
	else vNumStdErr = GetNumStdErr(&vNumStdErr, m_vP,m_cDist);

	//println("vNumStdErr", vNumStdErr);
	println("");
    println("Table 1. Estimated parameters.");

//	println("m_vP", m_vP);
    decl i;	
    decl aspar = ParNames(); //<1;2;3;4>; //GetParNames();
// For One comp without leverage and DCM(1,1) ==> use Analytical std errors

	decl mpar;
    decl cdfloss = GetcDfLoss();
	decl tval, tprob;

	if (m_vStruct[0]==1 && m_vStruct[1]==1 && m_vStruct[3]==0 && m_iTwoComp==0){
		decl cMaxpq = max(m_vStruct[0], m_vStruct[1]);
		decl vAnaStdErr = new matrix[1][m_cP];
		vAnaStdErr = GetAnaStdErr(&vAnaStdErr, m_vP,m_cDist,  m_cT, cMaxpq);
		tval = m_vP./vAnaStdErr';
		tprob = tailt(fabs(tval), m_cT-cdfloss);
    	mpar = m_vP ~ vAnaStdErr'~vNumStdErr' ~ tval ~ tprob ;
		println("%-15s","", "%-12s", "Coefficient", "%-14s", "Ana.Std.Error", "%-14s", "Num.Std.Error","%-9s", "t-stats*", "%-9s", "p-value*");
    	for (i = 0; i < rows(m_vP); ++i){			  
			println("%-15s", aspar[i], " ","%-11.4g", mpar[i][0], " ","%-13.4g",mpar[i][1], " ", "%-13.4g", mpar[i][2], " ","%-8.4g",mpar[i][3]," ","%-8.4f",mpar[i][4]);
			}
			
//		mpar = m_vP ~ vAnaStdErr'~vNumStdErr' ;
//		println("%-15s","", "%-12s", "Coefficient", "%-14s", "Ana.Std.Error", "%-14s", "Num.Std.Error");
//    	for (i = 0; i < rows(m_vP); ++i){			  
//			println("%-15s", aspar[i], " ","%-11.4g", mpar[i][0], " ","%-13.4g",mpar[i][1], " ", "%-13.4g", mpar[i][2]);
//			}


		println("-- * Indicates statistic based on analytical standard errors --");
	}
	else {
		tval = m_vP./vNumStdErr';
		tprob = tailt(fabs(tval), m_cT-cdfloss);
    	mpar = m_vP ~ vNumStdErr' ~ tval ~ tprob ;		
		println("%-15s","", "%-12s", "Coefficient", "%-11s", "Std.Error","%-8s", "t-value", "%-8s", "p-value");
    	for (i = 0; i < rows(m_vP); ++i){			  
			println("%-15s", aspar[i], " ","%-11.4g", mpar[i][0], " ","%-10.4g",mpar[i][1], " ","%-7.4g",mpar[i][2]," ","%-7.4f",mpar[i][3]);
			}
   }

//    	mpar = m_vP ~ vNumStdErr' ;		
//		println("%-15s","", "%-12s", "Coefficient", "%-11s", "Std.Error");
//    	for (i = 0; i < rows(m_vP); ++i){			  
//			println("%-15s", aspar[i], " ","%-11.4g", mpar[i][0], " ","%-10.4g",mpar[i][1]);
//			}
//   }

   
PrintLikStats();
}

DySco::DCSOutputParLaTex()
{
if (m_iModelStatus!= MS_ESTIMATED) {oxwarning("Need to estimate model first");exit(0);}
decl vNumStdErr= new matrix[1][m_cP];
//println(m_vP~m_cDist);
vNumStdErr=GetNumStdErr(&vNumStdErr, m_vP,m_cDist);
println("");
println("Table 1 in LaTex format. Estimated parameters.");
decl i;	
decl aspar = ParNames(); //<1;2;3;4>; //GetParNames();
decl mpar;
decl cdfloss = GetcDfLoss();
decl tval, tprob;

if (m_vStruct[0]==1 && m_vStruct[1]==1 && m_vStruct[3]==0 && m_iTwoComp==0 && m_cDist!=GB2){
	decl cMaxpq = max(m_vStruct[0], m_vStruct[1]);
	decl vAnaStdErr = new matrix[1][m_cP];
	vAnaStdErr = GetAnaStdErr(&vAnaStdErr, m_vP,m_cDist,  m_cT, cMaxpq);
	tval = m_vP./vAnaStdErr';
	tprob = tailt(fabs(tval), m_cT-cdfloss);
    mpar = m_vP ~ vAnaStdErr'~vNumStdErr' ~ tval ~ tprob ;
	println("\\begin{table} \\centering \\begin{tabular}{|l|ccccc|} \\hline");
	println("%-15s","Parameter &", "%-12s", "Coefficient &", "%-14s", "Ana.Std.Error &", "%-14s", "Num.Std.Error &","%-9s", "t-stats &", "%-9s", "p-value* \\\\");
    for (i = 0; i < rows(m_vP); ++i){			  
			println("%-15s", aspar[i], "& ","%-11.4g", mpar[i][0], " &","%-13.4g",mpar[i][1], " &", "%-13.4g", mpar[i][2], "& ","%-8.4g",mpar[i][3]," &","%-8.4f",mpar[i][4], "\\\\" );
			}
	}
	else {
		tval = m_vP./vNumStdErr';
		tprob = tailt(fabs(tval), m_cT-cdfloss);
    	mpar = m_vP ~ vNumStdErr' ~ tval ~ tprob ;		
		println("\\begin{table} \\centering \\begin{tabular}{|l|ccccc|} \\hline");
		println("%-15s","Parameter &", "%-12s", "Coefficient &", "%-11s", "Std.Error &","%-8s", "t-value &", "%-8s", "p-value \\\\");
    	for (i = 0; i < rows(m_vP); ++i){			  
			println("%-15s", aspar[i], " & ","%-11.4g", mpar[i][0], "& ","%-10.4g",mpar[i][1], "& ","%-7.4g",mpar[i][2]," &","%-7.4f",mpar[i][3],"\\\\");
			}
   }
PrintLikStatsLatex();
println("\\hline \\end{tabular} \\end{table}");
}

DySco::MEMOutputPar()
{
	println("");
    println("Table 1. Estimated parameters.");
    decl i;
    decl aspar;
	if (m_iModelClass==MEM) aspar = ParNames();
	//else if (m_iModelClass==LANNE) aspar =LanneParNames();
	decl vNumStdErr = new matrix[m_cP][1];
	vNumStdErr=GetNumStdErrMEM(&vNumStdErr, m_vP,m_cDist);
	
    decl mpar; 
    decl cdfloss = GetcDfLoss();
	decl tval, tprob;
//	tval = m_vP./vNumStdErr';
//	tprob = tailt(fabs(tval), m_cT-cdfloss);


if (m_vStruct[0]==1 && m_vStruct[1]==1 && m_vStruct[3]==0 && m_iTwoComp==0){
		decl cMaxpq = max(m_vStruct[0], m_vStruct[1]);
		decl vAnaStdErr = new matrix[1][m_cP];
		vAnaStdErr = GetAnaStdErr(&vAnaStdErr, m_vP,m_cDist,  m_cT, cMaxpq);
		tval = m_vP./vAnaStdErr';
		tprob = tailt(fabs(tval), m_cT-cdfloss);
    	mpar = m_vP ~ vAnaStdErr'~vNumStdErr' ~ tval ~ tprob ;
		println("%-15s","", "%-12s", "Coefficient", "%-14s", "Ana.Std.Error", "%-14s", "Num.Std.Error","%-9s", "t-stats*", "%-9s", "p-value*");
    	for (i = 0; i < rows(m_vP); ++i){			  
			println("%-15s", aspar[i], " ","%-11.4g", mpar[i][0], " ","%-13.4g",mpar[i][1], " ", "%-13.4g", mpar[i][2], " ","%-8.4g",mpar[i][3]," ","%-8.4f",mpar[i][4]);
			}
		println("-- * Indicates statistic based on analytical standard errors --");
	}
else {
		tval = m_vP./vNumStdErr';
		tprob = tailt(fabs(tval), m_cT-cdfloss);
    	mpar = m_vP ~ vNumStdErr' ~ tval ~ tprob ;		
		println("%-15s","", "%-12s", "Coefficient", "%-11s", "Std.Error","%-8s", "t-value", "%-8s", "p-value");
    	for (i = 0; i < rows(m_vP); ++i){			  
			println("%-15s", aspar[i], " ","%-11.4g", mpar[i][0], " ","%-10.4g",mpar[i][1], " ","%-7.4g",mpar[i][2]," ","%-7.4f",mpar[i][3]);
			}
   }

//
//
//	GetMEMAnaStdErr(vAnaStdErr, const vP, const iDist, const cT, const cMaxpq)
//	//println(m_vNumStdErr);
//	//println(m_mNumHess);
//	//println(m_mNumCov);
////	println(m_vP);
//	println(vNumStdErr);
////	println(tval');
////	println(tprob');
//    mpar = m_vP ~ vNumStdErr' ~ tval ~ tprob ;		
//	println("%-15s","", "%-12s", "Coefficient", "%-11s", "Std.Error","%-8s", "t-value", "%-8s", "p-value");
//    for (i = 0; i < rows(m_vP); ++i){			  
//		println("%-15s", aspar[i], " ","%-11.4g", mpar[i][0], " ","%-10.4g",mpar[i][1], " ","%-7.4g",mpar[i][2]," ","%-7.4f",mpar[i][3]);
//		}
PrintLikStats();
}

/////////////////////////////////////////////////////////////////////////
///////////////////////////// Graphical Output /////////////////////////
////////////////////////////////////////////////////////////////////////
DySco::GetEps()
{ // purpose: get series of epsilons for 1...T // Also gets Scores for 1..T
 // requires: LambdaFilter or MEMLambdaFilter
 GetPredict();
 m_vEps = m_mData[][0]./exp(m_vLambda);
}

DySco::DrawEps(const iPos)
{DrawTMatrix(iPos,m_vEps[0:m_iT2sel]',{"In Sample \\varepsilon_{t}"});
 DrawDensity(iPos+1, m_vEps[0:m_iT2sel]', {"In sample distribution of Eps"});
}

DySco::DrawEpsOoS(const iPos)
{// Purpose: draw out of sample Eps
  DrawTMatrix(iPos,m_vEps[m_iT2sel+1:]',{"Out of Sample \\varepsilon_{t}"});
}

DySco::GetPredict()
{//  Purpose: Take parameters, obtain lambda's obtain predicted observations for t=1...T
// // Both in and out of sample Predictions! // //
decl vLambda = new matrix[rows(m_mData[][0])][1];
if (m_iModelClass==SCORE){
if (m_vStruct[0]>0 && m_vStruct[1]>0){
//	println("m_vStruct" ,m_vStruct);
//	println("m_vP ",m_vP);
	decl vU = new matrix[rows(m_mData[][0])][1];
	if (m_iTwoComp == 0)		LambdaFilter(&vLambda, &vU, m_mData[][0], m_vP, m_vStruct, m_cDist);
	else 				  		LambdaFilter2comp(&vLambda, &vU, m_mData[][0], m_vP, m_vStruct, m_cDist);
	m_vU = vU; delete vU;
	}
else m_vLambda = m_vP[0];
}
else if (m_iModelClass==MEM){
	if (m_iTwoComp==0) 		MEMLambdaFilter(&vLambda, m_mData[][0], m_vP, m_vStruct);
	else 					MEMLambdaFilter2comp(&vLambda, m_mData[][0], m_vP, m_vStruct);

	// get scores
	decl vU;
	decl eps = m_mData[1:][0]./exp(vLambda[1:]);
	decl vShapePara = GetShapePara(m_vP, m_vStruct, m_cDist);
	
	if (m_cDist==EX) 				vU = -1+eps;
	else if (m_cDist==GA)			vU = -vShapePara[0]+eps;
	else if (m_cDist==WBL)			vU = -1+(eps).^(vShapePara[0]);
	else if (m_cDist==GGA)			vU = (eps).^(vShapePara[0]) - vShapePara[1];
	else if (m_cDist==LNORM)		vU = log(eps);
	else if (m_cDist==LLOG)			vU = -1+2*eps.^(vShapePara[0])./(1+eps.^(vShapePara[0]));
	else if (m_cDist==BURR)			vU = -1+(vShapePara[1]+1)*eps.^(vShapePara[0])./(1+eps.^(vShapePara[0]));
	else if (m_cDist==F)			vU = -0.5*vShapePara[0]+ 0.5*(vShapePara[0]+vShapePara[1])*vShapePara[0]*eps./(vShapePara[1]+vShapePara[0]*eps);
	else if (m_cDist==DAG)			vU = -vShapePara[1]+ (vShapePara[1]+1)*eps.^(vShapePara[0])./(1+eps.^(vShapePara[0]));
	else if (m_cDist==GB2)			vU = -vShapePara[1]+(vShapePara[1]+vShapePara[2])*eps.^(vShapePara[0])./(1+eps.^(vShapePara[0]));
	m_vU = vU; delete vU;	
	}

m_vLambda = vLambda;
delete vLambda;


if (m_cDist ==EX) 			m_vMu = exp(m_vLambda);
decl vShapePara = GetShapePara(m_vP, m_vStruct,m_cDist);
if (m_cDist == GA) 	m_vMu = exp(m_vLambda)*vShapePara[0];
else if (m_cDist == WBL) 	m_vMu = exp(m_vLambda).*gammafact(1+1/vShapePara);
else if (m_cDist == LNORM) 	{println("shape para ", vShapePara); m_vMu =
exp(m_vLambda+ 0.5*vShapePara); }
else if (m_cDist == BURR) {
							decl const_mu = gammafact(1+1/vShapePara[0])*gammafact(vShapePara[1]-1/vShapePara[0])/(gammafact(vShapePara[1]));
							m_vMu = const_mu*exp(m_vLambda);
						}
else if (m_cDist == DAG) 	m_vMu = exp(m_vLambda)*gammafact(vShapePara[1]+1/vShapePara[0])*gammafact(1-1/vShapePara[0])/gammafact(vShapePara[1]);
else if (m_cDist == LLOG) 	m_vMu = exp(m_vLambda)*M_PI/(sin(M_PI/vShapePara)*vShapePara);
else if (m_cDist == F) 		m_vMu = exp(m_vLambda)*vShapePara[1]/(vShapePara[1]-2);
else if (m_cDist == GGA) 	m_vMu = exp(m_vLambda)*gammafact(vShapePara[0]+1/vShapePara[1])/gammafact(vShapePara[0]);
else if (m_cDist == GB2) 	m_vMu = exp(m_vLambda)*Beta(vShapePara[1]+1/vShapePara[0],vShapePara[2]-1/vShapePara[0])/Beta(vShapePara[1],vShapePara[2]);
}

DySco::DrawPredict(const iPos)
{// Purpose: draw actual observations vs predicted observations
GetPredict();					
DrawTMatrix(iPos,m_mData[m_iT1sel:m_iT2sel][0]'|m_vMu[m_iT1sel:m_iT2sel]',{"Observations", "Model"});
//DrawTMatrix(iPos,m_mData[600:1000][0]'|m_vMu[600:1000]',{"Observations", "Model"});
decl s = sprint(GetDistriString(),"(", m_vStruct[0],",",m_vStruct[1],") model");
DrawTitle(iPos, s);
}

DySco::DrawPredictOoS(const iPos)
{// Purpose: draw actual observations vs predicted observations out of sample
GetPredict();
decl aux = 0;//strfind(m_asNames, "Range");
DrawTMatrix(iPos,m_mData[m_iT2sel+1:][0]'|m_vMu[m_iT2sel+1:]',{"Observations", "predictions"});
decl s = sprint("Out of Sample for ", GetDistriString(),"(", m_vStruct[0],",",m_vStruct[1],") model");
DrawTitle(iPos, s);
}

DySco::probgengamma(const x, const nu , const beta, const alpha)
{// Purpose: get cfd of generalized gamma distribution
 decl p;
 p = gammafunc((x./beta).^(alpha),nu)./gammafact(nu);
 return p;
}

DySco::probloglogistic(const x, const scale, const beta)
{
 decl p;
 p = 1-(1+(x./scale).^(beta)).^(-1);
 return p;
}

DySco::probburr(const x, const scale, const beta, const alpha)
{
 decl p;
 p = 1- (1+(x./scale).^(beta)).^(-alpha);
 return p;	 
}

DySco::probdag(const x, const scale, const nu, const xi)
{decl p;
 p =  (1+ (x./scale).^(-nu)).^(-xi);
 return p;
}

DySco::probgb2(const x,const scale, const nu, const xi,const varsig)
{decl p;
// wrong, requires hypergemoetric distributions, but not available in OxMetrics
 p = Beta(xi,varsig)^(-1)*betafunc((x./scale).^(nu), xi, varsig);
return p;
}

DySco::GetPit(const vEps, const vShapePara)
{// Purpose: get vector of cdf's required for PIT graphs
 decl x; // return value
 if (m_cDist == EX) 			x = probexp(vEps,1);
 else if (m_cDist == GA) 		x = probgamma(vEps, vShapePara,1);
 else if (m_cDist == WBL) 		x = probweibull(vEps, 1,vShapePara);
 else if (m_cDist == LNORM)		x = probmvn(log(vEps),vShapePara); // use standard normal instead		 problogn(
 else if (m_cDist == GGA) 		x = probgengamma(vEps, vShapePara[1], 1, vShapePara[0]);
 else if (m_cDist == LLOG) 		x = probloglogistic(vEps,1,vShapePara);
 else if (m_cDist == BURR)		x = probburr(vEps, 1, vShapePara[0], vShapePara[1]);
 else if (m_cDist == DAG)		x = probdag(vEps,1,vShapePara[0],vShapePara[1]);
 else if (m_cDist == F)			x = probf(vEps,vShapePara[0], vShapePara[1]);
 else if (m_cDist == GB2)       x = probgb2(vEps,1,vShapePara[0], vShapePara[1],vShapePara[2]);
 return x;
}

DySco::LannePit(const vY, const mu1, const mu2, const gamma1, const gamma2, const ProbPi)
{// purpose: for vector of observations and mu1, mu2, gamma1, gamma2, get PIT for gamma mixture model in Lanne
decl T = rows(vY);
decl z = zeros(T,1);
decl z1 = zeros(T,1);
decl z2 = zeros(T,1);

//z1 = ProbPi*mu1.^(-2*gamma1).*(gammafact(gamma1)-gammafunc(gamma1,mu1.*vY./gamma1))/gammafact(gamma1);
//z2 = (1-ProbPi)*mu2.^(-2*gamma2).*(gammafact(gamma2)-gammafunc(gamma2,mu2.*vY./gamma2))/gammafact(gamma2);
z1 = ProbPi.*(gammafact(gamma1)-gammafunc(gamma1,vY.*gamma1./mu1))/gammafact(gamma1);
z2 = (1-ProbPi).*(gammafact(gamma2)-gammafunc(gamma2,vY.*gamma2./mu2))/gammafact(gamma2);
z = z1+z2;
//println(z[0:10]~z1[0:10]~z2[0:10]);
//println(gammafunc(gamma1,vY[0:10].*gamma1./mu1[0:10]));
//println((gammafact(gamma1)-gammafunc(gamma1,vY[0:10].*gamma1./mu1[0:10]))~(gammafact(gamma1)-gammafunc(gamma1,vY[0:10].*gamma1./mu1[0:10]))/gammafact(gamma1) );
return z;
}
														  
//DySco::DrawLannePit(const iPos)
//{// gets PITs and draws them for LANNE model
//decl mu1 = LanneFilter(m_vPar[0:1+m_vStruct[0]], m_mData, m_vStruct[0]);
//decl mu2 = LanneFilter(m_vPar[2+m_vStruct[0]:3+m_vStruct[0]+m_vStruct[0]], m_mData, m_vStruct[0]);
//
//decl gamma1 = m_vPar[4+2*m_vStruct[0]];
//decl gamma2 = m_vPar[5+2*m_vStruct[0]];
//decl ProbPi = m_vPar[6+2*m_vStruct[0]];
//decl z = LannePit(m_mData, mu1,  mu2, gamma1,gamma2, ProbPi);
//z = sortbyc(z,0);
//decl Raux = range(1/m_cBins,1-1/m_cBins,1/m_cBins);
//DrawHistogram(iPos, countc(z,Raux)');
//decl s = sprint(" $y_{t}$ PITs Gamma Mixture MEM(",m_vStruct[0],",1) model" );
//DrawTitle(iPos,s);
//
////DrawTMatrix(iPos+1,(m_mData./mu1)',{"y/mu1"});
////DrawTMatrix(iPos+2,(m_mData./mu2)',{"y/mu2"});
////DrawTMatrix(iPos+3,gammafunc(gamma1,m_mData.*gamma1./mu1), {"gammafunc"});
////DrawTMatrix(iPos+4, gammafact(gamma1)-gammafunc(gamma1,m_mData.*gamma1./mu1) ,{"differnce gamma func"});
//
//ShowDrawWindow();
//}



DySco::GetScores()
{// purpose: get scores after estimation
// Transform m_vU's into m_vScore depending on distribution
decl vShapePara = GetShapePara(m_vP, m_vStruct, m_cDist);
if (m_cDist == EX)			m_vScore = m_vU+1;
else if (m_cDist ==GA) 		m_vScore = m_vU+vShapePara;
else if (m_cDist == WBL) 	m_vScore = m_vU +1;
else if (m_cDist ==LNORM) 	m_vScore = m_vU;
else if (m_cDist ==LLOG) 	m_vScore = 0.5*(m_vU+1);
else if (m_cDist == BURR)	m_vScore = (m_vU+1)./(vShapePara[1]+1);
else if (m_cDist == DAG) 	m_vScore = (m_vU+vShapePara[1])./(vShapePara[1]+1);
else if (m_cDist == F)		m_vScore = (2*m_vU + vShapePara[0])./(vShapePara[0]+vShapePara[1]);
else if (m_cDist == GGA) 	m_vScore = m_vU+vShapePara[1];
else if (m_cDist == GB2)	m_vScore = (m_vU+vShapePara[1])./(vShapePara[1]+vShapePara[2]);
}


DySco::GetPitScore(const vScore, const vShapePara)
{// Purpose: get vector of cdf's required for PIT graphs
 decl x; // return value
 if (m_cDist == EX) 			x = probexp(vScore,1);
 else if (m_cDist == GA) 		x = probgamma(vScore,vShapePara, 1);
 else if (m_cDist == WBL) 		x = probexp(vScore, 1);
 else if (m_cDist == LNORM)		x = probmvn(vScore,vShapePara); // use standard normal instead
 else if (m_cDist == GGA) 		x = probgamma(vScore, vShapePara[1],1);
 else if (m_cDist == LLOG) 		x = probbeta(vScore,1,1);
 else if (m_cDist == BURR)		x = probbeta(vScore, 1, vShapePara[1]);
 else if (m_cDist == DAG)		x = probbeta(vScore, vShapePara[1] ,1);
 else if (m_cDist == F)			x = probbeta(vScore, 0.5*vShapePara[0], 0.5*vShapePara[1]);
 else if (m_cDist == GB2)		x = probbeta(vScore, vShapePara[1], vShapePara[2]);
 return x;
}


DySco::GetPitConfUpper(const avF, const cUpper)
{
 avF[0] = probbinomial(cUpper, m_cT, 1/m_cBins)-0.99;
 return 1;
}

DySco::GetPitConfLower(const avF, const cLower)
{
 avF[0] =  probbinomial(cLower, m_cT, 1/m_cBins)-0.01;
 return 1;
}

DySco::SetNumberBins(const cBins)
{
return m_cBins = cBins;
}

DySco::DrawPitEps(const iPos, const vEps)
{// auxiliary functions for pits
//println("shape aprea", m_vShapePara);
decl vShapePara = GetShapePara(m_vP, m_vStruct, m_cDist);
decl x = GetPit(vEps, vShapePara);
//DrawTMatrix(9, x',{"PIT RAW"});
//DrawAcf(10, x', {"PIT ACF"}, fabs(log(m_cT)), 1,0,1,2,1);
//DrawDensity(11, x', {"In sample distribution of PIT"}, FALSE, TRUE);
//println("meanc PIT", meanc(x));
//x = sortbyc(x,0);
//decl Raux = range(1/m_cBins,1-1/m_cBins,1/m_cBins);
//DrawHistogram(iPos, countc(x,Raux)');
DrawDensity(iPos, x', {}, FALSE, TRUE, FALSE, FALSE, FALSE,m_cBins);
decl s = sprint(" $\\varepsilon_{t}$ PITs ",GetDistriString(),"(", m_vStruct[0],",",m_vStruct[0],") model" );
DrawTitle(iPos,s);

}


DySco::DrawPitEpsIS(const iPos)
{ GetEps();
decl vShapePara = GetShapePara(m_vP, m_vStruct, m_cDist);
decl x = GetPit(m_vEps[m_iT1sel:m_iT2sel], vShapePara);
DrawDensity(iPos, x', {}, FALSE, TRUE, FALSE, FALSE, FALSE,m_cBins);
decl s = sprint(" $\\varepsilon_{t}$ PITs ",GetDistriString(),"(", m_vStruct[0],",",m_vStruct[0],") model" );
DrawTitle(iPos,s);
}

DySco::DrawPitEpsOoS2(const iPos)
{ GetEps();
decl vShapePara = GetShapePara(m_vP, m_vStruct, m_cDist);
decl x = GetPit(m_vEps[m_iT2sel:], vShapePara);
DrawDensity(iPos, x', {}, FALSE, TRUE, FALSE, FALSE, FALSE,m_cBins);
decl s = sprint(" $\\varepsilon_{t}$ PITs ",GetDistriString(),"(", m_vStruct[0],",",m_vStruct[0],") model" );
DrawTitle(iPos,s);
}

DySco::DrawPitEpsOoS(const iPos, const vEps)
{// auxiliary functions for pits
decl vShapePara = GetShapePara(m_vP, m_vStruct, m_cDist);
decl x = GetPit(vEps, vShapePara);
DrawDensity(iPos, x', {}, FALSE, TRUE, FALSE, FALSE, FALSE,m_cBins);
//x = sortbyc(x,0);
//decl Raux = range(1/m_cBins,1-1/m_cBins,1/m_cBins);
//DrawHistogram(iPos, countc(x,Raux)');
decl s = sprint("Out of Sample $\\varepsilon_{t}$ PITs ",GetDistriString(),"(", m_vStruct[0],",",m_vStruct[1],") model" );
DrawTitle(iPos,s);
}

DySco::DrawPitScore(const iPos, const vScore)
{  decl vShapePara = GetShapePara(m_vP, m_vStruct, m_cDist);
decl x = GetPitScore(vScore, vShapePara);
DrawDensity(iPos, x', {}, FALSE, TRUE, FALSE, FALSE, FALSE,m_cBins);
//x = sortbyc(x,0);
//decl Raux = range(1/m_cBins,1-1/m_cBins,1/m_cBins);
//DrawHistogram(iPos, countc(x,Raux)');
decl s = sprint(" $u_{t}$ PITs ",GetDistriString(),"(", m_vStruct[0],",",m_vStruct[1],") model" );
DrawTitle(iPos,s);
}

DySco::DrawPitScoreOoS(const iPos, const vScore)
{ decl vShapePara = GetShapePara(m_vP, m_vStruct, m_cDist);
decl x = GetPitScore(vScore, vShapePara);
DrawDensity(iPos, x', {}, FALSE, TRUE, FALSE, FALSE, FALSE,m_cBins);
//x = sortbyc(x,0);
//decl Raux = range(1/m_cBins,1-1/m_cBins,1/m_cBins);
//DrawHistogram(iPos, countc(x,Raux)');
decl s = sprint("Out of Sample $u_{t}$ PITs ",GetDistriString(),"(", m_vStruct[0],",",m_vStruct[1],") model" );
DrawTitle(iPos,s);
}

DySco::LBQStatPitEps(const vEps, const vLag)
{ decl vShapePara = GetShapePara(m_vP, m_vStruct, m_cDist);
decl x = GetPit(vEps, vShapePara);
LBQTest(x, vLag);
}

DySco::QQPlot(const vX, const cDist, const vShapePara, const iPos)
{
decl vXSorted, vEmpQ, vTheoQ;
 vXSorted = sortbyc(vX,0);
 vEmpQ = quantilec(vXSorted, range(0,1,(1/m_cT))); //range(0,1,(1/m_cT))); // compute 0.01, 0.02. ... 1 quantiles

 if (cDist == EX)  			vTheoQ = probexp(vEmpQ,1);
 else if (cDist == GA)		vTheoQ = probgamma(vEmpQ, vShapePara[0],1);
 else if (cDist == WBL)		vTheoQ = probweibull(vEmpQ,1, vShapePara);
 else if (cDist == LNORM)	vTheoQ = probmvn(log(vEmpQ),vShapePara);
 else if (cDist == GGA)		vTheoQ = probgengamma(vEmpQ, vShapePara[1], 1, vShapePara[0]);
 else if (cDist == LLOG)	vTheoQ = probloglogistic(vEmpQ,1,vShapePara);
 else if (cDist == BURR)	vTheoQ = probburr(vEmpQ, 1, vShapePara[0], vShapePara[1]);
 else if (cDist == DAG)		vTheoQ = probdag(vEmpQ,1,vShapePara[0],vShapePara[1]);
 else if (cDist == F)		vTheoQ = probf(vEmpQ, vShapePara[0], vShapePara[1]);
 else if (cDist == GB2)		vTheoQ = probgb2(vEmpQ,1,vShapePara[0],vShapePara[1],vShapePara[2]);

 decl x = range(0,1,(1/m_cT));
 DrawXMatrix(iPos, vTheoQ' |x, {}, x,{});
}

DySco::QQPlotDrawEps(const iPos, const vEps)
{ decl vShapePara = GetShapePara(m_vP, m_vStruct, m_cDist);
 QQPlot(vEps, m_cDist, vShapePara, iPos);
 decl s = sprint("$\\varepsilon_{t}$ PP Plot ", GetDistriString(), "(",m_vStruct[0],",",m_vStruct[1],") model");
 DrawTitle(iPos,s);
 DrawLegend(iPos,0,0,TRUE);
}
DySco::QQPlotDrawEpsOoS(const iPos, const vEps)
{ decl vShapePara = GetShapePara(m_vP, m_vStruct, m_cDist);
 QQPlot(vEps, m_cDist, vShapePara, iPos);
 decl s = sprint("Out of Sample $\\varepsilon_{t}$ PP Plot ", GetDistriString(), "(",m_vStruct[0],",",m_vStruct[1],") model");
 DrawTitle(iPos,s);
 DrawLegend(iPos,0,0,TRUE);
}

DySco::QQPlotDrawScore(const iPos, const vU)
{decl vShapePara = GetShapePara(m_vP, m_vStruct, m_cDist);
 QQPlotScores(iPos, m_cDist, vShapePara, vU);
 decl s = sprint("$u_{t}$ PP Plot ", GetDistriString(), "(",m_vStruct[0],",",m_vStruct[1],") model");
 DrawTitle(iPos,s);
 DrawLegend(iPos,0,0,TRUE);
}

DySco::QQPlotDrawScoreOoS(const iPos, const vU)
{  decl vShapePara = GetShapePara(m_vP, m_vStruct, m_cDist);
 QQPlotScores(iPos, m_cDist, vShapePara, vU);
 decl s = sprint("Out of Sample $u_{t}$ PP Plot ", GetDistriString(), "(",m_vStruct[0],",",m_vStruct[1],") model");
 DrawTitle(iPos,s);
 DrawLegend(iPos,0,0,TRUE);
}

DySco::QQPlotScores(const iPos, const cDist, const vShapePara, const vScores)
{
 decl vXSorted, vEmpQ, vTheoQ;
 vXSorted = sortbyc(vScores,0);
 vEmpQ = quantilec(vXSorted, range(0,1,(1/m_cT))); //range(0,1,(1/m_cT))); // compute 0.01, 0.02. ... 1 quantiles

 if (cDist == EX)  			vTheoQ = probexp(vEmpQ,1);
 else if (cDist == GA)		vTheoQ = probgamma(vEmpQ, vShapePara[0],1);
 else if (cDist == WBL)		vTheoQ = probexp(vEmpQ,1);
 else if (cDist == LNORM)	vTheoQ = probmvn(vEmpQ,vShapePara);
 else if (cDist == LLOG)	vTheoQ = probbeta(vEmpQ,1,1);
 else if (cDist == BURR)	vTheoQ = probbeta(vEmpQ, 1, vShapePara[1]);
 else if (cDist == DAG) 	vTheoQ = probbeta(vEmpQ, vShapePara[1], 1);
 else if (cDist == F)		vTheoQ = probbeta(vEmpQ, 0.5*vShapePara[0], 0.5*vShapePara[1]);
 else if (cDist == GGA)	    vTheoQ = probgamma(vEmpQ,vShapePara[1],1);
 else if (cDist == GB2)		vTheoQ = probbeta(vEmpQ,vShapePara[1],vShapePara[2]);

 decl x = range(0,1,(1/m_cT));
 DrawXMatrix(iPos, vTheoQ' |x, {}, x,{});
}

DySco::QQPlotPits(const iPos,const vScore)
{
decl vShapePara = GetShapePara(m_vP, m_vStruct, m_cDist);

decl x = GetPitScore(vScore, vShapePara);
decl xGaussian = quann(x);
DrawQQ(iPos, xGaussian', {"PIT QQ Plot"},QQ_N_SE,0,1);
//ShowDrawWindow();
}

DySco::QQPlotPits2(const iPos)
{
decl vShapePara = GetShapePara(m_vP, m_vStruct, m_cDist);

decl x = GetPitScore(m_vScore, vShapePara);
decl xGaussian = quann(x);
DrawQQ(iPos, xGaussian', {"PIT QQ Plot"},QQ_N_SE,0,1);
ShowDrawWindow();
}


DySco::DrawEpsACF(const iPos)
{
 DrawAcf(iPos, m_vEps', {}, fabs(log(m_cT)), 1,0,1,2,1);
 decl s = sprint("$\\varepsilon_{t}$ ACF ", GetDistriString(), "(",m_vStruct[0],",",m_vStruct[1],") model");
 DrawTitle(iPos,s);
 DrawLegend(iPos,0,0,TRUE);
 //DrawAxis(iPos,0,ANCHOR_USER,minc(acf(m_vEps,fabs(log(m_cT)))), maxc(acf(m_vEps,fabs(log(m_cT)))),0,5,1,0);
 //DrawAxis(iPos,1,ANCHOR_USER,0,fabs(log(m_cT)),0,0.1,0.01,0);
// DrawAxisAuto(iPos,0);
}

DySco::DrawEpsACFOoS(const iPos)
{
 DrawAcf(iPos, m_vEps[m_iT2sel+1:]', {}, fabs(log(m_cT)), 1,0,1,2,1);
  decl s = sprint("Out of Sample $\\varepsilon_{t}$ ACF ", GetDistriString(), "(",m_vStruct[0],",",m_vStruct[1],") model");
 DrawTitle(iPos,s);
 DrawLegend(iPos,0,0,TRUE);
}

DySco::DrawScoreACF(const iPos)
{
 DrawAcf(iPos, m_vU[0:m_iT2sel]', {}, fabs(log(m_cT)), 1,0,1,2,1);
 decl s = sprint("$u_{t}$ ACF ", GetDistriString(), "(",m_vStruct[0],",",m_vStruct[1],") model");
 DrawTitle(iPos,s);
 DrawLegend(iPos,0,0,TRUE);
}

DySco::DrawScoreACFOoS(const iPos)
{
 DrawAcf(iPos, m_vU[m_iT2sel+1:]', {}, fabs(log(m_cT)), 1,0,1,2,1);
 decl s = sprint("Out of Sample $u_{t}$ ACF ", GetDistriString(), "(",m_vStruct[0],",",m_vStruct[1],") model");
 DrawTitle(iPos,s);
 DrawLegend(iPos,0,0,TRUE);
}


DySco::DrawEpsSqACF(const iPos)
{	
 DrawAcf(iPos, (m_vEps[0:m_iT2sel].^(2))', {"$\\varepsilon_{t}^{2}$"}, fabs(log(m_cT)), 1,0,1,2,1);
  decl s = sprint("$\\varepsilon_{t}^{2}$ ACF ", GetDistriString(), "(",m_vStruct[0],",",m_vStruct[1],") model");
 DrawTitle(iPos,s);
 DrawLegend(iPos,0,0,TRUE);
}

DySco::DrawEpsSqACFOoS(const iPos)
{	
 DrawAcf(iPos, (m_vEps[m_iT2sel+1:].^(2))', {"Out of Sample $\\varepsilon_{t}^{2}$"}, fabs(log(m_cT)), 1,0,1,2,1);
  decl s = sprint("Out of Sample $\\varepsilon_{t}^{2}$ ACF ", GetDistriString(), "(",m_vStruct[0],",",m_vStruct[1],") model");
 DrawTitle(iPos,s);
 DrawLegend(iPos,0,0,TRUE); 
}

DySco::DrawScoreSqACF(const iPos)
{
 DrawAcf(iPos, (m_vU[0:m_iT2sel].^(2))', {"$u_{t}^{2}$"}, fabs(log(m_cT)), 1,0,1,2,1);
   decl s = sprint("$u_{t}^{2}$ ACF ", GetDistriString(), "(",m_vStruct[0],",",m_vStruct[1],") model");
 DrawTitle(iPos,s);
 DrawLegend(iPos,0,0,TRUE);
}

DySco::DrawScoreSqACFOoS(const iPos)
{
 DrawAcf(iPos, (m_vU[m_iT2sel+1:].^(2))', {"Out of Sample $u_{t}^{2}$"}, fabs(log(m_cT)), 1,0,1,2,1);
 decl s = sprint("Out of Sample $u_{t}^{2}$ ACF ", GetDistriString(), "(",m_vStruct[0],",",m_vStruct[1],") model");
 DrawTitle(iPos,s);
 DrawLegend(iPos,0,0,TRUE);
}

DySco::LBQTest(const ma, const vLag)
{// LB Q test
 // Input:		ma: Tx1 series of inputs
 //				vLag: vector of lags to be tested

 if (vLag != <>){
 decl T = sizer(ma);
 decl Iter = sizer(vLag);
 decl vACF=<>,q,marg, vT, vLagAux;

format(200);
 println("\t \t Lag \t \t QStat \t  p-Value")	;
 for( decl ant=0; ant<Iter;++ant){
 	vT = T*ones(1,vLag[ant]);
	vLagAux = range(1,vLag[ant]);
	vACF = acf(ma, vLag[ant]);
	q = T*(T+2)*sumc(vACF[1:vLag[ant]].^(2)./(vT-vLagAux)')	;
	marg = tailchi(q, vLag[ant]);
	println(vLag[ant]~q~marg);
 }
 println("H0 : No serial correlation ==> Accept H0 when prob. is High [Q < Chisq(lag)]");
 return 1;
 }
 else oxwarning("LBQTest requires vector of lags as input"); 		  
}


DySco::LBQTest2(const vLag)
{// LB Q test
 // Input:		ma: Tx1 series of inputs  default is m_vEps[0:m_iT2sel]
 //				vLag: vector of lags to be tested
 decl  ma = m_vEps[0:m_iT2sel];
 if (vLag != <>){
 decl T = sizer(ma);
 decl Iter = sizer(vLag);
 decl vACF=<>,q,marg, vT, vLagAux;

format(200);
 println("\t \t Lag \t \t QStat \t  p-Value")	;
 for( decl ant=0; ant<Iter;++ant){
 	vT = T*ones(1,vLag[ant]);
	vLagAux = range(1,vLag[ant]);
	vACF = acf(ma, vLag[ant]);
	q = T*(T+2)*sumc(vACF[1:vLag[ant]].^(2)./(vT-vLagAux)')	;
	marg = tailchi(q, vLag[ant]);
	println(vLag[ant]~q~marg);
 }
 println("H0 : No serial correlation ==> Accept H0 when prob. is High [Q < Chisq(lag)]");
 return 1;
 }
 else oxwarning("LBQTest requires vector of lags as input"); 		  
}


DySco::LBQTestScores(const vLag)
{// LB Q test
 // Input:		ma: Tx1 series of inputs  default is m_vEps[0:m_iT2sel]
 //				vLag: vector of lags to be tested
 decl  ma = m_vScore[0:m_iT2sel];
 if (vLag != <>){
 decl T = sizer(ma);
 decl Iter = sizer(vLag);
 decl vACF=<>,q,marg, vT, vLagAux;

format(200);
 println("\t \t Lag \t \t QStat \t  p-Value")	;
 for( decl ant=0; ant<Iter;++ant){
 	vT = T*ones(1,vLag[ant]);
	vLagAux = range(1,vLag[ant]);
	vACF = acf(ma, vLag[ant]);
	q = T*(T+2)*sumc(vACF[1:vLag[ant]].^(2)./(vT-vLagAux)')	;
	marg = tailchi(q, vLag[ant]);
	println(vLag[ant]~q~marg);
 }
 println("H0 : No serial correlation ==> Accept H0 when prob. is High [Q < Chisq(lag)]");
 return 1;
 }
 else oxwarning("LBQTest requires vector of lags as input"); 		  
}



DySco::APGT(const cd, const ng, const np) 
/*   
**  Purpose: Computes adjusted Pearson Chi-square Goodness-of-fit test (Vlaar, 93) 
** 
**  Input  : cd = values of the cumulative distribution function 
**           ng = classification in "ng" groups  
**           np = model with "np" parameters 
** 
*/ 
{ 
	if (ng != <>) 
	{ 
		decl x2,v,num1,num,enu,t,i,npp,k,stat,chi1,chi2,mx; 
		t=rows(cd); 
	    println("Adjusted Pearson Chi-square Goodness-of-fit test"); 
		k=sizer(ng); 
		stat = <>; chi1 = <>;	chi2 = <>; 		 
		for (i = 0; i < k; i++)   
		{ 
	    	v = range(1/ng[i], 1, 1/ng[i]); 
	    	num1 = countc(cd[][0], v); 
			num  = num1[0:sizer(num1)-2]; 
			npp = limits(np|ng[i]-2)[0]; 
			enu = t/ng[i]; 
	    	x2 = ((num-enu).^2)/enu; 
			stat |= sumc(x2); 
			chi1 |= tailchi(sumc(x2),ng[i]-1); 
			chi2 |= tailchi(sumc(x2),ng[i]-1-npp); 
		} 
		mx = ng ~ stat ~ chi1 ~ chi2; 
	    print("%c",{"Cells","Statistic","P-Value(lag-1)","P-Value(lag-k-1)"}, 
			  "%cf", {"%6.0f", "%12.4f", "%17.6f", "%20.6f"}, mx); 
		println("Rem.: k = # estimated parameters");  
		println("---------------"); 
	}		   
	return 1; 
} 


DySco::FEM(const forc, const obs, const bOoS)  
/*   
**  Purpose: computes Forecast Error Measures 
** 
**  Format : FEM(const forc, const obs);
**
**	input  : forc = scale forecast
**			 obs  = observed series and variance
**			bOoS  = boolean: 0: In sample, 1: out of sample (add predictive likelihood)
** 
**  Output : MSE = Mean Squared Error 
**           MAE = Mean Absolute Error 
**           RMSE = Root Mean Squared Error 
**           MAPE = Mean Absolute Percentage Error 
**           AMAPE = Adjusted Mean Absolute Percentage Error 
**           THEIL = Theil Inequality Coefficient 
**           LL = logarithmic loss function    
**  
*/ 
{  
	decl i,MSE,MAE,MAPE,AMAPE,PCS,LL,RMSE,THEIL, maxc, maxc2,a,b; 
	decl res1, res2; 
	decl xp = new matrix[2][10];
	decl Tforc=rows(forc); 
	res1 = (obs[][0] - forc[][0]);   
	res2 = (obs[][0] + forc[][0]);
	
	MSE = (sumc(res1.^2)/rows(res1)); 
	MAE = (sumc(fabs(res1))/rows(res1)); 
	RMSE = sqrt((sumc(res1.^2)/rows(res1))); 
	MAPE = (sumc(fabs(res1./obs))/rows(res1)); 
	AMAPE = (sumc(fabs(res1./res2))/rows(res1));
	LL = (sumc(log((fabs(forc./obs))).^2)/rows(res1)); 
	THEIL = sqrt((sumc(res1.^2))*(1/(1+Tforc)))/ (sqrt((sumc(obs.^2))*(1/(1+Tforc))) + sqrt((sumc(forc.^2))*(1/(1+Tforc)))); 
   
	xp = MSE~MAE~RMSE~MAPE~AMAPE~THEIL~LL;
	
	if (bOoS==1){ // add predictive likelihood to evaluation measures
		decl OoSLik;  decl vLikContr= new matrix[1][m_cT-m_iT2sel];
		if (m_cDist == EX) 			EXLikVal(&OoSLik, m_vP,  m_iT2sel); 
		else if (m_cDist == GA) 	OoSLik = GALikVal(m_vP, m_iT2sel);
		else if (m_cDist == WBL)	OoSLik = WBLLikVal(m_vP, m_iT2sel);
		else if (m_cDist == LNORM)	OoSLik = LNORMLikVal(m_vP, m_iT2sel);
		else if (m_cDist == GGA) 	OoSLik = GGALikVal(m_vP, m_iT2sel);
		else if (m_cDist == LLOG)	OoSLik = LLOGLikVal(m_vP, m_iT2sel);
		else if (m_cDist == BURR)	OoSLik = BURRLikVal(m_vP, m_iT2sel);
		else if (m_cDist == DAG)	OoSLik = DAGLikVal(m_vP, m_iT2sel);
		else if (m_cDist == F) 		OoSLik = FLikVal(m_vP, m_iT2sel);
		else if (m_cDist == GB2)	OoSLik = GB2LikVal(m_vP,m_iT2sel);
	//println("OoSLik ", OoSLik);

	xp ~= OoSLik;
	println("%r", {"Mean Squared Error(MSE)","Mean Absolute Error(MAE)", 
	               "Root Mean Squared Error(RMSE)", 
	               "Mean Absolute Percentage Error(MAPE)", "Adjusted Mean Absolute Percentage Error(AMAPE)", 
	               "Theil Inequality Coefficient(TIC)", "Logarithmic Loss Function(LL)", "Predictive Likelihood"},  
				    xp');
	}
	else if(bOoS==0) {
		println("%r", {"Mean Squared Error(MSE)","Mean Absolute Error(MAE)", 
	               "Root Mean Squared Error(RMSE)", 
	               "Mean Absolute Percentage Error(MAPE)", "Adjusted Mean Absolute Percentage Error(AMAPE)", 
	               "Theil Inequality Coefficient(TIC)", "Logarithmic Loss Function(LL)"},  
				    xp');
	}

//	println("---------------");
//	print("Forecast Evaluation Measures"); 
//	println("%r", {"Mean Squared Error(MSE)","Mean Absolute Error(MAE)", 
//	               "Root Mean Squared Error(RMSE)", 
//	               "Mean Absolute Percentage Error(MAPE)", "Adjusted Mean Absolute Percentage Error(AMAPE)", 
//	               "Theil Inequality Coefficient(TIC)", "Logarithmic Loss Function(LL)", "Predictive Likelihood"},  
//				    xp');
return 1; 
}



DySco::FEM2()  
/*   
**  Purpose: computes Forecast Error Measures 
** 
**  Format : FEM(const forc, const obs);
**
**	input  : forc = scale forecast
**			 obs  = observed series and variance
**			bOoS  = boolean: 0: In sample, 1: out of sample (add predictive likelihood)
** 
**  Output : MSE = Mean Squared Error 
**           MAE = Mean Absolute Error 
**           RMSE = Root Mean Squared Error 
**           MAPE = Mean Absolute Percentage Error 
**           AMAPE = Adjusted Mean Absolute Percentage Error 
**           THEIL = Theil Inequality Coefficient 
**           LL = logarithmic loss function    
**  
*/ 
{
  	decl obs = m_mData[0:m_iT2sel];
	decl forc = m_vMu[0:m_iT2sel];
	decl bOoS = 0;
	decl i,MSE,MAE,MAPE,AMAPE,PCS,LL,RMSE,THEIL, maxc, maxc2,a,b; 
	decl res1, res2; 
	decl xp = new matrix[2][10];
	decl Tforc=rows(forc); 
	res1 = (obs[][0] - forc[][0]);   
	res2 = (obs[][0] + forc[][0]);
	
	MSE = (sumc(res1.^2)/rows(res1)); 
	MAE = (sumc(fabs(res1))/rows(res1)); 
	RMSE = sqrt((sumc(res1.^2)/rows(res1))); 
	MAPE = (sumc(fabs(res1./obs))/rows(res1)); 
	AMAPE = (sumc(fabs(res1./res2))/rows(res1));
	LL = (sumc(log((fabs(forc./obs))).^2)/rows(res1)); 
	THEIL = sqrt((sumc(res1.^2))*(1/(1+Tforc)))/ (sqrt((sumc(obs.^2))*(1/(1+Tforc))) + sqrt((sumc(forc.^2))*(1/(1+Tforc)))); 
   
	xp = MSE~MAE~RMSE~MAPE~AMAPE~THEIL~LL;
	
	if (bOoS==1){ // add predictive likelihood to evaluation measures
		decl OoSLik;  decl vLikContr= new matrix[1][m_cT-m_iT2sel];
		if (m_cDist == EX) 			EXLikVal(&OoSLik, m_vP,  m_iT2sel); 
		else if (m_cDist == GA) 	OoSLik = GALikVal(m_vP, m_iT2sel);
		else if (m_cDist == WBL)	OoSLik = WBLLikVal(m_vP, m_iT2sel);
		else if (m_cDist == LNORM)	OoSLik = LNORMLikVal(m_vP, m_iT2sel);
		else if (m_cDist == GGA) 	OoSLik = GGALikVal(m_vP, m_iT2sel);
		else if (m_cDist == LLOG)	OoSLik = LLOGLikVal(m_vP, m_iT2sel);
		else if (m_cDist == BURR)	OoSLik = BURRLikVal(m_vP, m_iT2sel);
		else if (m_cDist == DAG)	OoSLik = DAGLikVal(m_vP, m_iT2sel);
		else if (m_cDist == F) 		OoSLik = FLikVal(m_vP, m_iT2sel);
		else if (m_cDist == GB2)    OoSLik = GB2LikVal(m_vP,m_iT2sel);
	//println("OoSLik ", OoSLik);

	xp ~= OoSLik;
	println("%r", {"Mean Squared Error(MSE)","Mean Absolute Error(MAE)", 
	               "Root Mean Squared Error(RMSE)", 
	               "Mean Absolute Percentage Error(MAPE)", "Adjusted Mean Absolute Percentage Error(AMAPE)", 
	               "Theil Inequality Coefficient(TIC)", "Logarithmic Loss Function(LL)", "Predictive Likelihood"},  
				    xp');
	}
	else if(bOoS==0) {
		println("%r", {"Mean Squared Error(MSE)","Mean Absolute Error(MAE)", 
	               "Root Mean Squared Error(RMSE)", 
	               "Mean Absolute Percentage Error(MAPE)", "Adjusted Mean Absolute Percentage Error(AMAPE)", 
	               "Theil Inequality Coefficient(TIC)", "Logarithmic Loss Function(LL)"},  
				    xp');
	}

//	println("---------------");
//	print("Forecast Evaluation Measures"); 
//	println("%r", {"Mean Squared Error(MSE)","Mean Absolute Error(MAE)", 
//	               "Root Mean Squared Error(RMSE)", 
//	               "Mean Absolute Percentage Error(MAPE)", "Adjusted Mean Absolute Percentage Error(AMAPE)", 
//	               "Theil Inequality Coefficient(TIC)", "Logarithmic Loss Function(LL)", "Predictive Likelihood"},  
//				    xp');

println(m_mData~m_vEps~m_vScore);
return 1; 
}


DySco::PrintSampleStats()
{

decl aux, nxt=0;
//aux = 0;//strfind(m_asNames, "Range");
println("--- Sample vs Model Characteristics ---");
println("%-15s", "", "%-15s", "Mean", "%-15s", "Variance", "%-15s", "Skew", "%-12s", "Kurt");
decl aux1 = moments(m_mData[0:m_iT2sel][0],4);
decl aux2 = moments(m_vMu[0:m_iT2sel],4);
println("%-15s", "Observations", "%-15.2g", aux1[1], "%-15.2g", aux1[2], "%-12.4g", aux1[3], "%-12.4g",aux1[4] );
println("%-15s", "Predictions", "%-15.2g", aux2[1], "%-15.2g", aux2[2], "%-12.4g", aux2[3], "%-12.4g",aux2[4] );
if (m_iModelClass==SCORE){
	decl aux3 = moments(m_vU[0:m_iT2sel],4);
	println("%-15s", "Scores", "%-15.2g", aux3[1], "%-15.2g", aux3[2], "%-12.4g", aux3[3], "%-12.4g",aux3[4] );
	}
if(m_iT2sel+1!=m_cT){
decl aux4 = moments(m_mData[m_iT2sel+1:][0],4);
decl aux5 = moments(m_vMu[m_iT2sel+1:],4);
println("");
println("--- Out of Sample");
println("%-15s", "Observations", "%-15.2g", aux4[1], "%-15.2g", aux4[2], "%-12.4g", aux4[3], "%-12.4g",aux4[4] );
println("%-15s", "Forecasts", "%-15.2g", aux5[1], "%-15.2g", aux5[2], "%-12.4g", aux5[3], "%-12.4g",aux5[4] );
if (m_iModelClass==SCORE){
	decl aux6 = moments(m_vU[m_iT2sel+1:],4);
	println("%-15s", "Scores", "%-15.2g", aux6[1], "%-15.2g", aux6[2], "%-12.4g", aux6[3], "%-12.4g",aux6[4] );
	}
}
		
}

/////////////////////////////////////////////////////
////////////////// SIMULATION STUFF ////////////////
///////////////////////////////////////////////////

/// Note: works for DCS(1,1) model only at the moment
DySco::SetTruePar(const vBeta)
{    m_vTruePar=vBeta;
//    SetStartPar(vBeta);
}


DySco::GASimulate(mCoefVal, const vElapsedTime, const cT, const cRep, const vPtrue, const sFilename, const sFileNameData, const iPos)
{decl i=0,cfailed=0, dfunc=0,ir=0, time=0;
for (i = cfailed = 0;i<cRep;) { //; i++){
	println("BFGS iteration step: ", i);
	decl vY = new matrix[cT+100][1];
	decl cY,cLambda=vPtrue[0], cU=0, vPest;
	for(decl t=0; t<cT+100;t++){
	cLambda =vPtrue[0]*(1-vPtrue[1]) + vPtrue[1]*cLambda + vPtrue[2]*cU;
	cY  = rangamma(1,1,vPtrue[3],exp(-cLambda));
	vY[t][0]=cY;
	cU = -vPtrue[3]+cY/exp(cLambda);
	}
	m_mData = vY[100:];
	delete vY,cLambda,cY,cU;	
	vPest = vPtrue; // starting values	
	time = timer();
	ir = MaxBFGS(GALik, &vPest, &dfunc, 0,0);
	decl b = vPest[1]^(2)-2*vPest[1]*vPest[2]*vPest[3]+vPest[2]^(2)*(1+vPest[3])*vPest[3];
	if (b>=1){oxwarning("b>1"); println("b is:" ,b ~vPest');continue;}

	if (ir!=0&& ir!=1){
		++cfailed;
		continue;
	}
	vElapsedTime[0][i] = (timer()-time)/100;
	mCoefVal[0][i][] = vPest'~ir;
	++i;
	delete vPest, time, m_mData, ir;
}
Report(MaxBFGS,vElapsedTime[0],mCoefVal[0],vPtrue,cRep,cT,cfailed);
//Save result in database
decl dbase = new Database();  
dbase.Create(1,1,1,m_cRep,1);
dbase.Append(mCoefVal[0][][]~vElapsedTime[0], {"Const","AR", "MA", "GA", "time" });
dbase.SaveXls(sFilename);

Draw(iPos,sFilename,mCoefVal[0][][:3]);

}


DySco::WBLSimulate(const mCoefVal, const vElapsedTime, const cT, const cRep, const vPtrue, const sFilename, const sFileNameData, const iPos)
{
decl i=0,cfailed=0, dfunc=0,ir=0, time=0;
for (i = cfailed = 0;i<cRep;) { //; i++){
	println("BFGS iteration step: ", i);
	decl vY = new matrix[cT+100][1];
	decl cY,cLambda=vPtrue[0], cU=0, vPest;
	for(decl t=0; t<cT+100;t++){
	cLambda =vPtrue[0]*(1-vPtrue[1]) + vPtrue[1]*cLambda + vPtrue[2]*cU;
	cY  = ranweibull(1,1,exp(-cLambda*vPtrue[3]), vPtrue[3]);
	vY[t][0] = cY;
	cU = (cY/exp(cLambda))^(vPtrue[3])-1;
	}
	m_mData = vY[100:];
	delete vY,cLambda,cY,cU;	
	vPest = vPtrue; // starting values	
//	MaxControl(2000, 0,TRUE);
	time = timer();
	ir = MaxBFGS(WBLLik, &vPest, &dfunc, 0, 0);
	decl b = vPest[1]^(2)-2*vPest[1]*vPest[2]*vPest[3]+2*vPest[2]^(2)*vPest[3]^(2);
	if (b>=1){oxwarning("b>1"); println("b is:" ,b ~vPest');continue;}

	if (ir!=0){// && ir!=1){
		++cfailed;
		continue;
	}
	vElapsedTime[0][i] = (timer()-time)/100;
	mCoefVal[0][i][] = vPest'~ir;
	++i;
	delete vPest, time, m_mData, ir;
}
Report(MaxBFGS,vElapsedTime[0],mCoefVal[0],vPtrue,cRep,cT,cfailed);
//Save result in database
decl dbase = new Database();  
dbase.Create(1,1,1,m_cRep,1);
dbase.Append(mCoefVal[0][][]~vElapsedTime[0], {"Const","AR", "MA", "GA", "time" });
dbase.SaveXls(sFilename);

Draw(iPos,sFilename,mCoefVal[0][][:3]);

}

DySco::LLOGSimulate(const mCoefVal, const vElapsedTime, const cT, const cRep, const vPtrue, const sFilename, const sFileNameData, const iPos)
{
decl i=0,cfailed=0, dfunc=0,ir=0, time=0;
for (i = cfailed = 0;i<cRep;) { //; i++){
	println("BFGS iteration step: ", i);
	decl vY = new matrix[cT+100][1];
	decl cY,cLambda=vPtrue[0], cU=0, vPest;
	for(decl t=0; t<cT+100;t++){
	cLambda =vPtrue[0]*(1-vPtrue[1]) + vPtrue[1]*cLambda + vPtrue[2]*cU;
	decl y= ranu(1,1);
	cY = (y/(1-y))^(1/vPtrue[3])*exp(cLambda);
	vY[t][0] = cY;
	cU = -1+2*cY^(vPtrue[3])*exp(-vPtrue[3]*cLambda)/(1+cY^(vPtrue[3])*exp(-vPtrue[3]*cLambda));
	}
	m_mData = vY[100:];
	delete vY,cLambda,cY,cU;	
	vPest = vPtrue; // starting values	
//	MaxControl(2000, 0,TRUE);
	time = timer();
	ir = MaxBFGS(LLOGLik, &vPest, &dfunc, 0, 0);
	decl b = vPest[1]^(2)-2*vPest[1]*vPest[2]*vPest[3]/3 + 2*vPest[2]^(2)*vPest[3]^(2)/15;
	if (b>=1){oxwarning("b>1"); println("b is:" ,b ~vPest');continue;}

	if (ir!=0 && ir!=1){
		++cfailed;
		continue;
	}
	vElapsedTime[0][i] = (timer()-time)/100;
	mCoefVal[0][i][] = vPest'~ir;
	++i;
	delete vPest, time, m_mData, ir;
}
Report(MaxBFGS,vElapsedTime[0],mCoefVal[0],vPtrue,cRep,cT,cfailed);
//Save result in database
decl dbase = new Database();  
dbase.Create(1,1,1,m_cRep,1);
dbase.Append(mCoefVal[0][][]~vElapsedTime[0], {"Const","AR", "MA", "NU", "time" });
dbase.SaveXls(sFilename);
//data.SaveXls(sFileNameData);

Draw(iPos,sFilename,mCoefVal[0][][:3]);

}

DySco::BURRSimulate(const mCoefVal, const vElapsedTime, const cT, const cRep, const vPtrue, const sFilename, const sFileNameData, const iPos)
{
decl i,cfailed, dfunc=0,ir=0, time=0;
for (i = cfailed = 0;i<cRep;) { //; i++){
	println("BFGS iteration step: ", i);
	decl vY = new matrix[cT+100][1];
	decl cY,cLambda=vPtrue[0], cU=0, vPest;
	for(decl t=0; t<cT+100;t++){
	cLambda =vPtrue[0]*(1-vPtrue[1]) + vPtrue[1]*cLambda + vPtrue[2]*cU;
	decl y= ranu(1,1);
	cY = exp(cLambda)*((1-y)^(-1/vPtrue[4])-1)^(1/vPtrue[3]);
	vY[t][0] = cY;
	cU = -1+(1+vPtrue[4])*cY^(vPtrue[3])*exp(-vPtrue[3]*cLambda)/(1+cY^(vPtrue[3])*exp(-vPtrue[3]*cLambda));
	}
	m_mData = vY[100:];
	delete vY,cLambda,cY,cU;	
	vPest = vPtrue; // starting values	
//	MaxControl(2000, 0,TRUE);
	time = timer();
	ir = MaxBFGS(BURRLik, &vPest, &dfunc, 0, 0);
	decl b =  vPest[1]^(2)-2*vPest[1]*vPest[2]*vPest[3]*vPest[4]/(2+vPest[4])+2*vPest[2]*vPest[3]*vPest[4]*(1+vPest[4])^(2)*gammafact(2+vPest[4])/gammafact(5+vPest[4]);
	if (b>=1){oxwarning("b>1"); println("b is:" ,b ~vPest');continue;}

	if (ir!=0 && ir!=1){
		++cfailed;
		continue;
	}
	decl vNumScore;
	vElapsedTime[0][i] = (timer()-time)/100;
	mCoefVal[0][i][] = vPest'~ir;
	++i;
	delete vPest, time, m_mData, ir;
}
Report(MaxBFGS,vElapsedTime[0],mCoefVal[0],vPtrue,cRep,cT,cfailed);
//Save result in database
decl dbase = new Database();  
dbase.Create(1,1,1,m_cRep,1);
dbase.Append(mCoefVal[0][][]~vElapsedTime[0], {"Const","AR", "MA", "Beta","Alpha", "time" });
dbase.SaveXls(sFilename);
//data.SaveXls(sFileNameData);

Draw(iPos,sFilename,mCoefVal[0][][:3]);

}

DySco::FSimulate(const mCoefVal, const vElapsedTime, const cT, const cRep, const vPtrue, const sFilename, const sFileNameData, const iPos)
{decl i=0,cfailed=0, dfunc=0,ir=0, time=0;
for (i = cfailed = 0;i<cRep;) { //; i++){
	println("BFGS iteration step: ", i);
	decl vY = new matrix[cT+100][1];
	decl cY,cLambda=vPtrue[0], cU=0, vPest, eps=0;
	for(decl t=0; t<m_cT+100;t++){
	cLambda = vPtrue[0]*(1-vPtrue[1]) + vPtrue[1]*cLambda + vPtrue[2]*cU;
	eps =ranf(1, 1, vPtrue[3], vPtrue[4]);
	cY = eps*exp(cLambda);
	vY[t][0] = cY;	
	cU = -0.5*vPtrue[3]+ 0.5*(vPtrue[3]+vPtrue[4])*vPtrue[3]*eps/(vPtrue[4]+vPtrue[3]*eps);
	}
	m_mData = vY[100:];
	delete vY,cLambda,cY,cU;	
	vPest = vPtrue; // starting values	
//	MaxControl(2000, 0,TRUE);
	time = timer();
	ir = MaxBFGS(FLik, &vPest, &dfunc, 0, 0);
	decl b = vPest[1]^(2)-2*vPest[1]*vPest[2]*vPest[3]/3 + 2*vPest[2]^(2)*vPest[3]^(2)/15;
	if (b>=1){oxwarning("b>1"); println("b is:" ,b ~vPest');continue;}

	if (ir!=0 && ir!=1){
		++cfailed;
		continue;
	}
	vElapsedTime[0][i] = (timer()-time)/100;
	mCoefVal[0][i][] = vPest'~ir;
	++i;

	delete vPest, time, m_mData, ir;
}

Report(MaxBFGS,vElapsedTime[0],mCoefVal[0],vPtrue,cRep,cT,cfailed); //Save result in database
decl dbase = new Database();  
dbase.Create(1,1,1,m_cRep,1);
dbase.Append(mCoefVal[0][][]~vElapsedTime[0], {"Const","AR", "MA", "NU","NU2","time" });
dbase.SaveXls(sFilename);
Draw(iPos,sFilename,mCoefVal[0][][:3]);
}

DySco::Draw(const iPos,const sFileName, const mCoefVal)
{
DrawDensity(iPos, mCoefVal', {"Const", "AR", "MA", "GA"}, TRUE,FALSE, TRUE, FALSE, FALSE,0,2);
ShowDrawWindow();
SaveDrawWindow(sprint(sFileName,".eps"));
}



DySco::NormTest1(const vX)                 // function with one argument
{
    decl xs, ct, skew, kurt, test;  // variables must be declared

    ct = rows(vX);                                 // sample size
    xs = standardize(vX);

    skew = sumc(xs .^ 3) / ct;                 // sample skewness
    kurt = sumc(xs .^ 4) / ct;                 // sample kurtosis
    test = sqr(skew) * ct / 6 + sqr(kurt - 3) * ct / 24;

return test | tailchi(test, 2);  // [0][0]: test, [1][0]: p-value
}

DySco::Report(const sMethod, const cElapsedTime, const mCoefVal, const vP, const cRep, const cT,const cRepFailed)
{
println("-------------------- ");
println("%-35s","Estimated model:",sprint("Score Gamma"));
println("%-35s","Maximisiation method:",sprint(sMethod));
println("%-35s","Time:",timestr(today()));
println("%-65s", " averaged Elapsed Time: ", sprint(sumc(cElapsedTime)/cRep, " sec"));
println("%-35s", "cT:", cT);
println("%-35s", "cRep:", cRep);
println("%-35s", "True Para: ", vP');
decl mom = moments(mCoefVal[][:3], 4,0);
decl cMeanBias = diagonal(vP - mom[1][]);
println("%-35s", "Mean Bias:" , cMeanBias);
decl cRSME=zeros(4,1);
cRSME =sqrt(sumr((mCoefVal[][:3]'-vP).^(2))/cRep);
println("%-35s", "RMSE (cmp to True):", cRSME');
println("%-35s", "Skew Estimate:", mom[3][]);
println("%-35s", "Kurto Estimate:", mom[4][]);
decl vNormTest= zeros(4,2);
for (decl j=0;j<4;j++){
	vNormTest[j][] = NormTest1(mCoefVal[][j]);
}
println("%-35s", "Normality Test", vNormTest[][0]');
println("%-35s", "Normality Test p val", vNormTest[][1]');
println("%-35s", "Not Converged:", cRepFailed);
println("-------------------- ");
println("");
}
/*--------------------------------------------------------------------------
 * Interface Code
 *
 *       (C) Philipp Andres 2012
 *
 *--------------------------------------------------------------------------*/
DySco::SendSpecials()
{    	 return <>; // {"Constant"};
}

DySco::SendVarStatus()
{    
    decl aArr;			
			aArr = {{ "&Y: Endogenous", 'Y',STATUS_GROUP + STATUS_ENDOGENOUS , Y_VAR}};
			aArr ~= {{ "&L: Leverage Variable", 'L', STATUS_GROUP , X_VAR}};
			// Index variables
    		aArr ~= {
             		{ "&i: ID var", 'i', STATUS_GROUP4, ID_VAR}
					// ,{ "&y: year index", 'y', STATUS_GROUP+STATUS_ONEONLY, YEAR_VAR}
					// ,{ "&w: weight", 'w', STATUS_GROUP + STATUS_ONEONLY, WEIGHT_VAR}
            		};
    return aArr;
    return Modelbase::SendVarStatus();   
}
			
DySco::ReceiveData()		 
{
	Modelbase::ReceiveData(); 
}
 
DySco::ReceiveModel()			 
{
	DeSelect(); // clear current variables
	Modelbase::ReceiveModel();
	//Select(X_VAR, "OxPackGetData"("SelGroup", X_VAR)); 	
	
	decl asL;
    GetGroupNames(X_VAR,&asL);
   
} 
						
DySco::SendMenu(const sMenu)
{								 
//		 if (sMenu == "ModelClass"){
//		return {{ "&1: Dynamic Conditional Score \tAlt+S", "modelclass0",m_iModelClass == SCORE},
//			 { "&2: Engle & Russel 1998 \tAlt+M","modelclass1",	m_iModelClass == MEM}
//			 };
//    }
//
//
if (sMenu == "Model"){
return
	{  { "Formulate ...\tAlt+F", "OP_FORMULATE"},
	   {"log-ACD model (Engle & Russell) \tAlt+M", "OP_MEM"},
	   { "Model &Settings...\tAlt+S", "OP_SETTINGS"},
	    { "&Estimates...\tAlt+E", "OP_ESTIMATE"},
		0,
		{ "&Progress...", "OP_PROGRESS"},
		0,
		{ "Options...\tAlt+O", "OP_OPTIONS"},
		0	   
		,{"Print Output Table ", "outputtable"}
		,{"Print LaTeX table ", "outputlatex"}
		,0
		,{"BibTeX Citation \tAlt+C", "bibtexcite"}
		,0
		//,{"Simulate \tAlt+S", "OP_SIMULATE"}
	 };
  }	 	
	if (sMenu == "Test"){
		decl aArr = {};        
            aArr ~= {
					 {"In Sample:"}
					 ,{"&3: Graphical Analysis", "OP_TEST_GRAPHICS"}
					 ,{"&4: Test Summary", "OP_TEST_SUMMARY"}
                   };
        return aArr;
	}

	if (sMenu == "Help"){
		decl aArr = {{"Handbook", "OP_HELP"}};        
		return aArr;
	}
	
}

DySco::DoSimDlg()
{
 decl adlg;
 m_iDeterministic = 0;  
 adlg = {
 		{"Monte Carlo simulations"},
 		{"Similation Settings", CTL_GROUP,1},
		{"Choose model", CTL_SELECT, 0, "Gamma|Weibull|Log Logistic|Burr|F|Dagum", "m_iModel"},
		{"True parameter vector", CTL_MATRIX, <-4;0.96;0.05;2>,"m_vTruePar"},
		{"Model specifications", CTL_GROUP,1},
		{"Number of time periods, T", CTL_INT,1000,"m_cT"},
//		{"# AR terms, p", CTL_INT,1, "cAR"},
//		{"# lagged innovation terms, q", CTL_INT, 1, "cMA"},
//		{"Two Components (0 or 1)", CTL_INT,0, "m_cTwoComp"},
//		{"Include Leverage Effects", CTL_CHECK,0,"m_cLev"},
		{"Starting Values", CTL_GROUP,1},
		{"Use True Parameter vector", CTL_CHECK,1,"m_vParStart"},
		//{"User Specified:", CTL_MATRIX, zeros(4,1), "m_vParStart"},
		//{"",CTL_SUBGROUP,0},
		{"Monte Carlo settings", CTL_GROUP,1},
		{"Number of replications", CTL_INT,3,"m_cR"},
//		{"BFGS, exact score", CTL_RADIO, 0, "m_iAlgorithm"},
		{"BFGS, approximate score", CTL_RADIO},
//		{"Method of Scoring", CTL_RADIO},
		{"Burn-in period" , CTL_INT,100, "m_cBurnIn"},
		{"Save Estimates under", CTL_STRING, "MCEstimates.xls", "sEstimates"},
//		{"Show coefficient distribution", CTL_SELECT,1,"No|Yes", "sGraphs"},
		{"Graph Filename", CTL_STRING, "GraphsMC.eps", "sGraphs"}
		};
		decl asoptions, aVals;
		if ("OxPackDialog"("Simulation", adlg, &asoptions,&aVals))
			{
//			println("aVals", aVals);
			if (aVals[0]==0)
				SetDistribution(GA);
			else if (aVals[0]==1)
				SetDistribution(WBL);				
			else if (aVals[0]==2)
				SetDistribution(LLOG);
			else if (aVals[0]==3)
				SetDistribution(BURR);
			else if (aVals[0]==4)
				SetDistribution(F);
			else if (aVals[0]==5)
				SetDistribution(DAG);

			m_vTruePar = aVals[1];
			m_cT = aVals[2];
//			m_vStruct[0] = aVals[3];
//			m_vStruct[1] = aVals[4];
//			m_iTwoComp = aVals[5];

//			if (aVals[6]!=0 && aVals[6]!=1){
//				println("Specify either with (1) /without (0) leverage.\n Inserted value of: ", aVals[6], "is not admissable");
//				exit(0);
//				}
//			else m_vStruct[3]=aVals[6];
		   m_vStruct = <1;1;1;0;0;0>;
		   m_iT1sel = 0;
		   m_iT2sel = m_cT-1;
		   m_iTwoComp = 0;
			if(aVals[3]==1) m_vParStart  = m_vTruePar;
			m_cRep	= aVals[4];
			decl sEstimatesName = aVals[6];
//			decl iShowGraph = aVals[7];
			decl sGraphName = aVals[7];

	decl mCoefVal = new matrix[m_cRep][rows(m_vTruePar)];
	decl vElapsedTime = new matrix[m_cRep][1];

	println("----------------------------------------- ");
	println("--         Start Simulation	     					 -- ");
	println("----------------------------------------- ");
	
	if (m_cDist==GA)
		GASimulate(&mCoefVal,&vElapsedTime, m_cT,m_cRep, m_vTruePar,  sEstimatesName, sGraphName, 0);
	if (m_cDist==WBL)
		WBLSimulate(&mCoefVal,&vElapsedTime, m_cT,m_cRep, m_vTruePar,  sEstimatesName, sGraphName, 0);
	if (m_cDist==LLOG)
		LLOGSimulate(&mCoefVal,&vElapsedTime, m_cT,m_cRep, m_vTruePar,  sEstimatesName, sGraphName, 0);
	if (m_cDist==BURR)
		BURRSimulate(&mCoefVal,&vElapsedTime, m_cT,m_cRep, m_vTruePar,  sEstimatesName, sGraphName, 0);
	if (m_cDist==F)
	    FSimulate(&mCoefVal,&vElapsedTime, m_cT,m_cRep, m_vTruePar,  sEstimatesName, sGraphName, 0);

 return TRUE; 
	}
	else
	return FALSE;
}

DySco::DoSettingsDlg()
{
	decl adlg, asoptions, aValues;
	adlg =
	{
	{"Model settings",CTL_GROUP,1},
	{"Choose model", CTL_SELECT, 1, "Exponential|Gamma|Weibull|Log Normal|Log Logistic|Burr|Dagum|F|Generalized Gamma|GB2", "m_iModel"},
	{"# AR lags", CTL_INT, 1, "m_vStruct[0]"},
	{"# MA lags", CTL_INT, 1, "m_vStruct[1]"},
	{"Include Leverage Effects", CTL_CHECK, FALSE,"m_vStruct[3]"},
//	{"One/Two Components", CTL_SELECT, 0, "One Component|Two Components", "m_iTwoComp"},
	{"# short run AR coefficients", CTL_INT,0,"m_vStruct[4]"},
	{"# short run MA coefficients", CTL_INT,0,"m_vStruct[5]"}
	//}	
	};
	
	if ("OxPackDialog"("Model Settings", adlg, &asoptions,&aValues))
	{		if (aValues[0]==0)
				SetDistribution(EX);
			else if (aValues[0]==1)
				SetDistribution(GA);
		//	else if (aValues[0]==2)
			//	SetDistribution(IGA);
			else if (aValues[0]==2)
				SetDistribution(WBL);				
			else if (aValues[0]==3)
				SetDistribution(LNORM);
			else if (aValues[0]==4)
				SetDistribution(LLOG);
			else if (aValues[0]==5)
				SetDistribution(BURR);
			else if (aValues[0]==6)
				SetDistribution(DAG);
			else if (aValues[0]==7)
				SetDistribution(F);
			else if (aValues[0]==8)
				SetDistribution(GGA);
			else if (aValues[0]==9)
				SetDistribution(GB2);	
			//SetShapeClass();
			//println(aValues[1][0]);
			//println(isarray(aValues));
			decl vStruct; decl aux=new matrix[5][1];
			aux[0] = aValues[1][0];
			aux[1] = aValues[2];
			aux[2] = aValues[3];
			aux[3] = aValues[4];
			aux[4] = aValues[5];		

//			println(aux);
			SetModel(aux);
//			vStruct = <aValues[1]; aValues[2]>;
////			vStruct[1] = aValues[1];
////			vStruct[1] = aValues[2];
//			println(vStruct);
////			println(sumc(vStruct));
//			println(<aValues[1],aValues[2],aValues[3],aValues[4],aValues[5]>);
//
//			SetModel(<aValues[1];aValues[2];aValues[3];aValues[4];aValues[5]>);
//println(m_vStruct);
//		m_vStruct[0] 		= 	aValues[1];
//		m_vStruct[1]		= 	aValues[2];
//		m_vStruct[3] 		=	aValues[3];	
//		m_vStruct[4]	 	= 	aValues[4];
//		m_vStruct[5]		= 	aValues[5];
		SetTwoComp();
		//println("m_iTwoComp" , m_iTwoComp);
		
	 return TRUE; 
	}
	else
	return FALSE;	
}

DySco::DoEstimateDlg()
{
 decl adgl, asoptions, aValues;
 decl m_iRoutine=0;
 adgl =
 {{"Choose Sample Period", CTL_GROUP,1},
  {"Start Period", CTL_INT,0,"m_iT1est"},
  {"End Period", CTL_INT,sizer(m_mData[][0]), "m_iT2est"},
  {"Choose Maximisation Routine", CTL_GROUP,1}};
  
  if(m_iModelClass==SCORE){
  	adgl ~= {{"Method of Scoring", CTL_RADIO, 2, "m_cAlgorithm"},
  			{"BFGS, exact derivatives", CTL_RADIO},// 1, "m_cAlgorithm"},
  			{"BFGS, approximate derivatives", CTL_RADIO},
			{"Newton-Raphson", CTL_RADIO},
			{"BHHH", CTL_RADIO}};
	}
  else if (m_iModelClass==MEM) {
	adgl ~= {{"BFGS, exact derivatives", CTL_RADIO, 1, "m_cAlgorithm"},// 1, "m_cAlgorithm"},
  			{"BFGS, approximate derivatives", CTL_RADIO}};
  	}
  if(m_iModelClass==SCORE) adgl ~= {{"Newton-Raphson", CTL_RADIO}};
  if(m_iModelClass==SCORE) adgl ~= {{"BHHH", CTL_RADIO}};
  adgl ~= {{"Further Options", CTL_GROUP,1},
  //{"# Iterations", CTL_INT,1000,"MaxControl"},
  {"Starting Values", CTL_GROUP,1},
  {"Automatic", CTL_RADIO,0, "m_vStartValues"},
  {"User Defined", CTL_RADIO},
  {"User Defined Starting Values:", CTL_SUBGROUP,1},
  { "Values", CTL_MATRIX,	zeros(1,6), 1, 9, 1, 1, "m_iModelBased" },
  {"",CTL_SUBGROUP,0}
  };
if ("OxPackDialog"("Maximisation Options", adgl ,&asoptions,&aValues)){
//	decl aux1, aux2;
//	if(aValues[0]==1) aux1 = -1;
//	else aux1 = aValues[0];
//	if(aValues[1]==sizer(m_mData)) aux2 = -1;
//	else aux2 = aValues[1];
//
	SetInSamplePeriod(aValues[0], aValues[1]);
	//SetSelSample(aux1,1,aux2,1);
	//delete aux1, aux2;
	//SetcT();

	if (aValues[2]==0) SetRoutine(SCORING);
	else if (aValues[2]==1){
		SetRoutine(BFGS); SetScore(1);}
	else if (aValues[2]==2){
		SetRoutine(BFGS); SetScore(1);}
	else if (aValues[2] == 3) SetRoutine(NR);
	else if (aValues[2] ==4) SetRoutine(BHHH);

	//m_cMaxIter = aValues[3];

	if (aValues[3]==0)			m_iAutoStartValues = 1;
	else if (aValues[3]==1)		{m_iAutoStartValues=0; UserStartValues(aValues[4]);}

	//InitGlobals();
	InitData();		
  }
}

DySco::DoTestGraphicsDlg()
{
decl adgl, asoptions, aValues;
 adgl ={
  {"In Sample Graphs"},
  {"Epsilons", CTL_GROUP,1},	
  {"Predictions", CTL_CHECK,0, "redictions"},
  {"ACF Epsilons", CTL_CHECK,1,"acfeps"},
  {"ACF Epsilons^(2)",CTL_CHECK,0,"acfesp2"},
  {"PITs Epsilons", CTL_CHECK,1,"piteps"},
  {"QQ Plot", CTL_CHECK,1, "qqeps"},
  {"QQ Plot PIT", CTL_CHECK,1,"QQPits"}
  };
  if (m_iModelClass==SCORE) {
  adgl~= {{"Scores", CTL_GROUP,1},
  {"ACF Scores", CTL_CHECK,1, "acfscores"},
  {"ACF Scores^(2)",CTL_CHECK,0, "acfscores2"},
  {"PITs Scores", CTL_CHECK,1, "pitscores"},
  {"QQ Plot", CTL_CHECK,1, "QQscores"}
  };
  } 
  if ((m_iT2sel+1)!=m_cTAll){
  adgl ~= { {"Out of Sample Graphs"},
  	{"Epsilons", CTL_GROUP,1},	
  	{"Predictions", CTL_CHECK,0, "OoSredictions"},
  	{"ACF Epsilons", CTL_CHECK,1,"OoSacfeps"},
  	{"ACF Epsilons^(2)",CTL_CHECK,0,"OoSacfesp2"},
  	{"PITs Epsilons", CTL_CHECK,1,"OoSpiteps"},
  	{"QQ Plot", CTL_CHECK,1, "OoSqqeps"}
	};
	if (m_iModelClass==SCORE) {
	adgl~= {{"Scores", CTL_GROUP,1},
  	{"ACF Scores", CTL_CHECK,1, "OoSacfscores"},
  	{"ACF Scores^(2)",CTL_CHECK,0, "OoSacfscores2"},
  	{"PITs Scores", CTL_CHECK,1, "OoSpitscores"},
  	{"QQ Plot", CTL_CHECK,1, "OoSQQscores"}
  	};
};
};
	if ("OxPackDialog"("Graphical Diagnostics", adgl ,&asoptions,&aValues))
	{
	 GetEps();
	 //if (m_iModelClass==SCORE)
	 GetScores();
	 decl nxt = 0 ; // graph window count
		if (m_iModelClass==SCORE){
		if (aValues[0] == 1){
			DrawPredict(nxt);
			nxt += 1;
			}
		if (aValues[1] == 1){
			DrawEpsACF(nxt);
			nxt += 1;
			}
		if (aValues[2] == 1){
			DrawEpsSqACF(nxt);
			nxt += 1;
			}
		if (aValues[3] == 1){
			DrawPitEps(nxt, m_vEps[0:m_iT2sel]);
			nxt += 1;
			}
		if (aValues[4] == 1){
			QQPlotDrawEps(nxt,m_vEps[0:m_iT2sel]);
			nxt += 1;
			}

		if (aValues[5] == 1){
			QQPlotPits(nxt, m_vScore[0:m_iT2sel]);
			nxt += 1;
			}
			
		if (aValues[6] == 1){
			DrawScoreACF(nxt);
			nxt += 1;
			}
		if (aValues[7] == 1){
			DrawScoreSqACF(nxt);
			nxt += 1;
			}
		if (aValues[8] == 1){
			DrawPitScore(nxt, m_vScore[0:m_iT2sel]);
			nxt += 1;
			}
		if (aValues[9] == 1){
			QQPlotDrawScore(nxt, m_vScore[0:m_iT2sel]);
			nxt += 1;
			}
	

			
	    if (m_iT2sel+1!=m_cTAll){
		   	//GetEpsOoS();
			if (aValues[10]==1){
				DrawPredictOoS(nxt);
				nxt += 1;
			}
			if (aValues[11]==1){
				DrawEpsACFOoS(nxt);
				nxt += 1;
				}
			if (aValues[12]==1){
				DrawEpsSqACF(nxt);
				nxt += 1;
				}
			if (aValues[13]==1){
				DrawPitEpsOoS(nxt, m_vEps[m_iT2sel+1:]);
				nxt += 1;
				}
			if (aValues[14]==1){
				QQPlotDrawEpsOoS(nxt, m_vEps[m_iT2sel+1:]);
				nxt += 1;
				}
			if (aValues[15]==1){
				DrawScoreACFOoS(nxt);
				nxt += 1;
				}
			if (aValues[16]==1){
				DrawScoreSqACFOoS(nxt);
				nxt += 1;
				}
			if (aValues[17]==1){
				DrawPitScoreOoS(nxt, m_vScore[m_iT2sel+1:]);
				nxt += 1;
				}
			if (aValues[18]==1){
				QQPlotDrawScoreOoS(nxt, m_vScore[m_iT2sel+1:]);
				nxt += 1;
				}
			}
		 }
		 else if (m_iModelClass== MEM){
		 	if (aValues[0] == 1){
			DrawPredict(nxt);
			nxt += 1;
			}
		if (aValues[1] == 1){
			DrawEpsACF(nxt);
			nxt += 1;
			}
		if (aValues[2] == 1){
			DrawEpsSqACF(nxt);
			nxt += 1;
			}
		if (aValues[3] == 1){
			DrawPitEps(nxt, m_vEps[0:m_iT2sel]);
			nxt += 1;
			}
		if (aValues[4] == 1){
			QQPlotDrawEps(nxt,m_vEps[0:m_iT2sel]);
			nxt += 1;
			}

		if (aValues[5] == 1){
			QQPlotPits(nxt, m_vScore[0:m_iT2sel]);
			nxt += 1;
			}

	    if (m_iT2sel+1!=m_cTAll){
			if (aValues[6]==1){
				DrawPredictOoS(nxt);
				nxt += 1;
			}
			if (aValues[7]==1){
				DrawEpsACFOoS(nxt);
				nxt += 1;
				}
			if (aValues[8]==1){
				DrawEpsSqACF(nxt);
				nxt += 1;
				}
			if (aValues[9]==1){
				DrawPitEpsOoS(nxt, m_vEps[m_iT2sel+1:]);
				nxt += 1;
				}
			if (aValues[10]==1){
				QQPlotDrawEpsOoS(nxt,m_vEps[m_iT2sel+1:]);
				nxt += 1;
				}
			}			
		 }			
		ShowDrawWindow();
		}
}

//DySco::DoTestDlg()
//{																	 
//decl adgl, asoptions, aValues;
// adgl ={
//  {"In Sample Tests", CTL_GROUP,1},
//  {"Performance Tests", CTL_CHECK,1, "pertest"},
//  {"Sample Statistics", CTL_CHECK,1,"Sample"},
//  {"LBQ Test: Residuals", CTL_CHECK,1, "LBQ"},
//  { "Lags", CTL_MATRIX, <5;10;50>, 1, 9, 1, 1, "Lags" },
//  {"LBQ Test: Residual PITs", CTL_CHECK,1, "LBQ"},
//  { "Lags", CTL_MATRIX, <5;10;50>, 1, 9, 1, 1, "Lags" }	
//  };
//  if (m_iModelClass==SCORE){
//  adgl~= {{"LBQ Test: Scores", CTL_CHECK,1, "LBQ2"},
//  { "Lags", CTL_MATRIX, <5;10;50>, 1, 9, 1, 1, "Lags" }};
//  }
//   // {"Pearson's chi-squared test", CTL_CHECK,0, "Chi2"}
//
// if (m_iT2sel+1!=m_cTAll){
// adgl~= {
//  {"Out of Sample Tests", CTL_GROUP,1},
//  {"Performance Tests", CTL_CHECK,1, "pertest"},
//  {"LBQ Test: Residuals", CTL_CHECK,1, "LBQ"},
//  { "Lags", CTL_MATRIX, <5;10;50>, 1, 9, 1, 1, "Lags" },
//  {"LBQ Test: Residual PITs", CTL_CHECK,1, "LBQ"},
//  { "Lags", CTL_MATRIX, <5;10;50>, 1, 9, 1, 1, "Lags" }	
//  };
//  if (m_iModelClass==SCORE){
//  adgl ~={{"LBQ Test: Scores", CTL_CHECK,1, "LBQ2"},
//  { "Lags", CTL_MATRIX, <5;10;50>, 1, 9, 1, 1, "Lags2" }
//  }; }}
//   
//if ("OxPackDialog"("Diagnostic Tests", adgl ,&asoptions,&aValues)){
//	// println(aValues);
//	  GetPredict();
//	  GetEps();
//	 if (m_iModelClass==SCORE){
//	 if(aValues[0]==1){
//	  PrintSampleStats();
//	  println("");
//	 }
//	 if (aValues[1]==1){
//		println("------------------------------");
//		//decl aux = 0;//strfind(m_asNames, "Range");
//		println("-----------------------------------------");
//		println("---- In-Sample Evaluation Measures ----");
//		println("-----------------------------------------");
//		FEM(m_vMu[0:m_iT2sel],m_mData[0:m_iT2sel][0],0); //m_mData[0:m_iT2est][1]);
//		println("------------------------------");
//		}
//	 if (aValues[2] == 1){
//	 	println("");
//		println(" --- LBQ Test for Residuals --- ");
//		LBQTest(m_vEps[0:m_iT2sel], aValues[3]);
//		println("");
//		}
//	if (aValues[4] == 1){
//	 	println("");
//		println(" --- LBQ Test for Residual PITs --- ");
//		LBQStatPitEps(m_vEps[0:m_iT2sel], aValues[5]);
//		println("");
//		}
//
//	if (aValues[6] ==1){
//		println(" --- LBQ Test for Scores --- ");
//		LBQTest(m_vU[0:m_iT2sel], aValues[7]);
//		println("");
//		}
//	if(m_iT2sel+1!=m_cTAll){
//		println("");
//		println("-----------------------------------------");
//		println("------- Out of Sample Statistics ------- ");
//		println("-----------------------------------------");
////		if(aValues[6]==1){
////	  		PrintSampleStats();
////	  		println("");
////		 }
//		if (aValues[8]==1){
//			println("------------------------------");
//			println("Forecast Evaluation Measures");
//			//decl aux = 0; //strfind(m_asNames, "Range");
//			FEM(m_vMu[m_iT2sel+1:], m_mData[m_iT2sel+1:][0],1);
//			println("------------------------------");
//		}
//
//		if (aValues[9]==1){
//			println(" --- LBQ Test for Out of Sample Residuals --- ");
//			LBQTest(m_vEps[m_iT2sel+1:], aValues[10]);
//			println("");
//		}
//
//		if (aValues[11]==1){
//			println(" --- LBQ Test for Out of Sample Residual PITs --- ");
//			LBQStatPitEps(m_vEps[m_iT2sel+1:], aValues[12]);
//			println("");
//		}
//
//		
//		if(aValues[13]==1){
//			println(" --- LBQ Test for Out of Sample Scores --- ");
//			LBQTest(m_vU[m_iT2sel+1:], aValues[14]);
//			println("");
//		}
//	}
//	}
//	else if (m_iModelClass==MEM){
//	 if(aValues[0]==1){
//	  PrintSampleStats();
//	  println("");
//	 }
//	 if (aValues[1]==1){
//		println("------------------------------");
//		println("In-Sample Evaluation Measures");
//		//decl aux = 0; //strfind(m_asNames, "Range");
//		FEM(m_vMu[0:m_iT2sel],m_mData[0:m_iT2sel][0],0); //m_mData[0:m_iT2est][1]);
//		println("------------------------------");
//		}
//	 
//	 if (aValues[2] == 1){
//	 	println("");
//		println(" --- LBQ Test for Residuals --- ");
//		LBQTest(m_vEps[0:m_iT2sel], aValues[3]);
//		println("");
//		}
//
//	if (aValues[4] == 1){
//	 	println("");
//		println(" --- LBQ Test for Residual PITs --- ");
//		LBQStatPitEps(m_vEps[0:m_iT2sel], aValues[5]);
//		println("");
//		}
//
//
//	if(m_iT2sel+1!=m_cTAll){
//		println("");
//		println("Out of Sample Statistics");
////		if(aValues[4]==1){
////	  		PrintSampleStats();
////	  		println("");
////		 }
//		if (aValues[6]==1){
//			println("------------------------------");
//			println("Forecast Evaluation Measures");
//			//decl aux = 0;//strfind(m_asNames, "Range");
//			FEM(m_vMu[m_iT2sel+1:], m_mData[m_iT2sel+1:][0],1);
//			println("------------------------------");
//		}
//
//		if (aValues[7]==1){
//			println(" --- LBQ Test for Out of Sample Residuals --- ");
//			LBQTest(m_vEps[m_iT2sel+1:], aValues[8]);
//			println("");
//		}
//	
//		if (aValues[9] == 1){
//	 		println("");
//			println(" --- LBQ Test for Out of Sample Residual PITs --- ");
//			LBQStatPitEps(m_vEps[m_iT2sel+1:], aValues[10]);
//			println("");
//
//	   }}
//	   }
//	}
//}

DySco::DoTestDlg()
{																	 
decl adgl, asoptions, aValues;
 adgl ={
  {"In Sample Tests", CTL_GROUP,1},
  {"Performance Tests", CTL_CHECK,1, "pertest"},
  {"Sample Statistics", CTL_CHECK,1,"Sample"},
  {"LBQ Test: Residuals", CTL_CHECK,1, "LBQ"},
  { "Lags", CTL_MATRIX, <5;10;50>, 1, 9, 1, 1, "Lags" },
  {"LBQ Test: Residual PITs", CTL_CHECK,1, "LBQ"},
  { "Lags", CTL_MATRIX, <5;10;50>, 1, 9, 1, 1, "Lags" }	
  };
 // if (m_iModelClass==SCORE){
  adgl~= {{"LBQ Test: Scores", CTL_CHECK,1, "LBQ2"},
  { "Lags", CTL_MATRIX, <5;10;50>, 1, 9, 1, 1, "Lags" }};
 // }
   // {"Pearson's chi-squared test", CTL_CHECK,0, "Chi2"}

 if (m_iT2sel+1!=m_cTAll){
 adgl~= {
  {"Out of Sample Tests", CTL_GROUP,1},
  {"Performance Tests", CTL_CHECK,1, "pertest"},
  {"LBQ Test: Residuals", CTL_CHECK,1, "LBQ"},
  { "Lags", CTL_MATRIX, <5;10;50>, 1, 9, 1, 1, "Lags" },
  {"LBQ Test: Residual PITs", CTL_CHECK,1, "LBQ"},
  { "Lags", CTL_MATRIX, <5;10;50>, 1, 9, 1, 1, "Lags" }	
  };
 // if (m_iModelClass==SCORE){
  adgl ~={{"LBQ Test: Scores", CTL_CHECK,1, "LBQ2"},
  { "Lags", CTL_MATRIX, <5;10;50>, 1, 9, 1, 1, "Lags2" }
  }; //}
  }
   
if ("OxPackDialog"("Diagnostic Tests", adgl ,&asoptions,&aValues)){
	// println(aValues);
	  GetPredict();
	  GetEps();
//	 if (m_iModelClass==SCORE){
	 if(aValues[0]==1){
	  PrintSampleStats();
	  println("");
	 }
	 if (aValues[1]==1){
		println("------------------------------");
		//decl aux = 0;//strfind(m_asNames, "Range");
		println("-----------------------------------------");
		println("---- In-Sample Evaluation Measures ----");
		println("-----------------------------------------");
		FEM(m_vMu[0:m_iT2sel],m_mData[0:m_iT2sel][0],0); //m_mData[0:m_iT2est][1]);
		println("------------------------------");
		}
	 if (aValues[2] == 1){
	 	println("");
		println(" --- LBQ Test for Residuals --- ");
		LBQTest(m_vEps[0:m_iT2sel], aValues[3]);
		println("");
		}
	if (aValues[4] == 1){
	 	println("");
		println(" --- LBQ Test for Residual PITs --- ");
		LBQStatPitEps(m_vEps[0:m_iT2sel], aValues[5]);
		println("");
		}

	if (aValues[6] ==1){
		println(" --- LBQ Test for Scores --- ");
		LBQTest(m_vU[0:m_iT2sel], aValues[7]);
		println("");
		}
	if(m_iT2sel+1!=m_cTAll){
		println("");
		println("-----------------------------------------");
		println("------- Out of Sample Statistics ------- ");
		println("-----------------------------------------");
//		if(aValues[6]==1){
//	  		PrintSampleStats();
//	  		println("");
//		 }
		if (aValues[8]==1){
			println("------------------------------");
			println("Forecast Evaluation Measures");
			//decl aux = 0; //strfind(m_asNames, "Range");
			FEM(m_vMu[m_iT2sel+1:], m_mData[m_iT2sel+1:][0],1);
			println("------------------------------");
		}

		if (aValues[9]==1){
			println(" --- LBQ Test for Out of Sample Residuals --- ");
			LBQTest(m_vEps[m_iT2sel+1:], aValues[10]);
			println("");
		}

		if (aValues[11]==1){
			println(" --- LBQ Test for Out of Sample Residual PITs --- ");
			LBQStatPitEps(m_vEps[m_iT2sel+1:], aValues[12]);
			println("");
		}

		
		if(aValues[13]==1){
			println(" --- LBQ Test for Out of Sample Scores --- ");
			LBQTest(m_vU[m_iT2sel+1:], aValues[14]);
			println("");
		}
	}
//	}
//	else if (m_iModelClass==MEM){
//	 if(aValues[0]==1){
//	  PrintSampleStats();
//	  println("");
//	 }
//	 if (aValues[1]==1){
//		println("------------------------------");
//		println("In-Sample Evaluation Measures");
//		//decl aux = 0; //strfind(m_asNames, "Range");
//		FEM(m_vMu[0:m_iT2sel],m_mData[0:m_iT2sel][0],0); //m_mData[0:m_iT2est][1]);
//		println("------------------------------");
//		}
//	 
//	 if (aValues[2] == 1){
//	 	println("");
//		println(" --- LBQ Test for Residuals --- ");
//		LBQTest(m_vEps[0:m_iT2sel], aValues[3]);
//		println("");
//		}
//
//	if (aValues[4] == 1){
//	 	println("");
//		println(" --- LBQ Test for Residual PITs --- ");
//		LBQStatPitEps(m_vEps[0:m_iT2sel], aValues[5]);
//		println("");
//		}
//
//
//	if(m_iT2sel+1!=m_cTAll){
//		println("");
//		println("Out of Sample Statistics");
////		if(aValues[4]==1){
////	  		PrintSampleStats();
////	  		println("");
////		 }
//		if (aValues[6]==1){
//			println("------------------------------");
//			println("Forecast Evaluation Measures");
//			//decl aux = 0;//strfind(m_asNames, "Range");
//			FEM(m_vMu[m_iT2sel+1:], m_mData[m_iT2sel+1:][0],1);
//			println("------------------------------");
//		}
//
//		if (aValues[7]==1){
//			println(" --- LBQ Test for Out of Sample Residuals --- ");
//			LBQTest(m_vEps[m_iT2sel+1:], aValues[8]);
//			println("");
//		}
//	
//		if (aValues[9] == 1){
//	 		println("");
//			println(" --- LBQ Test for Out of Sample Residual PITs --- ");
//			LBQStatPitEps(m_vEps[m_iT2sel+1:], aValues[10]);
//			println("");
//
//	   }
//}
//	   }
	}
}		

DySco::ReceiveMenuChoice(const sMenu)
{
   decl ret;
   switch_single(sMenu)
   {
	case "OP_SIMULATE":
		if (DoSimDlg())
			// uncomment to move to Formulate Dlg after clicking ok in the simulate dlg
			//"OxPackSendMenuChoice"("OP_FORMULATE");
		return 1;
	
	case "OP_FORMULATE":
		SetModelClass(SCORE);
		if (DoFormulateDlg(-1))
			"OxPackSendMenuChoice"("OP_SETTINGS");
		return 1;

	case "OP_MEM":
		SetModelClass(MEM);
		if (DoFormulateDlg(-1))
			"OxPackSendMenuChoice"("OP_SETTINGS");
		return 1;
		
	case "OP_SETTINGS":
		if(DoSettingsDlg())		
			"OxPackSendMenuChoice"("OP_ESTIMATE");
		return 1;
		
	case "OP_ESTIMATE":{
	 if(DoEstimateDlg())
			"OxPackSendMenuChoice"("OP_ESTIMATE");
	 return 1;
	 }
	 
	case "OP_OPTIONS":{
		// process user actions for options dialog
		Modelbase::ReceiveMenuChoice(sMenu);
	}

	case "OP_TEST_GRAPHICS":{
		// call graphical diagnostics window
		if(DoTestGraphicsDlg())
			"OxPackSendMenuChoice"("OP_TESTGRAPHICS");
		return 1;
	}
	case "OP_TEST_SUMMARY":{
		if(DoTestDlg())
		"OxPackSendMenuChoice"("OP_TESTGRAPHICS");
		return 1;
	}

	case "outputtable":
		Output();
		return TRUE;
    case  "outputlatex":
		DCSOutputParLaTex();
		return TRUE;
    case "bibtexcite":
		Citation();
		return TRUE;
//	default:
	return Modelbase::ReceiveMenuChoice(sMenu);		
}
}





