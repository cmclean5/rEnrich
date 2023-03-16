#include "Headers.h"
#include "NetworkEnrichment.h"

//----------------
//NOTE:
// When editing this file, i.e. adding or removing functions, need to carry out the following tasks:
// 1) edit file NAMESPACE, adding/removing the function name in export(...)
// Next run following in R, inside the package:
// 2) $> cd BlockModelr/
// 3) $> R
// 4)  > library(Rcpp)
// 5)  > Rcpp::compileAttributes()
// running 5) will create new RcppExports.R in R/ and
// new RcppExports.cpp in src/
//----------------

// Global
NetworkEnrichment *enrD=0;

// [[Rcpp::export]]
void test(){
    //cout << "  Hello! " << endl;
    }

//' Reset calculation environment.
//' @noRd
// [[Rcpp::export]]
void reset(){ enrD=0; }

//' Clean up the calculation environment remove all results.
//' @noRd
// [[Rcpp::export]]
void erase(){ if(enrD !=0){ delete enrD; } }

//' Load data for calculation.
//'
//' This function reads group/cluster membership and element
//' annotation data for analysis. It does not return anything,
//' it just create required structures in memory.
//'
//' @param x group/cluster membership data with element name in the
//' first column and group id in the second
//' @param anno element annotation data with term id, term name and
//' the element name in columns one, two and three respectively.
//' @noRd
//'
// [[Rcpp::export]]
void load(Rcpp::DataFrame x,
          Rcpp::DataFrame anno ){

  int i,j,k, x_size, anno_size;

  Rcpp::StringVector V1, V2, V3;

  //Create Network Enrichment object and pass our input files to it.
  reset();

  int x_cols    = x.length();
  int x_rows    = x.nrows();

  int anno_cols = anno.length();
  int anno_rows = anno.nrows();

  if( (x_cols == 2 && x_rows > 0) &&
      (anno_cols == 3 && anno_rows > 0) ){


     //cout << "> loading dataframe..." << endl;
    //cout << "> x_cols: " << x_cols << " anno_cols: " << anno_cols << endl;

    // set size and fill the community membership dataset
    x_size = x_rows * x_cols;
    string *MEMBERSHIP = new string[x_size];

    // fill membership dataset
    V1 = x[0];
    V2 = x[1];

    for(k=0; k<x_size; k++){
      i = floor(k/x_cols);
      j = k % x_cols;

      Rcpp::String v1(V1[i]);
      Rcpp::String v2(V2[i]);

      if( j == 0 ){ MEMBERSHIP[(i*x_cols)+j] = v1.get_cstring(); }
      if( j == 1 ){ MEMBERSHIP[(i*x_cols)+j] = v2.get_cstring(); }

    }

    // set size and fill the annotation dataset
    anno_size = anno_rows * anno_cols;
    string *ANNO = new string[anno_size];

    V1 = anno[0];
    V2 = anno[1];
    V3 = anno[2];

    for(k=0; k<anno_size; k++){
      i = floor(k/anno_cols);
      j = k % anno_cols;

      Rcpp::String v1(V1[i]);
      Rcpp::String v2(V2[i]);
      Rcpp::String v3(V3[i]);

      if( j == 0 ){ ANNO[(i*anno_cols)+j] = v1.get_cstring(); }
      if( j == 1 ){ ANNO[(i*anno_cols)+j] = v2.get_cstring(); }
      if( j == 2 ){ ANNO[(i*anno_cols)+j] = v3.get_cstring(); }

    }

    enrD = new NetworkEnrichment(MEMBERSHIP, x_rows, x_cols,
                                 ANNO, anno_rows, anno_cols);

  }

  //cout << "> done!" << endl;

}

//' Run enrichment analysis.
//'
//' This function perform actual calculations.
//'
//' @param useChi2 =0,
//' @param useOneSided =0,
//' @param useTwoSided =1,
//' @param useMaxSS =0,
//' @param runPerm =0,
//' @param singlePerm =0,
//' @param useSeed =0,
//' @param setNOP =0,
//' @param pesudoCount =-1.0,
//' @param FDRmeth ="BY"
//' @noRd
//'
// [[Rcpp::export]]
void run( Rcpp::IntegerVector useChi2=0,
          Rcpp::IntegerVector useOneSided=0,
          Rcpp::IntegerVector useTwoSided=1,
          Rcpp::IntegerVector useMaxSS=0,
          Rcpp::IntegerVector runPerm=0,
          Rcpp::IntegerVector singlePerm=0,
          Rcpp::IntegerVector useSeed=0,
          Rcpp::IntegerVector setNOP=0,
          Rcpp::NumericVector pesudoCount=-1.0,
          Rcpp::String FDRmeth="BY" ){

  int i,cal1;

  Rcpp::StringVector FDRset(3);
  FDRset[0] = "BH"; FDRset[1] = "BY"; FDRset[2] = "BL";
  Rcpp::String OUTDIR, Ext;

  bool useFDR=false;

  OUTDIR="";
  Ext="";

  if( enrD != 0 ){

  //cout << "> running enrichment..." << endl;

  //
  //if( useRelDist ){ enrD->setRelDist( true ); }

  //
  if( useChi2[0] ){ enrD->setChi2( true ); }

  //calculate depletion
  if( useOneSided[0] ){ enrD->oneSided(); }

  //calculate two-sided
  if( useTwoSided[0] ){ enrD->twoSided(); }

   //calculate two-sided
  if( useMaxSS[0] ){ enrD->maxSS(); }

  // Set Random number seed
  if( useSeed[0] != 0 ){ enrD->seedOffSet(true, useSeed[0]); }

  // Set number of permutations
  if( setNOP[0] != 0 ){ enrD->setNoP( setNOP[0] ); }

  // set pseudo count
  if( pesudoCount[0] != -1 ){ enrD->setPesudoCount( pesudoCount[0] ); }

  // Set which FDR method to use: BY (default), BH, BL.
  for( i=0; i<FDRset.size(); i++ ){
    if( FDRset[i] == FDRmeth ){
      enrD->setFDRmethod( FDRmeth.get_cstring() );
      useFDR = true;
    }
  }

  cal1 = enrD->calculateOverlapinCommunities(runPerm[0],
                                             OUTDIR.get_cstring(),
                                             Ext.get_cstring(),
                                             useFDR,
                                             false,
                                             singlePerm[0]);
  }

  //cout << "> done." << endl;

}

//' Aggregate results into the final matrix.
//'
//' @param printTwoSided print two sided p-values
//' @param usePrintAlt print also alternative side, i.e. depletion
//' @param usePrintCnew print old community ids
//' @param usePrintID print annotation ID if \code{TRUE},
//' else annotation description
//' @param usePrintAn print annotation type size
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::CharacterMatrix getResults(Rcpp::IntegerVector printTwoSided=1,
                                 Rcpp::IntegerVector usePrintAlt=0,
                                 Rcpp::IntegerVector usePrintCnew=0,
                                 Rcpp::IntegerVector usePrintID=0,
                                 Rcpp::IntegerVector usePrintAn=1){

  /*
   printTwoSided[1] == print two sided pvalues
   usePrintalt[1]   == print also alternative side, i.e. depletion
   usePrintCnew[0]  == print old community ids
   usePrintID[1]    == print annotation ID, else annotation description
   usePrintAn[1]    == print annotation type size.
  */

  int N,M,F;
  int i,m,f,k,K;
  int SIGMASIZE,index;
  int Coffset,Foffset,Fsize;

  bool printTS=false;
  bool printAlt=false;
  bool printCnew=false;
  bool printID=false;
  bool printAn=false;

  if( printTwoSided[0]==1 )   { printTS=true; }
  if( usePrintAlt[0]==1 )     { printAlt=true; }
  if( usePrintCnew[0]==1 )    { printCnew=true; }
  if( usePrintID[0]==1 )      { printID=true; }
  if( usePrintAn[0]==1 )      { printAn=true; }

  Rcpp::CharacterMatrix RES;

  if( enrD != 0 ){

    //cout << "> get results..." << endl;

    int BUFFERSIZE = enrD->getBufferSize();
    char buffer    [BUFFERSIZE];
    char bufferOR  [BUFFERSIZE];

    //
    if( printAlt ){ enrD->setALT( false ); }

    //Whether to print annotation IDs, or description
    if( printCnew ){ enrD->setPrintCNEW( true ); }

    //Whether to print annotation IDs, or description
    if( printID ){ enrD->setPrintID( true ); }

    //Whether to print annotation Size
    if( printAn ){ enrD->setPrintAn( true ); }

    N         = enrD->getN();
    M         = enrD->getM();
    F         = enrD->getF();
    index     = enrD->getANNOindex();
    SIGMASIZE = enrD->getSigmaSize();

    int bonferroni [SIGMASIZE];
    for( i=0; i<SIGMASIZE; i++){
      bonferroni[i] = enrD->calBonferroni(i, (double)(M*F));
    }

    Coffset=2;
    if( printAlt ){ Foffset = 12; } else { Foffset = 9; }
    Fsize  = ((F*Foffset)+Coffset);

    RES    = Rcpp::CharacterMatrix(M,Fsize);

    // Set colnames of RES
    Rcpp::CharacterVector cn (Fsize);
    cn[0] = "C"; cn[1] = "Cn";
    for( f=0, k=2; f<F; f++){
      (printID ? cn[k++] = enrD->getAnnoID(index, f) :
        cn[k++] = enrD->getAnnoDes(index, f));
      if( printAlt ){
        cn[k++] = "actual overlap";
        cn[k++] = "expected overlap";
        cn[k++] = "OR";
        cn[k++] = "95% CI";
        cn[k++] = "p-value";
        cn[k++] = "p.adjusted";
        cn[k++] = "perm.p_value";
        cn[k++] = "BC";
        cn[k++] = "p-value(ALT)";
        cn[k++] = "p.adjusted(ALT)";
        cn[k++] = "BC(ALT)";
      } else {
        cn[k++] = "actual overlap";
        cn[k++] = "expected overlap";
        cn[k++] = "OR";
        cn[k++] = "95% CI";
        cn[k++] = "p-value";
        cn[k++] = "p.adjusted";
        cn[k++] = "perm.p_value";
        cn[k++] = "BC";
      }
    }

    colnames(RES) = cn;

    //--- loop over all communities
    for(m=0; m<M; m++){

      if( printCnew ){ RES(m,0) = (enrD->getCom0(m)-enrD->getKOFFSET()); }
      else           { RES(m,0) = enrD->getCom0(m); }

      RES(m,1) = enrD->getCom1(m);

      //--- loop over each annotation type
      for(f=0, k=2; f<F; f++ ){

        int ele = (m*F)+f;

        string starsEO = "";
        string starsET = "";
        string starsDO = "";
        string starsDT = "";

        for(i=0; i<SIGMASIZE; i++){
          if( enrD->getPvalue(ele)   <= bonferroni[i] ){ starsEO += "*";  }
          if( enrD->getPvalueT(ele)  <= bonferroni[i] ){ starsET += "*";  }
          if( enrD->getPvalueD(ele)  <= bonferroni[i] ){ starsDO += "*";  }
          if( enrD->getPvalueDT(ele) <= bonferroni[i] ){ starsDT += "*";  }
        }

        //--- Hypergeometric mean
        double mn = 0.0;
        mn  = double(enrD->getCom1(m)) * double(enrD->getAnnoK( index, f));
        mn /= double(N);

        //--- Odds Ratio
        double ors = 0;
        double orL = 0;
        double orU = 0;

        double a   = double(enrD->getOverlap(ele));
        double b   = double(enrD->getCom1(m) - enrD->getOverlap(ele));
        double c   = double(enrD->getAnnoK(index,f) - enrD->getOverlap(ele));
        double d   = double(N - enrD->getAnnoK(index,f) + enrD->getOverlap(ele) - enrD->getCom1(m));

        enrD->calculateOddsRatio( a, b, c, d, ors, orL, orU );

        sprintf(bufferOR,"[%.2f,%.2f]",orL, orU);
        //---

        //--- print Anno size
        sprintf(buffer,"%d",enrD->getAnnoK(index,f));

        (printAn ? RES(m,k++) = string(buffer) : RES(m,k++) = "");

        double PermFrac = enrD->getPseudoCount()+enrD->getPermute(ele);
        PermFrac /= enrD->getNoP();
        PermFrac *= 100;

        if( printAlt ){
          RES(m,k++) = enrD->getOverlap(ele);
          RES(m,k++) = mn;
          RES(m,k++) = ors;
          RES(m,k++) = bufferOR;
          if( printTS ){
            RES(m,k++) = enrD->getPvalueT(ele);
          } else {
            RES(m,k++) = enrD->getPvalue(ele);
          }
          if( printTS ){
            RES(m,k++) = enrD->getPadjustedT(ele);
          } else {
            RES(m,k++) = enrD->getPadjusted(ele);
          }
          RES(m,k++) = PermFrac;
          (printTS ? RES(m,k++) = starsET : RES(m,k++) = starsEO);
          if(printTS){
            RES(m,k++) = enrD->getPvalueDT(ele);
          } else {
            RES(m,k++) = enrD->getPvalueD(ele);
          }
          if( printTS ){
            RES(m,k++) = enrD->getPadjustedDT(ele);
          } else {
            RES(m,k++) = enrD->getPadjustedD(ele);
          }
          (printTS ? RES(m,k++) = starsDT : RES(m,k++) = starsDO);
        } else {
          RES(m,k++) = a;
          RES(m,k++) = mn;
          RES(m,k++) = ors;
          RES(m,k++) = bufferOR;
          if( printTS ){
            RES(m,k++) = enrD->getPvalueT(ele);
          } else {
            RES(m,k++) = enrD->getPvalue(ele);
          }
          if( printTS ){
            RES(m,k++) = enrD->getPadjustedT(ele);
          } else {
            RES(m,k++) = enrD->getPadjusted(ele);
          }
          RES(m,k++) = PermFrac;
          (printTS ? RES(m,k++) = starsET : RES(m,k++) = starsEO);
        }

        }
    }

    //cout << "> done." << endl;

  } else { RES = Rcpp::CharacterMatrix(0,0); }

  return RES;

}
