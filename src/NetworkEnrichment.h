#ifndef NETWORKENRICHMENT_H
#define NETWORKENRICHMENT_H

#include "Headers.h"
#include "buildSets.h"

class NetworkEnrichment : buildSets {
 public:
  NetworkEnrichment();
  NetworkEnrichment( string*, int, int, string*, int, int );
  NetworkEnrichment( string*, int, int, vector<string*>, vector<int>, vector<int> );
  ~NetworkEnrichment();

  //--------------------------------------------------------
  //  GET AND SET FUNCTIONS
  //--------------------------------------------------------
  void seedOffSet(bool=false, int=-1);
  
  void setNoP         ( int );
  void setPesudoCount ( double );

  void setPrintCNEW( bool );
  void setPrintID  ( bool );
  void setPrintAn  ( bool );
  void setALT      ( bool );

  void setExpectedOverlap( bool );
  void setFoldChange     ( bool );
  void setRelDist        ( bool );
  void setRCfisher       ( bool );
  void setChi2           ( bool );
  
  void oneSided();
  void twoSided();

  void maxSS();
  
  int   getKOffset();
  void  setKOffset( int = 1 );
  
  double getMINOVERLAP ( int );
  void   setMINOVERLAP ( int , int  );
  void   setANNOindex( int );
  void   setFDRmethod ( const char* );

  int    calculateOverlapinCommunities(bool,          const char*, const char*, bool = false, bool = false, bool = false );
  int    calculateOverlapinCommunities(int, int,      const char*, const char*, bool = true );
  int    calculateOverlapinNetwork    (int, int,      const char*, const char*, bool = true );
  int    calculateOverlapinNetwork    (int, int, int, const char*, const char*, bool = true );
  int    calculateOverlapinCommunities(int, bool,     const char*, const char*, bool = false, bool = false, bool = false );
  
  void calculateInteractionDistance(int, int, int, int, int, int, int, int &, double &);

  void calculateSampleSpace( int, int, vector<pairIntInt> & );
  void calculateOddsRatio( double, double, double, double, double &, double &, double & );  

  int getN();
  int getM();
  int getF();
  int getFsize(int);
  int getSigmaSize();
  double calBonferroni(int, double);
  double getPvalue(int);
  double getPvalueT(int);
  double getPvalueD(int);
  double getPvalueDT(int);
  double getPadjusted(int);
  double getPadjustedT(int);
  double getPadjustedD(int);
  double getPadjustedDT(int);
  double getOverlap(int);
  double getPchi2(int);
  double getPadjustedChi2(int);
  int    getAnnoK(int, int);
  string getAnnoID(int, int);
  string getAnnoName(int, int);
  string getAnnoDes(int, int);
  int getCom0(int);
  int getCom1(int);
  int getCom2(int);
  double getMuCab(int);
  int getKOFFSET();
  int getANNOindex();
  int getPermute(int);
  int getNoP();
  int getPseudoCount();
  int getBufferSize();
  
 private:
  void               freeMemory();
  void               setSeed(bool=false, int=0);
  unsigned long int  getSeed();
  bool _min( double, double );
  bool _max( double, double );

  void geneAssociations( int, int, int, int &, int & );
  void geneAssociations( int, int, int[] );
  void overlapinNetwork();
  void overlapinNetork( int, int );
  void overlapinNetork( int, int, int );
  void permutation( double );
  double prob_overlap( int, int, int, int );//Hypergeometric distribution
  double prob_overlap( int, int, int, int, int, int, int, int );//Hypergeometric distribution
  double prob_overlap( int, int, int, int, int );//intersection between three sets

  void overlapinComsHypergeometricTestRnd(bool = true);
  void overlapinComsHypergeometricTest   ();
  void overlapinComsHypergeometricTest   (int, int);
  void overlapinNetHypergeometricTest    (int, int);
  void overlapinNetHypergeometricTest    (int, int, int);  
  void CalculateFDR_BY( vector<pairDoubInt>, double, double &, double &, int &, vector<pairDoubInt> &);
  void CalculateFDR_BH( vector<pairDoubInt>, double, double &, double &, int &, vector<pairDoubInt> &);
  void CalculateFDR_BL( vector<pairDoubInt>, double, double &, double &, int &, vector<pairDoubInt> &);
  void calculateFDR( int=1, int=-1, int=-1, int=-1 );

  void printOverlapinCommunities   (const char *, const char *, bool, bool = false);
  void printOverlapinCommunitiesAlt(const char *, const char *, bool);
  void printOverlapinCommunities   (int, int, const char *, const char *, bool);
  void printOverlapinNetwork       (int, int, const char *, const char *, bool);
  void printOverlapinNetwork       (int, int, int, const char *, const char *, bool);
  
  void chi2_rxc       ( double &, double *, int, int, bool = true, bool = false );

  //set output file delineator(s), in print functions
  static const int DELSIZE = 1;
  char dels[DELSIZE]; //use for output file delineator(s
  
  
  //GSL random number and seed
  unsigned long int seed;
  unsigned long int seedOffset;
  gsl_rng *g;

vector<string> baseNAME; //base names for print function, can change.
  //vector<string> setNames;
  vector<string> FDRmethods;
  
  vector<bool> freedMemory;//flags to indicate if memory needs freeing 
  int ANNOindex;           //Selected annotation set index
  int NoP;                 //No: of permutation studies to run  
  int FDRtest;             //Select which FDR test to run
  int KOFFSET;             //Value added to community number  
  bool isOFFSET;           //Flag if we have offset community numbers

  bool printMeanMu;        //Flag expected overlap  [Default true]
  bool printFC;            //Flag print fold-change [Default false]
  bool printALT;           //Flag print alternative enrichment / depletion [Default true]

  bool printCnew;
  bool printID;            //Flag where we print annotation IDs or descriptions.
  bool printAn;            //Flag where we print annotation Size.

  bool printOneSided;      //print enrichment / depletion of one-sided
  bool printTwoSided;      //print enrichment / depletion of two-sided [Default] 

  bool calRelDist;

  bool useMaxSS;

  bool useRCfisher;
  bool useChi2;
  
  double FDR;
  double PV;
  double SIGMA;
  int    LEVEL;
  int    TESTS;

  double pesudocount; //to avoid p-values of zero, in the permutation test
  
  //global variable values
  static const int MAXEXPO     = 700; //see prob_overlap functions
  static const int SIGMASIZE   = 3;
  static const int OVERLAPSIZE = 4; 
  static const int BUFFERSIZE  = 250;
  
  double sigma[SIGMASIZE];
  double bonferroni[SIGMASIZE];
  double MINOVERLAP[OVERLAPSIZE];

  //arrays for Hypergeometric tests  
  int* comSIZE;
  int* geneCOM;
  int* annoSIZE;
  
  double* studies;
  double* overlap;
  double* muab;
  double* muCab;
  double* nab;

  double* p_values;
  double* padjusted;
  double* permute;
  
  double* p_valuesD;
  double* padjustedD;
  double* permuteD;

  double* p_valuesDT;
  double* padjustedDT;
  double* permuteDT;
  
  double* p_valuesT;
  double* padjustedT;
  double* permuteT;

  double* p_dist;
  double* padjustedRD;
  double* reldist;

  double* p_exfisher;
  double* padjustedEXF;
  double* p_chi2;
  double* padjustedCHI2;
  
 protected:

};

#endif
