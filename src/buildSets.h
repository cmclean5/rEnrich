//define guards, so headers are declare only once.
#ifndef BUILDSETS_H
#define BUILDSETS_H

#include "Headers.h"

//structs
//define useful list struct objects needed for annotation files
struct LISTst {
  char *annoID;
  char *annoDES;
  int   ID;
  int   K;
};

struct sortpairDoubInt{
  bool operator()(const std::pair<double,int> &l,
                  const std::pair<double,int> &r) const {
    return l.first < r.first;
    }
};


struct sortpairIntInt{
  bool operator()(const std::pair<int,int> &l,
                  const std::pair<int,int> &r) const {

    if( l.first < r.first ) return true;
    if( l.first > r.first ) return false;
    return l.second < r.second;    

    }
};


struct sortpairStrInt{
  bool operator()(const std::pair<std::string,int> &l,
                  const std::pair<std::string,int> &r) const {

    if( l.first.compare(r.first) < 0 ) return true; 
    if( l.first.compare(r.first) > 0 ) return false;
    return l.second < r.second;
  
    }
};


class buildSets {
   
public:
  buildSets();
  ~buildSets();

  void addSets(string*, int, int, string*, int, int ); 
  
private:
  
  //functions
  void    freeSpace();
  LISTst *createList( string*, int, int );
  void    freqofComslist( bool = false, int = 0 );
  LISTst *freqofAnnolist( LISTst*, int, int & );
  LISTst *removeDuplicateIDs( LISTst *, int &, LISTst *, int );
  
  //variables
  bool    freedClist; 
  bool    freedAlist;
  bool    freedANNOS;
  
  
 protected:
  void addKOffset( int = 0 );
  
  typedef std::pair<double,int> pairDoubInt;
  typedef std::pair<int,int>    pairIntInt;
  typedef std::pair<string,int> pairStrInt;

  typedef std::tuple<int,int,int> tripleInt;
    

  //for membership file
  int     Mmin_old;
  int     Mmax_old;
  int     M_old;
  int     Mmin;
  int     Mmax;
  int     Clines;
  int     Ccols;
  LISTst *Clist;
  vector<tripleInt> COMS; 

  //for annotation files
  vector<int>     Alines;
  int             Acols;
  vector<LISTst*> Alist;
  vector<LISTst*> ANNOS;

  //make these protected
  int N;
  int M;
  int F;
  vector<int> Fsize;
    
};

#endif
