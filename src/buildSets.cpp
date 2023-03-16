#include "buildSets.h"

buildSets::buildSets() {

  this->M        = 0;
  this->N        = 0;
  this->F        = 0;
  this->Mmin_old = 0;
  this->Mmax_old = 0;
  this->Mmin     = 0;
  this->Mmax     = 0; 
  
  this->Clines   = 0;
  this->Ccols    = 0;
  this->Clist    = 0;
  
  this->freedClist = false;
  this->freedAlist = false;
  this->freedANNOS = false;
  
}

buildSets::~buildSets() { freeSpace(); }


void buildSets::addSets( string* m,    int m_rows, int m_cols,
                         string* anno, int a_rows, int a_cols ){

  int i,Ff;
  
  this->M          = 0;
  this->N          = 0;
  this->F          = 0;
  this->M_old      = 0;
  this->Mmin_old   = 0;
  this->Mmax_old   = 0;
  this->Mmin       = 0;
  this->Mmax       = 0;
  
  this->freedClist = false;
  this->freedAlist = false;
  this->freedANNOS = false;
 
  
  //read-in membership data.frame
  Clines = m_rows;
  Ccols  = m_cols;  
  Clist  = createList( m, m_rows, m_cols );

  //read-in annotation data.frame 
  this->Alines.push_back(a_rows);
  this->Alist.push_back(createList( anno, a_rows, a_cols ));

  //calculate unqiue number of annotations types & the sizes.  
  freqofComslist();

  for( i=0; i<Alines.size(); i++ ){
    Ff = 0;
    this->ANNOS.push_back( freqofAnnolist(Alist[i], Alines[i], Ff) );
    this->Fsize.push_back( Ff );
  }
  
  //check for duplicate Entrez IDs in each annotation type, which may exist
  //in the input annotation files.
  for( i=0; i<Alines.size(); i++ ){
    Alist[i] = removeDuplicateIDs( Alist[i], Alines[i], ANNOS[i], Fsize[i] );    
  }

  N = Clines;
  M = COMS.size();
  F = Fsize[0];
    
  
}

void buildSets::addSets( string* m, int m_rows, int m_cols,
                         vector<string*> anno, vector<int> a_rows, vector<int> a_cols ){

  int i,Ff;
  
  this->M          = 0;
  this->N          = 0;
  this->F          = 0;
  this->M_old      = 0;
  this->Mmin_old   = 0;
  this->Mmax_old   = 0;
  this->Mmin       = 0;
  this->Mmax       = 0;
  
  this->freedClist = false;
  this->freedAlist = false;
  this->freedANNOS = false;
 
  //read-in membership data.frame
  Clines = m_rows;
  Ccols  = m_cols;  
  Clist  = createList( m, m_rows, m_cols );
  
  //read-in annotation data.frames
  for(i=0; i<anno.size(); i++){
    if( a_rows[i] > 0 ){
      this->Alines.push_back(a_rows[i]);
      this->Alist.push_back(createList( anno[i], a_rows[i], a_cols[i] ));
    }
  }

  //calculate unqiue number of annotations types & the sizes.  
  freqofComslist();

  for( i=0; i<Alines.size(); i++ ){
    Ff = 0;
    this->ANNOS.push_back( freqofAnnolist(Alist[i], Alines[i], Ff) );
    this->Fsize.push_back( Ff );
  }

  //check for duplicate Entrez IDs in each annotation type, which may exist
  //in the input annotation files.
  for( i=0; i<Alines.size(); i++ ){
    Alist[i] = removeDuplicateIDs( Alist[i], Alines[i], ANNOS[i], Fsize[i] );    
  }

  N = Clines;
  M = COMS.size();
  F = Fsize[0];
      
}



LISTst* buildSets::createList( std::string *data, int NROWS, int NCOLS ){
  
  int i,j,k,data_size;

  string v1,v2,v3;
  
  LISTst *LIST = (LISTst*)calloc(NROWS,sizeof(LISTst));

  data_size = NROWS*NCOLS;

  
  if( NCOLS == 3 ){
    //annotation data.frame
      
    for(k=0; k<data_size; k++){
      i = floor(k/NCOLS);
      j = k % NCOLS;
	
      if( j == 0 ){ 
        v1=data[(i*NCOLS)+j]; 
        LIST[i].annoID = (char*)malloc((strlen(v1.c_str())+1)*sizeof(char));
        strcpy(LIST[i].annoID,v1.c_str());
      }
      
      if( j == 1 ){ 
        v2=data[(i*NCOLS)+j];
        LIST[i].annoDES = (char*)malloc((strlen(v2.c_str())+1)*sizeof(char));
        strcpy(LIST[i].annoDES,v2.c_str()); 
      }

      if( j == 2 ){ 
        v3=data[(i*NCOLS)+j];
        LIST[i].ID = atoi(v3.c_str());
      }

    }

  }
    
  if( NCOLS == 2 ){
    // membership data.frame
    
    for(k=0; k<data_size; k++){
      i = floor(k/NCOLS);
      j = k % NCOLS;

      if( j == 0 ){
        v1=data[(i*NCOLS)+j]; 
        LIST[i].ID = atoi(v1.c_str());
      }

      if( j == 1 ){
        v2=data[(i*NCOLS)+j];
        LIST[i].K = atoi(v2.c_str());
      }

    }

  }
  
  return LIST;
  
}


void buildSets::freqofComslist( bool addOffset, int Koffset ){

  int i, j, counter, Knew;
  bool found;

  int *temp = (int*)calloc( Clines,sizeof(int) );
  
  Mmin=Clist[0].K;
  Mmax=Clist[0].K;

  for( i=0; i<Clines; i++ ){

    temp[i] = -1;
    
    if( Clist[i].K > Mmax ) Mmax = Clist[i].K;
    if( Clist[i].K < Mmin ) Mmin = Clist[i].K;
    
  }

  //--- reorder communities from 1 to (Max-Mmin)
  Knew    = 1;
  counter = Mmin;
  
  while( counter <= Mmax ){

    COMS.push_back( tripleInt(0, 0, counter) );
    
    found=false;
    for( i=0; i<Clines; i++ ){
      if( Clist[i].K == counter ){
	temp[i] = Knew;
	found=true;
      }
    }

    if(found) Knew++;

    counter++;

  }
  //----
  
  
  //--- unique list of communities and sizes
  for( i=0; i<Clines; i++ ){
    for( j=0; j<COMS.size(); j++ ){
      if( std::get<2>(COMS[j]) == Clist[i].K ){
	std::get<1>(COMS[j]) = std::get<1>(COMS[j]) + 1;
	std::get<0>(COMS[j]) = temp[i];
	continue;
      }	
    }
  }
  //----

  //--- save the rescaled community numbers
  for(i=0; i<Clines; i++){
    Clist[i].K = temp[i];
  }
  
  //delete temp array
  if( temp != 0 ){ free(temp); }

  
  //-- remove communities with zero size
  for (i=0; i < COMS.size(); i++){
    if( std::get<1>(COMS[i]) == 0 ){
      COMS.erase(COMS.begin() + i);
      i--;
    }
  }
 
  //Store number of communities
  M = COMS.size();  

}


LISTst* buildSets::freqofAnnolist( LISTst *_ALIST, int _ALINES, int &_F ){

  int i,j;

  vector<pairStrInt> ANNOIDS; 

  for( i=0; i<_ALINES; i++ ){
    if(_ALIST[i].annoID != NULL){//std::string can't handle NULL, undefinded behaviour
      ANNOIDS.push_back(pairStrInt(string(_ALIST[i].annoID),0));
    }
  }

  //---unique list of anno IDs & names
  sort(ANNOIDS.begin(),ANNOIDS.end(),sortpairStrInt());
  ANNOIDS.erase(unique(ANNOIDS.begin(),ANNOIDS.end()),ANNOIDS.end());

  _F = ANNOIDS.size();  
  LISTst *LIST = (LISTst*)calloc(_F,sizeof(LISTst));

  for( i=0; i<ANNOIDS.size(); i++ ){
    LIST[i].annoID = (char*)malloc((strlen(ANNOIDS[i].first.c_str())+1)*sizeof(char));
    strcpy(LIST[i].annoID, ANNOIDS[i].first.c_str());
    LIST[i].K  = 0;
    LIST[i].ID = 0;
  }

  for( i=0; i<_ALINES; i++ ){
    if(_ALIST[i].annoID != NULL){//std::string can't handle NULL, undefinded behaviour
      for( j=0; j<ANNOIDS.size(); j++ ){
	if( ANNOIDS[j].first.compare(_ALIST[i].annoID) == 0 ){
	  LIST[j].annoDES = (char*)malloc((strlen(_ALIST[i].annoDES)+1)*sizeof(char));
	  strcpy(LIST[j].annoDES, _ALIST[i].annoDES);
	  //LIST[j].K++;
	  //ANNOIDS[j].second++; continue;
	  continue;
	}	
      }

    }
  }

  vector<pairStrInt>().swap(ANNOIDS);  //free space of vector
  
  return LIST;
  
}

LISTst* buildSets::removeDuplicateIDs( LISTst *_ALIST, int &_ALINES, LISTst *_ANNOS, int _F ){

  int i,j;

  vector<pairStrInt> ANNOIDS; 

  for( i=0; i<_ALINES; i++ ){
    if(_ALIST[i].annoID != NULL){//std::string can't handle NULL, undefinded behaviour
      ANNOIDS.push_back(pairStrInt(string(_ALIST[i].annoID),_ALIST[i].ID));
    }
  }

  //---unique list of anno IDs & gene IDs
  sort(ANNOIDS.begin(),ANNOIDS.end(),sortpairStrInt());
  ANNOIDS.erase(unique(ANNOIDS.begin(),ANNOIDS.end()),ANNOIDS.end());
  
  //---delete _ALIST
  for(i=0; i<_ALINES; i++){
    if(_ALIST[i].annoID  != NULL) free(_ALIST[i].annoID);
    if(_ALIST[i].annoDES != NULL) free(_ALIST[i].annoDES);
  }
  free(_ALIST);
  _ALINES = 0;
  //---
  
  //---create new _ALIST struct
  _ALINES          = ANNOIDS.size();
  LISTst *ALISTtmp = (LISTst*)calloc(_ALINES,sizeof(LISTst));
  //---
  
  //loop through the set of annotation types, and fill _ALIST
  for( i=0; i<_ALINES; i++ ){

    for( j=0; j<_F; j++ ){

      if( ANNOIDS[i].first.compare(_ANNOS[j].annoID) == 0 ){

     ALISTtmp[i].annoID = (char*)malloc((strlen(_ANNOS[j].annoID)+1)*sizeof(char));
     strcpy(ALISTtmp[i].annoID, _ANNOS[j].annoID);
	
     ALISTtmp[i].annoDES = (char*)malloc((strlen(_ANNOS[j].annoDES)+1)*sizeof(char));
     strcpy(ALISTtmp[i].annoDES, _ANNOS[j].annoDES);
	
     ALISTtmp[i].ID = ANNOIDS[i].second; 	
     _ANNOS[j].K++;
	
      }

    }
  }

  vector<pairStrInt>().swap(ANNOIDS);  //free space of vector

  return(ALISTtmp);
  
}


void buildSets::freeSpace(){
  
 int i,j;

 if( Clist != 0 && freedClist == false ){ 
   for(i=0; i<Clines; i++){ 
     if(Clist[i].annoID != NULL) free(Clist[i].annoID);
     if(Clist[i].annoDES != NULL)free(Clist[i].annoDES);      
   }
   free(Clist);
   freedClist = true;
 }

 
 if( freedAlist == false ){
   for(i=0; i<Alines.size(); i++){
     for(j=0; j<Alines[i]; j++){
       if(Alist[i][j].annoID  != NULL) free(Alist[i][j].annoID);
       if(Alist[i][j].annoDES != NULL) free(Alist[i][j].annoDES);
     }
     free(Alist[i]);
   }
   freedAlist = true;
 }

 if( freedANNOS == false ){
   for(i=0; i<Fsize.size(); i++){
     for(j=0; j<Fsize[i]; j++){
       if(ANNOS[i][j].annoID  != NULL) free(ANNOS[i][j].annoID);
       if(ANNOS[i][j].annoDES != NULL) free(ANNOS[i][j].annoDES);
     }
     free(ANNOS[i]);
   }
   freedANNOS = true;
 }
  
 
}


void buildSets::addKOffset( int addOffset ){

  M_old    = M;
  Mmin_old = Mmin;
  Mmax_old = Mmax;

  COMS.clear();
  
  if( addOffset < 0 ){ addOffset = 0; }
  
  //consider K = 0
  //freqofComslist( true, addOffset );

  /*
  cout << "------------------- "<< endl;
  cout << "Before resetting"    << endl;
  cout << "M    : " << M_old    << endl;
  cout << "K min: " << Mmin_old << endl;
  cout << "K max: " << Mmax_old << endl;

  cout << "After resetting" << endl;
  cout << "M    : " << M    << endl;
  cout << "K min: " << Mmin << endl;
  cout << "K max: " << Mmax << endl;
  cout << "------------------- "<< endl;
  */
  
}

