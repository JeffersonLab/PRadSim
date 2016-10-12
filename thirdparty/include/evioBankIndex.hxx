//  evioBankIndex.hxx
//
// creates bank index for serialized event, indexes all banks including container banks
// eventually need to switch to std::tuple instead of custom struct
//
//
//  Author:  Elliott Wolin, JLab, 15-aug-2012



#ifndef _evioBankIndex_hxx
#define _evioBankIndex_hxx


#include "evioException.hxx"
#include "evioUtil.hxx"


namespace evio {

using namespace std;
using namespace evio;


// holds bank index info
typedef struct {
  int containerType;            // bank container type
  int contentType;              // bank content type
  int depth;                    // depth in hierarchy
  const uint32_t *bankPointer;  // pointer to first word of bank
  int bankLength;               // length of bank in 32-bit words
  const void *data;             // pointer to first word of data of bank
  int dataLength;               // length of data item array units of item
} bankIndex;
  

// compares tagNums to each other, first by tag, then by num
struct tagNumComp {
  bool operator() (const tagNum &lhs, const tagNum &rhs) const {
    if(lhs.first<rhs.first) {
      return(true);
    } else if(lhs.first>rhs.first) {
      return(false);
    } else {
      return(lhs.second<rhs.second);
    }    
  }
};


// bank index map and range
typedef multimap<tagNum,bankIndex,tagNumComp> bankIndexMap;
typedef pair< bankIndexMap::const_iterator, bankIndexMap::const_iterator > bankIndexRange;




//-----------------------------------------------------------------------
//-----------------------------------------------------------------------


/**
 * Creates bank index for serialized event.
 *
 * Note that a given tag/num may appear more than once in event and map.
 */
class evioBankIndex {

public:
  evioBankIndex(int maxDepth=0);
  evioBankIndex(const uint32_t *buffer, int maxDepth=0);
  virtual ~evioBankIndex();


public:
  bool parseBuffer(const uint32_t *buffer, int maxDepth);
  bool tagNumExists(const tagNum& tn) const;
  int tagNumCount(const tagNum& tn) const;
  bankIndexRange getRange(const tagNum& tn) const;
  bankIndex getBankIndex(const tagNum &tn) const throw(evioException);
  int getMaxDepth();


public:
  bankIndexMap tagNumMap;     /**<Holds index to one or more banks having tag/num.*/

private:
  int maxDepth;


public:
  /**
   * Returns length and pointer to data, NULL if container bank, bad tagNum or wrong data type.
   *
   * @param tn tagNum
   * @param pLen Pointer to int to receive data length, set to 0 upon error
   * @return Pointer to data, NULL on error
   */
  template <typename T> const T* getData(const tagNum &tn, int *pLen) throw (evioException) {

    bankIndexMap::const_iterator iter = tagNumMap.find(tn);

    if((iter!=tagNumMap.end()) && ((((*iter).second).contentType)==evioUtil<T>::evioContentType())) {
      *pLen=((*iter).second).dataLength;
      return(static_cast<const T*>(((*iter).second).data));
    } else {
      *pLen=0;
      return(NULL);
    }
  }


  /**
   * Returns length and pointer to data, assumes valid bankIndex
   *
   * @param bi bankIndex
   * @param pLen Pointer to int to receive data length, set to 0 for bad type
   * @return Pointer to data, NULL on bad type
   */
  template <typename T> const T* getData(const bankIndex &bi, int *pLen) throw (evioException) {

    if(bi.contentType==evioUtil<T>::evioContentType()) {
      *pLen=(bi.dataLength);
      return(static_cast<const T*>(bi.data));
    } else {
      *pLen=0;
      return(NULL);
    }
  }

};

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

} // namespace evio


#endif
  
