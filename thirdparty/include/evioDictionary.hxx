// evioDictionary.hxx
//
// parses XML and creates a pair of name<-->tagNum maps
//
// probably not thread-safe due to static expat handlers
//
//
// basic dictionary element:
//    <dictEntry name="" tag="" num="">
//
// This simply makes entries in the two maps. It can occur anywhere in the XML.
//
//
//
// hierarchical dictionary elements:
//    <bank name="" tag="" num="">
//      <bank name="" tag="" num="">
//         <leaf name="" tag="" num=""/>
//         <leaf name="" tag="" num=""/>
//         <leaf name="" tag="" num=""/>
//      </bank>
//    </bank>
//
// Here the name entered in the map reflects the position of the bank or leaf in
//   the full hierarchy.  The full name is a concatenation of the hierarchy of names
//   with a separator character between them (e.g. '.' or '/')
//
//
//
//  Author:  Elliott Wolin, JLab, 13-apr-2012



#ifndef _evioDictionary_hxx
#define _evioDictionary_hxx


#include "evioTypedefs.hxx"
#include "evioException.hxx"

#include <expat.h>



namespace evio {

using namespace std;
using namespace evio;


// xml tag holding straight dictionary entries
const string dictEntryTag = "dictentry";



/**
 * Parses XML dictionary string and holds two maps, one for each lookup direction.
 */
class evioDictionary {

public:
  evioDictionary();
  evioDictionary(const string &dictXML, const string &sep=".");
  evioDictionary(ifstream &dictIFS, const string &sep=".");
  virtual ~evioDictionary();


public:
  bool parseDictionary(const string &dictionaryXML);
  tagNum getTagNum(const string &name) const throw(evioException);
  string getName(tagNum tn) const throw(evioException);
  string getDictionaryXML(void) const;
  void setSeparator(const string &sep);
  string getSeparator(void) const;
  string toString(void) const throw(evioException);


private:
  static void startElementHandler(void *userData, const char *xmlname, const char **atts);
  static void endElementHandler(void *userData, const char *xmlname);
  string dictionaryXML;
  string separator;
  string parentPrefix;
  bool parentIsLeaf;


public:
  map<tagNum,string> getNameMap;     /**<Gets node name given tag/num.*/
  map<string,tagNum> getTagNumMap;   /**<Gets tag/num given node name.*/
};


}

#endif
  
