// evioUtil.hxx

//  Author:  Elliott Wolin, JLab, 26-apr-2012


//  must do:
//   check for duplicate map entries in dictionary, evioBankIndex, etc.
//   append mode and random access i/o
//   update word doc or docx (conflict?)


//  should do:
//   dictionary checksum
//   std::shared_ptr
//   API for creating and reading composite banks and data
//   multi-threading:  defaultToStringconfig?
//   who does dictionary belong to?  tree?  toStringConfig?


//  would like to do:
//   improve java evio tree query and manipulation compatibility


// not sure:
//   can string array contain NULL string?
//   change how evioToStringConfig and default config used
//   evioParseBankAt(int index_of_bank_in_event), i.e. random access to banks in events?


// probably not:
//   cMsg channel
//   ET channel
//   XML channel and xml2evio
//   decompress/compress on input/output (gzip, bzip, etc.)
//   pipes, named pipes on input
//   scheme for exception type codes
//   convert vectors to arrays



#ifndef _evioUtil_hxx
#define _evioUtil_hxx


#include <list>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <functional>
#include <memory>
#include <utility>
#include <cstring>
#include <typeinfo>
#include <expat.h>

#ifdef vxworks
#include <iostream.h>
#else
#include <iostream>
#include <iomanip>
#endif

#include "evio.h"
#include "evioTypedefs.hxx"
#include "evioException.hxx"
#include "evioDictionary.hxx"
#include "evioChannel.hxx"



#ifdef sun
#define __FUNCTION__ "unknown"
#endif



/** @mainpage  EVIO Event I/O Package.
 * @author Elliott Wolin.
 * @version 2.0.
 * @date 5-Feb-2007.
 *
 * @section intro Introduction
 * The EVIO package is an object-oriented extension of the original EVIO C event I/O utility.  
 *
 * The base utility reads and writes EVIO format events in an event buffer to and from disk.  
 * Events are blocked into standard block sizes, and endian swapping is performed if necessary.
 *
 * This package maps the event in the event buffer into an in-memory tree-like bank hierarchy.  
 * Event trees can be queried, modified, or created from scratch.  In-memory trees can be automatically 
 * serialized into a buffer and written to disk.
 * 
 * @warning Internally EVIO uses only the unambiguous types int8_t, uint8_t, ... int64_t, uint64_t.  
 * The long and long long data types are not supported, as their interpretion varies among different compilers
 * and architectures.  The unambiguous types are compatible with char, short, and int wherever I have checked.
 * But only the 64-bit types int64_t and uint64_t work consistently across architectures.
 *
 * @warning The safest route, of course, is to use the unambiguous types exclusively.  
 * But the example programs use char, short, int, and int64_t and so far work fine.  No guarantees 
 * for the future, though...
 */


/** All evio symbols reside in the evio namespace.*/
namespace evio {

using namespace std;
using namespace evio;


// evio classes
class evioStreamParserHandler;
class evioStreamParser;
class evioDOMTree;
class evioDOMNode;
class evioDOMContainerNode;
template<typename T> class evioDOMLeafNode;
class evioSerializable;
template <typename T> class evioUtil;
class evioToStringConfig;



//-----------------------------------------------------------------------------
//--------------------------- evio Classes ------------------------------------
//-----------------------------------------------------------------------------


/**
 * General utilities.
 */
class evioUtilities {

public:
  static uint32_t *appendToBuffer(const uint32_t *buffer, ContainerType bufferType, const uint32_t *structure, ContainerType structureType)
    throw(evioException);
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Configuration options for toString() method.
 *   max_depth:          depth to convert to, 0 means no limit, default 0.
 *   no_data:            true to not dump data of leaf nodes, default false.
 *   indentSize          indent size for increasing bank level depth, default 3
 *   toStringDictionary  dictionary to use, overrides default dictionary
 */
class evioToStringConfig {

public:
  evioToStringConfig();
  evioToStringConfig(const evioDictionary *dictionary);
  evioToStringConfig(const evioDictionary &dictionary);
  virtual ~evioToStringConfig();


private:
  evioToStringConfig(const string &s) {};  // eliminates auto-coercion of string to dictionary


public:
  bool xtod;                     /**<True print unsigned values as decimal.*/
  bool noData;                   /**<True do not dump leaf node data.*/
  int maxDepth;                  /**<Max depth to dump.*/
  int indentSize;                /**<Indent size per unit of depth.*/
  bool verbose;                  /**<Turn on verbose mode.*/
  vector<uint16_t> bankOk;       /**<Vector of bank tags to dump.*/
  vector<uint16_t> noBank;       /**<Vector of bank tags to skip.*/
  vector<string> bankNameOk;     /**<Vector of bank names to dump.*/
  vector<string> noBankName;     /**<Vector of bank names to skip.*/

  const evioDictionary *toStringDictionary;   /**<Dictionary to use.*/


protected:
  void init(void);


public:
  virtual void setDictionary(const evioDictionary *dict)  {toStringDictionary=dict;}
  virtual void setDictionary(const evioDictionary &dict)  {toStringDictionary=&dict;}
  virtual const evioDictionary *getDictionary(void) const {return(toStringDictionary);}
  bool skipNode(const evioDOMNodeP pNode) const;
};



// default toString config
static evioToStringConfig defaultToStringConfig;



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Interface defines node and leaf handlers for use with evioStreamParser.
 * Separate handlers defined for container nodes and leaf nodes.
 */
class evioStreamParserHandler {

public:
  virtual void *containerNodeHandler(int bankLength, int containerType, int contentType, uint16_t tag, uint8_t num, 
                                     int depth, const uint32_t *bankPointer, int payloadLength, const uint32_t *payload, void *userArg) = 0;
  virtual void *leafNodeHandler(int bankLength, int containerType, int contentType, uint16_t tag, uint8_t num, 
                                int depth, const uint32_t *bankPointer, int dataLength, const void *data, void *userArg) = 0;
  virtual ~evioStreamParserHandler(void) {};
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Stream parser dispatches to evioStreamParserHandler handlers when node or leaf reached.
 */
class evioStreamParser {

public:
  void *parse(const uint32_t *buf, evioStreamParserHandler &handler, void *userArg)
    throw(evioException);
  virtual ~evioStreamParser(void) {};

  
private:
  void *parseBank(const uint32_t *buf, int bankType, int depth, 
                  evioStreamParserHandler &handler, void *userArg) throw(evioException);

};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/** 
 * Virtual class represents an evio node in memory, concrete sub-classes evioDOMContainerNode and evioDOMLeafNode
 *   are hidden from users.
 * Users work with nodes via this class, and create nodes via factory method createEvioDOMNode.
 * Factory model ensures nodes are created on heap.  
 */
class evioDOMNode {


  friend class evioDOMTree;    /**<Allows evioDOMTree class to manipulate nodes.*/


protected:
  evioDOMNode(evioDOMNodeP parent, uint16_t tag, uint8_t num, int contentType) throw(evioException);
  evioDOMNode(evioDOMNodeP parent, const string &name, const evioDictionary *dictionary, int contentType) throw(evioException);


public:
  virtual ~evioDOMNode(void);


private:
  evioDOMNode(const evioDOMNode &node) throw(evioException);
  bool operator=(const evioDOMNode &node) const {return(false);}


// public factory methods for node creation
public:
  static evioDOMNodeP createEvioDOMNode(uint16_t tag, uint8_t num, ContainerType cType=BANK) throw(evioException);
  template <typename T> static evioDOMNodeP createEvioDOMNode(uint16_t tag, uint8_t num) throw(evioException);
  template <typename T> static evioDOMNodeP createEvioDOMNode(uint16_t tag, uint8_t num, const vector<T> &tVec) throw(evioException);
  template <typename T> static evioDOMNodeP createEvioDOMNode(uint16_t tag, uint8_t num, const T* t, int len) throw(evioException);
  static evioDOMNodeP createEvioDOMNode(uint16_t tag, uint8_t num, const evioSerializable &o, ContainerType cType=BANK) 
    throw(evioException);
  static evioDOMNodeP createEvioDOMNode(uint16_t tag, uint8_t num, void (*f)(evioDOMNodeP c, void *userArg), void *userArg, 
                                        ContainerType cType=BANK) throw(evioException);
  template <typename T> static evioDOMNodeP createEvioDOMNode(uint16_t tag, uint8_t num, T *t,
                                                              void *userArg, ContainerType cType=BANK) throw(evioException);
  template <typename T> static evioDOMNodeP createEvioDOMNode(uint16_t tag, uint8_t num, T *t, 
                                                              void* T::*mfp(evioDOMNodeP c, void *userArg),
                                                              void *userArg, ContainerType cType=BANK) throw(evioException);

  static evioDOMNodeP createEvioDOMNode(uint16_t tag, uint8_t num, uint16_t formatTag, const string &formatString,
                                        uint16_t dataTag, uint8_t dataNum, const vector<uint32_t> &tVec) throw(evioException);

  static evioDOMNodeP createEvioDOMNode(uint16_t tag, uint8_t num, uint16_t formatTag, const string &formatString,
                                        uint16_t dataTag, uint8_t dataNum, const uint32_t* t, int len) throw(evioException);

  static evioDOMNodeP createUnknownEvioDOMNode(uint16_t tag, uint8_t num, const vector<uint32_t> &tVec) throw(evioException);
  static evioDOMNodeP createUnknownEvioDOMNode(uint16_t tag, uint8_t num, const uint32_t *t, int len) throw(evioException);


  static evioDOMNodeP createEvioDOMNode(const string &name, const evioDictionary *dictionary, ContainerType cType=BANK) throw(evioException);
  template <typename T> static evioDOMNodeP createEvioDOMNode(const string &name, const evioDictionary *dictionary) throw(evioException);
  template <typename T> static evioDOMNodeP createEvioDOMNode(const string &name, const evioDictionary *dictionary, const vector<T> &tVec)
    throw(evioException);
  template <typename T> static evioDOMNodeP createEvioDOMNode(const string &name, const evioDictionary *dictionary, const T* t, int len)
    throw(evioException);
  static evioDOMNodeP createEvioDOMNode(const string &name, const evioDictionary *dictionary, const evioSerializable &o, ContainerType cType=BANK) 
    throw(evioException);
  static evioDOMNodeP createEvioDOMNode(const string &name, const evioDictionary *dictionary, void (*f)(evioDOMNodeP c, void *userArg),
                                        void *userArg, ContainerType cType=BANK) throw(evioException);
  template <typename T> static evioDOMNodeP createEvioDOMNode(const string &name, const evioDictionary *dictionary, T *t,
                                                              void *userArg, ContainerType cType=BANK) throw(evioException);
  template <typename T> static evioDOMNodeP createEvioDOMNode(const string &name, const evioDictionary *dictionary, T *t, 
                                                              void* T::*mfp(evioDOMNodeP c, void *userArg),
                                                              void *userArg, ContainerType cType=BANK) throw(evioException);
  static evioDOMNodeP createEvioDOMNode(const string &name, const evioDictionary *dictionary, uint16_t formatTag, const string &formatString,
                                        uint16_t dataTag, uint8_t dataNum, const vector<uint32_t> &tVec) throw(evioException);
  
  static evioDOMNodeP createEvioDOMNode(const string &name, const evioDictionary *dictionary, uint16_t formatTag, const string &formatString,
                                        uint16_t dataTag, uint8_t dataNum, const uint32_t* t, int len) throw(evioException);
  static evioDOMNodeP createUnknownEvioDOMNode(const string &name, const evioDictionary *dictionary, const vector<uint32_t> &tVec)
    throw(evioException);
  static evioDOMNodeP createUnknownEvioDOMNode(const string &name, const evioDictionary *dictionary, const uint32_t *t, int len)
    throw(evioException);


public:
  virtual void addNode(evioDOMNodeP node) throw(evioException);
  void append(const string &s) throw(evioException);
  void append(const char *s) throw(evioException);
  void append(char *s) throw(evioException);
  void append(const char **ca, int len) throw(evioException);
  void append(char **ca, int len) throw(evioException);

  template <typename T> void append(T tVal) throw(evioException);
  template <typename T> void append(const vector<T> &tVec) throw(evioException);
  template <typename T> void append(const T* tBuf, int len) throw(evioException);
  template <typename T> void replace(const vector<T> &tVec) throw(evioException);
  template <typename T> void replace(const T* tBuf, int len) throw(evioException);


public:
  virtual evioDOMNodeP cut(void) throw(evioException);
  virtual void cutAndDelete(void) throw(evioException);
  virtual evioDOMNodeP move(evioDOMNodeP newParent) throw(evioException);


public:
  virtual bool operator==(uint16_t tag) const;
  virtual bool operator!=(uint16_t tag) const;
  bool operator==(tagNum tnPair) const;
  bool operator!=(tagNum tnPair) const;


public:
  evioDOMNode& operator<<(evioDOMNodeP node) throw(evioException);
  evioDOMNode& operator<<(const string &s) throw(evioException);
  evioDOMNode& operator<<(const char *s) throw(evioException);
  evioDOMNode& operator<<(char *s) throw(evioException);
  template <typename T> evioDOMNode& operator<<(T tVal) throw(evioException);
  template <typename T> evioDOMNode& operator<<(const vector<T> &tVec) throw(evioException);


public:
  evioDOMNodeList *getChildList(void) throw(evioException);
  evioDOMNodeListP getChildren(void) throw(evioException);
  template <class Predicate> evioDOMNodeListP getChildren(Predicate pred) throw(evioException);
  template <typename T> vector<T> *getVector(void) throw(evioException);


public:
  virtual string toString(void) const;
  virtual string getHeader(int depth, const evioToStringConfig *config = &defaultToStringConfig) const = 0;
  virtual string getBody(int depth, const evioToStringConfig *config = &defaultToStringConfig) const = 0;
  virtual string getFooter(int depth, const evioToStringConfig *config = &defaultToStringConfig) const = 0;
  virtual int getSize(void) const = 0;

  const evioDOMNodeP getParent(void) const;
  int getContentType(void) const;
  const evioDOMTree *getParentTree(void) const;
  bool isContainer(void) const;
  bool isLeaf(void) const;


protected:
  static string getIndent(int depth, int size = 3);


protected:
  evioDOMNodeP parent;           /**<Pointer to node parent.*/
  evioDOMTree *parentTree;       /**<Pointer to parent tree if this node is the root.*/
  int contentType;               /**<Content type.*/


public:
  uint16_t tag;            /**<The node tag, max 16-bits depending on container type.*/
  uint8_t num;             /**<The node num, max 8 bits, used by BANK and String container types (2-word header).*/
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Sub-class of evioDOMNode represents an evio container node.
 * Only accessible to users via pointer to evioDOMNode object.
 */
class evioDOMContainerNode : public evioDOMNode {

  friend class evioDOMNode;    /**<Allows evioDOMNode to use private subclass methods.*/


protected:
  evioDOMContainerNode(evioDOMNodeP parent, uint16_t tag, uint8_t num, ContainerType cType) throw(evioException);
  virtual ~evioDOMContainerNode(void);

  evioDOMContainerNode(const evioDOMContainerNode &cNode) throw(evioException);
  bool operator=(const evioDOMContainerNode &node);


public:
  virtual string getHeader(int depth, const evioToStringConfig *config = &defaultToStringConfig) const;
  virtual string getBody(int depth, const evioToStringConfig *config = &defaultToStringConfig) const;
  virtual string getFooter(int depth, const evioToStringConfig *config = &defaultToStringConfig) const;
  virtual int getSize(void) const;


public:
  evioDOMNodeList childList;   /**<STL List of pointers to children.*/

};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Sub-class of evioDOMNode represents an evio leaf node.
 * Only accessible to users via pointer to evioDOMNode object.
 */
template <typename T> class evioDOMLeafNode : public evioDOMNode {

  friend class evioDOMNode;     /**<Allows evioDOMNode to use private subclass methods.*/


protected:
  evioDOMLeafNode(evioDOMNodeP par, uint16_t tag, uint8_t num) throw(evioException);
  evioDOMLeafNode(evioDOMNodeP par, uint16_t tag, uint8_t num, const vector<T> &v) throw(evioException);
  evioDOMLeafNode(evioDOMNodeP par, uint16_t tag, uint8_t num, const T *p, int ndata) throw(evioException);
  virtual ~evioDOMLeafNode(void);

  evioDOMLeafNode(const evioDOMLeafNode<T> &lNode) throw(evioException);
  bool operator=(const evioDOMLeafNode<T> &lNode);


public:
  virtual string getHeader(int depth, const evioToStringConfig *config = &defaultToStringConfig) const;
  virtual string getBody(int depth, const evioToStringConfig *config = &defaultToStringConfig) const;
  virtual string getFooter(int depth, const evioToStringConfig *config = &defaultToStringConfig) const;
  virtual int getSize(void) const;


public:
  vector<T> data;    /**<Vector<T> of node data.*/
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Sub-class of evioDOMLeafNode<uint32_t> represents composite leaf node.
 * Only accessible to users via pointer to evioDOMNode object.
 */
class evioCompositeDOMLeafNode : public evioDOMLeafNode<uint32_t> {


  friend class evioDOMNode;     /**<Allows evioDOMNode to use private subclass methods.*/


protected:
  evioCompositeDOMLeafNode(evioDOMNodeP par, uint16_t tag, uint8_t num, uint16_t formatTag, const string &formatString, 
                           uint16_t dataTag, uint8_t dataNum, const vector<uint32_t> &v) throw(evioException);
  evioCompositeDOMLeafNode(evioDOMNodeP par, uint16_t tag, uint8_t num, uint16_t formatTag, const string &formatString, 
                           uint16_t dataTag, uint8_t dataNum, const uint32_t *p, int ndata) throw(evioException);
  ~evioCompositeDOMLeafNode(void);

  evioCompositeDOMLeafNode(const evioCompositeDOMLeafNode &lNode) throw(evioException);
  bool operator=(const evioCompositeDOMLeafNode &lNode);


public:
  virtual string getBody(int depth, const evioToStringConfig *config = &defaultToStringConfig) const;
  virtual int getSize(void) const;


public:
  uint16_t formatTag;      /**<Tag to use for the internal format bank.*/
  string formatString;     /**<The format string.*/
  uint16_t dataTag;        /**<Tag to use for the internal data bank.*/
  uint8_t dataNum;         /**<Num to use for the internal data bank.*/
};



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Represents an evio tree/event in memory.
 * Tree root is an evioDOMNode.
 */
class evioDOMTree : public evioStreamParserHandler, public evioChannelBufferizable {


public:
  evioDOMTree(void) throw(evioException);
  evioDOMTree(evioDictionary *dictionary) throw(evioException);
  evioDOMTree(const evioChannel &channel, const string &name = "evio") throw(evioException);
  evioDOMTree(const evioChannel *channel, const string &name = "evio") throw(evioException);
  evioDOMTree(const uint32_t *buf, const string &name = "evio") throw(evioException);
  evioDOMTree(evioDOMNodeP node, const string &name = "evio") throw(evioException);
  evioDOMTree(uint16_t tag, uint8_t num, ContainerType cType=BANK, const string &name = "evio") throw(evioException);
  evioDOMTree(const string &bankName, ContainerType cType=BANK, const string &name = "evio") throw(evioException);
  evioDOMTree(tagNum tn, ContainerType cType=BANK, const string &name = "evio") throw(evioException);
  evioDOMTree(const string &bankName, evioDictionary *dictionary, ContainerType cType=BANK, const string &name = "evio") throw(evioException);
  virtual ~evioDOMTree(void);


private:
  evioDOMTree(const evioDOMTree &tree) throw(evioException);
  bool operator=(const evioDOMTree &tree);


public:
  void clear(void) throw(evioException);
  void addBank(evioDOMNodeP node) throw(evioException);
  template <typename T> void addBank(uint16_t tag, uint8_t num, const vector<T> &dataVec) throw(evioException);
  template <typename T> void addBank(uint16_t tag, uint8_t num, const T* dataBuf, int dataLen) throw(evioException);
  template <typename T> void addBank(tagNum tn, const vector<T> &dataVec) throw(evioException);
  template <typename T> void addBank(tagNum tn, const T* dataBuf, int dataLen) throw(evioException);
  template <typename T> void addBank(const string &name, const vector<T> &dataVec) throw(evioException);
  template <typename T> void addBank(const string &name, const T* dataBuf, int dataLen) throw(evioException);

  void addBank(uint16_t tag, uint8_t num, uint16_t formatTag, const string &formatString, 
               uint16_t dataTag, uint8_t dataNum, const vector<uint32_t> &dataVec) throw(evioException);
  void addBank(uint16_t tag, uint8_t num, uint16_t formatTag, const string &formatString, 
               uint16_t dataTag, uint8_t dataNum, const uint32_t *t, int len) throw(evioException);
  void addBank(tagNum tn, uint16_t formatTag, const string &formatString, 
               uint16_t dataTag, uint8_t dataNum, const vector<uint32_t> &dataVec) throw(evioException);
  void addBank(tagNum tn, uint16_t formatTag, const string &formatString, 
               uint16_t dataTag, uint8_t dataNum, const uint32_t *t, int len) throw(evioException);
  void addBank(const string &name, uint16_t formatTag, const string &formatString, 
               uint16_t dataTag, uint8_t dataNum, const vector<uint32_t> &dataVec) throw(evioException);
  void addBank(const string &name, uint16_t formatTag, const string &formatString, 
               uint16_t dataTag, uint8_t dataNum, const uint32_t *t, int len) throw(evioException);


  evioDOMNodeP createNode(const string &name, ContainerType cType=BANK) const throw(evioException);
  template <typename T> evioDOMNodeP createNode(const string &name, const vector<T> &tVec) const throw(evioException);
  template <typename T> evioDOMNodeP createNode(const string &name, const T* t, int len) const throw(evioException);
  evioDOMNodeP createNode(const string &name, const evioSerializable &o, ContainerType cType=BANK) const throw(evioException);
  evioDOMNodeP createNode(const string &name, void (*f)(evioDOMNodeP c, void *userArg), void *userArg, ContainerType cType=BANK) const
    throw(evioException);
  template <typename T> evioDOMNodeP createNode(const string &name, T *t, void *userArg, ContainerType cType=BANK) const throw(evioException);
  template <typename T> evioDOMNodeP createNode(const string &name, T *t, 
                                                void* T::*mfp(evioDOMNodeP c, void *userArg),
                                                void *userArg, ContainerType cType=BANK) const throw(evioException);
  evioDOMNodeP createNode(const string &name, uint16_t formatTag, const string &formatString, 
                          uint16_t dataTag, uint8_t dataNum, const vector<uint32_t> &dataVec) const throw(evioException);
  evioDOMNodeP createNode(const string &name, uint16_t formatTag, const string &formatString, 
                          uint16_t dataTag, uint8_t dataNum, const uint32_t *t, int len) const throw(evioException);



public:
  evioDOMTree& operator<<(evioDOMNodeP node) throw(evioException);


public:
  int getSerializedLength(void) const throw(evioException);
  int toEVIOBuffer(uint32_t *buf, int size) const throw(evioException);


public:
  evioDOMNodeListP getNodeList(void) throw(evioException);
  evioDOMNodeListP getNodeList(const string &name) throw(evioException);
  template <class Predicate> evioDOMNodeListP getNodeList(Predicate pred) throw(evioException);
  template <class Predicate> evioDOMNodeP getFirstNode(Predicate pred) throw(evioException);
  template <typename T> vector<T> *getVectorUnique(void) throw(evioException);
  template <typename T, class Predicate> vector<T> *getVectorUnique(Predicate pred) throw(evioException);


public:
  string toString(void) const;
  string toString(const evioToStringConfig *config) const;
  string toString(const evioToStringConfig &config) const;

  const evioDictionary *getDictionary(void) const;
  void setDictionary(const evioDictionary *dict);
  void setDictionary(const evioDictionary &dict);


private:
  evioDOMNodeP parse(const uint32_t *buf) throw(evioException);
  int getSerializedLength(const evioDOMNodeP pNode) const throw(evioException);
  int toEVIOBuffer(uint32_t *buf, const evioDOMNodeP pNode, int size) const throw(evioException);
  void toOstream(ostream &os, const evioDOMNodeP node, int depth, const evioToStringConfig *config = &defaultToStringConfig) const 
    throw(evioException);
  template <class Predicate> evioDOMNodeList *addToNodeList(evioDOMNodeP pNode, evioDOMNodeList *pList, Predicate pred) throw(evioException);
  template <class Predicate> evioDOMNodeP findFirstNode(evioDOMNodeP pNode, Predicate pred) throw(evioException);


private:
  void *containerNodeHandler(int bankLength, int constainerType, int contentType, uint16_t tag, uint8_t num, 
                             int depth, const uint32_t *bankPointer, int payloadLength, const uint32_t *payload, void *userArg);
  void *leafNodeHandler(int bankLength, int containerType, int contentType, uint16_t tag, uint8_t num, 
                        int depth, const uint32_t *bankPointer, int dataLength, const void *data, void *userArg);
  

public:
  evioDOMNodeP root;                 /**<Pointer to root node of tree.*/
  string name;                       /**<Name of tree.*/
  const evioDictionary *dictionary;  /**<Dictionary to use for this tree.*/

};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Interface for object serialization.
 * Just defines serialize method.
 */
class evioSerializable {

public:
  virtual void serialize(evioDOMNodeP node) const throw(evioException) = 0;
  virtual ~evioSerializable(void) {};
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------



// include templates
#include "evioUtilTemplates.hxx"



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


}  // namespace evio


#endif
