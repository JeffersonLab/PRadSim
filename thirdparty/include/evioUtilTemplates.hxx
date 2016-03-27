// evioUtilTemplates.hxx

//  all the evio templates, some prototypes in evioUtil.hxx

//  ejw, 5-dec-2006

#include <iostream>
#include <iomanip>

#include "evio.h"
#include "evioTypedefs.hxx"
#include "evioException.hxx"
#include "evioDictionary.hxx"
#include "evioChannel.hxx"


#ifndef _evioUtilTemplates_hxx
#define _evioUtilTemplates_hxx


template <typename T> class typeIs;


//--------------------------------------------------------------
//----------------------- local class  -------------------------
//--------------------------------------------------------------


/**
 * Templated utility class has method that returns content type based on typename T.
 * Complete specializations supplied for all defined types.
 * Cumbersome, but this is the only way to do it on solaris...ejw, 9-jan-2007.
 * Note that the long data type is not supported due to ambiguities on different architectures.
 */
template <typename T> class evioUtil {

public:
  static int evioContentType(void) throw(evioException) {
    throw(evioException(0,"?evioUtil<T>::evioContentType...unsupported data type: "+string(typeid(T).name()),
                        __FILE__,__FUNCTION__,__LINE__));
    return(0x0);
  }
};

template <> class evioUtil<uint32_t>     {public: static int evioContentType(void)  throw(evioException) {return(0x1);}};
template <> class evioUtil<float>        {public: static int evioContentType(void)  throw(evioException) {return(0x2);}};
template <> class evioUtil<string>       {public: static int evioContentType(void)  throw(evioException) {return(0x3);}};
template <> class evioUtil<string&>      {public: static int evioContentType(void)  throw(evioException) {return(0x3);}};
template <> class evioUtil<int16_t>      {public: static int evioContentType(void)  throw(evioException) {return(0x4);}};
template <> class evioUtil<uint16_t>     {public: static int evioContentType(void)  throw(evioException) {return(0x5);}};
template <> class evioUtil<int8_t>       {public: static int evioContentType(void)  throw(evioException) {return(0x6);}};
template <> class evioUtil<uint8_t>      {public: static int evioContentType(void)  throw(evioException) {return(0x7);}};
template <> class evioUtil<double>       {public: static int evioContentType(void)  throw(evioException) {return(0x8);}};
template <> class evioUtil<int64_t>      {public: static int evioContentType(void)  throw(evioException) {return(0x9);}};
template <> class evioUtil<uint64_t>     {public: static int evioContentType(void)  throw(evioException) {return(0xa);}};
template <> class evioUtil<int32_t>      {public: static int evioContentType(void)  throw(evioException) {return(0xb);}};


//-----------------------------------------------------------------------------
//--------------------- evioDOMNode templated methods -------------------------
//-----------------------------------------------------------------------------


/** 
 * Static factory method to create empty leaf node of type T.
 * @param tag Node tag
 * @param num Node num
 * @return Pointer to new node
 */
template <typename T> evioDOMNodeP evioDOMNode::createEvioDOMNode(uint16_t tag, uint8_t num) throw(evioException) {
  return(new evioDOMLeafNode<T>(NULL,tag,num));
}


//-----------------------------------------------------------------------------


/** 
 * Static factory method to create empty leaf node of type T.
 * @param name Node name
 * @param dictionary Dictionary to use
 * @return Pointer to new node
 */
template <typename T> evioDOMNodeP evioDOMNode::createEvioDOMNode(const string &name, const evioDictionary *dictionary) throw(evioException) {

  if(dictionary!=NULL) {
    tagNum tn = dictionary->getTagNum(name);
    return(new evioDOMLeafNode<T>(NULL,tn.first,tn.second));
  } else {
    throw(evioException(0,"?evioDOMNode constructor...NULL dictionary for bank name: " + name,__FILE__,__FUNCTION__,__LINE__));
  }

}


//-----------------------------------------------------------------------------


/** 
 * Static factory method to create leaf node of type T.
 * @param tag Node tag
 * @param num Node num
 * @param tVec vector<T> of data
 * @return Pointer to new node
 */
template <typename T> evioDOMNodeP evioDOMNode::createEvioDOMNode(uint16_t tag, uint8_t num, const vector<T> &tVec)
  throw(evioException) {
  return(new evioDOMLeafNode<T>(NULL,tag,num,tVec));
}


//-----------------------------------------------------------------------------


/** 
 * Static factory method to create leaf node of type T.
 * @param name Node name
 * @param dictionary Dictionary to use
 * @param tVec vector<T> of data
 * @return Pointer to new node
 */
template <typename T> evioDOMNodeP evioDOMNode::createEvioDOMNode(const string &name, const evioDictionary *dictionary, const vector<T> &tVec)
  throw(evioException) {

  if(dictionary!=NULL) {
    tagNum tn = dictionary->getTagNum(name);
    return(new evioDOMLeafNode<T>(NULL,tn.first,tn.second,tVec));
  } else {
    throw(evioException(0,"?evioDOMNode constructor...NULL dictionary for bank name: " + name,__FILE__,__FUNCTION__,__LINE__));
  }

}


//-----------------------------------------------------------------------------


/** 
 * Static factory method to create leaf node of type T.
 * @param tag Node tag
 * @param num Node num
 * @param t Pointer to array containg data of type T
 * @param len Length of array
 * @return Pointer to new node
 */
template <typename T> evioDOMNodeP evioDOMNode::createEvioDOMNode(uint16_t tag, uint8_t num, const T* t, int len)
  throw(evioException) {
  return(new evioDOMLeafNode<T>(NULL,tag,num,t,len));
}


//-----------------------------------------------------------------------------


/** 
 * Static factory method to create leaf node of type T.
 * @param name Node name
 * @param dictionary Dictionary to use
 * @param t Pointer to array containg data of type T
 * @param len Length of array
 * @return Pointer to new node
 */
template <typename T> evioDOMNodeP evioDOMNode::createEvioDOMNode(const string &name, const evioDictionary *dictionary, const T* t, int len)
  throw(evioException) {

  if(dictionary!=NULL) {
    tagNum tn = dictionary->getTagNum(name);
    return(new evioDOMLeafNode<T>(NULL,tn.first,tn.second,t,len));
  } else {
    throw(evioException(0,"?evioDOMNode constructor...NULL dictionary for bank name: " + name,__FILE__,__FUNCTION__,__LINE__));
  }

}


//-----------------------------------------------------------------------------


/** 
 * Static factory method to create container node using serialize method to fill container.
 * @param name Node name
 * @param dictionary Dictionary to use
 * @param t Pointer to object having serialize method
 * @param userArg User arg passed to serialize method
 * @param cType Container node content type
 * @return Pointer to new node
 */
template <typename T> evioDOMNodeP evioDOMNode::createEvioDOMNode(const string &name, const evioDictionary *dictionary, T *t, 
                                                                  void *userArg, ContainerType cType) throw(evioException) {

  if(dictionary!=NULL) {
    tagNum tn = dictionary->getTagNum(name);
    evioDOMContainerNode *c = new evioDOMContainerNode(NULL,tn.first,tn.second,cType);
    t->serialize(c,userArg);
    return(c);
  } else {
    throw(evioException(0,"?evioDOMNode constructor...NULL dictionary for bank name: " + name,__FILE__,__FUNCTION__,__LINE__));
  }

}


//-----------------------------------------------------------------------------


/** 
 * Static factory method to create container node filled via object method pointer.
 * @param tag Node tag
 * @param num Node num
 * @param t Pointer to object
 * @param mfp Pointer to object method used to fill container
 * @param userArg User arg passed to object method
 * @param cType Container node content type
 * @return Pointer to new node
 */
template <typename T> evioDOMNodeP evioDOMNode::createEvioDOMNode(uint16_t tag, uint8_t num, T *t, 
                                                                  void* T::*mfp(evioDOMNodeP c, void *userArg),
                                                                  void *userArg, ContainerType cType) throw(evioException) {
  evioDOMContainerNode *c = new evioDOMContainerNode(NULL,tag,num,cType);
  t->mfp(c,userArg);
  return(c);
}


//-----------------------------------------------------------------------------


/** 
 * Static factory method to create container node filled via object method pointer.
 * @param name Node name
 * @param dictionary Dictionary to use
 * @param t Pointer to object
 * @param mfp Pointer to object method used to fill container
 * @param userArg User arg passed to object method
 * @param cType Container node content type
 * @return Pointer to new node
 */
template <typename T> evioDOMNodeP evioDOMNode::createEvioDOMNode(const string &name, const evioDictionary *dictionary, T *t, 
                                                                  void* T::*mfp(evioDOMNodeP c, void *userArg),
                                                                  void *userArg, ContainerType cType) throw(evioException) {

  if(dictionary!=NULL) {
    tagNum tn = dictionary->getTagNum(name);
    evioDOMContainerNode *c = new evioDOMContainerNode(NULL,tn.first,tn.second,cType);
    t->mfp(c,userArg);
    return(c);
  } else {
    throw(evioException(0,"?evioDOMNode constructor...NULL dictionary for bank name: " + name,__FILE__,__FUNCTION__,__LINE__));
  }

}


//-----------------------------------------------------------------------------


/** 
 * Appends one data element to leaf node.
 * Must be done this way because C++ forbids templated virtual functions...ejw, dec-2006
 * @param tVal Data to be added.
 */

template <typename T> void evioDOMNode::append(T tVal) throw(evioException) {

  static int errCount = 0;


  // first try for exact match between leaf node type and data value
  if(contentType==evioUtil<T>::evioContentType()) {
    evioDOMLeafNode<T> *l = static_cast< evioDOMLeafNode<T>* >(this);
    l->data.push_back(tVal);
    return;
  }
  
  // no match, try all possibilities
  if((++errCount%2)==100) cout << "evioDOMNode::append...type mismatch, slowly trying every possibility, count = " 
                               << errCount << endl;
  

  {
    if(contentType==evioUtil<int32_t>::evioContentType()) {
      evioDOMLeafNode<int32_t> *l = static_cast<evioDOMLeafNode<int32_t>*>(this);
      l->data.push_back((int32_t)tVal);
      return;
    }
  }

  {
    if(contentType==evioUtil<uint32_t>::evioContentType()) {
      evioDOMLeafNode<uint32_t> *l = static_cast<evioDOMLeafNode<uint32_t>*>(this);
      l->data.push_back((uint32_t)tVal);
      return;
    }
  }

  {
    if(contentType==evioUtil<double>::evioContentType()) {
      evioDOMLeafNode<double> *l = static_cast<evioDOMLeafNode<double>*>(this);
      l->data.push_back((double)tVal);
      return;
    }
  }

  {
    if(contentType==evioUtil<float>::evioContentType()) {
      evioDOMLeafNode<float> *l = static_cast<evioDOMLeafNode<float>*>(this);
      l->data.push_back((float)tVal);
      return;
    }
  }

  {
    if(contentType==evioUtil<int16_t>::evioContentType()) {
      evioDOMLeafNode<int16_t> *l = static_cast<evioDOMLeafNode<int16_t>*>(this);
      l->data.push_back((int16_t)tVal);
      return;
    }
  }

  {
    if(contentType==evioUtil<uint16_t>::evioContentType()) {
      evioDOMLeafNode<uint16_t> *l = static_cast<evioDOMLeafNode<uint16_t>*>(this);
      l->data.push_back((uint16_t)tVal);
      return;
    }
  }

  {
    if(contentType==evioUtil<int64_t>::evioContentType()) {
      evioDOMLeafNode<int64_t> *l = static_cast<evioDOMLeafNode<int64_t>*>(this);
      l->data.push_back((int64_t)tVal);
      return;
    }
  }

  {
    if(contentType==evioUtil<uint64_t>::evioContentType()) {
      evioDOMLeafNode<uint64_t> *l = static_cast<evioDOMLeafNode<uint64_t>*>(this);
      l->data.push_back((uint64_t)tVal);
      return;
    }
  }

  {
    if(contentType==evioUtil<uint8_t>::evioContentType()) {
      evioDOMLeafNode<int8_t> *l = static_cast<evioDOMLeafNode<int8_t>*>(this);
      l->data.push_back((int8_t)tVal);
      return;
    }
  }

  {
    if(contentType==evioUtil<uint8_t>::evioContentType()) {
      evioDOMLeafNode<uint8_t> *l = static_cast<evioDOMLeafNode<uint8_t>*>(this);
      l->data.push_back((uint8_t)tVal);
      return;
    }
  }


  // no match, must not be a leaf node
  throw(evioException(0,"?evioDOMNode::append...unable to append",__FILE__,__FUNCTION__,__LINE__));
}


//-----------------------------------------------------------------------------


/** 
 * Appends vector of data to leaf node.
 * Must be done this way because C++ forbids templated virtual functions...ejw, dec-2006
 * @param tVec vector<T> of data to add to leaf node
 */
template <typename T> void evioDOMNode::append(const vector<T> &tVec) throw(evioException) {
  if(contentType!=evioUtil<T>::evioContentType())
    throw(evioException(0,"?evioDOMNode::append...not appropriate node",__FILE__,__FUNCTION__,__LINE__));
  evioDOMLeafNode<T> *l = static_cast<evioDOMLeafNode<T>*>(this);
  l->data.insert(l->data.end(), tVec.begin(), tVec.end());
}


//-----------------------------------------------------------------------------


/** 
 * Appends array of data to leaf node.
 * Must be done this way because C++ forbids templated virtual functions...ejw, dec-2006
 * @param tBuf Buffer of data of type T
 * @param len Length of buffer
 */
template <typename T> void evioDOMNode::append(const T* tBuf, int len) throw(evioException) {
  if(contentType!=evioUtil<T>::evioContentType())
    throw(evioException(0,"?evioDOMNode::append...not appropriate node",__FILE__,__FUNCTION__,__LINE__));
  evioDOMLeafNode<T> *l = static_cast<evioDOMLeafNode<T>*>(this);
  l->data.insert(l->data.end(), tBuf, tBuf + len);
}


//-----------------------------------------------------------------------------


/** 
 * Replaces data in leaf node.
 * Must be done this way because C++ forbids templated virtual functions...ejw, dec-2006
 * @param tVec vector<T> of data
 */
template <typename T> void evioDOMNode::replace(const vector<T> &tVec) throw(evioException) {
  if(contentType!=evioUtil<T>::evioContentType())
    throw(evioException(0,"?evioDOMNode::replace...not correctl leaf node",__FILE__,__FUNCTION__,__LINE__));
  evioDOMLeafNode<T> *l = static_cast<evioDOMLeafNode<T>*>(this);
  l->data = tVec;
}


//-----------------------------------------------------------------------------


/** 
 * Replaces data in leaf node.
 * Must be done this way because C++ forbids templated virtual functions...ejw, dec-2006
 * @param tBuf Buffer of data of type T
 * @param len Length of buffer
 */
template <typename T> void evioDOMNode::replace(const T* tBuf, int len) throw(evioException) {
  if(contentType!=evioUtil<T>::evioContentType())
    throw(evioException(0,"?evioDOMNode::replace...not appropriate node",__FILE__,__FUNCTION__,__LINE__));
  evioDOMLeafNode<T> *l = static_cast<evioDOMLeafNode<T>*>(this);
  l->data.assign(tBuf, tBuf + len);
}


//-----------------------------------------------------------------------------


/** 
 * Returns pointer to leaf node vector<T> of data
 * Must be done this way because C++ forbids templated virtual functions...ejw, dec-2006
 * @return Pointer to vector<T>
 */
template <typename T> vector<T> *evioDOMNode::getVector(void) throw(evioException) {
  if(contentType!=evioUtil<T>::evioContentType())return(NULL);
  evioDOMLeafNode<T> *l = static_cast<evioDOMLeafNode<T>*>(this);
  return(&l->data);
}


//-----------------------------------------------------------------------------


/** 
 * Appends single data value to leaf node
 * @param tVal Data to be added
 * @return Reference to this
 */
template <typename T> evioDOMNode& evioDOMNode::operator<<(T tVal) throw(evioException) {
  append(tVal);
  return(*this);
}


//-----------------------------------------------------------------------------


/** 
 * Appends vector of data to leaf node.
 * @param tVec vector<T> of data to add to leaf node
 * @return Reference to this
 */
template <typename T> evioDOMNode& evioDOMNode::operator<<(const vector<T> &tVec) throw(evioException) {
  append(tVec);
  return(*this);
}


//-----------------------------------------------------------------------------


/** 
 * Returns list of children satisfying predicate
 * @param pred Predicate
 * @return List of children satisfying predicate
 */
template <class Predicate> evioDOMNodeListP evioDOMNode::getChildren(Predicate pred) throw(evioException) {
  evioDOMNodeListP l = getChildren();
  if(l.get()==NULL)return(l);
  l->remove_if(not1(pred));
  return(l);
}


//-----------------------------------------------------------------------------
//--------------------- evioDOMLeafNode templated methods ---------------------
//-----------------------------------------------------------------------------


/** 
 * Leaf node constructor used internally.
 * @param par Parent node
 * @param tag Node tag
 * @param num Node num
 */
template <typename T> evioDOMLeafNode<T>::evioDOMLeafNode(evioDOMNodeP par, uint16_t tag, uint8_t num)
  throw(evioException) : evioDOMNode(par,tag,num,evioUtil<T>::evioContentType()) {
}


//-----------------------------------------------------------------------------


/** 
 * Leaf node constructor used internally.
 * @param par Parent node
 * @param tag Node tag
 * @param num Node num
 * @param v vector<T> of data
 */
template <typename T> evioDOMLeafNode<T>::evioDOMLeafNode(evioDOMNodeP par, uint16_t tag, uint8_t num, const vector<T> &v)
  throw(evioException) : evioDOMNode(par,tag,num,evioUtil<T>::evioContentType()), data(v) {
}


//-----------------------------------------------------------------------------


/** 
 * Leaf node constructor used internally.
 * @param par Parent node
 * @param tag Node tag
 * @param num Node num
 * @param p Pointer to array containg data of type T
 * @param ndata Length of array
 */
template <typename T> evioDOMLeafNode<T>::evioDOMLeafNode(evioDOMNodeP par, uint16_t tag, uint8_t num, const T* p, int ndata) 
  throw(evioException) : evioDOMNode(par,tag,num,evioUtil<T>::evioContentType()), data(p, p + ndata) {
}


//-----------------------------------------------------------------------------


/*
 * Destructor
 */
template <typename T> evioDOMLeafNode<T>::~evioDOMLeafNode(void) {
}


//-----------------------------------------------------------------------------


/**
 * Returns XML string containing header needed by toString
 * @param depth Current depth
 * @return XML string
 */
template <typename T> string evioDOMLeafNode<T>::getHeader(int depth, const evioToStringConfig *config) const {

  ostringstream os;
  string indent = ((config==NULL)?getIndent(depth):getIndent(depth,config->indentSize));


  // get node name
  string name;
  if((config!=NULL)&&(config->toStringDictionary!=NULL)) {
    map<tagNum,string>::const_iterator iter = config->toStringDictionary->getNameMap.find(tagNum(tag,num));
    if(iter!=config->toStringDictionary->getNameMap.end()) name=(*iter).second;
  }
  if(name.size()<=0) name = evGetTypename(parent==NULL?BANK:parent->getContentType());


  os << indent
     <<  "<" << name << " content=\"" << evGetTypename(contentType) << "\""
     << " data_type=\"" << hex << showbase << contentType << noshowbase << dec
     << "\" tag=\"" << tag;

  if((parent==NULL)||((parent->getContentType()==0xe)||(parent->getContentType()==0x10)))
    os << dec << "\" num=\"" << (int)num;

  if((config!=NULL)&&(config->verbose)) {
    os << dec << "\" nwords=\"" << getSize();
  }


  os << "\">" << endl;


  return(os.str());
}


//-----------------------------------------------------------------------------


/**
 * Returns XML string containing body needed by toString
 * @param depth Current depth
 * @return XML string
 */
template <typename T> string evioDOMLeafNode<T>::getBody(int depth, const evioToStringConfig *config) const {

  ostringstream os;
  string indent = ((config==NULL)?getIndent(depth):getIndent(depth,config->indentSize));
  string indent2 = indent + "       ";
  string spaces = "     ";

  int wid,swid;
  switch (contentType) {

  case 0x0:
  case 0x1:
  case 0x2:
  case 0xb:
    wid=5;
    swid=10;
    break;
  case 0x4:
  case 0x5:
    wid=8;
    swid=6;
    break;
  case 0x6:
  case 0x7:
    wid=8;
    swid=4;
    break;
  case 0x8:
  case 0x9:
  case 0xa:
    wid=2;
    swid=28;
    break;
  default:
    wid=1;
    swid=30;
    break;
  }



  // dump data...odd what has to be done for 1-byte types 0x6,0x7 due to bugs in ostream operator <<
  int16_t k;
  typename vector<T>::const_iterator iter;
  for(iter=data.begin(); iter!=data.end();) {

    os << indent2;
    for(int j=0; (j<wid)&&(iter!=data.end()); j++) {
      switch (contentType) {

      case 0x0:
      case 0x1:
      case 0x5:
      case 0xa:
        os.width(swid);
        if((config!=NULL)&&(config->xtod)) {
          os << dec << noshowbase << *iter << spaces;
        } else {
          os << hex << showbase << *iter << noshowbase << dec << spaces;
        }
        break;
      case 0x2:
        os.precision(6);
        os.width(swid);
        os << showpoint << *iter << noshowpoint << spaces;
        break;
      case 0x3:
        os << "<![CDATA[" << *iter << "]]>";
        break;
      case 0x6:
        k = (*((int16_t*)(&(*iter)))) & 0xff;
        if((k&0x80)!=0)k|=0xff00;
        os.width(swid);
        os << k << spaces;
        break;
      case 0x7:
        os.width(swid);
        if((config!=NULL)&&(config->xtod)) {
          os << dec << noshowbase << ((*(int*)&(*iter))&0xff) << spaces;
        } else {
          os << hex << showbase << ((*(int*)&(*iter))&0xff) << noshowbase << dec << spaces;
        }
        break;
      case 0x8:
        os.width(swid);
        os.precision(20);
        os << scientific << *iter << fixed << spaces;
        break;
      default:
        os.width(swid);
        os << *iter << spaces;
        break;
      }
      iter++;
    }
    os << dec << endl;

  }


  return(os.str());
}


//-----------------------------------------------------------------------------


/**
 * Returns XML string containing footer needed by toString
 * @param depth Current depth
 * @return XML string
 */
template <typename T> string evioDOMLeafNode<T>::getFooter(int depth, const evioToStringConfig *config) const {
  ostringstream os;

  // get node name
  string name;
  if((config!=NULL)&&(config->toStringDictionary!=NULL)) {
    map<tagNum,string>::const_iterator iter = config->toStringDictionary->getNameMap.find(tagNum(tag,num));
    if(iter!=config->toStringDictionary->getNameMap.end()) name=(*iter).second;
  }
  if(name.size()<=0) name = evGetTypename(parent==NULL?BANK:parent->getContentType());

  os << ((config==NULL)?getIndent(depth):getIndent(depth,config->indentSize)) << "</" << name << ">" << endl;
  return(os.str());
}



//-----------------------------------------------------------------------------


/**
 * Returns numnber of data elements
 * @return number of data elements
 */
template <typename T> int evioDOMLeafNode<T>::getSize(void) const {
  return(data.size());
}



//-----------------------------------------------------------------------------
//----------------------- evioDOMTree templated methods -----------------------
//-----------------------------------------------------------------------------


/**
 * Returns list of nodes in tree satisfying predicate.
 * @param pred Function object true if node meets predicate criteria
 * @return Pointer to node list (actually auto_ptr<>)
 */
template <class Predicate> evioDOMNodeListP evioDOMTree::getNodeList(Predicate pred) throw(evioException) {
  evioDOMNodeList *pList = addToNodeList(root, new evioDOMNodeList(), pred);
  return(evioDOMNodeListP(pList));
}  


//-----------------------------------------------------------------------------


/**
 * Returns pointer to first node in tree satisfying predicate.
 * Consider that the order of the search is undefined, so if there is
 *   more than one bank satisfying predicate this method may not return
 *   the same thing each time.  Currently it does a depth-first search, 
 *   but this could change.
 * I.e. to be safe only use this method when you are sure there is only 
 *   one bank that satisfies predicate.
 * @param pred Function object true if node meets predicate criteria
 * @return Pointer to node
 */
template <class Predicate> evioDOMNodeP evioDOMTree::getFirstNode(Predicate pred) throw(evioException) {
  return(findFirstNode(root,pred));
}  


//-----------------------------------------------------------------------------


/**
 * Returns pointer to first node in tree satisfying predicate.
 * @param pred Function object true if node meets predicate criteria
 * @return Pointer to node
 */
template <class Predicate> evioDOMNodeP evioDOMTree::findFirstNode(evioDOMNodeP pNode, Predicate pred) throw(evioException) {


  // check this node
  if(pred(pNode))return(pNode);
  if(pNode->isLeaf())return(NULL);
  
 
  // check children
  if(pNode->isContainer()) {
    evioDOMContainerNode *c = static_cast<evioDOMContainerNode*>(pNode);
    evioDOMNodeList::iterator iter;
    for(iter=c->childList.begin(); iter!=c->childList.end(); iter++) {
      evioDOMNodeP p=findFirstNode(*iter,pred);
      if(p!=NULL)return(p);
    }
  }   

  return(NULL);
}  


//-----------------------------------------------------------------------------


/**
 * Returns vector<T> from single node in tree containing vector<T>.
 * Throws exception if more than one node contains vector<T>.
 * @return Pointer to vector, NULL if no node contains vector<T>
 */
template <typename T> vector<T> *evioDOMTree::getVectorUnique(void) throw(evioException) {

  evioDOMNodeListP l = getNodeList(typeIs<T>());
  int s = l->size();

  if(s==0) {
    return(NULL);
  } else if (s==1) {
    return((l->front())->getVector<T>());
  } else {
    throw(evioException(0,"?evioDOMTree::getVectorUnique...more than one node found",__FILE__,__FUNCTION__,__LINE__));
  }
}  


//-----------------------------------------------------------------------------



/**
 * Returns vector<T> from single node in tree that is of type T AND satisfies predicate.
 * Throws exception if more than one node satisfies predicate.
 * @param pred Function object true if node satisfies predicate
 * @return Pointer to vector<T>, NULL if no node containing vector<T> satisfies predicate
 */
template <typename T, class Predicate> vector<T> *evioDOMTree::getVectorUnique(Predicate pred) throw(evioException) {
  evioDOMNodeListP l = getNodeList(pred);

  int c = count_if(l->begin(),l->end(),typeIs<T>());
  
  if(c==0) {
    return(NULL);
  } else if (c==1) {
    evioDOMNodeList::iterator iter = find_if(l->begin(),l->end(),typeIs<T>());
    return((*iter)->getVector<T>());
  } else {
    throw(evioException(0,"?evioDOMTree::getVectorUnique...more than one node found",__FILE__,__FUNCTION__,__LINE__));
  }

}  


//-----------------------------------------------------------------------------


/**
 * Adds node to node list, used internally by getNodeList.
 * @param pNode Node to check agains predicate
 * @param pList Current node list
 * @param pred true if node meets predicate criteria
 * @return Pointer to node list
 */
template <class Predicate> evioDOMNodeList *evioDOMTree::addToNodeList(evioDOMNodeP pNode, evioDOMNodeList *pList, Predicate pred)
  throw(evioException) {

  if(pNode==NULL)return(pList);


  // add this node to list
  if(pred(pNode))pList->push_back(pNode);
  
  
  // add children to list
  if(pNode->isContainer()) {
    evioDOMContainerNode *c = static_cast<evioDOMContainerNode*>(pNode);
    evioDOMNodeList::iterator iter;
    for(iter=c->childList.begin(); iter!=c->childList.end(); iter++) {
      addToNodeList(*iter,pList,pred);
    }
  }


  // return the list
  return(pList);
}


//-----------------------------------------------------------------------------


/**
 * Creates leaf node and adds it to tree root node.
 * @param tag Node tag
 * @param num Node num
 * @param dataVec vector<T> of data
 */
template <typename T> void evioDOMTree::addBank(uint16_t tag, uint8_t num, const vector<T> &dataVec) throw(evioException) {

  if(root==NULL) {
    root = evioDOMNode::createEvioDOMNode(tag,num,dataVec);
    root->parentTree=this;

  } else {
    if(!root->isContainer())throw(evioException(0,"?evioDOMTree::addBank...root not a container node",__FILE__,__FUNCTION__,__LINE__));
    evioDOMContainerNode* c = static_cast<evioDOMContainerNode*>(root);
    evioDOMNodeP node = evioDOMNode::createEvioDOMNode(tag,num,dataVec);
    c->childList.push_back(node);
    node->parent=root;
  }

  return;
}


//-----------------------------------------------------------------------------


/** 
 * Creates leaf node and adds it to tree root node.
 * @param tag Node tag
 * @param num Node num
 * @param dataBuf Pointer to array containg data of type T
 * @param dataLen Length of array
 */
template <typename T> void evioDOMTree::addBank(uint16_t tag, uint8_t num, const T* dataBuf, int dataLen)
  throw(evioException) {

  if(root==NULL) {
    root = evioDOMNode::createEvioDOMNode(tag,num,dataBuf,dataLen);
    root->parentTree=this;

  } else {
    if(!root->isContainer())throw(evioException(0,"?evioDOMTree::addBank...root not a container node",__FILE__,__FUNCTION__,__LINE__));
    evioDOMContainerNode* c = static_cast<evioDOMContainerNode*>(root);
    evioDOMNodeP node = evioDOMNode::createEvioDOMNode(tag,num,dataBuf,dataLen);
    c->childList.push_back(node);
    node->parent=root;
  }

  return;
}


//-----------------------------------------------------------------------------


/**
 * Creates leaf node and adds it to tree root node.
 * @param tn Leaf node tagNum
 * @param dataVec vector<T> of data
 */
template <typename T> void evioDOMTree::addBank(tagNum tn, const vector<T> &dataVec) throw(evioException) {
  addBank(tn.first,tn.second,dataVec);
  return;
}


//-----------------------------------------------------------------------------


/** 
 * Creates leaf node and adds it to tree root node.  @param tn Leaf
 * node tagNum @param dataBuf Pointer to array containg data of type T
 * @param dataLen Length of array
 */
template <typename T> void evioDOMTree::addBank(tagNum tn, const T* dataBuf, int dataLen) throw(evioException) {
  addBank(tn.first,tn.second,dataBuf,dataLen);
  return;
}


//-----------------------------------------------------------------------------


/**
 * Creates leaf node and adds it to tree root node.
 * @param name Leaf node name
 * @param dataVec vector<T> of data
 */
template <typename T> void evioDOMTree::addBank(const string &name, const vector<T> &dataVec) throw(evioException) {

  if(dictionary!=NULL) {
    tagNum tn = dictionary->getTagNum(name);
    addBank(tn,dataVec);
  } else {
    throw(evioException(0,"?evioDOMTree::addBank...no dictionary",__FILE__,__FUNCTION__,__LINE__));
  }

  return;
}


//-----------------------------------------------------------------------------


/** 
 * Creates leaf node and adds it to tree root node.
 * @param name Node name
 * @param dataBuf Pointer to array containg data of type T
 * @param dataLen Length of array
 */
template <typename T> void evioDOMTree::addBank(const string &name, const T* dataBuf, int dataLen) throw(evioException) {

  if(dictionary!=NULL) {
    tagNum tn = dictionary->getTagNum(name);
    addBank(tn.first,tn.second,dataBuf,dataLen);
  } else {
    throw(evioException(0,"?evioDOMTree::addBank...no dictionary",__FILE__,__FUNCTION__,__LINE__));
  }

  return;
}


//-----------------------------------------------------------------------------


/** 
 * Creates leaf node.
 * @param name Node name
 * @param tVec Vector of values
 * @return Pointer to new node
 */
template <typename T> evioDOMNodeP evioDOMTree::createNode(const string &name, const vector<T> &tVec) const throw(evioException) {
  return(evioDOMNode::createEvioDOMNode(name,dictionary,tVec));
}


//-----------------------------------------------------------------------------


/** 
 * Creates leaf node.
 * @param name Node name
 * @param t Pointer to array of values
 * @param len Length of array
 * @return Pointer to new node
 */
template <typename T> evioDOMNodeP evioDOMTree::createNode(const string &name, const T* t, int len) const throw(evioException) {
  return(evioDOMNode::createEvioDOMNode(name,dictionary,t,len));
}


//-----------------------------------------------------------------------------


/** 
 * Creates container node.
 * @param name Node name
 * @param t Pointer to object having serialize method
 * @param userArg User arg passed to serialize method
 * @parme cType Container node type
 * @return Pointer to new node
 */
template <typename T> evioDOMNodeP evioDOMTree::createNode(const string &name, T *t, void *userArg, ContainerType cType) const 
  throw(evioException) {
  return(evioDOMNode::createEvioDOMNode(name,dictionary,t,userArg,cType));
}


//-----------------------------------------------------------------------------


/** 
 * Creates container node.
 * @param name Node name
 * @param t Pointer to object having serialize method
 * @param userArg User arg passed to serialize method
 * @parme cType Container node type
 * @return Pointer to new node
 */
template <typename T> evioDOMNodeP evioDOMTree::createNode(const string &name, T *t, 
                                              void* T::*mfp(evioDOMNodeP c, void *userArg),
                                              void *userArg, ContainerType cType) const throw(evioException) {
  return(evioDOMNode::createEvioDOMNode(name,dictionary,t,mfp,userArg,cType));
}



//-----------------------------------------------------------------------------
//--------------------- Misc Function Objects ---------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Boolean function object compares to content type for typename T.
 * Note that this function object cannot deal with unknown or composite node types since they do not have a distinguishing type.
 */
template <typename T> class typeIs : public unary_function<const evioDOMNodeP,bool> {

public:
  typeIs(void) : type(evioUtil<T>::evioContentType()) {}

  bool operator()(const evioDOMNodeP node) const {
    int nodeType = node->getContentType();
    if((nodeType==0x0)||(nodeType==0xf))nodeType=0x1;   // coerce for unknown and composite types
    return(nodeType==type);
  }

private:
 int type;
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Boolean function object compares on type.
 */
class typeEquals : public unary_function<const evioDOMNodeP,bool> {

public:
  typeEquals(int aType) : type(aType) {}
  bool operator()(const evioDOMNodeP node) const {return(node->getContentType()==type);}
private:
  int type;
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Boolean function object compares on tag.
 */
class tagEquals : public unary_function<const evioDOMNodeP,bool> {

public:
  tagEquals(uint16_t aTag) : tag(aTag) {}

  tagEquals(const string &name, evioDictionary *dictionary) : tag(0) {
    if(dictionary!=NULL) {
      tagNum tn = dictionary->getTagNum(name); 
      tag=tn.first;
    } else {
      cerr << "?tagEquals...NULL dictionary, using tag=0" << endl;
    }
  }

  tagEquals(const string &name, evioDictionary &dictionary) : tag(0) {
    tagNum tn = dictionary.getTagNum(name); 
    tag=tn.first;
  }

  bool operator()(const evioDOMNodeP node) const {return(node->tag==tag);}

private:
  uint16_t tag;
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Boolean function object compares on num.
 */
class numEquals : public unary_function<const evioDOMNodeP,bool> {

public:
  numEquals(uint8_t aNum) : num(aNum) {}

  numEquals(const string &name, evioDictionary *dictionary) : num(0) {
    if(dictionary!=NULL) {
      tagNum tn = dictionary->getTagNum(name); 
      num=tn.second;
    } else {
      cerr << "?numEquals...NULL dictionary, using num=0" << endl;
    }
  }

  numEquals(const string &name, evioDictionary &dictionary) : num(0) {
    tagNum tn = dictionary.getTagNum(name); 
    num=tn.second;
  }

  bool operator()(const evioDOMNodeP node) const {return(node->num==num);}

private:
  uint8_t num;
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Boolean function object compares on tag and num.
 */
class tagNumEquals : public unary_function<const evioDOMNodeP, bool> {

public:
  tagNumEquals(uint16_t aTag, uint8_t aNum) : tag(aTag), num(aNum) {}

  tagNumEquals(tagNum tn) : tag(tn.first), num(tn.second) {}

  tagNumEquals(const string &name, evioDictionary *dictionary) : tag(0), num(0) {
    if(dictionary!=NULL) {
      tagNum tn = dictionary->getTagNum(name); 
      tag=tn.first;
      num=tn.second;
    } else {
      cerr << "?tagNumEquals...NULL dictionary, using tag=0,num=0" << endl;
    }
  }

  tagNumEquals(const string &name, evioDictionary &dictionary) : tag(0), num(0) {
    tagNum tn = dictionary.getTagNum(name); 
    tag=tn.first;
    num=tn.second;
  }

  tagNumEquals(const string &name, const evioDOMTree *tree) : tag(0), num(0) {
    if((tree!=NULL)&&(tree->dictionary!=NULL)) {
      tagNum tn = tree->dictionary->getTagNum(name); 
      tag=tn.first;
      num=tn.second;
    } else {
      cerr << "?tagNumEquals...NULL tree or dictionary, using tag=0,num=0" << endl;
    }
  }

  tagNumEquals(const string &name, const evioDOMTree &tree) : tag(0), num(0) {
    if(tree.dictionary!=NULL) {
      tagNum tn = tree.dictionary->getTagNum(name); 
      tag=tn.first;
      num=tn.second;
    } else {
      cerr << "?tagNumEquals...NULL tree, using tag=0,num=0" << endl;
    }
  }


  bool operator()(const evioDOMNodeP node) const {return((node->tag==tag)&&(node->num==num));}

private:
  uint16_t tag;
  uint8_t num;
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Boolean function object compares on parent content type.
 */
class parentTypeEquals : public unary_function<const evioDOMNodeP, bool> {

public:
  parentTypeEquals(int aType) : type(aType) {}
  bool operator()(const evioDOMNodeP node) const {return((node->getParent()==NULL)?false:(node->getParent()->getContentType()==type));}
private:
  int type;
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Boolean function object compares on parent tag.
 */
class parentTagEquals : public unary_function<const evioDOMNodeP, bool> {

public:
  parentTagEquals(uint16_t aTag) : tag(aTag) {}
  bool operator()(const evioDOMNodeP node) const {return((node->getParent()==NULL)?false:(node->getParent()->tag==tag));}
private:
  uint16_t tag;
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Boolean function object compares on parent num.
 */
class parentNumEquals : public unary_function<const evioDOMNodeP, bool> {

public:
  parentNumEquals(uint8_t aNum) : num(aNum) {}
  bool operator()(const evioDOMNodeP node) const {return((node->getParent()==NULL)?false:(node->getParent()->num==num));}
private:
  uint8_t num;
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Boolean function object compares on parent tag and num.
 */
class parentTagNumEquals : public unary_function<const evioDOMNodeP, bool> {

public:
  parentTagNumEquals(uint16_t aTag, uint8_t aNum) : tag(aTag), num(aNum) {}
  parentTagNumEquals(tagNum tn) : tag(tn.first), num(tn.second) {}
  bool operator()(const evioDOMNodeP node) const {
    return((node->getParent()==NULL)?false:((node->getParent()->tag==tag)&&(node->getParent()->num==num)));}
private:
  uint16_t tag;
  uint8_t num;
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Boolean function object compares on parent tag and num.
 */
class parentNameEquals : public unary_function<const evioDOMNodeP, bool> {

public:
  parentNameEquals(const string &name, const evioDictionary *dictionary) : tag(0), num(0) {
    if(dictionary!=NULL) {
      tagNum tn = dictionary->getTagNum(name); 
      tag=tn.first;
      num=tn.second;
    } else {
      cerr << "?parentNameEquals...NULL dictionary, using tag=0,num=0" << endl;
    }
  }

  parentNameEquals(const string &name, const evioDictionary &dictionary) : tag(0), num(0) {
    tagNum tn = dictionary.getTagNum(name); 
    tag=tn.first;
    num=tn.second;
  }

  parentNameEquals(const string &name, const evioDOMTree *tree) : tag(0), num(0) {
    if((tree!=NULL)&&(tree->dictionary!=NULL)) {
      tagNum tn = tree->dictionary->getTagNum(name); 
      tag=tn.first;
      num=tn.second;
    } else {
      cerr << "?parentNameEquals...NULL tree or dictionary, using tag=0,num=0" << endl;
    }
  }

  parentNameEquals(const string &name, const evioDOMTree &tree) : tag(0), num(0) {
    if(tree.dictionary!=NULL) {
      tagNum tn = tree.dictionary->getTagNum(name); 
      tag=tn.first;
      num=tn.second;
    } else {
      cerr << "?parentNameEquals...NULL dictionary, using tag=0,num=0" << endl;
    }
  }

  bool operator()(const evioDOMNodeP node) const {
    return((node->getParent()==NULL)?false:((node->getParent()->tag==tag)&&(node->getParent()->num==num)));}

private:
  uint16_t tag;
  uint8_t num;
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Boolean function object true if container node.
 */
class isContainer : public unary_function<const evioDOMNodeP,bool> {

public:
  isContainer(void) {}
  bool operator()(const evioDOMNodeP node) const {return(node->isContainer());}
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Boolean function object true if leaf node.
 */
class isLeaf : public unary_function<const evioDOMNodeP,bool> {

public:
  isLeaf(void) {}
  bool operator()(const evioDOMNodeP node) const {return(node->isLeaf());}
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Function object streams node->toString() to cout.
 */
class toCout: public unary_function<const evioDOMNodeP,void> {

public:
  toCout(void) {}
  void operator()(const evioDOMNodeP node) const {cout << node->toString() << endl;}
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


#endif
