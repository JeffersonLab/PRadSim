// evioChannel.hxx

//  Author:  Elliott Wolin, JLab, 18-feb-2010


#ifndef _evioChannel_hxx
#define _evioChannel_hxx


#include <stdint.h>
#include "evioException.hxx"



using namespace std;


namespace evio {


class evioDOMTree;
class evioDictionary;



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Pure virtual class defines single method that implements serialization of 
 *  object to evio buffer.  
 **/
class evioChannelBufferizable {

public:
  virtual int toEVIOBuffer(uint32_t *buf, int size) const throw(evioException) = 0;
  virtual ~evioChannelBufferizable() {}
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Defines EVIO I/O channel functionality.
 * Sub-class gets channel-specific info from constructor and implements 
 *   evioChannel methods.
 *
 * Some default methods throw unsupported feature exception.
 *
 * Users should be able to program to this interface.
 */
class evioChannel {

public:
  evioChannel(void) : dictionary(NULL) {}
  evioChannel(evioDictionary *dict) : dictionary(dict) {}
  virtual ~evioChannel(void) {};

  virtual void open(void) throw(evioException) = 0;

  virtual bool read(void) throw(evioException) = 0;
  virtual bool read(uint32_t *myBuf, int length) throw(evioException) = 0;
  virtual bool readAlloc(uint32_t **buffer, uint32_t *bufLen) throw(evioException) = 0;
  virtual bool readNoCopy(void) throw(evioException) = 0;
  virtual bool readRandom(uint32_t eventNumber) throw(evioException) {
    throw(evioException(0,"?evioChannel::readRandom...unsupported method"));
  }

  virtual void write(void) throw(evioException) = 0;
  virtual void write(const uint32_t* myBuf) throw(evioException) = 0;
  virtual void write(const evioChannel &channel) throw(evioException) = 0;
  virtual void write(const evioChannel *channel) throw(evioException) = 0;
  virtual void write(const evioChannelBufferizable &o) throw(evioException) = 0;
  virtual void write(const evioChannelBufferizable *o) throw(evioException) = 0;

  virtual void close(void) throw(evioException) = 0;
  virtual int ioctl(const string &request, void *argp) throw(evioException) = 0;


protected:
  evioDictionary *dictionary;


public:
  virtual const evioDictionary *getDictionary(void) const {return(dictionary);}

  virtual const uint32_t *getBuffer(void) const throw(evioException) = 0;
  virtual int getBufSize(void) const = 0;
  virtual const uint32_t *getNoCopyBuffer(void) const throw(evioException) = 0;

  virtual const uint32_t *getRandomBuffer(void) const throw(evioException) {
    return(NULL);
    //    throw(evioException(0,"?evioChannel::getRandomBuffer...unsupported method"));
  }
  virtual void getRandomAccessTable(uint32_t *** const table, uint32_t *len) const throw(evioException) {
    throw(evioException(0,"?evioChannel::getRandomAccessTable...unsupported method"));
  }

};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


} // namespace evio


#endif
