// evioBufferChannel.hxx

//  Author:  Elliott Wolin, JLab, 12-Apr-2012


#ifndef _evioBufferChannel_hxx
#define _evioBufferChannel_hxx


#include <iostream>
#include <stdint.h>
#include "evioChannel.hxx"
#include "evioUtil.hxx"
#include "evio.h"


using namespace std;


namespace evio {


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Implements evioChannel functionality for I/O to and from user-supplied evio buffer.
 */
class evioBufferChannel : public evioChannel {

public:
  evioBufferChannel(uint32_t *streamBuf, int bufLen, const string &mode = "r", int size=100000) throw(evioException);
  evioBufferChannel(uint32_t *streamBuf, int bufLen, evioDictionary *dict, const string &mode = "r", int size=100000) throw(evioException);
  virtual ~evioBufferChannel(void);


  void open(void) throw(evioException);

  bool read(void) throw(evioException);
  bool read(uint32_t *myEventBuf, int length) throw(evioException);
  bool readAlloc(uint32_t **buffer, uint32_t *bufLen) throw(evioException);
  bool readNoCopy(void) throw(evioException);
  bool readRandom(uint32_t bufferNumber) throw(evioException);

  void write(void) throw(evioException);
  void write(const uint32_t *myEventBuf) throw(evioException);
  void write(const evioChannel &channel) throw(evioException);
  void write(const evioChannel *channel) throw(evioException);
  void write(const evioChannelBufferizable &o) throw(evioException);
  void write(const evioChannelBufferizable *o) throw(evioException);

  void close(void) throw(evioException);
  int ioctl(const string &request, void *argp) throw(evioException);


  const uint32_t *getBuffer(void) const throw(evioException);
  int getBufSize(void) const;
  const uint32_t *getNoCopyBuffer(void) const throw(evioException);

  const uint32_t *getRandomBuffer(void) const throw(evioException);
  void getRandomAccessTable(uint32_t *** const table, uint32_t *len) const throw(evioException);


  const uint32_t *getStreamBuffer(void) const throw(evioException);
  int getStreamBufSize(void) const;
  string getMode(void) const;
  uint32_t getEVIOBufferLength(void) const throw(evioException);
  string getBufferXMLDictionary(void) const;


private:
  uint32_t *streamBuf;           /**<Pointer to user-supplied stream i/o buffer.*/
  int streamBufSize;             /**<Size of user-supplied stream buffer.*/
  string mode;                   /**<Open mode, "r" or "ra" or "w" or "a".*/
  int handle;                    /**<Internal evio handle.*/
  uint32_t *buf;                 /**<Pointer to internal event buffer.*/
  int bufSize;                   /**<Size of internal buffer.*/
  const uint32_t *noCopyBuf;     /**<Pointer to no copy buffer.*/
  const uint32_t *randomBuf;     /**<Pointer to random read buffer.*/
  string bufferXMLDictionary;    /**<XML dictionary in buffer.*/
  bool createdBufferDictionary;  /**<true if dictionary created from buffer.*/
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

} // namespace evio


#endif
