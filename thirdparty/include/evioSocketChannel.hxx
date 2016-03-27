// evioSocketChannel.hxx

//  Author:  Elliott Wolin, JLab, 12-Apr-2012


#ifndef _evioSocketChannel_hxx
#define _evioSocketChannel_hxx


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
 * Implements evioChannel functionality for I/O to and from socket.
 */
class evioSocketChannel : public evioChannel {

public:
  evioSocketChannel(int socFd, const string &mode = "r", int size=100000) throw(evioException);
  evioSocketChannel(int socFd, evioDictionary *dict, const string &mode = "r", int size=100000) throw(evioException);
  virtual ~evioSocketChannel(void);

  void open(void) throw(evioException);

  bool read(void) throw(evioException);
  bool read(uint32_t *myEventBuf, int length) throw(evioException);
  bool readAlloc(uint32_t **buffer, uint32_t *bufLen) throw(evioException);
  bool readNoCopy(void) throw(evioException);

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

  string getMode(void) const;
  int getSocketFD(void) const throw(evioException);
  string getSocketXMLDictionary(void) const;


private:
  int sockFD;                     /**<Socket file descriptor.*/
  string mode;                    /**<Open mode, "r" or "w".*/
  int handle;                     /**<Internal evio handle.*/
  uint32_t *buf;                  /**<Pointer to internal event socket.*/
  int bufSize;                    /**<Size of internal socket.*/
  const uint32_t *noCopyBuf;      /**<Pointer to no copy buffer.*/
  string socketXMLDictionary;     /**<XML dictionary in socket.*/
  bool createdSocketDictionary;   /**<true if created dictionary from socket.*/
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

} // namespace evio


#endif
