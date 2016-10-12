// evioFileChannel.hxx

//  Author:  Elliott Wolin, JLab, 18-feb-2010


#ifndef _evioFileChannel_hxx
#define _evioFileChannel_hxx


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
 * Implements evioChannel functionality for I/O to and from files.
 * Basically a wrapper around the original evio C library.
 */
class evioFileChannel : public evioChannel {

public:
  evioFileChannel(const string &fileName, const string &mode = "r", int size = 1000000) throw(evioException);
  evioFileChannel(const string &fileName, evioDictionary *dict,
                  const string &mode = "r", int size = 1000000) throw(evioException);
  evioFileChannel(const string &fileName, evioDictionary *dict, const uint32_t *firstEvent,
                  const string &mode = "w", int size = 1000000) throw(evioException);
  virtual ~evioFileChannel(void);


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

  string getFileName(void) const;
  string getMode(void) const;
  string getFileXMLDictionary(void) const;


private:
  string filename;            /**<Name of evio file.*/
  string mode;                /**<Open mode, "r" for read, "ra" for random access read,
                                * "w" for write, "a" for append, "s" for splitting while writing.*/
  int handle;                 /**<Internal evio handle.*/
  uint32_t *buf;              /**<Pointer to internal event buffer.*/
  int bufSize;                /**<Size of internal event buffer.*/
  const uint32_t *firstEvent; /**<Pointer first event buffer.*/
  const uint32_t *noCopyBuf;  /**<Pointer to no copy event buffer.*/
  const uint32_t *randomBuf;  /**<Pointer to random read buffer.*/
  string fileXMLDictionary;   /**<XML dictionary in file.*/
  bool createdFileDictionary; /**<true if internally created new dictionary from file.*/
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

} // namespace evio


#endif
