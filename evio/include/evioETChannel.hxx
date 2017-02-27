// evioETChannel.hxx

//  Author:  Elliott Wolin, JLab, 27-Nov-2012


#ifndef _evioETChannel_hxx
#define _evioETChannel_hxx


#include <iostream>
#include <stdint.h>
#include <evioChannel.hxx>
#include <evioUtil.hxx>
#include <evio.h>

extern "C" {
#include "et.h"
}


using namespace std;


namespace evio {


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Implements evioChannel functionality for I/O to and from ET system
 */
class evioETChannel : public evioChannel {

public:
  evioETChannel(et_sys_id et_system_id, et_att_id et_attach_id, const string &mode = "r", int chunk=1, int et_mode=ET_SLEEP)
    throw(evioException);
  evioETChannel(et_sys_id et_system_id, et_att_id et_attach_id, evioDictionary *dict, const string &mode = "r", 
                int chunk=1, int et_mode=ET_SLEEP) throw(evioException);
  virtual ~evioETChannel(void);


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
  int getChunkSize(void) const;
  string getBufferXMLDictionary(void) const;


private:
  et_sys_id et_system_id;        /**<ET system id.*/
  et_att_id et_attach_id;        /**<ET attach id.*/

  et_event **ETBuffers;          /**<Array to hold ET buffer pointers.*/
  string mode;                   /**<Open mode, "r" or "rw" or "w".*/
  int chunk;                     /**<Number of ET buffers to fetch at one go.*/
  int et_mode;                   /**<ET new/get mode.*/
  string bufferXMLDictionary;    /**<XML dictionary in buffer.*/
  bool createdBufferDictionary;  /**<true if dictionary created from buffer.*/
  
  int etBufReceived;             /**<Number of ET buffers received.*/
  int etBufUsed;                 /**<Number of ET buffers used.*/
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

} // namespace evio


#endif
