      integer function strlen(st)
      integer i
      character st*(*)
      i = len(st)
      do while (st(i:i) .eq. ' ')
        i = i - 1
      enddo
      strlen = i
      return
      end
