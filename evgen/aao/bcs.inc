c      include file BCS
c
c_begin_doc
c  RCS ID string
c  $Id: bcs.inc,v 1.1 2000/05/18 23:07:48 lcsmith Exp $
c
c  Documentation for include file BCS
c
c  Purpose: Standard common for BOS banks
c  --------
c
c  Used in routines: all
c  ----------------
c
c  Notes:
c  ------
c
c  Author:   Dieter Cords      Created:  Mon Jun 12 17:50:51 EDT 1995
c  -------
c
c  Major revisions:
c  JAM - changed the parameter to 700,000
c  ----------------
c
c_end_doc
c
c
c  Remember NO NESTED INCLUDES
c
c  COMMON BLOCK NAME should be the same as include file name!
c
c_begin_var
c  variables for include file BCS:
c  -----------------------------------
c  RCS information: 
      CHARACTER*(132)  CIBCS
      DATA CIBCS/
     1'$Id: bcs.inc,v 1.1 2000/05/18 23:07:48 lcsmith Exp $'
     2/
c----6----------------------------------------------------------------72
c  
c  Variables contained in BCS common block
c  --------------------------------------------
c
c_end_var
c
      INTEGER Nbcs
      PARAMETER (Nbcs = 700000)
c
      INTEGER IW(Nbcs), JUNK(5)
      REAL RW(Nbcs)
      INTEGER*2 IW16(2*Nbcs)
      BYTE IW08(4*Nbcs)
      EQUIVALENCE (IW(1),RW(1),IW08(1),IW16(1))
c 
c IWCOL and IWROW provide a quick way to find the number of columns
c and rows for a particular bank
c
c
      COMMON/BCS/JUNK,IW
      INTEGER IWROW(Nbcs), IWCOL(Nbcs), IWNAME(Nbcs), IWNR(Nbcs), 
     1        IWNXT(Nbcs)
      EQUIVALENCE (IW(1), IWCOL(6))
      EQUIVALENCE (IW(1), IWROW(5))
      EQUIVALENCE (IW(1), IWNAME(4))
      EQUIVALENCE (IW(1), IWNR(3))
      EQUIVALENCE (IW(1), IWNXT(2))
c
c
      save /BCS/
c----6----------------------------------------------------------------72




