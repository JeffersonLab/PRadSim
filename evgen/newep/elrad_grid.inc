      integer nt1,nt2,nvv,nph
      parameter(nt1=210)
      parameter(nt2=140)
      parameter(nvv=240)
      parameter(nph=180)
      common/grid/grt1(nt1),grt2(nt2),grv(nvv,0:7),grv0(nvv),grphi(nph)
      common/dist/distsit(0:nt1+nt2),distart(0:nt1+nt2),
     .            distsiv(0:nvv),distarv(0:nvv),
     .            distsip(0:nph),distarp(0:nph)
      real*8 grt1,grt2,grv,grv0,grphi,distsit,distart,distsiv,distarv,
     .       distsip,distarp
