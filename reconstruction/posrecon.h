#ifndef posrecon_h
#define posrecon_h 1

typedef struct{
  double x[1728], y[1728];
}CellPosition;

CellPosition Position;
double Energy[1728];

double GetXCoord(int index) {
  double xcoord = 0.;

  if(index < 1152) {
    int copyNo = index;
    if(index >= 560) copyNo = index + 2;
    if(index >= 592) copyNo = index + 4;
    int xindex = copyNo % 34;
    xcoord = (double)(xindex - 16)*2.*1.025 - 1.025;
  }
  if(index >= 1152) {
    int copyNo = index - 1152;
    int GroupIndex = copyNo/144;
    int BlockIndex = copyNo%144;
    int fNoBlocksX = 24;
    int fNoBlocksY = 6;
    double fHalfsizeX = 1.91;
    double fHalfsizeY = 1.91;
    double  XOffset = fNoBlocksX*fHalfsizeX - 34.85;
   double  YOffset = fNoBlocksY*fHalfsizeY + 34.85;
  
    int  YIndex = BlockIndex/24;
    int  XIndex = BlockIndex%24;
 
    if(GroupIndex == 0)  xcoord = -(fNoBlocksX - 1)*fHalfsizeX + XIndex*2.*fHalfsizeX + XOffset;
    if(GroupIndex == 2)  xcoord = -(fNoBlocksX - 1)*fHalfsizeX + XIndex*2.*fHalfsizeX - XOffset;
    if(GroupIndex == 1)  xcoord = (fNoBlocksY - 1)*fHalfsizeY - YIndex*2.*fHalfsizeY - YOffset;
    if(GroupIndex == 3)  xcoord = (fNoBlocksY - 1)*fHalfsizeY - YIndex*2.*fHalfsizeY + YOffset;
  }  
  return xcoord;
}

double GetYCoord(int index) {
  double ycoord = 0.;
  if(index < 1152) {
    int copyNo = index;
    if(index >= 560) copyNo = index + 2;
    if(index >= 592) copyNo = index + 4;
    int yindex = copyNo/34;
    ycoord = (double)(16 - yindex)*2.*1.025 + 1.025;
  }

  if(index >= 1152) {
    int copyNo = index - 1152;
    int GroupIndex = copyNo/144;
    int BlockIndex = copyNo%144;
    int fNoBlocksX = 24;
    int fNoBlocksY = 6;

    double fHalfsizeX = 1.91;
    double fHalfsizeY = 1.91;
    double  XOffset = fNoBlocksX*fHalfsizeX - 34.85;
    double  YOffset = fNoBlocksY*fHalfsizeY + 34.85;

    int  YIndex = BlockIndex/24;
    int  XIndex = BlockIndex%24;

    if(GroupIndex == 0)  ycoord = (fNoBlocksY - 1)*fHalfsizeY - YIndex*2.*fHalfsizeY + YOffset;
    if(GroupIndex == 2)  ycoord = (fNoBlocksY - 1)*fHalfsizeY - YIndex*2.*fHalfsizeY - YOffset;
    if(GroupIndex == 1)  ycoord = -(fNoBlocksX - 1)*fHalfsizeX + XIndex*2.*fHalfsizeX + XOffset;
    if(GroupIndex == 3)  ycoord = -(fNoBlocksX - 1)*fHalfsizeX + XIndex*2.*fHalfsizeX - XOffset;
  }
  return ycoord;
}

double GetDistance(int index1, int index2) {
  double D;
  double x1, x2, y1, y2;
  x1 = GetXCoord(index1);
  y1 = GetYCoord(index1);
  x2 = GetXCoord(index2);
  y2 = GetYCoord(index2);

  D = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
  return D;
}

int GetXIndex(int index) {
  int copyNo = index;
 
  int xindex = copyNo % 34;

  return xindex;
}

int GetYIndex(int index) {
  int copyNo = index;
 
  int yindex = copyNo/34;  

  return yindex;
}




#endif
