#include <stdio.h>
#include<string.h>
#include<stdlib.h>
#include"matrixCalclate.h"

int main(){
    MTX mtx1,mtx2;
    //mtx1=GenerateMTX(3,3);
    mtx1=FscanMTX("mtxData.txt");
    ShowMTX(mtx1);
    //InputMTX(&mtx1);
    //printf("%f",);
    //ShowMTX(CofactorExpandMTX(0,1,mtx1));
    
    //printf("mtx=%f",DeterminantMTX(mtx1));
    mtx2=InverseMTX(mtx1);
    //printf("test\n");
    //FprintMTX(mtx2,"mtxTest2.txt");
    ShowMTX(mtx2);
}

