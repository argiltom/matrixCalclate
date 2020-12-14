#include <stdio.h>
#include<string.h>
#include<stdlib.h>
#include"matrixCalclate.h"

int main(){
    MTX mtx1=GenerateMTX(3,3);
    MTX mtx2=GenerateMTX(3,3);
    InputMTX(&mtx1);
    //InputMTX(&mtx2);
    //printf("%f",);
    ShowMTX(CofactorExpandMTX(0,1,mtx1));
    mtx2=InverseMTX(mtx1);
    ShowMTX(InverseMTX(mtx1));
    ShowMTX(MultiplyMTX(mtx1,mtx2));
}
