#include <stdio.h>
#include<string.h>
#include<stdlib.h>
#include"matrixCalclate.h"

int main(){
    MTX mtx1=GenerateMTX(2,2);
    MTX mtx2=GenerateMTX(2,2);
    MTX mtx3=GenerateMTX(2,1);
    InputMTX(&mtx1);
    InputMTX(&mtx3);
    //InputMTX(&mtx2);
    //printf("%f",);
    //ShowMTX(CofactorExpandMTX(0,1,mtx1));
    mtx2=InverseMTX(mtx1);
    ShowMTX(MultiplyMTX(mtx2,mtx3));
}
