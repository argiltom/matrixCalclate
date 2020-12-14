#include <stdio.h>
#include<string.h>
#include<stdlib.h>
#include"matrixCalclate.h"

int main(){
    MTX mtx1=GenerateMTX(2,3);
    MTX mtx2=GenerateMTX(3,2);
    InputMTX(&mtx1);
    InputMTX(&mtx2);
    ShowMTX(MultiplyMTX(mtx1,mtx2));
}
