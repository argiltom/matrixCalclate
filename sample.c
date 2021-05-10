#include <stdio.h>
#include<string.h>
#include<stdlib.h>
#include"matrixCalclate.h"
//サンプルコード
//コンパイルコマンド: gcc sample.c -lm
int main(){
    MTX mtx1,mtx2;
    //行列をファイルから入力する
    
    char* stringFileName="InputMatrix.txt";//ここにファイル名を入力

    mtx1=FscanMTX(stringFileName);
    printf("%s",stringFileName);
    printf("の行列↓\n");
    ShowMTX(mtx1);
    printf("この行列の逆行列は:\n");
    mtx2=InverseMTX(mtx1);
    ShowMTX(mtx2);

    
    FprintMTX(mtx2,"outputResult.txt");
    printf("outputResult.txtというファイル名で逆行列の演算結果をカレントディレクトリに出力しました");

    printf("元の行列とその逆行列の積↓\n");
    ShowMTX(MultiplyMTX(mtx1,mtx2));
    return 0;
}