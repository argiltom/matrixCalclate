#ifndef matrixCalclate
#define matrixCalclate

#include<math.h>
#include<stdlib.h>
//これは,行列演算用の自作ライブラリです。
typedef struct {
    //行  行数
    int h;//行
    //列  列数
    int w;//列
    //行列の要素
    double **element;
} MTX;
//行列を生成し,全要素0で返す
MTX GenerateMTX(int h,int w){
    MTX ret;
    int i,j;
    ret.h=h;
    ret.w=w;
    ret.element=(double**)malloc(sizeof(double*)*h);
    for(i=0;i<h;i++){
        ret.element[i]=(double*)malloc(sizeof(double)*w);
    }
    for(i=0;i<h;i++){
        for(j=0;j<w;j++){
            ret.element[i][j]=0;
        }
    }
    return ret;
}

MTX MultiplyMTX(MTX mtx1,MTX mtx2){
    MTX ret;
    if(mtx1.w!=mtx2.h){
        printf("第一引数の列数と,第二引数の行数が一致していません!\n");
        return ret;
    }
    ret=GenerateMTX(mtx1.h,mtx2.w);
    int i,j,k;
    for(i=0;i<ret.h;i++){
        for(j=0;j<ret.w;j++){
            for(k=0;k<mtx1.w;k++){//mtx2.hでも良い
                ret.element[i][j]+=mtx1.element[i][k]*mtx2.element[k][j];
            }
        }
    }
    return ret;
}

void ShowMTX(MTX mtx){
    int i,j;
    for(i=0;i<mtx.h;i++){
        for(j=0;j<mtx.w;j++){
            printf(" %f",mtx.element[i][j]);
        }
        printf("\n");
    }
}

void InputMTX(MTX *mtx){
    int i,j;
    for(i=0;i<mtx->h;i++){
        for(j=0;j<mtx->w;j++){
            printf("mtx[%d][%d]=",i,j);
            scanf("%lf",&(mtx->element[i][j]));
        }
        printf("\n");
    }
}


//複素数演算
typedef struct complex_number{
    double Re;
    double Im;
} Complex_Num;
Complex_Num Complex_add(Complex_Num a,Complex_Num b){
    Complex_Num ans;
    ans.Re=a.Re+b.Re;
    ans.Im=a.Im+b.Im;
    return ans;
}
Complex_Num Complex_sub(Complex_Num a,Complex_Num b){//a-b
    Complex_Num ans;
    ans.Re=a.Re-b.Re;
    ans.Im=a.Im-b.Im;
    return ans;
}
Complex_Num Complex_conjugate(Complex_Num a){
    Complex_Num ans;
    ans.Re=a.Re;
    ans.Im=a.Im*-1;
    return ans;
}
Complex_Num Complex_multiply(Complex_Num a,Complex_Num b){
    Complex_Num ans;
    ans.Re=a.Re*b.Re-a.Im*b.Im;
    ans.Im=a.Im*b.Re+a.Re*b.Im;
    return ans;
}
Complex_Num Complex_divid(Complex_Num a,Complex_Num b){//a/b
    Complex_Num ans;
    ans.Re=(a.Re*b.Re+a.Im*b.Im)/(b.Re*b.Re+b.Im*b.Im);
    ans.Im=(a.Im*b.Re-a.Re*b.Im)/(b.Re*b.Re+b.Im*b.Im);
    return ans;
}
double Complex_norm(Complex_Num a){
    return sqrt(a.Re*a.Re+a.Im*a.Im);
}
Complex_Num W(int m,int N){//回転子 //引数は分子　分母の順
    Complex_Num ans;
    ans.Re=cos(2*M_PI*m/N);
    ans.Im=-sin(2*M_PI*m/N);
    return ans;
}

int block(int n,Complex_Num *input,Complex_Num *output){//n=2^mで固定
    int i;
    //F(i)=f[i%(N/2)]+W(i,N)*f[i%(N/2)+N/2]

    //イメージ//output[i]=input[i%(n/2)]+W(i,n)*input[i%(n/2)+n/2];
    for(i=0;i<n;i++){
    output[i]=Complex_add(input[i%(n/2)],Complex_multiply(W(i,n),input[i%(n/2)+n/2]));
    }
    return n;
}
//ノート参照
int inverse_block(int n,Complex_Num *input,Complex_Num *output){//n=2^mで固定 //動作確認済み
    int i,A;
    //F(i)=f[i%(N/2)]+W(i,N)*f[i%(N/2)+N/2]
    Complex_Num one;
        one.Re=1;
        one.Im=0;
    //イメージ//output[i]=input[i%(n/2)]+W(i,n)*input[i%(n/2)+n/2];
    for(i=0;i<n;i++){
        A=i%(n/2);
        
        if(i<n/2){
            //イメージ//output[bit_inverse[i]]=(input[iN2]-W(iN2,n)*input[iN2+n/2]/W(iN2+n/2,n))/(1-W(iN2,n)/W(iN2+n/2,n));
            //miss//output[i]=Complex_divid( Complex_sub( input[A], Complex_multiply( Complex_divid(W(A,n),W(A+n/2)),input[A+n/2] ) ),(Complex_sub(one,Complex_divid(W(A,n),W(A+n/2)))) );
           output[i]=Complex_divid(Complex_sub(Complex_multiply(W(A+n/2,n),input[A]),Complex_multiply(W(A,n),input[A+n/2])),Complex_sub(W(A+n/2,n),W(A,n)));
           // printf("[i]=%d\n",bit_inverse[i]);
        }
        if(i>=n/2){
            //イメージ//output[bit_inverse[i]]=(input[iN2]-input[iN2+n/2])/(W(iN2,n)-W(iN2+n/2,n));
            //miss//output[i]=Complex_divid( Complex_sub(input[A],input[A+n/2]),Complex_sub(W(A,n),W(A+n/2,n)) );
            output[i]=Complex_divid(Complex_sub(input[A],input[A+n/2]),Complex_sub(W(A,n),W(A+n/2,n)));
            // printf("[i]=%d\n",bit_inverse[i]);
        }
    }
    return n;
}



void FFT(int n,Complex_Num *input,Complex_Num *output){ 
    //並び替え
    int m=0,step,i,bit_inverse[n],j,bit;//mにはlog(2,n)が格納される //bitはbit演算用のtemp変数
    

    while(n!=pow(2,m)){//2^mを検出
        m++;
        if(m>n){
            printf("入力された配列数は2^mを満たしていません 処理を中断します");
            return;
        }
    }
    
    //インプットとアウトプットの間の中間配列の配列[step総数][配列n]を確保　//[あるステップ][配列nの要素]とできる　//ステップ毎の配列要素が残る
    Complex_Num **temp;
    temp=(Complex_Num **)malloc(sizeof(Complex_Num*)*m);//step==mの処理はoutputへの代入だけだからm+1個分は要らない
    for(i=0;i<m;i++){
        temp[i]=(Complex_Num *)malloc(sizeof(Complex_Num)*n);
    }

    for(i=0;i<n;i++){//初期化
        bit_inverse[i]=0;
    }


    for(step=0;step<=m;step++){
        if(step==0){//input配列の並び替えを行い　並び替えた結果をtemp[0][i]に代入
            //配列番号をビット列に直し,逆から読んであげた値をtemp[0][i]に代入すればいい
            bit=pow(2,m-1);//100000 のような、　文字数mのバイナリデータを作成　j分左にずらしながら&をとり,1が入力数値(i)と重なるところがあるのなら,1をｊ分、右にずらして加算する
            for(i=0;i<n;i++){
                for(j=0;j<m;j++){
                    if(i&(bit>>j)) bit_inverse[i]+=1<<j;//左からjつ目を見て,それが1なら,　1をjつ右にずらして加算すれば線対照的にbitを表現できる
                }
                //printf("bit_inverse[%d]=%d\n",i,bit_inverse[i]);//テスト
            }

            for(i=0;i<n;i++){
                temp[0][i]=input[bit_inverse[i]];
            }
            
        }
        else if(step!=m){//ブロック関数を呼び,temp[step][i]に temp[step-1][i]のblock()した結果を代入
            for(i=0;i<n;i=i+pow(2,step)){//2^stepの長さに区切って、それぞれのブロックごとの先頭ポインタを飛ばし、blook関数で処理をする
            block(pow(2,step),&temp[step-1][i],&temp[step][i]);
            }
        }
        else if(step==m){//output関数にtemp[step-1][i]をブロックにかけた結果を代入する
            for(i=0;i<n;i=i+pow(2,step)){//2^stepの長さに区切って、それぞれのブロックごとの先頭ポインタを飛ばし、blook関数で処理をする
            block(pow(2,step),&temp[step-1][i],output);
            }
        }
    }
    //free
    for(i=0;i<m;i++){
        free(temp[i]);
    }
    free(temp);
    

   
}

void iFFT(int n,Complex_Num *input,Complex_Num *output){ //動作確認済み7/4、iFFT完成セリ！！！
    //並び替え
    int i,A,bit_inverse[n],j,bit,step,m=0;
    while(n!=pow(2,m)){//2^mを検出
        m++;
        if(m>n){
            printf("入力された配列数は2^mを満たしていません 処理を中断します");
            return;
        }
    }
    for(i=0;i<n;i++){//初期化
        bit_inverse[i]=0;
    }

    bit=pow(2,m-1);//100000 のような、　文字数mのバイナリデータを作成　j分左にずらしながら&をとり,1が入力数値(i)と重なるところがあるのなら,1をｊ分、右にずらして加算する
            for(i=0;i<n;i++){
                for(j=0;j<m;j++){
                    if(i&(bit>>j)) bit_inverse[i]+=1<<j;//左からjつ目を見て,それが1なら,　1をjつ右にずらして加算すれば線対照的にbitを表現できる
                }
                //printf("bit_inverse[%d]=%d\n",i,bit_inverse[i]);
            }

    Complex_Num **temp;
    temp=(Complex_Num **)malloc(sizeof(Complex_Num*)*(m+1));//step==mの処理はoutputへの代入だけだからm+1個分は要らない
    for(i=0;i<=m;i++){//ここだ！　セグメントエラーの原因！！
        temp[i]=(Complex_Num *)malloc(sizeof(Complex_Num)*n);
    }

    //メインの処理,
    for(i=0;i<n;i++){
            temp[m][i]=input[i];
        }
    for(step=m;step>=0;step--){
        if(step!=0){//インバースブロック関数を呼び
            for(i=0;i<n;i=i+pow(2,step)){//2^stepの長さに区切って、それぞれのブロックごとの先頭ポインタを飛ばし、blook関数で処理をする
            inverse_block(pow(2,step),&temp[step][i],&temp[step-1][i]);
            }
        }
        else if(step==0){//結果を入力する
            for(i=0;i<n;i++){
            output[bit_inverse[i]]=temp[0][i];
            }
        }
    }
    //free
    for(i=0;i<m+1;i++){
        free(temp[i]);
    }
    free(temp);



}

void DFT(int n,Complex_Num *input,Complex_Num *output){
    int i,k;
    Complex_Num temp;
    for(i=0;i<n;i++){
        output[i].Re=0;
        output[i].Im=0;
    }
    for(i=0;i<n;i++){
        for(k=0;k<n;k++){
            temp.Re=cos(2*M_PI*k*i/n);
            temp.Im=-sin(2*M_PI*k*i/n);
            output[i]=Complex_add(output[i],Complex_multiply(input[k],temp));
        }
    }
}

void iDFT(int n,Complex_Num *input,Complex_Num *output){
    int i,k;
    Complex_Num temp;
    for(i=0;i<n;i++){
        output[i].Re=0;
        output[i].Im=0;
    }
    for(i=0;i<n;i++){
        for(k=0;k<n;k++){
            temp.Re=cos(2*M_PI*k*i/n);
            temp.Im=sin(2*M_PI*k*i/n);
            output[i]=Complex_add(output[i],Complex_multiply(input[k],temp));
        }
        output[i].Re=output[i].Re/n;
        output[i].Im=output[i].Im/n;
    }
}

int DCT(int n,double *input,double *output){//input とoutputは配列
    int i,k;
    //outputをイニシャライズ
    for(i=0;i<n;i++){
        output[i]=0;
    }
    

    for(k=0;k<n;k++){ //n-1にしたのが間違えだった・・・反省！6/25 
        for(i=0;i<n;i++){
            output[k]+=input[i]*cos( ((2*i+1)*k*M_PI)/(2*n) );
        }
    }
    return 0;
}

int iDCT(int n,double *input,double *output){//input とoutputは配列
    int i,k;
    //outputをイニシャライズ
    for(i=0;i<n;i++){
        output[i]=0;
    }


    for(i=0;i<n;i++){ //ループの含有関係をDCTと逆にした　この関数を編集する際はインクリメント変数に要注意！！！！
        for(k=0;k<n;k++){

            if(k==0)output[i]+=(0.5*input[k])*2/n; //元の式のnはここではiとなる

            if(k>0)output[i]+=input[k]*cos( ((2*i+1)*k*M_PI)/(2*n) )*2/n;
        }
    }
    return 0;
}

void dbl_dim_DCT(int w,int h,double **input,double **output){//横縦の順番で引数を受け取る　配列も[x][y]で考える //二次元DCT
    int i,j;
    double temp[w][h],Ttemp[h][w],Tresult[h][w];//転置行列を定義する
    //縦方向へ変換
    for(i=0;i<w;i++){
        DCT(h,input[i],temp[i]);
    }
    //転置する
    for(i=0;i<w;i++){
        for(j=0;j<h;j++){
            Ttemp[j][i]=temp[i][j];
        }
    }
    //横方向へ変換
    for(j=0;j<h;j++){
        DCT(w,Ttemp[j],Tresult[j]);
    }
    //転地して結果に代入
    for(i=0;i<w;i++){
        for(j=0;j<h;j++){
            output[i][j]=Tresult[j][i];
        }
    }
}



void dbl_dim_iDCT(int w,int h,double **input,double **output){//横縦の順番で引数を受け取る　配列も[x][y]で考える //二次元逆DCT
    int i,j;
    double temp[w][h],Ttemp[h][w],Tresult[h][w];//転置行列を定義する
    //縦方向へ変換
    for(i=0;i<w;i++){
        iDCT(h,input[i],temp[i]);
    }
    //転置する
    for(i=0;i<w;i++){
        for(j=0;j<h;j++){
            Ttemp[j][i]=temp[i][j];
        }
    }
    //横方向へ変換
    for(j=0;j<h;j++){
        iDCT(w,Ttemp[j],Tresult[j]);
    }
    //転地して結果に代入
    for(i=0;i<w;i++){
        for(j=0;j<h;j++){
            output[i][j]=Tresult[j][i];
        }
    }
}

void dbl_dim_DFT(int w,int h,Complex_Num **input,Complex_Num **output){//横縦の順番で引数を受け取る　配列も[x][y]で考える //二次元DCT
    int i,j;
    Complex_Num temp[w][h],Ttemp[h][w],Tresult[h][w];//転置行列を定義する
    //縦方向へ変換
    for(i=0;i<w;i++){
        DFT(h,input[i],temp[i]);
    }
    //転置する
    for(i=0;i<w;i++){
        for(j=0;j<h;j++){
            Ttemp[j][i]=temp[i][j];
        }
    }
    //横方向へ変換
    for(j=0;j<h;j++){
        DFT(w,Ttemp[j],Tresult[j]);
    }
    //転地して結果に代入
    for(i=0;i<w;i++){
        for(j=0;j<h;j++){
            output[i][j]=Tresult[j][i];
        }
    }
}

void dbl_dim_iDFT(int w,int h,Complex_Num **input,Complex_Num **output){//横縦の順番で引数を受け取る　配列も[x][y]で考える //二次元DCT
    int i,j;
    Complex_Num temp[w][h],Ttemp[h][w],Tresult[h][w];//転置行列を定義する
    //縦方向へ変換
    for(i=0;i<w;i++){
        iDFT(h,input[i],temp[i]);
    }
    //転置する
    for(i=0;i<w;i++){
        for(j=0;j<h;j++){
            Ttemp[j][i]=temp[i][j];
        }
    }
    //横方向へ変換
    for(j=0;j<h;j++){
        iDFT(w,Ttemp[j],Tresult[j]);
    }
    //転地して結果に代入
    for(i=0;i<w;i++){
        for(j=0;j<h;j++){
            output[i][j]=Tresult[j][i];
        }
    }
}

void dbl_dim_FFT(int w,int h,Complex_Num **input,Complex_Num **output){//横縦の順番で引数を受け取る　配列も[x][y]で考える //二次元DCT
    int i,j;
    Complex_Num temp[w][h],Ttemp[h][w],Tresult[h][w];//転置行列を定義する
    //縦方向へ変換
    for(i=0;i<w;i++){
        FFT(h,input[i],temp[i]);
    }
    //転置する
    for(i=0;i<w;i++){
        for(j=0;j<h;j++){
            Ttemp[j][i]=temp[i][j];
        }
    }
    //横方向へ変換
    for(j=0;j<h;j++){
        FFT(w,Ttemp[j],Tresult[j]);
    }
    //転地して結果に代入
    for(i=0;i<w;i++){
        for(j=0;j<h;j++){
            output[i][j]=Tresult[j][i];
        }
    }
}

void dbl_dim_iFFT(int w,int h,Complex_Num **input,Complex_Num **output){//動作確認//横縦の順番で引数を受け取る　配列も[x][y]で考える //二次元DCT
    int i,j;
    Complex_Num temp[w][h],Ttemp[h][w],Tresult[h][w];//転置行列を定義する
    //縦方向へ変換
    for(i=0;i<w;i++){
        iFFT(h,input[i],temp[i]);
    }
    //転置する
    for(i=0;i<w;i++){
        for(j=0;j<h;j++){
            Ttemp[j][i]=temp[i][j];
        }
    }
    //横方向へ変換
    for(j=0;j<h;j++){
        iFFT(w,Ttemp[j],Tresult[j]);
    }
    //転地して結果に代入
    for(i=0;i<w;i++){
        for(j=0;j<h;j++){
            output[i][j]=Tresult[j][i];
        }
    }
}



#endif
