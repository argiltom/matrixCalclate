#ifndef matrixCalclate
#define matrixCalclate

#include<math.h>
#include<stdlib.h>
#include<stdio.h>
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
//行列にスカラーを掛ける
MTX MultiplyScalarMTX(double scalar,MTX mtx){
    MTX ret;
    ret=GenerateMTX(mtx.h,mtx.w);
    int i,j,k;
    for(i=0;i<ret.h;i++){
        for(j=0;j<ret.w;j++){
            ret.element[i][j]=mtx.element[i][j]*scalar;
        }
    }
    return ret;
}


void ShowMTX(MTX mtx){
    int i,j;
    printf("H=%d W=%d\n",mtx.h,mtx.w);
    for(i=0;i<mtx.h;i++){
        for(j=0;j<mtx.w;j++){
            printf(" %f",mtx.element[i][j]);
        }
        printf("\n");
    }
}
//第一引数:MTX 第二引数:char*
void FprintMTX(MTX mtx,char* fileName){
    int i,j;
    FILE *fp;
    fp=fopen(fileName,"w");
    if(fp==NULL) return;
    //行数と列数を書き込む1スペース挟んで
    fprintf(fp,"%d %d\n",mtx.h,mtx.w);
    //行列を書き込む
    for(i=0;i<mtx.h;i++){
        for(j=0;j<mtx.w;j++){
            fprintf(fp,"%f ",mtx.element[i][j]);
        }
        fprintf(fp,"\n");
    }
}

MTX FscanMTX(char* fileName){
    int i,j,H,W;
    FILE *fp;
    MTX ret;
    fp=fopen(fileName,"r");
    if(fp==NULL){
        printf("該当するファイルが存在しません!");
        exit(-1);
    }
    fscanf(fp,"%d %d",&H,&W);
    ret=GenerateMTX(H,W);
    for(i=0;i<ret.h;i++){
        for(j=0;j<ret.w;j++){
            fscanf(fp,"%lf",&ret.element[i][j]);
        }
    }
    return ret;
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
MTX AddMTX(MTX mtx1,MTX mtx2){
    MTX ret;
    if(mtx1.h!=mtx2.h||mtx1.w!=mtx2.w){
        printf("第一引数の行列の形と,第二引数の行列の形が一致していません!\n");
        return ret;
    }
    ret=GenerateMTX(mtx1.h,mtx1.w);
    int i,j,k;
    for(i=0;i<ret.h;i++){
        for(j=0;j<ret.w;j++){
            ret.element[i][j]=mtx1.element[i][j]+mtx2.element[i][j];
        }
    }
    return ret;
}

double Determinant(int n,double **a){//再帰関数//double hyouji(int n,double a[][n])//左から順に変数が宣言されている？・・・
    int i,j;//この関数は常に0行目について余因子展開を行う。　//bは世因子展開後のやつ
    
    if(n==1){//一次正方行列だった時
        return a[0][0];
    }
    else if(n==2){//二次正方行列だった時
        return a[0][0]*a[1][1]-a[1][0]*a[0][1];
    }

    else if(n==3){//三次正方行列だった時
        return a[0][0]*a[1][1]*a[2][2]+a[0][1]*a[1][2]*a[2][0]+a[0][2]*a[1][0]*a[2][1]-(a[0][2]*a[1][1]*a[2][0]+a[0][0]*a[1][2]*a[2][1]+a[0][1]*a[1][0]*a[2][2]);
    }
    else{//四次正方行列以上だった時
        int v=0,w=0,r,seifuhanntei;//vは行wは列を担当
        double ans=0;
        double **b;//何度も使いまわす。
        b=(double **)malloc(sizeof(double*)*(n-1));
        for(i=0;i<n-1;i++){//ここn-1じゃないと領域侵犯が発生する
            b[i]=(double *)malloc(sizeof(double)*(n-1));
        }
        //0行目について（a[0][n]について）余因子展開を行う
        
        
        for(r=0;r<n;r++){
            //a[0][r]での展開
            v=0;
            for(i=0;i<n;i++){//座標(0,r)で余因子展開　後でa[0][r]の値を乗算する。
                if(i==0) continue;
                w=0;
                for(j=0;j<n;j++){
                    if(j==r) continue;
                    b[v][w]=a[i][j];
                    w++;
                }
                v++;
            }
            //hyouji(n-1,b);//世因子展開が上手くいっているかを判定
            seifuhanntei=1;
            for(i=0;i<r;i++) seifuhanntei=seifuhanntei*(-1);//余因子展開の際の正負判定
            ans=ans+a[0][r]*Determinant(n-1,b)*seifuhanntei;
        }
        for(i=0;i<n-1;i++){//ここn-1じゃないと領域侵犯が発生する
            free(b[i]);
        }
        free(b);
        return ans;
    }
}

double DeterminantMTX(MTX mtx){
    if(mtx.h!=mtx.w){
        printf("行列式は正方行列でしか定義できません,行と列の数が違います！ 実行できません\n");
        exit(-1);
    }
    return Determinant(mtx.h,mtx.element);
}
//余因子行列を求める yとxが余因子展開の爆心地
MTX CofactorExpandMTX(int y,int x,MTX mtx){
    MTX ret=GenerateMTX(mtx.h-1,mtx.w-1);
    int i,j,k=0,l=0;
    for(i=0;i<mtx.h;i++){
        l=0;
        if(i==y) continue;
        for(j=0;j<mtx.w;j++){
            if(j==x) continue;
            ret.element[k][l]=mtx.element[i][j];
            l++;
        }
        k++;
    }
    return ret;
}

//軽量化版 返り値の値は第三引数の指す先へ返される 余因子行列を求める(逆行列用) yとxが余因子展開の爆心地
void CofactorExpandMTX2(int y,int x,MTX mtx,MTX *ret){
    int i,j,k=0,l=0;
    for(i=0;i<mtx.h;i++){
        l=0;
        if(i==y) continue;
        for(j=0;j<mtx.w;j++){
            if(j==x) continue;
            ret[0].element[k][l]=mtx.element[i][j];
            l++;
        }
        k++;
    }
}


//逆行列を求める
MTX InverseMTX(MTX mtx){
    MTX ret=GenerateMTX(mtx.h,mtx.w);
    MTX cofacter=GenerateMTX(mtx.h-1,mtx.w-1);//何度も使いまわす
    if(mtx.h!=mtx.w){
        printf("逆行列を求めるためには正方行列であることが前提です、実行できません\n");
        exit(-1);
    }
    double det=DeterminantMTX(mtx);
    if(det==0){
        printf("行列式が0の為逆行列は存在しません 全要素0の行列を返します");
        return ret;
    }
    int i,j,k;
    for(i=0;i<mtx.h;i++){
        for(j=0;j<mtx.w;j++){//P56より a*ij=(-1)^(i+j)*|Aji|
            CofactorExpandMTX2(j,i,mtx,&cofacter);//こっちの方が領域確保の回数が少なく処理が早い
            ret.element[i][j]=DeterminantMTX(cofacter)*pow(-1,i+j)/det;
            //ret.element[i][j]=DeterminantMTX(CofactorExpandMTX(j,i,mtx))*pow(-1,i+j)/det;
            if(ret.element[i][j]==0) ret.element[i][j]=0;
        }
    }
    return ret;
}
//転置行列を返す
MTX TransposeMTX(MTX mtx){
    MTX ret=GenerateMTX(mtx.w,mtx.h);
    int i,j;
    for(i=0;i<mtx.h;i++){
        for(j=0;j<mtx.w;j++){
            ret.element[j][i]=mtx.element[i][j];
        }
    }
    return ret;
}

//以下応用演算

//一般逆行列を求める
MTX GeneralInverseMTX(MTX mtx){
    MTX tMtx=TransposeMTX(mtx);
    MTX ret=MultiplyMTX(InverseMTX(MultiplyMTX(tMtx,mtx)),tMtx);
    return ret;
}
//マハラノビス距離を返す input1:クラスターの平均点(MTX(Vec)) input2:クラスターの分散共分散行列(MTX) input3:入力点(MTX(Vec))
double MahalanobisDistance(MTX ave,MTX covariance,MTX input){
    double retVal;
    MTX deviation=AddMTX(input,MultiplyScalarMTX(-1,ave));
    MTX mid = MultiplyMTX(MultiplyMTX(TransposeMTX(deviation),InverseMTX(covariance)),deviation);
    retVal=sqrt(mid.element[0][0]);
    return retVal;
}
//二次元回転行列を生成,第一引数(n)=2なら,2×2の回転行列を返し,n=3なら,同時座標表現の3×3の行列を返す 角度は第二引数で受け取り,単位はラジアン
MTX Rotation_2D_GenerateMTX(int n,double rad){
    MTX ret;
    if(n==2){
        ret=GenerateMTX(2,2);
        ret.element[0][0]=cos(rad);
        ret.element[0][1]=-sin(rad);
        ret.element[1][0]=sin(rad);
        ret.element[1][1]=cos(rad);
    }
    else if(n==3){//同時座標表現
        ret=GenerateMTX(3,3);
        ret.element[0][0]=cos(rad);
        ret.element[0][1]=-sin(rad);
        ret.element[0][2]=0;
        ret.element[1][0]=sin(rad);
        ret.element[1][1]=cos(rad);
        ret.element[1][2]=0;
        ret.element[2][0]=0;
        ret.element[2][1]=0;
        ret.element[2][2]=1;
    }else{
        printf("引数が不正です！！");
    }
    return ret;
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
