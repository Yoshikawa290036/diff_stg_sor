# 物質の移流拡散

![xyc9000000_Pe10e2_5](https://user-images.githubusercontent.com/102783602/193714485-bf9a058a-624e-42e4-b460-49f9f2e82f7e.png)

球座標系における、物質濃度の移流拡散現象を数値シミュレーションするコード。

## 支配方程式

代表流速を $1$，球の半径を $1$ とした．
ストークス流中の球形液滴周りの流れ関数 $\psi$ は，以下のように書き表される．
$$\psi = \dfrac{r^2 \sin^2{\theta}}{2} \left(1-\dfrac{2\lambda}{r}+ \dfrac{2\lambda -1}{r^3}\right),$$
ここで、 $\lambda$ は液滴内部・液滴外部の粘性比 $\kappa$ を用いて以下のように書き表される．
$$\lambda=\dfrac{3\kappa+2}{4\kappa+4}.$$

速度の半径方向成分 $u_r$ ，周方向成分 $u_\theta$ は，以下のように書き表される．
 $$ u_r=\dfrac{1}{r^2\sin{\theta}}\dfrac{\partial \psi}{\partial \theta}=\cos{\theta}\left(1-\dfrac{2\lambda}{r}+\dfrac{2\lambda-1}{r^3}\right),$$
 $$ u_\theta=-\dfrac{1}{r\sin{\theta}}\dfrac{\partial \psi}{\partial r}=\sin{\theta}\left(-1+\dfrac{\lambda}{r}+\dfrac{2\lambda-1}{2r^3}\right).$$


球座標系における，無次元化された物質濃度 $C$ の輸送方程式は，以下のように書き表される．
$$ \dfrac{\partial C}{\partial t}+u_r \dfrac{\partial C}{\partial r}+\dfrac{u_\theta}{r}\dfrac{\partial C}{\partial \theta}=\dfrac{2}{Pe}\left\{\dfrac{1}{r^2}\dfrac{\partial }{\partial r}\left(r^2 \dfrac{\partial C}{\partial r}\right) + \dfrac{1}{r^2\sin{\theta}}\dfrac{\partial }{\partial \theta}\left(\sin{\theta} \dfrac{\partial C}{\partial \theta}\right)\right\},$$
ここで，ペクレ数 $Pe$ の代表長さは球の直径である．無次元化された物質濃度 $C$ は，球表面で $1$ ，球から十分遠方で $0$ とした．
計算領域の最大半径 $R_{max}$ として，境界条件は以下のように書き表される．
$$ C=1 \;\;\;\;\; \text{at} \;\;\;\;\; r=1, $$
$$ C =0 \;\;\;\;\; \text{at} \;\;\;\;\; r=R_{max},$$
$$\dfrac{\partial C}{\partial \theta} =0 \;\;\;\;\; \text{at} \;\;\;\;\; \theta=0\;\text{and}\;\pi.$$

対流項は陽的に、拡散項は陰的(SOR法)に解く。


## 実行

#### Linux環境下でコンパイル
```
$ make
```
作成された a.out を実行する。

例えば inpというファイルを作成し、ペクレ数を入力する。実行時にそれを読み込ます。
```
$ ./a.out < inp > stdout &
```
stdout というファイルが作成され、標準出力が格納される。

stdout の中身を確認したいとき
```
$ tail -f stdout
```
