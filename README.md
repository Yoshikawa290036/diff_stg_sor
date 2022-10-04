# 物質の移流拡散

![xyc9000000_Pe10e2_5](https://user-images.githubusercontent.com/102783602/193714485-bf9a058a-624e-42e4-b460-49f9f2e82f7e.png)

スタッガード格子で定義

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
