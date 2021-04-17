# Subroutines

UMAT を用いた Abaqus のサブルーチンを開発するリポジトリ

## Folder tree

<pre>
├─Datas                 各種パラメータを格納するフォルダ
├─Images                matplotlibで作成した画像や各種画像の保存用フォルダ
├─src_fortran           fortran77ファイルを格納したフォルダ
└─src_python            pythonファイルを格納したフォルダ
</pre>

### 上記以外のファイル

- [.devcontainer.json](.devcontainer.json)  
  vscode 拡張機能である remote-container の設定用ファイル。コンテナの設定を記述。
- [setup.cfg](setup.cfg)  
  python 環境の設定を記述するファイル。現状 formatter と linter の設定だけ記述している。

## Python 環境構築

1. Docker,docker-compose をインストール
1. vscode の拡張機能`ms-vscode-remote.remote-containers`をインストール
1. 画面左下の><を押し、`Remote-Containers: Open Folder in Container`を選択

初回はコンテナの生成に時間がかかる。右下に表示されるプログレスバーが終わったら環境構築完了

## src_python フォルダ

- [draw.py](src_python/draw.py)  
  降伏局面を描く処理を記述している。matplotlib の contour というメソッドを利用。
- [identify_hill48_params.py](src_python/identify_hill48_params.py),
  [identify_yld2004_params.py](src_python/identify_yld2004_params.py)  
  各降伏関数で利用する異方性パラメータを計算し、また各方向の降伏応力、r 値を計算する処理を記述。各方向の計算は[Barlat 氏の yld2004-18p に関する論文]("https://www.sciencedirect.com/science/article/abs/pii/S0749641904001160)の appendix を参考にしている。
- その他のファイル  
  一応残しているだけで大した意味はないです。

## src_fortran フォルダ

- [BaseSubroutine.for](src_fortran/BaseSubroutine.for)  
  subroutine を作成する上で base となる subroutine
- [HardeningRules.for](src_fortran/HardeningRules.for)  
  硬化則についてその微分も併せて記述したファイル。subroutine を作成する際適宜コピーしてください。
- [Lib.for](stc_fortran/Lib.for)  
  数値計算を行う上で便利な関数を用意したもの
- [test.for](src_fortran/test.for)  
  fortran77 の挙動を確認するために用意したもの。これをコンパイルするための Makefile も用意している。使い方は
  1. gfortran をインストール
  1. test.for を適当に記述
  1. `make`コマンドを打つ
  1. エラーが発生しなければ Debug フォルダ下に実行用ファイルが作成される
- [TestSubroutine.for](src_fortran/TestSubroutine.for)  
  Abaqus で subroutine が動かせるか確かめるために用意したもの。使う必要ないです。
- [YieldFunctions.for](src_fortran/YieldFunctions.for)  
  降伏関数とその一次微分、二次微分を記述したもの。subroutine を作成する際適宜コピーしてください。
- [Yld2004Yld2004Swift.for](src_fortran/Yld2004Yld2004Swift.for)  
  Yld2004-18P を非関連流れ則に基づいて実装した subroutine

## ABAQUS documentation

http://130.149.89.49:2080/v2016/index.html
