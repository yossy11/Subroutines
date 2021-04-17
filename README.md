# Subroutines

UMAT を用いた Abaqus のサブルーチンを開発するリポジトリ

## Folder tree

<pre>
├─Datas                 各種パラメータを格納するフォルダ
├─Images                matplotlibで作成した画像や各種画像の保存用
├─src_fortran           fortran77ファイル
└─src_python            pythonファイル
</pre>

## Python 環境構築

1. Docker,docker-compose をインストール
1. vscode の拡張機能`ms-vscode-remote.remote-containers`をインストール
1. 画面左下の><を押し、`Remote-Containers: Open Folder in Container`を選択

初回はコンテナの生成に時間がかかる。右下に表示されるプログレスバーが終わったら環境構築完了

## Python ファイル一覧

- [draw.py]("./blob/master/src_python/draw.py")  
  降伏局面を描く処理を記述している。matplotlib の contour というメソッドを利用。
- identify_hill48_params.py  
  identify_yld2004_params.py  
  各降伏関数で利用する異方性パラメータを計算し、また各方向の降伏応力、r 値を計算する処理を記述。計算は[Barlat 氏の yld2004-18p に関する論文]("https://www.sciencedirect.com/science/article/abs/pii/S0749641904001160)の appendix を参考している。
- その他のファイル  
  一応残しているだけで大した意味はないです。

## ABAQUS documentation

http://130.149.89.49:2080/v2016/index.html
