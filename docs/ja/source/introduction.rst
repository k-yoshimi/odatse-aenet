はじめに
========

odatse-aenet とは
------------------

odatse-aenet は、機械学習ポテンシャル `AENET <http://ann.atomistic.net/>`_ (Atom-centered Neural Network) を
`ODAT-SE <https://github.com/issp-center-dev/ODAT-SE>`_ (Open Data Analysis Tool for Science and Engineering) フレームワークと統合するためのソルバーモジュールです。

ODAT-SE が提供する最適化アルゴリズムを用いて、AENET で構築した機械学習ポテンシャルによるパラメータ探索を実行できます。

利用可能なアルゴリズム
~~~~~~~~~~~~~~~~~~~~~~

ODAT-SE フレームワークでは、以下のアルゴリズムが利用可能です:

- Nelder-Mead 法 (``minsearch``)
- グリッド探索法 (``mapper``)
- ベイズ最適化 (``bayes``)
- レプリカ交換モンテカルロ法 (``exchange``)
- ポピュレーションアニーリングモンテカルロ法 (``pamc``)

本パッケージの構成
~~~~~~~~~~~~~~~~~~

- **AenetSolver**: ODAT-SE のソルバークラスとして AENET predict.x を実行し、エネルギー計算を行います。

開発者
------

- **谷田秀哉** （鳥取大学大学院 持続性社会創成科学研究科） --- 初期コード開発（修士論文）
- **星健夫** （核融合科学研究所） --- ODAT-SE 共同開発、研究指導
- **吉見一慶** （東京大学物性研究所） --- コード整理・パッケージ化

ライセンス
----------

本ソフトウェアは Mozilla Public License 2.0 (MPL-2.0) のもとで公開されています。

引用
----

ODAT-SE を利用した研究成果を公表する際は、以下の文献を引用してください:

  Y. Motoyama, K. Yoshimi, I. Mochizuki, H. Iwamoto, H. Ichinose, and T. Hoshi,
  "Data-analysis software framework 2DMAT and its application to experimental measurements for two-dimensional material structures",
  Computer Physics Communications **280**, 108465 (2022).
