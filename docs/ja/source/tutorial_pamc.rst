チュートリアル: AENET + PAMC による N2 二量体の最適化
=====================================================

本チュートリアルでは、AENET で構築した機械学習ポテンシャルと ODAT-SE の
ポピュレーションアニーリングモンテカルロ法（PAMC）を組み合わせて、
N2 二量体の最安定構造を探索する方法を説明します。

サンプルファイルは ``sample/aenet_pamc/`` にあります。

前提条件
--------

- AENET の ``predict.x`` がインストール済みであること
- 学習済みの ANN ポテンシャルファイル（``N.5t-5t.ann``）が準備済みであること
  （:doc:`tutorial_training` を参照）

ファイル構成
------------

.. code-block:: text

   sample/aenet_pamc/
   ├── do.sh                  # 実行スクリプト
   ├── input.toml             # ODAT-SE 設定ファイル
   ├── predict.in             # AENET 予測設定
   ├── template.xsf           # 構造テンプレート
   ├── plot_histogram.py      # 結果のヒストグラムプロット
   ├── summarize_results.py   # 温度ごとの結果集計
   └── batch_plot.py          # ヒストグラムの一括プロット

input.toml の説明
------------------

.. code-block:: toml

   [base]
   dimension = 1
   output_dir = "output"

パラメータ空間の次元を 1（N-N 結合距離のみ）に設定します。

.. code-block:: toml

   [solver]
   name = "aenet"

   [solver.config]
   aenet_exec_file = "predict.x"
   aenet_ann_potential = "N.5t-5t.ann"

   [solver.param]
   string_list = ["value_01"]

ソルバーとして AENET を使用し、``predict.x`` のパスと ANN ポテンシャルファイルを指定します。
``string_list`` はテンプレート内のプレースホルダーです。

.. code-block:: toml

   [algorithm]
   name = "pamc"
   label_list = ["z"]

   [algorithm.param]
   min_list = [0.5]
   max_list = [2.0]
   unit_list = [0.015]

   [algorithm.pamc]
   numsteps_per_T = 40
   downTsteps_resampling = 2
   Tnum = 41
   Tmin = 1
   Tmax = 1000
   replica_per_prop = 10
   Tlogspace = true
   resampling_bool = true
   use_so_for_resampling = true

PAMC アルゴリズムの設定:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - パラメータ
     - 説明
   * - ``min_list`` / ``max_list``
     - 探索範囲: N-N 距離 0.5〜2.0 Å
   * - ``unit_list``
     - モンテカルロステップの単位変化量: 0.015 Å
   * - ``Tnum``
     - 温度点数: 41
   * - ``Tmin`` / ``Tmax``
     - 温度範囲: 1〜1000 K
   * - ``Tlogspace``
     - 対数スケールで温度を分割
   * - ``numsteps_per_T``
     - 各温度でのモンテカルロステップ数: 40
   * - ``replica_per_prop``
     - 提案あたりのレプリカ数: 10
   * - ``resampling_bool``
     - リサンプリングの有効化

テンプレートファイル
--------------------

``template.xsf`` は N2 二量体の構造テンプレートです:

.. code-block:: text

   ATOMS
   N             0.0000000000        0.0000000000        0.0000000000
   N             0.0000000000        0.0000000000        value_01

``value_01`` が PAMC によって 0.5〜2.0 Å の範囲で探索される N-N 結合距離です。

実行方法
--------

.. code-block:: bash

   cd sample/aenet_pamc
   sh do.sh

または直接実行:

.. code-block:: bash

   odatse-aenet input.toml

MPI 並列実行:

.. code-block:: bash

   mpiexec -np 2 odatse-aenet input.toml

出力
----

計算結果は ``output/`` ディレクトリに生成されます:

- ``output/0/`` , ``output/1/`` , ...: 各 MPI ランクの出力
- 各ランク内に ``ColorMap.txt``, ``BestSteps.txt`` 等が生成
- 作業ディレクトリ ``Log{step}_{iset}/`` にソルバーの入出力ファイル

計算結果
--------

PAMC の結果として、各温度における N-N 結合距離のヒストグラムが得られます。

低温（T → 1 K）では、分布がエネルギー最小点（約 1.1 Å）付近にシャープに集中します。
高温（T → 1000 K）では、広い範囲にわたる平坦な分布となります。

.. note::

   実行後、出力データからヒストグラムを作成することで、
   温度ごとの N-N 結合距離の確率分布を可視化できます。

結果の解析
~~~~~~~~~~

出力ファイルから最適パラメータを確認するには:

.. code-block:: bash

   # 最も低いR-factorを持つステップを確認
   sort -k2 -n output/0/ColorMap.txt | head -5

これにより、エネルギーが最小となる N-N 結合距離（約 1.1 Å）が最適解として得られます。
この値は N2 分子の実験的な結合距離（1.098 Å）とよく一致します。

解析ツール
~~~~~~~~~~

PAMC の出力データを解析・可視化するためのスクリプトが用意されています。

**1. 結果の集計（summarize_results.py）**

各温度ステップごとに ``weight.txt`` から重みデータを集計し、
正規化/非正規化の重みファイルを出力します:

.. code-block:: bash

   python3 summarize_results.py -i input.toml -d my_data --idnum 0

**2. ヒストグラムプロット（plot_histogram.py）**

集計済みファイルから、指定した変数の重み付き 1D ヒストグラムを作成します:

.. code-block:: bash

   python3 plot_histogram.py <集計ファイル> <温度ステップ番号> tau x1

**3. 一括プロット（batch_plot.py）**

全温度ステップのヒストグラムを一括生成します:

.. code-block:: bash

   python3 batch_plot.py plot_histogram.py <集計ディレクトリ> --variable x1
